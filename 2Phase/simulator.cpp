#include "DiscreteProblem.hpp"
#include "Intel_Pardiso.hpp"
#include "NewtonSolver.hpp"
#include "R_Precond_Solver.hpp"
#include "Intel_Prec_GMRES.hpp"
#include "Intel_ILU0.hpp"
#include "fastl/containers/pod_vector_unbounded.hpp"
#include <fstream>
#include <dlib/matrix.h>

//typedef GENSOL::R_Precond_Solver< GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 > >  LINEARSOLVER;
//typedef GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 >  LINEARSOLVER;
typedef GENSOL::Intel_Pardiso                       LINEARSOLVER;
typedef GENSOL::NewtonSolver< DiscreteProblem, LINEARSOLVER > STDN;
typedef GENSOL::CSR_Matrix<double, int>            Jacobian;
typedef DiscreteProblem::Well					  vWell;

//typedef NewtonSolver< DiscreteProblem, LINEARSOLVER > STDN;

class simulator
{
	typedef DiscreteProblem::wellrate   w_rate;
	typedef DiscreteProblem::rate   rate;
private:
	std::size_t NX, NY, NZ;
	double      LX, LY, LZ;
	const std::size_t MAX_NLNITER = 10;
	const double      DT_INIT = 10.0;
	const double      DT_CUT = 0.5;
	const double      DT_GROW = 2.0;
	const double      DT_MAX = 35.0;
	const double OPT_NLNITER = 5;
	std::vector<double> times, pro_time, pro_dt;
	double Qo, Qw;
	dlib::matrix<double> m;
	std::vector<unsigned> ob_si;
	std::vector<double> vPORO;
	std::vector<double> vKX;
	std::vector<double> vKY;
	std::vector<double> vKZ;
	unsigned N_tc;
	unsigned Max_nz;
	unsigned Num_var;
	dlib::matrix<double> up;
	dlib::matrix<double> low;
	typedef struct
	{
		dlib::matrix<double> BHP;
		dlib::matrix<double> qt;
		dlib::matrix<double> qw;
		dlib::matrix<double> qo;
		dlib::matrix<double> wct;
	}CS_gradient;
	typedef struct
	{
		double BHP, qt, qw,qo,wct;
	}CS_value;



	void
		dump_solution(std::string filename, const DiscreteProblem::StateVector & v)
	{
			std::ofstream strm(filename);
			for (std::size_t i = 0; i < v.size(); ++i)
				strm << v[i].Po.value() << "\t"
				<< v[i].Sw.value() << std::endl;;
		}

	void
		dump_field(const char * filename,
		const std::vector<double> &phi,
		const std::vector<double> &kx,
		const std::vector<double> &ky,
		const std::vector<double> &kz)
	{
			std::ofstream strm(filename);
			for (std::size_t i = 0; i < phi.size(); ++i)
				strm << phi[i] << "\t"
				<< kx[i] << "\t"
				<< ky[i] << "\t"
				<< kz[i] << std::endl;
			strm.close();
		}

	void read_reservoir()
	{
		std::ifstream strm("./input/reservoir.dat");
		strm >> NX;
		strm >> NY;
		strm >> NZ;
		strm >> LX;
		strm >> LY;
		strm >> LZ;
	}
	void read_schedule(std::vector<double> &times)
	{
		double tmp;
		std::ifstream strm("./input/schedule.dat");
		while (strm>>tmp)
		{
			times.push_back(tmp);
		}
	}
	void m_to_vector(dlib::matrix<double> &si,
		std::vector<double> &phi,
		std::vector<double> &kx,
		std::vector<double> &ky,
		std::vector<double> &kz)
	{
		size_t n = si.size();
		size_t nc = n / 4;
		phi.resize(NX*NY*NZ);
		kx.resize(NX*NY*NZ);
		ky.resize(NX*NY*NZ);
		kz.resize(NX*NY*NZ);
		for (unsigned i = 0; i < nc; i++)
		{
			kx[i] = exp(si(i, 0));
			ky[i] = exp(si(i, 0));
			kz[i] = exp(si(i, 0));
			phi[i] = si(i + 3*nc, 0);
		}
	}

	void load_schedule(std::vector<vWell> &wells, dlib::matrix<double> &s,
		std::vector<double> &times, std::vector<double> &schedule, unsigned Np)
	{
		//check_schedule_size(wells.size(), times.size(), s.size());
		unsigned Nj = wells.size() - Np;
		for (unsigned i = 0; i < Np;i++)
		{
			wells[i].Pbh.resize(times.size());
			for (unsigned j = 0; j < times.size();j++)
			{
				wells[i].Pbh[j] = s(i*times.size()+j);
			}
		}
		for (unsigned i = Np; i < wells.size(); i++)
		{
			wells[i].QINJ.resize(times.size());
			for (unsigned j = 0; j < times.size(); j++)
			{
				wells[i].QINJ[j] = s(i*times.size() + j);
			}
		}
		schedule.resize(times.size());
		for (unsigned j = 0; j < times.size(); j++)
		{
			schedule[j] = times[j];
		}
	}

	void scaling(const dlib::matrix<double> &x)
	{
		dlib::matrix<double> & xn = const_cast<dlib::matrix<double>&> (x);
		xn = dlib::pointwise_multiply((x - low), 1 / (up - low));
	}
	void scaling_back(const dlib::matrix<double> &x)
	{
		dlib::matrix<double> & xn = const_cast<dlib::matrix<double>&> (x);
		xn = dlib::pointwise_multiply(x, (up - low)) + low;
	}


	void
		read_from_file_unsigned(std::string filename, dlib::matrix<unsigned> &v)
	{
			double tmp;
			size_t n = 0;
			std::ifstream strm(filename);
			while (strm >> tmp)
			{
				n++;
			}
			v.set_size(n, 1);
			strm.close();
			strm.open(filename);
			for (unsigned i = 0; i < n; i++)
			{
				strm >> v(i, 0);
			}
		}



	void
		read_from_file(std::string filename, dlib::matrix<double> &v)
	{
			double tmp;
			size_t n = 0;
			std::ifstream strm(filename);
			while (strm >> tmp)
			{
				n++;
			}
			v.set_size(n, 1);
			strm.close();
			strm.open(filename);
			for (unsigned i = 0; i < n; i++)
			{
				strm >> v(i, 0);
			}
		}


	double comput_NPV(std::vector<bool> isproducer,std::vector<double> &pro_time, std::vector<double> &pro_dt,
		std::vector<std::vector<rate>> &well_prod, std::vector<std::vector<CS_gradient>> &all_grad,dlib::matrix<double> &NPV_grad)
	{
		double ro,  rw,  rinj, b;
		std::ifstream input("./input/economic.dat");
		input >> ro;
		input >> rw;
		input >> rinj;
		input >> b;
		std::ofstream output("./output/NPV.out");
		std::ofstream output_rate("./output/rate.dat");
		double total = 0;
		NPV_grad.set_size(Num_var, 1);
		NPV_grad = 0;
		unsigned nw = well_prod[0].size();
		for (unsigned t = 0; t < pro_time.size();t++)
		{
			double tmp = 0;
			dlib::matrix<double> t_grad(Num_var,0);
			for (unsigned w = 0; w < nw;w++)
			{
				if (isproducer[w])
				{
					tmp += (ro*well_prod[t][w].qo - rw*well_prod[t][w].qw);
					t_grad += (ro*all_grad[t][w].qo - rw*all_grad[t][w].qw);
					output_rate << well_prod[t][w].qo << "   " << well_prod[t][w].qw << "   ";
				}
				else
				{
					tmp -= (rinj*well_prod[t][w].qw);
					t_grad -= (rinj*all_grad[t][w].qw);
					output_rate << well_prod[t][w].qw <<"   ";
				}
			}
			output_rate << std::endl;
			total += tmp*pro_dt[t] / pow((1 + b), pro_time[t] / 365.0);
			NPV_grad += t_grad*pro_dt[t] / pow((1 + b), pro_time[t] / 365.0);
		}
		output << total;
		output.close();
		return total;
	}


	double comput_NPV(std::vector<bool> isproducer, std::vector<double> &pro_time, std::vector<double> &pro_dt,
		std::vector<std::vector<rate>> &well_prod)
	{
		double ro, rw, rinj, b;
		std::ifstream input("./input/economic.dat");
		input >> ro;
		input >> rw;
		input >> rinj;
		input >> b;
		std::ofstream output("./output/NPV.out");
		std::ofstream output_rate("./output/rate.dat");
		double total = 0;
		unsigned nw = well_prod[0].size();
		for (unsigned t = 0; t < pro_time.size(); t++)
		{
			double tmp = 0;
			dlib::matrix<double> t_grad(Num_var, 0);
			for (unsigned w = 0; w < nw; w++)
			{
				if (isproducer[w])
				{
					tmp += (ro*well_prod[t][w].qo - rw*well_prod[t][w].qw);
					output_rate << well_prod[t][w].qo << "   " << well_prod[t][w].qw << "   ";
				}
				else
				{
					tmp -= (rinj*well_prod[t][w].qw);
					output_rate << well_prod[t][w].qw << "   ";
				}
			}
			output_rate << std::endl;
			total += tmp*pro_dt[t] / pow((1 + b), pro_time[t] / 365.0);
		}
		output << total;
		output.close();
		return total;
	}





	void vector_to_d(dlib::matrix<double> &d,
		std::vector<std::vector<double>> &qo,
		std::vector<std::vector<double>> &qw,
		std::vector<std::vector<double>> &pwf,
		std::size_t Nob)
	{
		std::size_t Np = qo.size(), Nj = pwf.size();
		std::size_t Nd = (2 * Np + Nj)*Nob;
		d.set_size(Nd, 1);
		for (unsigned i = 0; i < Np; i++)
		{
			for (unsigned j = 0; j < Nob; j++)
			{
				d(i*Nob + j, 0) = qo[i][j];
				d(Nob*(Np + i) + j, 0) = qw[i][j];
			}
		}
		for (unsigned i = 0; i < Nj; i++)
		{
			for (unsigned j = 0; j < Nob; j++)
			{
				//std::cout << 2 * Nob*Np + Nob*i + j << std::endl;
				d(2 * Nob*Np + Nob*i + j, 0) = pwf[i][j];
			}
		}
	}


	void dump_prod_rate(std::vector<double> &pro_time, 
		std::vector<std::vector<w_rate>> &well_prod)
	{
		std::ofstream output("./output/well_oil_rate.dat");
		std::ofstream output1("./output/well_water_rate.dat");
		//std::ofstream output2("./output/well_inject_rate.dat");
		unsigned nw = well_prod[0].size();
		for (unsigned t = 0; t < pro_time.size(); t++)
		{
			output << pro_time[t] << " ";
			for (unsigned w = 0; w < nw; w++)
			{
				output << well_prod[t][w].qo << " ";
			}
			output << std::endl;
		}
		for (unsigned t = 0; t < pro_time.size(); t++)
		{
			output1 << pro_time[t] << " ";
			for (unsigned w = 0; w < nw; w++)
			{
				output1 << well_prod[t][w].qw << " ";
			}
			output1 << std::endl;
		}
		//for (unsigned t = 0; t < pro_time.size(); t++)
		//{
		//	output2 << pro_time[t] << " ";
		//	for (unsigned w = 0; w < nw; w++)
		//	{
		//		output2 << well_prod[t][w].qinj << " ";
		//	}
		//	output2 << std::endl;
		//}
	}


public:
	std::vector<std::vector<CS_gradient>> CS_grad;
	std::vector<std::vector<CS_value>> CS_val;

	void sum_qo_qw(std::vector<double> &pro_time, std::vector<double> &pro_dt,
		std::vector<std::vector<w_rate>> &well_prod)
	{
		unsigned nw = well_prod[0].size();
		Qo = 0, Qw = 0;
		for (unsigned t = 0; t < pro_time.size(); t++)
		{
			double tmpo=0,tmpw = 0;
			for (unsigned w = 0; w < nw; w++)
			{
				tmpo+= well_prod[t][w].qo.value();
				tmpw+= well_prod[t][w].qw.value();
			}
			Qo+=tmpo*pro_dt[t] ;
			Qw+=tmpw*pro_dt[t];
		}
	}

	void output_classifer(std::vector<double> &pro_time, std::vector<double> &pro_dt,std::vector<double> bbt,
		std::vector<std::vector<w_rate>> &well_prod,int model_num)
	{
		unsigned np = bbt.size();
		std::vector<double> qo(np,0), qw(np,0);
		//std::cout << np <<std::endl;
		for (unsigned w = 0; w < np; w++)
		{
			for (unsigned t = 0; t < pro_time.size(); t++)
			{
				qo[w] += well_prod[t][w].qo.value()*pro_dt[t];
				qw[w] += well_prod[t][w].qw.value()*pro_dt[t];
			}
		}
		std::string filename;
		filename = "./output/classifer/model" + std::to_string(model_num) + ".out";
		std::ofstream out_claasifer(filename);
		for (unsigned i = 0; i < np;i++)
		{
			out_claasifer << qo[i] << std::endl << qw[i] << std::endl << bbt[i] << std::endl;
		}
	}
	void record_pressure(const DiscreteProblem::StateVector & v, std::vector<unsigned> &idx, std::vector<std::vector<double>> &pwf)
	{
		for (unsigned i = 0; i < idx.size(); i++)
		{
			pwf[i].push_back(v[idx[i]].Po.value());
		}
	}

	std::vector<double> transpose_prod(DiscreteProblem::Sparse_matrix &J_m,std::vector<double> &ajv,unsigned nwc)
	{
		std::vector<double> result(nwc,0);
		for (unsigned i = 0; i < J_m.val.size();i++)
		{
			result[J_m.col[i]]+= J_m.val[i] * ajv[J_m.row[i]];
		}
		return result;
	}
	std::vector<double> transpose_prod(DiscreteProblem::Sparse_matrix &J_m, std::vector<double> &ajv, unsigned nwc,unsigned nrh)
	{
		std::vector<double> result(nwc*nrh, 0);
		unsigned ny = ajv.size() / nrh;
		for (unsigned j = 0; j <nrh ;j++)
		{
			for (unsigned i = 0; i < J_m.val.size(); i++)
			{
				//std::cout << J_m.col[i] + j*nwc << "  " << J_m.row[i] + j*nwc <<"  "<< J_m.val[i] << std::endl;
				result[J_m.col[i] + j*nwc] += J_m.val[i] * ajv[J_m.row[i] + j*ny];
			}
		}

		return result;
	}

	void compute_well_cs(const std::vector<DiscreteProblem::CS> &cons,
		std::vector<unsigned> &w_cs,unsigned &nrh)
	{
		nrh = 0;
		w_cs.resize(cons.size());
		//count the number of constrain for each well
		for (unsigned w = 0; w < cons.size(); w++)
		{
			if (cons[w].BHP_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].qt_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].qo_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].qw_i)
			{
				w_cs[w]++;
				nrh++;
			}
			if (cons[w].wct_i)
			{
				w_cs[w]++;
				nrh++;
			}
		}
	}
	void conglo_v(std::vector<double> &ob, const std::vector<DiscreteProblem::CS> &cons, 
		const unsigned nv,const unsigned &nrh)
	{
		unsigned id=0;
		double temp;
		//////////////////////////////////////////////////////////////////////////////
		ob.resize(nrh*nv);
		for (unsigned w = 0; w < cons.size(); w++)
		{
			if (cons[w].BHP_i)
			{
				std::vector<double> x(nv,0);
				cons[w].BHP.extract(temp,x.data(),nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].qt_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].qt.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].qo_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].qo.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].qw_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].qw.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
			if (cons[w].wct_i)
			{
				std::vector<double> x(nv, 0);
				cons[w].wct.extract(temp, x.data(), nv);
				for (unsigned i = 0; i < nv; i++)	ob[i + id*nv] = x[i];
				id++;
			}
		}
	}
	void assign_gradient( std::vector<DiscreteProblem::CS> &cons,
		std::vector<CS_gradient> &grad, std::vector<double> &sum,unsigned nv)
	{
		unsigned id = 0;
		grad.resize(cons.size());
		for (unsigned w = 0; w < cons.size(); w++)
		{
			if (cons[w].BHP_i)
			{
				grad[w].BHP.set_size(nv,1);
				grad[w].BHP = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].BHP(i,0) = sum[i+id*nv];
				id++;
			}
			if (cons[w].qt_i)
			{
				grad[w].qt.set_size(nv,1);
				grad[w].qt = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].qt(i,0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].qo_i)
			{
				grad[w].qo.set_size(nv,1);
				grad[w].qo = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].qo(i,0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].qw_i)
			{
				grad[w].qw.set_size(nv,1);
				grad[w].qw = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].qw(i,0) = sum[i + id*nv];
				id++;
			}
			if (cons[w].wct_i)
			{
				grad[w].wct.set_size(nv,1);
				grad[w].wct = 0;
				for (unsigned i = 0; i < nv; i++) grad[w].wct(i,0) = sum[i + id*nv];
				id++;
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Backward_CS_gradient(std::vector<unsigned> & time,
		std::vector<CSR> & al_j,
		std::vector<DiscreteProblem::Sparse_matrix> & ol_j,
		std::vector<DiscreteProblem::Sparse_matrix> & al_jm,
		std::vector<DiscreteProblem::CS> &cons,
		std::vector<DiscreteProblem::CS> &cons_m,
		std::vector<CS_gradient> &grad)
	{
		unsigned nrh;
		unsigned nos = time.back() - 1;
		std::vector<unsigned> w_cs;
		compute_well_cs(cons, w_cs, nrh);
		LINEARSOLVER lnsolve(N_tc, Max_nz,nrh);
		std::size_t n = al_j[0].N(), nv = Num_var;
		std::vector<double> ajv(n*nrh, 0), tmp(nv*nrh), tmp1(n*nrh), sum(nv*nrh, 0);
		std::vector<double> oby, obm;
		// frist time to compute ajoint vector backward;

		conglo_v(oby, cons, n, nrh);
		conglo_v(obm, cons_m, nv, nrh);

		//for (unsigned i = 0; i < oby.size();i++)
		//{
		//	std::cout<<i<<"   " << oby[i] << std::endl;
		//}
		//frist run;
		lnsolve.solvet(al_j[nos], ajv, oby);
		tmp = transpose_prod(al_jm[nos], ajv, nv,nrh);
		////////////////////////// can simplify//////////////////////////////////////
		for (unsigned i = 0; i < nv*nrh; i++)
		{
			sum[i] += tmp[i];
		}
		//////////////////////////////////////////////////////////////////////////////
		//compute ajoint vector for each simulation time
		for (unsigned j = nos; j >0; j--)
		{
			tmp1 = transpose_prod(ol_j[j - 1], ajv, n,nrh);
			lnsolve.solvet(al_j[j - 1], ajv, tmp1);
			tmp = transpose_prod(al_jm[j - 1], ajv, nv,nrh);
			for (unsigned i = 0; i < nv*nrh; i++)
			{
				sum[i] += tmp[i];
			}
		}
		//summation of ajoint vecotor product of Fm match at time
		//////////////////////////////////////////////////////////////////////////////
		for (unsigned i = 0; i < nv*nrh; i++)
		{
			sum[i] -= obm[i];
		}
		//////////////////////////////////////////////////////////////////////////////
		assign_gradient(cons, grad, sum, nv);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	~simulator()
	{}

	simulator()
	{
		read_reservoir();
		read_schedule(times);

	}

	std::vector<int> chop_vector(std::size_t n, std::vector<unsigned> &time)
	{
		std::vector<int> choped(n + 1);

		for (unsigned i = 0; i <= n; i++)
		{
			choped[i] = time[i];
		}
		return choped;
	}

	void output_cs_gradient(unsigned N_var,unsigned N_cs,dlib::matrix<double> &NPV_grad,std::vector<std::vector<CS_gradient>> &all_grad)
	{

		dlib::matrix<double> cs_gradient(N_cs,N_var);
		cs_gradient = 0;
		unsigned i = 0;
		for (unsigned w = 0; w < all_grad[0].size(); w++)
		{
			for (unsigned t = 0; t < all_grad.size(); t++)
			{
				if (all_grad[t][w].BHP.nr()!=0)
				{
					set_rowm(cs_gradient, i) = trans(all_grad[t][w].BHP);
					i++;
				}
				if (all_grad[t][w].qt.nr() != 0)
				{
					set_rowm(cs_gradient, i) = trans(all_grad[t][w].qt);
					i++;
				}
				if (all_grad[t][w].qo.nr() != 0)
				{
					set_rowm(cs_gradient, i) = trans(all_grad[t][w].qo);
					i++;
				}
				if (all_grad[t][w].qw.nr() != 0)
				{
					set_rowm(cs_gradient, i) = trans(all_grad[t][w].qw);
					i++;
				}
				if (all_grad[t][w].wct.nr() != 0)
				{
					set_rowm(cs_gradient, i) = all_grad[t][w].wct;
					i++;
				}
			}
		}
		std::ofstream output("output/cs_gradient.dat");
		std::ofstream outputNPV("output/NPV_gradient.dat");
		output << cs_gradient << std::endl;
		outputNPV << NPV_grad << std::endl;
		output.close();
	}

	void initial_properties(dlib::matrix<double> &m){
		m_to_vector(m, vPORO, vKX, vKY, vKZ); //assign field properties
	}
	void get_cs_val(std::vector<std::vector<DiscreteProblem::CS>> &cons)
	{
		CS_val.resize(cons.size());
		for (unsigned i = 0; i < cons.size();i++){
			CS_val[i].resize(cons[i].size());
			for (unsigned j = 0; j < cons[i].size(); j++){
				CS_val[i][j].BHP = cons[i][j].BHP.value();
				CS_val[i][j].qt = cons[i][j].qt.value();
				CS_val[i][j].qo = cons[i][j].qo.value();
				CS_val[i][j].qw = cons[i][j].qw.value();
				CS_val[i][j].wct = cons[i][j].wct.value();
			}
		}
	}


	double run(const dlib::matrix<double> &s)
	{
		double NPV_v;
		scaling_back(s);
		std::vector<std::vector<rate>> well_prod;
		DiscreteProblem model(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
		model.loadwell_data(times, s);
		N_tc = model.N_tc;
		Max_nz = model.max_num_nnz();
		LINEARSOLVER lnsolver(model.N_tc, model.max_num_nnz());
		STDN	newton(model, MAX_NLNITER, 1);
		DiscreteProblem::StateVector uOld, uNew, aOld, aNew;
		model.initialize_state(uOld);
		double DT = DT_INIT;
		double time = 0.0;
		unsigned j = 0, ob_t = 0;
//////////////////////////////////////////////////////////////////////////////////////
		for (std::size_t i = 0; i < times.size(); i++)     //begin of simulation
		{
			while (time < times[i])
			{
				uNew = uOld;
				model.sum_t = time;
				STDN::report_t stdsmry = newton.solve_timestep(uNew, uOld, DT, model, lnsolver);
				if (stdsmry.is_converged)
				{
					uOld = uNew;
					time += DT;
					if (stdsmry.niter < OPT_NLNITER) DT *= DT_GROW;
					DT = std::min(DT, (times[i] - time));
					DT == 0 ? DT = DT_INIT : DT = DT;			//cutting time step to fit schedule;
					j++;
				}
				else
				{
					DT *= DT_CUT;
					//std::cout << "FAILED " << std::endl;
				}

			}
			model.extract_obj_der(i);	//extract derivetive for objective function get production rate
			pro_time.push_back(time);
			well_prod.push_back(model.well_rate);
			pro_dt.push_back(DT);
		}
		get_cs_val(model.H_constrain);
		NPV_v = comput_NPV(model.is_producer, pro_time, pro_dt, well_prod);
		scaling(s);
		return -NPV_v;
	}

	const dlib::matrix<double> grad(const dlib::matrix<double>& s)
	{
		Num_var = s.nr();
		CS_grad.clear();
		scaling_back(s);
		DiscreteProblem model(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
		DiscreteProblem model_m(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
		model.loadwell_data(times, s);
		model_m.loadwell_independent(times, s);
		N_tc = model.N_tc;
		Max_nz = model.max_num_nnz();
		LINEARSOLVER lnsolver(model.N_tc, model.max_num_nnz());
		STDN	newton(model, MAX_NLNITER, 1);
		DiscreteProblem::StateVector uOld, uNew, aOld, aNew;
		model.initialize_state(uOld);
		aOld.resize(uOld.size());
		aNew.resize(uOld.size());
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::vector<std::vector<rate>> well_prod;
		CSR		Jaco(model.N_tc, model.max_num_nnz());
		DiscreteProblem::Sparse_matrix		o_Jaco;
		std::vector<DiscreteProblem::Sparse_matrix> F_v;
		std::vector<DiscreteProblem::Sparse_matrix> F_yold;
		std::vector<CSR> F_y;
		DiscreteProblem::Sparse_matrix	j_v;
		std::vector<CS_gradient> grad;
		dlib::matrix<double> NPV_grad;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double DT = DT_INIT;
		double time = 0.0;
		unsigned j = 0, ob_t = 0;
		//begin of simulation
		for (std::size_t i = 0; i < times.size(); i++){
			while (time < times[i]){
				uNew = uOld;
				model.sum_t = time;
				model_m.sum_t = time;
				STDN::report_t stdsmry = newton.solve_timestep_new(uNew, uOld, DT, model, lnsolver, Jaco, o_Jaco);  //one simulation step
				if (stdsmry.is_converged){
					//////////////////////////////////////////////////////record jacobian//////////////////////////////////////////////////////////////////////
					model_m.initialize_state(uOld, uNew, aOld, aNew);// assign state to model_m;
					model_m.extract_F_der(aOld, aNew, DT, j_v);
					F_v.push_back(j_v);
					F_y.push_back(Jaco);
					F_yold.push_back(o_Jaco);
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					uOld = uNew;
					time += DT;
					if (stdsmry.niter < OPT_NLNITER) DT *= DT_GROW;
					DT = std::min(DT, (times[i] - time));
					DT == 0 ? DT = DT_INIT : DT = DT;			//cutting time step to fit schedule;
					j++;
				}
				else{
					DT *= DT_CUT;
				}
			}
			model.extract_obj_der(i);	//extract derivetive for objective function
			model_m.extract_obj_der(i);
			pro_time.push_back(time);
			well_prod.push_back(model.well_rate);
			pro_dt.push_back(DT);
			ob_si.push_back(j);
			Backward_CS_gradient(ob_si, F_y, F_yold, F_v, model.H_constrain[ob_t], model_m.H_constrain[ob_t], grad);
			ob_t++;
			CS_grad.push_back(grad);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
		get_cs_val(model.H_constrain);
		double NPV_value = comput_NPV(model.is_producer, pro_time, pro_dt, well_prod, CS_grad, NPV_grad);
		pro_time.clear();
		pro_dt.clear();
		ob_si.clear();
		scaling(s);
		return NPV_grad;
	}
	





	double run(dlib::matrix<double> &m, dlib::matrix<double> &s, dlib::matrix<double> &NPV_grad, double &NPV_v, bool output_i)
	{
		Num_var = s.nr();
		////////////////////////////////////////////////////////////////////////////////////////////////////
		m_to_vector(m, vPORO, vKX, vKY, vKZ); //assign field properties
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		DiscreteProblem model(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
		DiscreteProblem model_m(NX, NY, NZ, LX, LY, LZ, vKX, vKY, vKZ, vPORO);
		model.loadwell_data(times, s);
		model_m.loadwell_independent(times, s);
		N_tc = model.N_tc;
		Max_nz = model.max_num_nnz();
		LINEARSOLVER lnsolver(model.N_tc, model.max_num_nnz());
		STDN	newton(model, MAX_NLNITER, 1);
		DiscreteProblem::StateVector uOld, uNew, aOld, aNew;
		model.initialize_state(uOld);
		aOld.resize(uOld.size());
		aNew.resize(uOld.size());
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::vector<std::vector<rate>> well_prod;
		CSR		Jaco(model.N_tc, model.max_num_nnz());
		DiscreteProblem::Sparse_matrix		o_Jaco;
		std::vector<DiscreteProblem::Sparse_matrix> F_v;
		std::vector<DiscreteProblem::Sparse_matrix> F_yold;
		std::vector<CSR> F_y;
		DiscreteProblem::Sparse_matrix	j_v;
		std::vector<std::vector<CS_gradient>> all_grad;
		std::vector<CS_gradient> grad;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double DT = DT_INIT;
		double time = 0.0;
		unsigned j = 0, ob_t = 0;
		for (std::size_t i = 0; i < times.size(); i++)     //begin of simulation
		{
			while (time < times[i])
			{
				uNew = uOld;
				model.sum_t = time;
				model_m.sum_t = time;
				STDN::report_t stdsmry = newton.solve_timestep_new(uNew, uOld, DT, model, lnsolver, Jaco, o_Jaco);  //one simulation step
				if (stdsmry.is_converged)
				{
					//////////////////////////////////////////////////////record jacobian//////////////////////////////////////////////////////////////////////
					model_m.initialize_state(uOld, uNew, aOld, aNew);// assign state to model_m;
					model_m.extract_F_der(aOld, aNew, DT, j_v);
					F_v.push_back(j_v);
					F_y.push_back(Jaco);
					F_yold.push_back(o_Jaco);
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					uOld = uNew;
					time += DT;
					if (stdsmry.niter < OPT_NLNITER) DT *= DT_GROW;
					DT = std::min(DT, (times[i] - time));
					DT == 0 ? DT = DT_INIT : DT = DT;			//cutting time step to fit schedule;
					j++;
				}
				else
				{
					DT *= DT_CUT;
				}
			}
			model.extract_obj_der(i);	//extract derivetive for objective function
			model_m.extract_obj_der(i);
			pro_time.push_back(time);
			well_prod.push_back(model.well_rate);
			pro_dt.push_back(DT);
			ob_si.push_back(j);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			Backward_CS_gradient(ob_si, F_y, F_yold, F_v, model.H_constrain[ob_t], model_m.H_constrain[ob_t], grad);
			ob_t++;
			all_grad.push_back(grad);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (output_i)
			{
				std::string filename = "./output/time_field/results";
				std::string fileend = ".out";
				int t_out = int(time);
				filename = filename + std::to_string(t_out) + fileend;
				dump_solution(filename, uNew);
			}
		}
		NPV_v = comput_NPV(model.is_producer, pro_time, pro_dt, well_prod, all_grad, NPV_grad);
		if (output_i)
		{
			dump_solution("./output/results.out", uNew);
			dump_field("./output/field.out", vPORO, vKX, vKY, vKZ);
			output_cs_gradient(Num_var, model.Num_cs, NPV_grad, all_grad);
			//dump_prod_rate(pro_time, well_prod);
			system("pause");
		}
		pro_time.clear();
		pro_dt.clear();
		ob_si.clear();
		return NPV_v;
	}




};