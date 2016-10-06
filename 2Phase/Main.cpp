#include <dlib/matrix.h>
#include <dlib/optimization.h>
#include <fstream>
#include "simulator.cpp"
#include <cmath>
#include <random>
#include <time.h>
#include "optimize.h"
#include <algorithm> 

unsigned ne, nw, np, nc,nm, MAX_step;
double step,larg_step; 
dlib::matrix<double> m;
dlib::matrix<double> up(nc*nw, 1);
dlib::matrix<double> low(nc*nw, 1);

void
read_from_file(const char * filename, dlib::matrix<double> &v)
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
read_well_schedule(dlib::matrix<double> &v)
{
	std::ifstream strm("input/wellschedule.dat");
	unsigned nwc;
	strm >> nwc;
	v.set_size(nwc,0);
	for (unsigned i = 0; i < nwc;i++)
	{
		strm >> v(i);
	}

}
void read_parameter(unsigned &ne, unsigned &nw, unsigned &np, unsigned &nc,
	dlib::matrix<double> &up, dlib::matrix<double> &low)
{
	ne = 0; nw = 0; nc = 0;
	
	std::ifstream strm("./input/parameter.dat");
	double qmax, qmin, pmax, pmin;
	strm >> ne;
	strm >> nw;
	strm >> np;
	strm >> nc;
	strm >> qmin;
	strm >> qmax;
	strm >> pmin;
	strm >> pmax;
	strm >> MAX_step;
	strm >> step;
	strm >> larg_step;
	up.set_size(nc*nw,1);
	low.set_size(nc*nw, 1);
	dlib::set_subm(up, dlib::rectangle(0, 0, nc*np - 1, 0))=pmax;
	dlib::set_subm(up, dlib::rectangle(nc*np, 0, nc*nw - 1, 0)) = qmax;

	dlib::set_subm(low, dlib::rectangle(0, 0, nc*np - 1, 0)) = pmin;
	dlib::set_subm(low, dlib::rectangle(nc*np, 0, nc*nw - 1, 0)) = qmin;
}

void read_parameter(unsigned &ne, unsigned &nw, unsigned &np, unsigned &nc,unsigned &nm,
	dlib::matrix<double> &up, dlib::matrix<double> &low)
{
	ne = 0; nw = 0; nc = 0;

	std::ifstream strm("./input/parameter.dat");
	double qmax, qmin, pmax, pmin;
	strm >> ne;
	strm >> nw;
	strm >> np;
	strm >> nc;
	strm >> qmin;
	strm >> qmax;
	strm >> pmin;
	strm >> pmax;
	strm >> MAX_step;
	strm >> step;
	strm >> larg_step;
	strm >> nm;
	up.set_size(nc*nw, 1);
	low.set_size(nc*nw, 1);
	dlib::set_subm(up, dlib::rectangle(0, 0, nc*np - 1, 0)) = pmax;
	dlib::set_subm(up, dlib::rectangle(nc*np, 0, nc*nw - 1, 0)) = qmax;

	dlib::set_subm(low, dlib::rectangle(0, 0, nc*np - 1, 0)) = pmin;
	dlib::set_subm(low, dlib::rectangle(nc*np, 0, nc*nw - 1, 0)) = qmin;
}
double spherecov(double delt, double h)
{
	if (h<1)
	{
		return delt*delt*(1 - 3.0 / 2.0 * h + pow(h, 3) / 2.0);
	}
	else
	{
		return 0;
	}
}

dlib::matrix<double> generate_Cx(unsigned nw, unsigned nc, double control_dis,double delta)
{
	double tmp (0),rho(0);
	dlib::matrix<double> Cx(nw*nc,nw*nc);
	Cx = 0;
	for (unsigned n = 0; n < nc;n++)
	{
		for (unsigned c = 0; c < nc;c++)
		{
			tmp = abs(double(n - c)) / control_dis;
			rho = spherecov(delta,tmp);
			for (unsigned w = 0; w < nw; w++)
			{
				Cx(n+w*nc, c+w*nc) = rho;
				Cx(c + w*nc, n + w*nc) = rho;
			}
		}
	}

	//output << Cx << std::endl;
	return Cx;
}

void scaling(const dlib::matrix<double> &x)
{
	dlib::matrix<double> & xn = const_cast<dlib::matrix<double>&> (x);
	xn = dlib::pointwise_multiply((x - low),1/(up - low));
}
void scaling_back(const dlib::matrix<double> &x)
{
	dlib::matrix<double> & xn = const_cast<dlib::matrix<double>&> (x);
	xn = dlib::pointwise_multiply(x, (up - low)) + low;
}
void scaling_back(dlib::matrix<double> &x)
{
	for (unsigned i = 0; i<nc*nw; i++)
	{
		if (x(i, 0)>1)
		{
			x(i, 0) = 1;
		}
		else if (x(i, 0)<0)
		{
			x(i, 0) = 0;
		}
	}
	x = dlib::pointwise_multiply(x, (up - low)) + low;
}

void scaling_uc(dlib::matrix<double> &x)
{
	for (unsigned i = 0; i < nw*nc; i++)
	{
		if (x(i, 0)<=low(i, 0))
		{
			x(i,0) = -7.0;
		}
		else if (x(i, 0)>=up(i,0))
		{
			x(i, 0) = 7.0;
		}
		else
		{
			x(i, 0) = log((x(i, 0) - low(i, 0)) / (up(i, 0) - x(i, 0)));
		}
	}
	//dlib::pointwise_multiply((x(i, 0) - low(i, 0)), 1 / (up(i, 0) - x(i, 0)
}

void scaling_uc_back(dlib::matrix<double> &x)
{
	for (unsigned i = 0; i < nw*nc;i++){
		//std::cout << x(i, 0) << std::endl;
		if (x(i,0)<0){
			x(i, 0) = (up(i, 0)*exp(x(i, 0)) + low(i, 0)) / (1 + exp(x(i, 0)));
		}
		else{
			x(i, 0) = (low(i, 0)*exp(-x(i, 0)) + up(i, 0)) / (1 + exp(-x(i, 0)));
		}
	}
}

dlib::matrix<double> CS, CS_val, CS_range, val, lamda, sign_val;
std::vector<unsigned> cs_locx, cs_locy,cs_locw;
std::vector<bool> CS_sign;

simulator sim_model;


void load_constrain()
{
	std::ifstream input("input/constrain_inform.dat");
	unsigned ncs;
	bool is_larger;
	double val(0);
	input >> ncs;
	CS_val.set_size(ncs, 1);
	CS_range.set_size(ncs,1);
	cs_locx.resize(ncs);
	cs_locy.resize(ncs);
	cs_locw.resize(ncs);
	CS_sign.resize(ncs);
	for (unsigned i = 0; i < ncs;i++)	{
		input >> cs_locw[i];
		input >> cs_locx[i];
		input >> cs_locy[i];
		input >> is_larger;
		CS_sign[i]=is_larger;
		input >> val;
		if (is_larger){
			val = -val;
		}
		CS_val(i,0)=val;
		input >> CS_range(i, 0);
	}
}

void assemble_cs_gradient(std::vector<dlib::matrix<double>> &CS_gradient)
{
	unsigned ncs = CS_sign.size();
	CS_gradient.resize(ncs);
	for (unsigned i = 0; i < ncs;i++)	{
		if (cs_locy[i]==0){
			if (CS_sign[i]){
				set_colm(CS_gradient[i], 0) = sim_model.CS_grad[cs_locx[i]][cs_locw[i]].BHP;
			}
			else{
				set_colm(CS_gradient[i], 0) = -sim_model.CS_grad[cs_locx[i]][cs_locw[i]].BHP;
			}
		}
		else if (cs_locy[i] == 1){
			if (CS_sign[i]){
				set_colm(CS_gradient[i], 0) = sim_model.CS_grad[cs_locx[i]][cs_locw[i]].qt;
			}
			else{
				set_colm(CS_gradient[i], 0) = -sim_model.CS_grad[cs_locx[i]][cs_locw[i]].qt;
			}
		}
		else if (cs_locy[i] == 2){
			if (CS_sign[i]){
				set_colm(CS_gradient[i], 0) = sim_model.CS_grad[cs_locx[i]][cs_locw[i]].qo;
			}
			else{
				set_colm(CS_gradient[i], 0) = -sim_model.CS_grad[cs_locx[i]][cs_locw[i]].qo;
			}
		}
		else if (cs_locy[i] == 3) {
			if (CS_sign[i]){
				set_colm(CS_gradient[i], 0) = sim_model.CS_grad[cs_locx[i]][cs_locw[i]].qw;
			}
			else{
				set_colm(CS_gradient[i], 0) = -sim_model.CS_grad[cs_locx[i]][cs_locw[i]].qw;
			}
		}
		else if (cs_locy[i] == 4) {
			if (CS_sign[i]){
				set_colm(CS_gradient[i], 0) = sim_model.CS_grad[cs_locx[i]][cs_locw[i]].wct;
			}
			else{
				set_colm(CS_gradient[i], 0) = -sim_model.CS_grad[cs_locx[i]][cs_locw[i]].wct;
			}
		}
		else{
			std::cout << "Error in constrain type" << std::endl;
		}
	}
}

void assemble_cs(dlib::matrix<double> &cs_v)
{
	unsigned ncs = CS_sign.size();
	cs_v.set_size(ncs,1);
	for (unsigned i = 0; i < ncs; i++)	{
		if (cs_locy[i] == 0){
			if (CS_sign[i]){
				cs_v(i) = sim_model.CS_val[cs_locx[i]][cs_locw[i]].BHP;
			}
			else{
				cs_v(i) = -sim_model.CS_val[cs_locx[i]][cs_locw[i]].BHP;
			}
		}
		else if (cs_locy[i] == 1){
			if (CS_sign[i]){
				cs_v(i) = sim_model.CS_val[cs_locx[i]][cs_locw[i]].qt;
			}
			else{
				cs_v(i) = -sim_model.CS_val[cs_locx[i]][cs_locw[i]].qt;
			}
		}
		else if (cs_locy[i] == 2){
			if (CS_sign[i]){
				cs_v(i) = sim_model.CS_val[cs_locx[i]][cs_locw[i]].qo;
			}
			else{
				cs_v(i) = -sim_model.CS_val[cs_locx[i]][cs_locw[i]].qo;
			}
		}
		else if (cs_locy[i] == 3) {
			if (CS_sign[i]){
				cs_v(i) = sim_model.CS_val[cs_locx[i]][cs_locw[i]].qw;
			}
			else{
				cs_v(i) = -sim_model.CS_val[cs_locx[i]][cs_locw[i]].qw;
			}
		}
		else if (cs_locy[i] == 4) {
			if (CS_sign[i]){
				cs_v(i) = sim_model.CS_val[cs_locx[i]][cs_locw[i]].wct;
			}
			else{
				cs_v(i) = -sim_model.CS_val[cs_locx[i]][cs_locw[i]].wct;
			}
		}
		else{
			std::cout << "Error in constrain type" << std::endl;
		}
	}
}


///////////////////////////////////////////////////////Augmented largrangian function///////////////////////////////////////////////////////////////

double miu, alpha_bar = 0.1;
double omega_f_final(1e-7);
double omega_x_final(1e-5);
double alpha = 0.1;
double beta = 0.9;
double eta0;
double w = 0.5;
double eta_final(0.001);


double sign_function(const double &lamda, const double &miu, const double &val)
{
	double res(0);
	if (val<0) //constraint is voliated
	{
		res = (1 / (2 * miu))*val*val - miu*val;
	}
	else//constraint is satisfied
	{
		res = -0.5*lamda*lamda*miu;
	}
	return res;
}

double compute_constrain()
{
	dlib::matrix<double> cs_v;
	assemble_cs(cs_v);
	double sum = 0;
	for (unsigned i = 0; i < CS_val.nr(); i++)
	{
		val(i) = ( cs_v(i)+CS_val(i)) / abs(CS_range(i));
		sign_val(i) = sign_function(lamda(i), miu, val(i));
	}

	return dlib::sum(sign_val);
}

void initial_parameter(unsigned &num_CS, double &eta, double &omega_f,
	double &omega_x, const dlib::matrix<double> &x)
{
	num_CS = CS_val.nr();
	lamda.set_size(num_CS, 1); lamda = 0;
	//initial miu;
	val.set_size(num_CS, 1);
	sign_val.set_size(num_CS, 1);
	double obj_v = sim_model.run(x);
	compute_constrain();
	double cx = dlib::sum(val);
	miu = 10 * cx / obj_v;
	//compute initial miu;
	omega_f = std::min(miu, 0.1);
	omega_x = std::min(miu, 0.1);
	eta0 = pow(std::min(miu, 0.1), 0.1);
	eta = eta0;
}


double constrain_funct(const dlib::matrix<double> &x)
{
	double result(0);
	result += sim_model.run(x);
	result += compute_constrain();
	return result;
}

const dlib::matrix<double> constrain_grad(const dlib::matrix<double> &x)
{
	dlib::matrix<double> grad(x.nr(), 1);
	std::vector<dlib::matrix<double>> CS_gradient;
	grad = 0;
	grad += sim_model.grad(x);
	compute_constrain();
	assemble_cs_gradient(CS_gradient);
	for (unsigned i = 0; i < val.nr(); i++)
	{
		if (val(i)<0)
		{
			grad += (val(i) / miu - lamda(i))*CS_gradient[i];
		}
	}
	return grad;
}

void update_larange_lamda(const unsigned &num_CS, double & eta, double &omega_f, double &omega_x)
{
	for (unsigned i = 0; i < num_CS; i++)
	{
		lamda(i, 0) = std::max((lamda(i, 0) - val(i, 0) / miu), 0.0);
	}
	omega_f = std::max(omega_f*std::min(pow(miu, beta), w), omega_f_final);
	omega_x = std::max(omega_x*std::min(pow(miu, beta), w), omega_x_final);
	eta = std::max(eta*std::min(pow(miu, beta), w), eta_final);; //correct update
}

void update_larange_miu(double & eta, double &omega_f, double &omega_x)
{
	miu = 0.1*miu;
	//???????????????????????????????????????????????????? need to be fixed
	omega_f = std::max(omega_f*std::min(pow(miu, alpha), w), omega_f_final);
	omega_x = std::max(omega_x*std::min(pow(miu, alpha), w), omega_x_final);
	//????????????????????????????????????????????????????
	eta = std::max(eta*std::min(pow(miu, alpha), w), eta_final);
}

double compute_validation()
{
	double nv(0), result(0);
	for (unsigned i = 0; i<val.nr(); i++){
		if (val(i)<0){
			nv++;
			result += val(i)*val(i);
		}
	}
	if (nv != 0){
		result /= nv;
		return pow(result, 0.5);
	}
	else{
		return 0;
	}
}

dlib::matrix<double> larangian_funct(dlib::matrix<double> &x, bool output_i)
{
	scaling(x);
	unsigned num_CS;
	double eta, omega_f, omega_x, delta_cv, alpha;
	double gradient_norm;
	//clear previouse output/////////////////////////////////////////////////////
	std::ofstream output0("output/obj.dat");
	std::ofstream output("output/iterpoint.dat");
	output.close();
	output0.close();
	///////////////////////////////////////////////////////////////////////////////
	std::ofstream output1("output/violation.dat");
	std::ofstream output2("output/penalty_factor.dat");
	std::ofstream output3("output/outer_point.dat");
	optimization opt;
	initial_parameter(num_CS, eta, omega_f, omega_x, x);
	while (1){
		opt.find_min_new(dlib::bfgs_search_strategy(),
			dlib::objective_delta_stop_strategy(omega_f).be_verbose(),
			dlib::x_delta_stop_strategy(omega_x).be_verbose(),
			constrain_funct,
			constrain_grad, x, output_i);
		delta_cv = compute_validation();
		if (delta_cv<eta){
			if (eta = eta_final
				&& omega_f_final == omega_f
				&& omega_x_final == omega_x){
				return x;
			}
			else{
				update_larange_lamda(num_CS, eta, omega_f, omega_x);
			}
		}
		else{
			update_larange_miu(eta, omega_f, omega_x);
		}
		if (output_i){
			output1 << delta_cv / miu << std::endl;
			output2 << miu << std::endl;
			output3 << trans(x) << std::endl;
		}
	}
	scaling_back(x);
}




///////////////////////////////////////////////////////Augmented largrangian function end///////////////////////////////////////////////////////////////


//void F_D(simulator &model, double peturb_scale, dlib::matrix<double> &x)
//{
//	unsigned id=15;
//	double original, peturbed,NPV,peturbed_size;
//	dlib::matrix<double> m,sens,x_p;
//	read_from_file("input/m_pri.dat", m);
//	dlib::matrix<double> grad(x.size(),1);
//	original=model.run(m, x, sens, NPV, 0);
//	std::cout << sens << std::endl;
//	std::cout << "Press key to compute Finite Difference gradient" <<original<< std::endl;
//	system("pause");
//
//	for (unsigned i = 0; i < x.size();i++)
//	{
//		x_p = x;
//		if (i<x.size()/2)
//		{
//			peturbed_size = peturb_scale*(x(i)-4000);
//			x_p(i, 0) += peturbed_size;
//		}
//		else
//		{
//			std::cout << i << std::endl;
//			peturbed_size = peturb_scale*x(i);
//			x_p(i, 0) += peturbed_size;
//		}
//		peturbed=model.run(m, x_p, sens, NPV, 0);
//		grad(i, 0) = (peturbed - original) / peturbed_size;
//		std::cout<<"Peturbted:"<<peturbed << "   Gradient " << i << " is " << grad(i, 0) <<"  adjoint gradient:"<<sens(i,0) <<std::endl;
//	}
//	dlib::matrix<double> error = (grad - sens);
//	std::cout << "Adjoint" << std::endl;
//	for (unsigned i = 0; i < error.size(); i++)
//	{
//		std::cout << i << "  " << sens(i) << std::endl;
//	}
//	std::cout << "F_D " << std::endl;
//	for (unsigned i = 0; i < error.size(); i++)
//	{
//		std::cout << i << "  " << grad(i) << std::endl;
//	}
//	std::cout << "error " << std::endl;
//	for (unsigned i = 0; i < error.size();i++)
//	{
//		error(i) = error(i) / (grad(i)+10e-5);
//		std::cout << i<< "  " << error(i)*100 << std::endl;
//	}
//	std::ofstream out("output/validation/FD.dat"), out1("output/validation/adjoint.dat"), out2("output/validation/err.dat");
//	out << grad << std::endl;
//	out1 <<sens << std::endl;
//	out2 <<  error<< std::endl;
//	std::cout << error << std::endl;
//	system("pause");
//}


int main()
{
	double NPV_orig, NPV, NPV_max=0,original;
	dlib::matrix<double> x, x_t, Cx, Cxf, Cxi;
	dlib::matrix<double> x_max(nc*nw, 1),sens;
	std::vector<double> NPV_v;
	//initial data

	load_constrain();
	simulator model;

	read_from_file("input/m_pri.dat", m);
	read_from_file("input/wellschedule.dat", x);

	sim_model.initial_properties(m);
	//std::cout << model.run(m, x, sens, NPV, 0) << std::endl;
	std::cout << sim_model.run(x) << std::endl;
	//original = model.run(m, x, sens, NPV, 0);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Normal Optimization
	//read_parameter(ne, nw, np, nc ,up,low);
	//simulator model(m, x, NPV_orig, 1);
	//system("pause");
	//NPV_max = NPV_orig;
	//NPV = NPV_orig;
	//Cx = generate_Cx(nw, nc, 10.0, 0.01);
	//dlib::matrix<double> L = chol(Cx);
	//optimize_0_1_scale(MAX_step, step, x, NPV, Cx, L, NPV_max, x_max);
	//std::ofstream outputmax("output/max_x.dat");
	//outputmax << x_max << std::endl;
	//simulator model(m, x_max, NPV_orig, 1);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Robust Optimization


	//read_parameter(ne, nw, np, nc, nm, up, low);
	//NPV_v.resize(nm);
	//compute_model(nm, m, x, NPV_orig,NPV_v ,0);

	//std::ofstream output("output/oig_NPV.out");
	//for (int i = 0; i < nm; i++)
	//{
	//	output << NPV_v[i] << std::endl;;
	//}
	//std::cout << "The original expected NPV is "<<NPV_orig << std::endl;
	//system("pause");

	////std::ifstream input("input/orig_NPV.dat");
	////for (int i = 0; i < nm; i++)
	////{
	////	input >> NPV_v[i];
	////}


	//NPV_max = NPV_orig;
	//NPV = NPV_max;
	//std::cout << "The original expected NPV is " << NPV_max << std::endl;
	//Cx = generate_Cx(nw, nc, 10.0, 0.01);
	//dlib::matrix<double> L = chol(Cx);
	//robust_optimize_0_1_scale(MAX_step, step, x, NPV, Cx, L,NPV_v, NPV_max, x_max);


	//std::ofstream output1("output/final_NPV.out");
	//for (int i = 0; i < nm; i++)
	//{
	//	output1 << NPV_v[i] << std::endl;;
	//}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// selected ensemble run
	//read_parameter(ne, nw, np, nc, nm, up, low);
	//NPV_v.resize(nm);
	//compute_model(nm, m, x, NPV_orig, NPV_v, 0);
	//std::vector<double> NPV_n;
	//std::vector<unsigned> rank;
	//std::ofstream output("output/oig_NPV.out");
	//for (int i = 0; i < nm; i++)
	//{
	//	output << NPV_v[i] << std::endl;
	//}
	//generate_cluster(NPV_v, NPV_n, rank);
	//std::cout << "The original expected NPV is " << NPV_orig << std::endl;
	//////////**************************************************************************************************************
	//std::ofstream outrank("output/NPV_rank.out");
	//for (unsigned i = 0; i < rank.size();i++)
	//{
	//	outrank << rank[i] << std::endl;
	//}

	//for (unsigned i = 0; i < NPV_n.size(); i++)
	//{
	//	NPV_max += NPV_n[i];
	//}
	//NPV_max = NPV_max / double(NPV_n.size());
	//system("pause");

	//NPV = NPV_orig;
	//Cx = generate_Cx(nw, nc, 10.0, 0.01);
	//dlib::matrix<double> L = chol(Cx);
	//robust_optimize_0_1_scale_p(MAX_step, step, x, NPV, rank, Cx, L, NPV_n, NPV_max, x_max);

	//std::ofstream output1("output/final_NPV.out");
	//for (int i = 0; i < nm; i++)
	//{
	//	output1 << NPV_v[i] << std::endl;;
	//}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Multi-Objective Optimization
	//double w1=0.5, w2=0.5;
	//read_parameter(ne, nw, np, nc, nm, up, low);
	//NPV_v.resize(nm);
	//model_multi_obj(nm, m, x, NPV_orig, NPV_v,w1,w2, 0);
	//std::ofstream output("output/oig_NPV.out");
	//for (int i = 0; i < nm; i++)
	//{
	//	output << NPV_v[i] << std::endl;;
	//}
	//std::cout << "The original expected NPV is " << NPV_orig << std::endl;
	//system("pause");

	//NPV_max = NPV_orig;
	//NPV = NPV_orig;
	//Cx = generate_Cx(nw, nc, 10.0, 0.01);
	//dlib::matrix<double> L = chol(Cx);
	//multi_object_optimize_0_1_scale(MAX_step, step, x, NPV, Cx, L, NPV_v, NPV_max, x_max,w1,w2);

	//std::ofstream output1("output/final_NPV.out");
	//for (int i = 0; i < nm; i++)
	//{
	//	output1 << NPV_v[i] << std::endl;;
	//}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//clustering run


	//std::vector<double> Qo, Qw;

	//std::vector<double> NPV_n;
	//std::vector<unsigned> rank;
	//read_parameter(ne, nw, np, nc, nm, up, low);
	//NPV_v.resize(nm);

	//std::ifstream input("input/cluster.dat");
	//unsigned trank;
	//while (input>>trank)
	//{
	//	rank.push_back(trank);
	//	std::cout << trank << std::endl;
	//}
	//NPV_n.resize(rank.size());
	//compute_model_p(m, x, NPV_max, NPV_n, rank, 0);
	//system("pause");

	//NPV_max = 0;
	//for (unsigned i = 0; i < NPV_n.size(); i++)
	//{
	//	NPV_max += NPV_n[i];
	//}
	//NPV_max = NPV_max / double(NPV_n.size());

	//NPV = NPV_max;
	//Cx = generate_Cx(nw, nc, 10.0, 0.01);
	//dlib::matrix<double> L = chol(Cx);
	//robust_optimize_0_1_scale_p(MAX_step, step, x, NPV, rank, Cx, L, NPV_n, NPV_max, x_max);
	//compute_model(nm, m, x_max, NPV_orig, NPV_v, 0);
	//std::ofstream output1("output/final_NPV.out");
	//for (int i = 0; i < nm; i++)
	//{
	//	output1 << NPV_v[i] << std::endl;;
	//}

}