#pragma once
#include <dlib/matrix.h>
#include <fstream>
#include "simulator.cpp"
#include <cmath>
#include <random>
#include <time.h>


class optimize
{
public:
	optimize()
	{
		output1.open("output/NPV_seq.out");
		output2.open("output/step.out");
		read_parameter();
		NPV_v.resize(nm);
		read_from_file("input/m_pri.dat", m);
		read_from_file("input/hwellschedule.dat", x);

	}


private:
	unsigned ne, nw, np, nc, nm, MAX_step;
	double step, larg_step;
	dlib::matrix<double> m;
	dlib::matrix<double> up;
	dlib::matrix<double> low;
	std::ofstream output1;
	std::ofstream output2;
	std::ofstream output;
	std::vector<double> NPV_v;
	dlib::matrix<double> x, x_t, Cx, Cxf;
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
			v.set_size(nwc, 0);
			for (unsigned i = 0; i < nwc; i++)
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
		up.set_size(nc*nw, 1);
		low.set_size(nc*nw, 1);
		dlib::set_subm(up, dlib::rectangle(0, 0, nc*np - 1, 0)) = pmax;
		dlib::set_subm(up, dlib::rectangle(nc*np, 0, nc*nw - 1, 0)) = qmax;

		dlib::set_subm(low, dlib::rectangle(0, 0, nc*np - 1, 0)) = pmin;
		dlib::set_subm(low, dlib::rectangle(nc*np, 0, nc*nw - 1, 0)) = qmin;
	}

	void read_parameter()
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

	dlib::matrix<double> generate_Cx(unsigned nw, unsigned nc, double control_dis, double delta)
	{
		double tmp(0), rho(0);
		dlib::matrix<double> Cx(nw*nc, nw*nc);
		Cx = 0;
		for (unsigned n = 0; n < nc; n++)
		{
			for (unsigned c = 0; c < nc; c++)
			{
				tmp = abs(double(n - c)) / control_dis;
				rho = spherecov(delta, tmp);
				for (unsigned w = 0; w < nw; w++)
				{
					Cx(n + w*nc, c + w*nc) = rho;
					Cx(c + w*nc, n + w*nc) = rho;
				}
			}
		}
		return Cx;
	}
	void scaling(dlib::matrix<double> &x)
	{
		x = dlib::pointwise_multiply((x - low), 1 / (up - low));
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

	dlib::matrix<double> generate_sample(dlib::matrix<double> &x, dlib::matrix<double> &L)
	{
		dlib::matrix<double> z;
		std::mt19937 eng;
		eng.seed((unsigned)time(NULL));
		std::normal_distribution<double> distribution(0.0, 1);
		z.set_size(nw*nc, 1);
		for (int i = 0; i < nw*nc; ++i)
			z(i, 0) = distribution(eng);
		for (unsigned i = np*nc; i < nw*nc; i++)
		{
			if (x(i, 0) == 0 && z(i, 0)<0)
			{
				z(i, 0) = -z(i, 0);
			}
		}
		dlib::matrix<double> xp = x + L*z;
		return xp;
	}

	void compute_Cxf_robust_p(dlib::matrix<double> &x, std::vector<double> &NPV_m,
		dlib::matrix<double> &Cxf, dlib::matrix<double> &L, std::vector<unsigned> &rank)
	{
		dlib::matrix<double> xp(nc*nw, 1);
		double f_n;

		Cxf.set_size(nc*nw, 1);
		Cxf = 0;
		unsigned m_n;
		for (unsigned j = 0; j < rank.size(); j++)
		{
			m_n = rank[j];
			for (unsigned i = 0; i < ne; i++)
			{
				xp = generate_sample(x, L);
				scaling_back(xp);
				simulator model(m, xp, f_n, 0, m_n);
				scaling(xp);
				Cxf += (xp - x)*(f_n - NPV_m[j]);
				std::cout << "The " << i + j*ne + 1 << "th peturbation. " << std::endl;
			}
		}
		Cxf /= double(ne*rank.size());
	}

	void compute_Cxf_robust(dlib::matrix<double> &x, std::vector<double> &NPV_m,
		dlib::matrix<double> &Cxf, dlib::matrix<double> &L)
	{
		dlib::matrix<double> xp(nc*nw, 1);
		double f_n;
		Cxf.set_size(nc*nw, 1);
		Cxf = 0;
		for (unsigned j = 0; j < nm; j++)
		{
			for (unsigned i = 0; i < ne; i++)
			{
				xp = generate_sample(x, L);
				scaling_back(xp);
				simulator model(m, xp, f_n, 0, j);
				scaling(xp);
				Cxf += (xp - x)*(f_n - NPV_m[j]);
				std::cout << "The " << i + j*ne + 1 << "th peturbation. " << std::endl;
			}
		}
		Cxf /= double(ne*nm);
	}

	void compute_Cxf(dlib::matrix<double> &x, double f,
		dlib::matrix<double> &Cxf, dlib::matrix<double> &L)
	{
		std::vector<dlib::matrix<double>> xp(ne);
		std::vector<double> f_n(ne);

		Cxf.set_size(nc*nw, 1);
		Cxf = 0;
		for (unsigned i = 0; i < ne; i++)
		{
			xp[i] = generate_sample(x, L);
			scaling_back(xp[i]);
			simulator model(m, xp[i], f_n[i], 0);
			scaling(xp[i]);
			Cxf += (xp[i] - x)*(f_n[i] - f);
			std::cout << "The " << i + 1 << "th peturbation. " << std::endl;
		}
		Cxf /= double(ne);
	}
	dlib::matrix<double> compute_gradient(dlib::matrix<double> &Cx,
		dlib::matrix<double> &Cxf)
	{
		dlib::matrix<double> gradient = Cx*Cxf;
		unsigned t;
		double norm = norm_vector(gradient, t);
		gradient = gradient / norm;
		return gradient;
	}

	double norm_vector(dlib::matrix<double> &x, unsigned &t)
	{
		double max = 0;
		for (unsigned i = 0; i < x.size(); i++)
		{
			if (max < x(i))
			{
				max = x(i);
				t = i;
			}
		}
		return max;
	}

	dlib::matrix<double> enopt_gradient(dlib::matrix<double> &x, double f,
		dlib::matrix<double> &Cx, dlib::matrix<double> &L)
	{
		//////////////////////////////generate peturbation//////////
		std::vector<dlib::matrix<double>> xp(ne);
		dlib::matrix<double> x_p, x_a(nw*nc, 1), f_p(ne, 1);
		dlib::matrix<double> Cxf(nc*nw, 1);
		Cxf = 0;
		double fn, fa;
		x_a = 0;
		for (unsigned i = 0; i < ne; i++)
		{
			xp[i] = generate_sample(x, L);
			x_a += xp[i];
			scaling_back(xp[i]);
			simulator model(m, xp[i], f_p(i, 0), 0);
			scaling(xp[i]);
			std::cout << "The " << i + 1 << "th peturbation. " << std::endl;
		}
		fa = sum(f_p) / ne;
		x_a /= ne;
		for (unsigned i = 0; i < ne; i++)
		{
			Cxf += (xp[i] - x_a)*(f_p(i, 0) - fa);
		}
		Cxf = Cxf / double(ne - 1);
		dlib::matrix<double> grad = compute_gradient(Cx, Cxf);
		////////////////////////////////////////////////////////////
		//svd(xp, U, S, V);
		//for (unsigned i = 0; i < ne;i++)
		//{
		//	S(i, i) = 1 / S(i,i);
		//}
		//grad = V*S*trans(U)*f_p;
		//unsigned t;
		//double norm = norm_vector(grad, t);
		//grad = grad / norm;
		////////////////////////////////////////////////////////////
		return grad;
	}

	bool one_step(double &step, dlib::matrix<double> &x,
		double &f, dlib::matrix<double> &Cx, dlib::matrix<double> &L, unsigned &peturb_num)
	{
		dlib::matrix<double> Cxf, x_t(nc*nw, 1), x_max(nc*nw, 1), gradient(nc*nw, 1);
		double f_t = 0, f_m = 0;
		unsigned np = 0;
		gradient = enopt_gradient(x, f, Cx, L);
		while (f_t <= f)
		{
			x_t = x + step*gradient;
			for (unsigned i = 0; i<nc*nw; i++)
			{
				if (x_t(i, 0)>1)
				{
					x_t(i, 0) = 1;
				}
				else if (x_t(i, 0)<0)
				{
					x_t(i, 0) = 0;
				}
			}
			scaling_back(x_t);
			simulator model(m, x_t, f_t, 0);
			scaling(x_t);
			np++;
			if (f_t>f_m)  //record largest NPV;
			{
				f_m = f_t;
				x_max = x_t;
			}
			if (np>5)
			{
				x = x_max;
				f = f_m;
				step = std::min(4 * step, larg_step);
				return false;
			}
			step *= 0.5;
		}
		x = x_t;
		f = f_t;
		//std::ofstream update("output/x_updata.dat");
		//update << 2 * step*gradient << std::endl;
		// Prevent not to be trapped at boundary
		step = std::min(4 * step, larg_step);
		return true;
	}
	void optimize_0_1_scale(unsigned MAX_step, double &step, dlib::matrix<double> &x,
		double &NPV, dlib::matrix<double> &Cx, dlib::matrix<double> &L,
		double &NPV_max, dlib::matrix<double> &x_max)
	{
		unsigned peturb_num = 0;
		scaling(x);
		for (unsigned i = 0; i < MAX_step; i++)
		{

			one_step(step, x, NPV, Cx, L, peturb_num);
			//NPV_all.push_back(NPV);
			if (NPV>NPV_max)
			{
				NPV_max = NPV;
				x_max = x;
			}
			scaling_back(x);
			output.open("output/optimal_rate.dat");
			output << x << std::endl;
			output.close();
			output1 << NPV << std::endl;
			output2 << step << std::endl;
			scaling(x);
		}
	}
	void compute_model(int nm, dlib::matrix<double> &m, 
		dlib::matrix<double> &s, double &NPV_expected, 
		std::vector<double> &NPV_m, bool output_i)
	{
		double NPV_v = 0;
		NPV_expected = 0;
		for (int i = 0; i < nm; i++)
		{
			simulator model(m, s, NPV_v, output_i, i);
			NPV_m[i] = NPV_v;
			NPV_expected += NPV_v;
			std::cout << "The NPV of " << i + 1 << "th simulation run is " << NPV_v << std::endl;
		}
		NPV_expected = NPV_expected / double(nm);
	}


	bool one_step_robust(double &step, dlib::matrix<double> &x,
		double &f, dlib::matrix<double> &Cx, dlib::matrix<double> &L, std::vector<double> &NPV_m)
	{
		dlib::matrix<double> Cxf, x_t(nc*nw, 1), x_max(nc*nw, 1), gradient(nc*nw, 1);
		double f_t = 0, f_m = 0, orig_step(step);
		unsigned np = 0, nrsp = 0;
		compute_Cxf_robust(x, NPV_m, Cxf, L);
		gradient = compute_gradient(Cx, Cxf);
		while (f_t <= f)
		{
			x_t = x + step*gradient;
			for (unsigned i = 0; i<nc*nw; i++)
			{
				if (x_t(i, 0)>1)
				{
					x_t(i, 0) = 1;
				}
				else if (x_t(i, 0)<0)
				{
					x_t(i, 0) = 0;
				}
			}
			scaling_back(x_t);
			compute_model(nm, m, x_t, f_t, NPV_m, 0);
			scaling(x_t);
			np++;
			step *= 0.5;
			if (f_t>f_m)  //record largest NPV;
			{
				f_m = f_t;
				x_max = x_t;
			}
			if (np>5)
			{
				x = x_max;
				f = f_m;
				step = std::min(8 * step, larg_step);
				return false;
			}
		}
		x = x_t;
		f = f_t;
		step = std::min(4 * step, larg_step);
		return true;
	}

	void robust_optimize_0_1_scale(unsigned MAX_step, double &step, dlib::matrix<double> &x,
		double &NPV, dlib::matrix<double> &Cx, dlib::matrix<double> &L, std::vector<double> &NPV_m,
		double &NPV_max, dlib::matrix<double> &x_max)
	{
		scaling(x);
		for (unsigned i = 0; i < MAX_step; i++)
		{
			one_step_robust(step, x, NPV, Cx, L, NPV_m);
			//NPV_all.push_back(NPV);
			if (NPV>NPV_max)
			{
				NPV_max = NPV;
				x_max = x;
			}
			scaling_back(x);
			output.open("output/optimal_rate.dat");
			output << x << std::endl;
			output.close();
			output1 << NPV << std::endl;
			output2 << step << std::endl;
			scaling(x);
		}
	}


public:





private:




};