#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include "dlib/optimization/optimization_abstract.h"
#include "dlib/optimization/optimization_search_strategies.h"
#include "dlib/optimization/optimization_stop_strategies.h"
#include "dlib/optimization/optimization_line_search.h"

using namespace dlib;

class optimization
{
public:
	optimization()
	{};

	template <
		typename search_strategy_type,
		typename stop_strategy_type,
		typename stop_strategy_type_1,
		typename funct,
		typename funct_der,
		typename T
	>
	double find_min_new(
	search_strategy_type search_strategy,
	stop_strategy_type stop_strategy,
	stop_strategy_type_1 stop_strategy_1,
	funct& f,
	funct_der& der,
	T& x,
	bool & is_output
	)
	{

		T g, s;
		double f_value = f(x);
		g = der(x);
		double last_alpha = 0.1;
		while (stop_strategy.should_continue_search(x, f_value, g) 
			|| stop_strategy_1.should_continue_search(x, f_value, g))
		{
			s = search_strategy.get_next_direction(x, f_value, g);
			double alpha = backtracking_line_search_new(
				f,x, s,
				f_value,
				dot(g, s), // compute gradient for the line search
				last_alpha,
				search_strategy.get_wolfe_rho(),
				search_strategy.get_max_line_search_iterations());
			if (alpha == last_alpha)
				last_alpha = std::min(last_alpha * 10, 1.0);
			else
				last_alpha = alpha;
			// Take the search step indicated by the above line search
			x += alpha*s;
			g = der(x);
			if (is_output)
			{
				std::cout << (stop_strategy.should_continue_search(x, f_value, g)
					|| stop_strategy_1.should_continue_search(x, f_value, g)) << std::endl;
				std::cout << "x: " << trans(x) << std::endl;
				std::cout << "f_value: " << f_value << std::endl;
				std::cout << "-----------------------------------------------------------" << std::endl;
			}
		}

		output1.close();
		output2.close();
		return f_value;
	}

	template <
		typename search_strategy_type,
		typename stop_strategy_type,
		typename funct,
		typename funct_der,
		typename T
	>
	double find_min(
	search_strategy_type search_strategy,
	stop_strategy_type stop_strategy,
	funct& f,
	funct_der& der,
	T& x
	)
	{
		T g, s;
		double f_value = f(x);
		g = der(x);
		double last_alpha = 1;
		while (stop_strategy.should_continue_search(x, f_value, g))
		{
			s = search_strategy.get_next_direction(x, f_value, g);

			double alpha = backtracking_line_search(
				make_line_search_function(f, x, s, f_value),
				f_value,
				dot(g, s), // compute gradient for the line search
				last_alpha,
				search_strategy.get_wolfe_rho(),
				search_strategy.get_max_line_search_iterations());

			if (alpha == last_alpha)
				last_alpha = std::min(last_alpha * 10, 1.0);
			else
				last_alpha = alpha;
			// Take the search step indicated by the above line search
			x += alpha*s;
			g = der(x);
		}
		return f_value;
	}





	template <typename funct>
	double backtracking_line_search_new(
		funct f,
		dlib::matrix<double> x,
		dlib::matrix<double> grad,
		double &f0,
		double d0,
		double alpha,
		double rho,
		unsigned long max_iter
		)
	{
		if ((d0 > 0 && alpha > 0) ||
			(d0 < 0 && alpha < 0))
		{
			alpha *= -1;
		}
		bool have_prev_alpha = false;
		double prev_alpha = 0;
		double prev_val = 0;
		unsigned long iter = 0;
		while (true)
		{
			++iter;
			const double val = f(alpha*grad + x);
			if (val <= f0 + alpha*rho*d0 || iter >= max_iter)
			{
				f0 = val;
				return alpha;
			}
			else
			{
				// Interpolate a new alpha.  We also make sure the step by which we
				// reduce alpha is not super small.
				double step;
				if (!have_prev_alpha)
				{
					if (d0 < 0)
						step = alpha*put_in_range(0.1, 0.9, poly_min_extrap(f0, d0, val));
					else
						step = alpha*put_in_range(0.1, 0.9, poly_min_extrap(f0, -d0, val));
					have_prev_alpha = true;
				}
				else
				{
					if (d0 < 0)
						step = put_in_range(0.1*alpha, 0.9*alpha, poly_min_extrap(f0, d0, alpha, val, prev_alpha, prev_val));
					else
						step = put_in_range(0.1*alpha, 0.9*alpha, -poly_min_extrap(f0, -d0, -alpha, val, -prev_alpha, prev_val));
				}

				prev_alpha = alpha;
				prev_val = val;

				alpha = step;
			}
		}
	}

	// ----------------------------------------------------------------------------------------


private:
	std::ofstream output1;
	std::ofstream output2;

};