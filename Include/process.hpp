#ifndef PROCESS_HPP
#define PROCESS_HPP

#include <iostream>
#include <functional>
#include <algorithm>
#include <cmath>

#include "rand_var.hpp"

/*************************************************************************************
	*
	* General purpose functions and typedefs
	*
	*************************************************************************************/

/*
	* FUNC eq_spaced_times<n>
	*
	* Return a list of equally spaced times between 0 and T.
	* /!\ (n-1 intervals and n timestamps)
	*
	*/
template<unsigned n>
inline Array<n> eq_spaced_times(double T)
{
				double h = T / (n-1);

				Array<n> times;
				for(unsigned i=0; i<=(n-1); ++i)
								times.at(i) = i*h;

				return times;
}


typedef std::pair<double,double> State;

template<unsigned t_dim>
using Process = std::array<State,t_dim>;


/*
	* FUNC proc_values<t_dim>
	*
	* Take a process (as a 'std::array<State,t_dim>') and return
	* only the values (as an 'Array<t_dim>')
	*
	*/
template<size_t t_dim>
Array<t_dim> proc_values(const Process<t_dim>& proc)
{
				Array<t_dim> a;
				for(int i=0; i<t_dim; ++i)
								a.at(i) = proc.at(i).second;

				return a;
}


/*
	* FUNC proc_times<t_dim>
	*
	* Take a process (as a 'std::array<State,t_dim>') and return
	* only the times (as an 'Array<t_dim>')
	*
	*/
template<size_t t_dim>
Array<t_dim> proc_times(const Process<t_dim>& proc)
{
				Array<t_dim> a;
				for(int i=0; i<t_dim; ++i)
								a.at(i) = proc.at(i).first;

				return a;
}


/*************************************************************************************
	*
	* STRUCT StdBrownian<t_dim>
	*
	* Standard Brownian Motion at 't_dim' given times
	*
	*************************************************************************************/
template<unsigned t_dim>
struct StdBrownian{

				typedef Process<t_dim> result_type;
				static constexpr unsigned dim_alea =  ((t_dim % 2)==0) ? t_dim: t_dim+1;

				StdBrownian() = delete;
				StdBrownian(const Array<t_dim>& times):
								G_(makeFill<t_dim>(0.),make_stds_from_times(times)),
								times_(times) {}

				result_type operator()(const Array<dim_alea>& pt)
				{
								G_val = G_(pt);

								val_tp.at(0) = State(times_.at(0), G_val.at(0));
								for(i=1; i<t_dim; ++i)
												val_tp.at(i) = State(times_.at(i), val_tp.at(i-1).second + G_val.at(i));

								return val_tp;
				}

protected:
				Gaussian_Ind<t_dim> G_;
				const Array<t_dim> times_;

				Array<t_dim> G_val;
				int i;
				result_type val_tp;


				/*
					* Compute square roots of the time increment, to use as standard deviations for
					* Gaussian variables
					*/
				static Array<t_dim> make_stds_from_times(const Array<t_dim>& times)
				{
								Array<t_dim> a;

								a.at(0) = times.at(0);
								for(int i=1; i<t_dim; ++i)
												a.at(i) = std::sqrt(times.at(i) - times.at(i-1));

								return a;
				}
};


/*************************************************************************************
	*
	* STRUCT Black_Scholes<t_dim>
	*
	* Black_Scholes process at 't_dim' given times
	*
	*************************************************************************************/
template<unsigned t_dim>
struct Black_Scholes{

				typedef Process<t_dim> result_type;
				static constexpr unsigned dim_alea = ((t_dim % 2)==0) ? t_dim: t_dim+1;

				Black_Scholes() = delete;
				Black_Scholes(double r,double sigma,double x0,const Array<t_dim>& times):
								B_(times),
								r_(r), sigma_(sigma), x0_(x0)
				{
								bs_func = [x0,r,sigma] (State s){
												return State(s.first, x0 * std::exp( (r - 0.5*std::pow(sigma,2))*s.first + sigma*s.second ));
								};
				}

				result_type operator()(const Array<dim_alea>& pt)
				{
								val_tp = B_(pt);
								std::transform(val_tp.begin(),val_tp.end(),val_tp.begin(),bs_func);
								return val_tp;
				}

protected:
				StdBrownian<t_dim> B_;
				const double r_;
				const double sigma_;
				const double x0_;

				std::function<State(State)> bs_func;
				result_type val_tp;
};






#endif // PROCESS_HPP
