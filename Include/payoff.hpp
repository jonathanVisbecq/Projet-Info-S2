/*************************************************************************************
	*
	* Collection of payoffs (functionals) to be used in MC simulations
	*
	*************************************************************************************/
#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include "stl_headers.hpp"
#include "process.hpp"


/*************************************************************************************
	* STRUCT Last_Value<t_dim>
	*
	* Return last value of the process
	*
	*************************************************************************************/
template<unsigned t_dim>
struct Last_Value{

				Last_Value() {}

				double operator()(const Process<t_dim>& proc)
				{
								return proc.back().second;
				}
};



/*************************************************************************************
	* STRUCT Asian_Call<t_dim>
	*
	* Return asian call payoff on a given process
	*
	*************************************************************************************/
template<unsigned t_dim>
struct Asian_Call{

				Asian_Call() = delete;
				Asian_Call(double K, double r = 0., double T = 0.):
								K_(K), r_(r), T_(T) {}

				double operator()(const Process<t_dim>& proc)
				{
								s = 0.;
								for(i=0; i<t_dim; ++i)
												s += proc.at(i).second;

								s /= t_dim;

								return (s > K_) ? std::exp(-r_*T_)*(s-K_) : 0.;
				}


protected:
				double K_;
				double r_;
				double T_;

				double s;
				int i;

};

/*************************************************************************************
	* STRUCT Asian_Put<t_dim>
	*
	* Return asian put payoff on a given process
	*
	*************************************************************************************/
template<unsigned t_dim>
struct Asian_Put{

				Asian_Put() = delete;
				Asian_Put(double K, double r = 0., double T = 1.):
								K_(K), r_(r), T_(T) {}

				double operator()(const Process<t_dim>& proc)
				{
								s = 0.;
								for(i=0; i<t_dim; ++i)
												s += proc.at(i).second;

								s /= t_dim;

								return (K_ > s) ? std::exp(-r_*T_)*(K_ - s) : 0.;
				}


protected:
				double K_;
				double r_;
				double T_;

				double s;
				int i;

};


























#endif // PAYOFF_HPP
