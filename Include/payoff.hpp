/*************************************************************************************
	*
	* Collection of payoffs (functionals) to be used in MC simulations
	*
	* Constraints on payoffs:
	*		- template 't_dim' for the number of discrete times in the process
	*  - operator
	*
	*
	*************************************************************************************/
#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <functional>

#include "rand_var.hpp"
#include "process.hpp"


/*
	* STRUCT Composed_Dist<Dist>
	*
	* Allow composition of a distribution by a real function
	*
	* Implement constraints on Distribution type.
	*
	*/
template<typename Dist>
struct Composed_Dist{

				typedef double result_type;
				typedef std::function<double(typename Dist::result_type)> func_type;
				static constexpr unsigned dim_alea = Dist::dim_alea;

				Composed_Dist(const Dist& dist, const func_type& func):
								dist_(dist), func_(func) {}

				double operator()(const Array<dim_alea>& pt)
				{
								return func_(dist_(pt));
				}

protected:
				Dist dist_;
				func_type func_;
};


template<typename Dist>
Composed_Dist<Dist>
compose_dist(const Dist& dist,
													const typename Composed_Dist<Dist>::func_type& func)
{
				return Composed_Dist<Dist>(dist,func);
}





/*
	* STRUCT Last_Value<t_dim>
	*
	* Return last value of the process
	*
	*/
template<unsigned t_dim>
struct Last_Value{

				Last_Value() {}

				double operator()(const Process<t_dim>& proc)
				{
								return proc.back().second;
				}
};



/*
	* STRUCT Asian_Call<t_dim>
	*
	* Return asian call payoff on a given process
	*
	*/
template<unsigned t_dim>
struct Asian_Call{

				Asian_Call() = delete;
				Asian_Call(double K): K_(K) {}

				double operator()(const Process<t_dim>& proc)
				{
								for(i=0; i<t_dim; ++i)
												s += proc.at(i).second;

								s /= t_dim;

								return (s > K_) ? (s-K_): 0;
				}


protected:
				const double K_;

				double s;
				int i;

};


/*
	* STRUCT Basket_Call<t_dim>
	*
	* Return basket call payoff on a given process
	*
	* ------------------------------------------------
	* TODO: implement multidimensional processes
	* ------------------------------------------------
	*
	*/

/*
template<unsigned t_dim>
struct Basket_Call{

				Basket_Call() = delete;
				Basket_Call(const Array<t_dim>& alpha, double K): alpha_(alpha), K_(K) {}

				double operator()(const Process<t_dim>& proc){
								prod = dot(alpha_, proc_values(proc));
								return (prod > K_) ? (prod - K_): 0;
				}

protected:
				const Array<t_dim> alpha_;
				const double K_;

				double prod;
};
*/
























#endif // PAYOFF_HPP
