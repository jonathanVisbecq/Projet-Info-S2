/*************************************************************************************
	* Distribution types, implemented with the following constraints:
	*
	*						- typedef 'result_type'
	*						- static constexpr 'dim_alea' : the number of quasi/pseudo random numbers
	*								needed for the 'operator(const Array<Dist::dim_alea>&)'
	*						- operator(const Array<Dist::dim_alea>&) yielding 'result_type'
	*
	*************************************************************************************/
#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP


#include "uniform.hpp"
#include "gaussian_ind.hpp"

/*************************************************************************************
	* Randomized QMC implemented as distributions
*************************************************************************************/
#include "random_start_halton.hpp"
#include "shifted_qmc.hpp"

/*************************************************************************************
	* Processes implemented as distributions
	*************************************************************************************/
#include "process.hpp"


/*************************************************************************************
	* STRUCT Composed_Dist<Dist>
	* Distribution
	*
	* Allow composition of a distribution by a real function (argument's type  should be
	* the distribution's result_type
	*
	*************************************************************************************/
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

/*---------------------------------------------------------------------------------------
	* FUNC Composed_Dist<Dist>
	*
	* Make a distribution by composing a function and a distribution
	*--------------------------------------------------------------------------------------*/
template<typename Dist>
Composed_Dist<Dist>
compose_dist(const Dist& dist,
													const typename Composed_Dist<Dist>::func_type& func)
{
				return Composed_Dist<Dist>(dist,func);
}












#endif // DISTRIBUTION_HPP
