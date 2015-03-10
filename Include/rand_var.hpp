#ifndef RAND_VAR_HPP
#define RAND_VAR_HPP

#define _USE_MATH_DEFINES
#include <cmath>


#include <functional>
#include <algorithm>
#include <iostream>

#include "my_array.hpp"



/*
	* STRUCT Rand_Var
	*
	*
	* Constraints on type 'Dist':
	*						- typedef 'result_type'
	*						- operator(Array<dim>) yielding 'Result_Type'
	*
	* Constraint on type 'Equi_Gen':
	*						- operator() yielding Array<dim>
	*
	*/
template<typename Dist, typename Generator>
struct Rand_Var{

				typedef typename Dist::result_type result_type;

				Rand_Var(const Dist& dist, const Generator& gen):
								dist_(dist), gen_(gen), value_(dist_(gen_())) {}

				result_type operator()(){ return value_ = dist_(gen_()); }
				result_type current() const { return value_; }

protected:
				Dist dist_;
				Generator gen_;
				result_type value_;
};



/*
	* FUNC make_rvar<Dist,Generator>
	*
	* Bind a generator (of real in a unit hypercube> to a distribution.
	*
	*/
template<typename Dist, typename Generator>
Rand_Var<Dist,Generator>
make_rvar(const Dist& dist, const Generator& gen){

				return Rand_Var<Dist,Generator>(dist,gen);
}




/*
	* STRUCT Uniform
	*
	* Wrapper for uniform distribution
	*
	*/
struct Uniform{

				typedef double result_type;

				Uniform(double min, double max):
								min_(min), max_(max) {}

				result_type operator()(const Array<1>& pt) const{

								return min_ + (max_ - min_)*pt.at(0);
				}

protected:
				const double min_;
				const double max_;
};


/*
	* STRUCT Gaussian_Ind<dim>
	*
	* Gaussian distribution for multidimensionnal gaussian whose coordinates are
	* mutually independant.
	*
	* Use Box-Muller method to be compatible with quasi-MC.
	*
	*/
template<unsigned dim>
struct Gaussian_Ind{

				typedef Array<dim> result_type;
				static constexpr unsigned dim_pt(){ return ((dim % 2)==0) ? dim: dim+1; }

				Gaussian_Ind(Array<dim> mean, Array<dim> std):
								mean_(mean), std_(std) {}


				result_type operator()(const Array<dim_pt()>& pt){

								for(i = 0; i<dim; ++i){
												if(i%2 == 0){
																R = std::sqrt(- 2 * std::log(pt.at(i)));
																theta = 2 * M_PI * pt.at(i+1);
																val_tp.at(i) = R * std::sin(theta);
												}else{
																val_tp.at(i) = R * std::cos(theta);
												}
								}

								return (mean_ + val_tp) * std_;
				}

protected:
				const Array<dim> mean_;
				const Array<dim> std_;

				double R,theta,val1,val2;
				Array<dim> val_tp;
				int i;
};




#endif // RAND_VAR_HPP
