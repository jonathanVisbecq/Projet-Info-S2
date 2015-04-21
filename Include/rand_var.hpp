#ifndef RAND_VAR_HPP
#define RAND_VAR_HPP

#define _USE_MATH_DEFINES
#include <cmath>


#include <functional>
#include <algorithm>
#include <iostream>
#include <type_traits>

#include "my_array.hpp"



/*
	* STRUCT Rand_Var
	*
	* Constraints on type 'Dist':
	*						- typedef 'result_type'
	*						- static constexpr 'dim_alea' : the number of quasi/pseudo
	*						  random numbers needed to yield
	*						- operator(const Array<Dist::dim_alea>&) yielding 'result_type'
	*
	* Constraint on type 'Generator':
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
	* Bind a generator (of real in a unit hypercube) to a distribution.
	*
	*/
template<typename Dist, typename Generator>
Rand_Var<Dist,Generator>
make_rvar(const Dist& dist, const Generator& gen)
{
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
				static constexpr unsigned dim_alea = 1;

				Uniform(double min, double max):
								min_(min), max_(max) {}

				result_type operator()(const Array<dim_alea>& pt) const
				{
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
				static constexpr unsigned dim_alea = ((dim % 2)==0) ? dim: dim+1;

				Gaussian_Ind(Array<dim> mean, Array<dim> std):
								mean_(mean), std_(std) {}


				result_type operator()(const Array<dim_alea>& pt){

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


/*
	* STRUCT Shifted_QMC<Dist,QMC_Generator>
	*
	* Shifted quasi-monte carlo
	*
	*/
template<typename Dist,typename QMC_Generator>
struct Shifted_QMC{

				static_assert(std::is_same<double,typename Dist::result_type>::value,
																		"Distribution's result_type typedef should be double for shifted QMC");

				typedef double result_type;
				static constexpr unsigned dim_alea = Dist::dim_alea;

				Shifted_QMC() = delete;
				Shifted_QMC(unsigned N,const Dist& dist,QMC_Generator& gen):
								N_(N), dist_(dist), QMCk_(N)
				{
								for(int k=0; k<N_; ++k)
												QMCk_.at(k) = gen();
				}

				double operator()(const Array<dim_alea>& pt)
				{
								double x = 0;
								for(int k=0; k<N_; ++k)
								{
												// Add both points and take the fractional part coordinate-wise
												for(int l=0; l<dim_alea; ++l)
												{
																a.at(l) = pt.at(l) + QMCk_.at(k).at(l);
																a.at(l) -= (a.at(l)>=1) ? 1 : 0;
												}

												x += dist_(a);
								}

								return x / N_;
				}


protected:
				unsigned N_;
				Dist dist_;
				std::vector<Array<dim_alea> > QMCk_;

				Array<dim_alea> a;
};


template<typename Dist,typename Generator>
Shifted_QMC<Dist,Generator>
make_shifted_qmc(unsigned N,
																	const Dist& dist,
																	Generator& gen)
{
				return Shifted_QMC<Dist,Generator>(N,dist,gen);
}




















#endif // RAND_VAR_HPP
