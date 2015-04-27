#ifndef GAUSSIAN_IND_HPP
#define GAUSSIAN_IND_HPP

#include "stl_headers.hpp"
#include "my_array.hpp"

/*************************************************************************************
	* STRUCT Gaussian_Ind<dim>
	* Distribution
	*
	* Gaussian distribution for multidimensionnal gaussian whose coordinates are
	* mutually independant.
	*
	* Use Box-Muller method to be compatible with quasi-MC.
	*
	*************************************************************************************/
template<unsigned dim>
struct Gaussian_Ind{

				typedef Array<dim> result_type;
				static constexpr unsigned dim_alea = ((dim % 2)==0) ? dim: dim+1;

				Gaussian_Ind(const Array<dim>& mean = makeFill<dim>(0.),
																	const Array<dim>& std = makeFill<dim>(1.)):
								mean_(mean), std_(std) {}

				// Use more pt than necessary but easier to use this way
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
				Array<dim> mean_;
				Array<dim> std_;

				double R,theta,val1,val2;
				Array<dim> val_tp;
				int i;
};

#endif // GAUSSIAN_IND_HPP
