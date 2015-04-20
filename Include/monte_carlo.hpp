#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include <chrono>
#include <functional>
#include <iostream>

#include "linear_estimator.hpp"



/*
	* STRUCT Monte_Carlo<Generator>
	*
	* Class to represent Monte Carlo estimation.
	*
	* Constraints on type 'Dist':
	*						- typedef 'result_type'
	*						- operator(Generator&) yielding 'result_type'
	*
	*/
template <typename Dist>
struct Monte_Carlo : public Linear_Estimator {

				typedef typename Dist::result_type result_type;
				typedef std::function<double(const result_type&)> Func_Type;

				Monte_Carlo(const Dist& dist, const Func_Type& func):
								dist_(dist), func_(func) {}


				template<typename Generator>
				double operator()(Generator& gen, unsigned M){

								reinit();
								double x = 0;

								auto time_start = std::chrono::steady_clock::now();
								for(unsigned m=0; m<M; ++m) {
												x = func_(dist_(gen()));
												sum_ += x;
												sum_of_squares_ += x*x;
								}

								sample_size_ += M;
								time_span_ = std::chrono::steady_clock::now() - time_start;

								return mean_est();
				}


protected:
				Dist dist_;
				Func_Type func_;
};


/*
	* FUNC make_mc<Dist>
	*
	* Make a Monte_Carlo<Dist> object from a distribution and a function
	*
	*/
template<typename Dist>
Monte_Carlo<Dist>
make_mc(const Dist& dist, const typename Monte_Carlo<Dist>::Func_Type& func){
				return Monte_Carlo<Dist>(dist,func);
}





















#endif // MONTE_CARLO_HPP
