#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include <chrono>
#include <functional>
#include <iostream>
#include <type_traits>

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

				static_assert(std::is_same<double,typename Dist::result_type>::value,
																		"Distribution's result_type typedef should be double for Monte Carlo estimation");

				Monte_Carlo(const Dist& dist): dist_(dist) {}


				template<typename Generator>
				double operator()(Generator& gen, unsigned M)
				{
								reinit();
								double x = 0;

								auto time_start = std::chrono::steady_clock::now();
								for(unsigned m=0; m<M; ++m) {
												auto a = gen();
												x = dist_(a);
												sum_ += x;
												sum_of_squares_ += x*x;
								}

								sample_size_ += M;
								time_span_ = std::chrono::steady_clock::now() - time_start;

								return mean_est();
				}


protected:
				Dist dist_;
};


/*
	* FUNC make_mc<Dist>
	*
	* Make a Monte_Carlo<Dist> object from a distribution and a function
	*
	*/
template<typename Dist>
Monte_Carlo<Dist>
make_mc(const Dist& dist)
{
				return Monte_Carlo<Dist>(dist);
}





















#endif // MONTE_CARLO_HPP
