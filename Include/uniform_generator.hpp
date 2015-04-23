#ifndef UNIFORM_HPP
#define UNIFORM_HPP

#include <random>
#include <chrono>

#include "my_array.hpp"

inline std::mt19937_64& generator(){

				static double seed = std::chrono::system_clock::now().time_since_epoch().count();
//				static double seed = 0;
				static std::mt19937_64 generator = std::mt19937_64(seed);
				return generator;
}



/*************************************************************************************
	* STRUCT Uniform_Gen<dim>
	*
	* Multidimensional uniform random variable over [0,1]^dim
	*
	*************************************************************************************/
template<unsigned dim>
struct Uniform_Gen{

				typedef Array<dim> result_type;

				Uniform_Gen(): U_(0,1), pt_(){}

				result_type operator()()
				{
								for(auto it=pt_.begin(); it!=pt_.end(); ++it)
												*it = U_(generator());

								return pt_;
				}


protected:
				std::uniform_real_distribution<double> U_;
				result_type pt_;
};

#endif // UNIFORM_HPP
