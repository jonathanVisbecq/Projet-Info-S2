#ifndef UNIFORM_GENERATOR_HPP
#define UNIFORM_GENERATOR_HPP

#include "stl_headers.hpp"
#include "my_array.hpp"

/*---------------------------------------------------------------------------------------
	* FUNC generator
	*
	* Seed and create a Mersenne Twister 19937 generator (implementation from STL) only once
	* and return a refernece on the seeded generator.
	*--------------------------------------------------------------------------------------*/
inline std::mt19937_64& generator(){

//				static double seed = 0;
				static double seed = std::chrono::system_clock::now().time_since_epoch().count();
				static std::mt19937_64 generator = std::mt19937_64(seed);
				return generator;
}



/*************************************************************************************
	* STRUCT Uniform_Gen<dim>
	* Generator
	*
	* Multidimensional uniform random variable over [0,1]^dim.
	*
	* Thanks to the use of function generator() all instances of this class behave like
	* independant uniform random variables.
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


/*************************************************************************************
	* STRUCT Uniform_Gen_Fixed<dim>
	* Generator
	*
	* Multidimensional uniform random variable over [0,1]^dim.
	*
	* Here every instance wil generate the same sequence of pseudo-random numbers
	*
	*************************************************************************************/
template<unsigned dim>
struct Uniform_Gen_Fixed{

				typedef Array<dim> result_type;

				Uniform_Gen_Fixed(double seed = 0.):
								U_(0,1), pt_()
				{
								gen_ = std::mt19937_64(seed);
				}

				result_type operator()()
				{
								for(auto it=pt_.begin(); it!=pt_.end(); ++it)
												*it = U_(gen_);

								return pt_;
				}


protected:
				std::uniform_real_distribution<double> U_;
				std::mt19937_64 gen_;
				result_type pt_;
};

#endif // UNIFORM_GENERATOR_HPP
