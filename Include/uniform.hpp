#ifndef UNIFORM_HPP
#define UNIFORM_HPP

#include "my_array.hpp"

/*************************************************************************************
	* STRUCT Uniform
	* Distribution
	*
	* Uniform distribution over an interval
	*
	*************************************************************************************/
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

#endif // UNIFORM_HPP
