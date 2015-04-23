/*************************************************************************************
	* Generator types, implemented with the following constraints:
	*
	*						- template (unsigned) dim: dimension of pseudo/quasi-random vectors generated
	*						- operator() yielding Array<dim>
	*
	*************************************************************************************/
#ifndef GENERATOR_HPP
#define GENERATOR_HPP

// Pseudo-random generators
#include "uniform_generator.hpp"

// QMC generators
#include "kakutani.hpp"
#include "tore.hpp"
#include "faure.hpp"

#endif // GENERATOR_HPP
