#ifndef SCRIPTS_HPP
#define SCRIPTS_HPP

#include <iostream>
#include <fstream>

#include "rand_var.hpp"
#include "uniform_generator.hpp"
#include "tore.hpp"
#include "monte_carlo.hpp"



/*
	* First test of Tu04
	*
	* Compare confidence intervals and computation times of pure MC with L = M*n simulations
	* and RQMC using SQRT sequences with M simulations and n quasi-random numbers.
	*
	*/
struct Unit_Triangle{

				typedef double result_type;
				static constexpr unsigned dim_alea = 2;

				result_type operator()(const Array<dim_alea>& pt) { return (pt.at(1)>pt.at(0)) ? 1 : 0; }
};


void compare_ciAndTime_MCvsRQMC();


/*
	* Second test of Tu04
	*
	* Analyse convergence rate of second and third normalized moments with L = M*n simulations
	* for MC and M simulations using n quasi-random numbers each for RQMC
	*
	*
	*/


#endif // SCRIPTS_HPP
