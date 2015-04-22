#include "scripts.hpp"

void compare_ciAndTime_MCvsRQMC(std::ostream &stream)
{
				int M = 1000;
				constexpr int N_max = 500;

				Array<N_max> mc_time, mc_ci, rqmc_time, rqmc_ci = makeFill<N_max>(0.);


				Unit_Triangle test_func;
				SQRT<2> sqrt_gen;
				Uniform_Gen<2> u_gen;

				// Monte Carlo object
				auto MC = make_mc(test_func);

				std::cout << "Starting simulations..." << std::endl;
				for(int n=1; n<N_max; ++n)
				{
								// Monte Carlo object out of a RQMC variable
								auto RQMC = make_mc(make_shifted_qmc(n,test_func,sqrt_gen));

								// Launch Monte Carlo simulations
								MC(u_gen,n*M);
								RQMC(u_gen,M);

								// Record half-confidence intervals and computation time
								mc_ci.at(n-1) = MC.ci();
								rqmc_ci.at(n-1) = RQMC.ci();
								mc_time.at(n-1) = MC.time();
								rqmc_time.at(n-1) = RQMC.time();
				}


				std::cout << "Writing to stream..." << std::endl;
				for(int i=1; i<N_max; ++i)
				{
							stream << i*M << " "       << mc_ci.at(i-1)   << " " << rqmc_ci.at(i-1)   <<
																								" "	      << mc_time.at(i-1) << " " << rqmc_time.at(i-1) <<
																								std::endl << std::endl;
				}
}
