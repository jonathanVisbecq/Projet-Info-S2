#include "scripts.hpp"

void compare_ciAndTime_MCvsRQMC()
{
				int M = 1000;
				constexpr int N_max = 500;

				// Container to store half confidence interval and computation time
				// For MC
				Array<N_max> mc_time, mc_ci = makeFill<N_max>(0.);
				// For shifted QMC
				Array<N_max> sqmc_time, sqmc_ci = makeFill<N_max>(0.);
				// For random-start halton
				Array<N_max> rsHalt_time, rsHal_ci = makeFill<N_max>(0.);


				// Payoff function, QMC generator for randomized QMC and pseudo-random number generator
				Unit_Triangle test_func;
				SQRT<2> sqrt_gen;
				Uniform_Gen<2> u_gen;

				// Monte Carlo object
				auto MC = make_mc(test_func);

				std::cout << "Starting simulations..." << std::endl;
				for(int n=1; n<N_max; ++n)
				{
								// Shifted QMC
								auto SQMC = make_mc(make_shifted_qmc(n,test_func,sqrt_gen));
								// Random-start Halton
								auto RSHALT = make_mc(make_randStart_halton(n,test_func));

								// Launch Monte Carlo simulations
								MC(u_gen,n*M);
								SQMC(u_gen,M);
								RSHALT(u_gen,M);

								// Record half-confidence intervals and computation time
								mc_ci.at(n-1) = MC.ci();
								sqmc_ci.at(n-1) = SQMC.ci();
								rsHal_ci.at(n-1) = RSHALT.ci();
								mc_time.at(n-1) = MC.time();
								sqmc_time.at(n-1) = SQMC.time();
								rsHalt_time.at(n-1) = RSHALT.time();
				}

				std::cout << "Writing to stream..." << std::endl;

				std::ofstream stream("../Data/compare_ciAndTime_MCvsRQMC.dat");
				for(int i=1; i<N_max; ++i)
				{
							stream << i*M << " " << mc_ci.at(i-1)   << " " << sqmc_ci.at(i-1)   << " " << rsHal_ci.at(i-1)    <<
																								" "	<< mc_time.at(i-1) << " " << sqmc_time.at(i-1) << " " << rsHalt_time.at(i-1) <<
																								std::endl << std::endl;
				}
				stream.close();
}
