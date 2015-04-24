#include "scripts.hpp"

void compare_ciAndTime_MCvsRQMC()
{
				int M = 1000;
				constexpr int N_max = 1000;

				// Container to store half confidence interval and computation time
				// For MC
				Array<N_max> mc_time, mc_ci = makeFill<N_max>(0.);
				// For shifted QMC
				Array<N_max> sqmc_time, sqmc_ci = makeFill<N_max>(0.);
				// For random-start halton
				Array<N_max> rsHalt_time, rsHalt_ci = makeFill<N_max>(0.);


				// Payoff function, QMC generator for randomized QMC and pseudo-random number generator
				Unit_Triangle test_func;
				SQRT<2> sqrt_gen;
				Uniform_Gen<2> u_gen;


				std::cout << "Starting simulations..." << std::endl;
				for(int n=1; n<=N_max; ++n)
				{
								// Monte Carlo object
								auto MC = make_mc(test_func);
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
								rsHalt_ci.at(n-1) = RSHALT.ci();
								mc_time.at(n-1) = MC.time();
								sqmc_time.at(n-1) = SQMC.time();
								rsHalt_time.at(n-1) = RSHALT.time();
				}

				std::cout << "Writing to stream..." << std::endl;

				std::ofstream stream("../Data/compare_ciAndTime_MCvsRQMC.dat");
				for(int n=1; n<=N_max; ++n)
				{
							stream << n*M << " " << mc_ci.at(n-1)   << " " << sqmc_ci.at(n-1)   << " " << rsHalt_ci.at(n-1)
														<< " "	<< mc_time.at(n-1) << " " << sqmc_time.at(n-1) << " " << rsHalt_time.at(n-1)
														<< std::endl;
				}
				stream.close();
}

void compare_moments()
{
				constexpr unsigned M = 1000;
				constexpr unsigned N_max = 1000;

				std::cout << "Function 1..." << std::endl;
				compare_moments<M,N_max,Func_1>("func1");
				std::cout << "Function 2..." << std::endl;
				compare_moments<M,N_max,Func_2>("func2");
				std::cout << "Function 3..." << std::endl;
				compare_moments<M,N_max,Func_3>("func3");
				std::cout << "Function 4..." << std::endl;
				compare_moments<M,N_max,Func_4>("func4");
				std::cout << "Function 5..." << std::endl;
				compare_moments<M,N_max,Func_5>("func5");
}













