#include <iostream>
#include <fstream>
#include <chrono>
#include <functional>
#include <random>
#include <cmath>

#include "rand_var.hpp"
#include "process.hpp"
#include "monte_carlo.hpp"
#include "kakutani.hpp"
#include "faure.hpp"
#include "tore.hpp"
#include "uniform_generator.hpp"
#include "payoff.hpp"
#include "scripts.hpp"

void test0()
{
				Uniform_Gen<2> u_gen;
				Array<2> a = u_gen();

				auto G = make_rvar(Gaussian_Ind<2>(),Halton_Fast<2>(a));

				for(int i=0; i<1000; ++i)
				{
								G();
								std::cout << G.current().at(0) << " " << G.current().at(1) << std::endl;
				}
}

void test1()
{
				Array<20> times = eq_spaced_times<20>(1);

				auto B = make_rvar(StdBrownian<20>(times),Uniform_Gen<20>());

				for(int i=0; i<100; ++i){
								auto val = B();
								for(auto it=val.begin(); it!=val.end(); ++it)
												std::cout  << it->first << "    " << it->second << std::endl;

								std::cout << std::endl;
				}
}

void test2()
{
				Array<20> times = eq_spaced_times<20>(1);

				auto B = make_rvar(StdBrownian<20>(times),SQRT<20>());

				for(int i=0; i<100; ++i){
								auto val = B();
								for(auto it=val.begin(); it!=val.end(); ++it)
												std::cout  << it->first << "    " << it->second << std::endl;

								std::cout << std::endl;
				}
}

void test3()
{
				Array<20> times = eq_spaced_times<20>(1);

				auto BS = make_rvar(Black_Scholes<20>(0.4,1,1,times),Halton<20>());

				for(int i=0; i<100; ++i){
								auto val = BS();
								for(auto it=val.begin(); it!=val.end(); ++it)
												std::cout  << it->first << "    " << it->second << std::endl;

								std::cout << std::endl;
				}
}

void testMC1()
{
				Array<20> times = eq_spaced_times<20>(1);

				auto dist = compose_dist(Black_Scholes<20>(0.6,1,1,times),
																													Last_Value<20>());
				auto MC = make_mc(dist);

				constexpr unsigned M = 250;
				constexpr unsigned N = 500;

				// Uniform pseudo-random numbers generator
				std::cout << "Uniform generator and Nb=" << M*N << std::endl;
				Uniform_Gen<20> u_gen;
				MC(u_gen,M*N);
				std::cout << MC << std::endl << std::endl;

				// Halton QMC
				std::cout << "Halton generator and M=" << M*N << std::endl;
				Halton_Fast<20> hal_gen(u_gen());
				MC(hal_gen,M*N);
				std::cout << MC.mean_est() << std::endl <<
																	MC.time() << std::endl << std::endl;

				// Faure QMC
				std::cout << "Faure generator and Nb=" << M*N << std::endl;
				Faure<20> fau_gen;
				MC(fau_gen,M*N);
				std::cout << MC.mean_est() << std::endl <<
																	MC.time() << std::endl << std::endl;


				// SQRT QMC
				std::cout << "SQRT generator and M=" << M*N << std::endl;
				SQRT<20> sqrt_gen;
				MC(sqrt_gen,M*N);
				std::cout << MC.mean_est() << std::endl  <<
																	MC.time() << std::endl << std::endl;

				// Shifted SQRT QMC
				std::cout << "Shifted SQRT generator; M=" << M << " ; N=" << N << std::endl;
				auto sqmc = make_shifted_qmc(N,dist,sqrt_gen);
				auto MC_shifted = make_mc(sqmc);

				MC_shifted(u_gen,M);
				std::cout << MC_shifted << std::endl << std::endl;

				// Random start Halton
				std::cout << "Random start Halton generator; M=" << M << " ; N=" << N << std::endl;
				auto rdHalt = make_randStart_halton(N,dist);
				auto MC_rdStartHalt = make_mc(rdHalt);

				MC_rdStartHalt(u_gen,M);
				std::cout << MC_rdStartHalt << std::endl << std::endl;

				// Benchmark
				std::cout << "Real expectation is " << std::exp(0.6) << std::endl;

}



int main(){

//				test0();

//				test1();

//				test2();

//				test3();

//				testMC1();

				compare_ciAndTime_MCvsRQMC();




    return 0;
}

















