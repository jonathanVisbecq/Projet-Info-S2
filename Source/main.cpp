#include <iostream>
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



void test1(){
				Array<20> times = eq_spaced_times<20>(1);

				auto B = make_rvar(StdBrownian<20>(times),Uniform_Gen<20>());

				for(int i=0; i<100; ++i){
								auto val = B();
								for(auto it=val.begin(); it!=val.end(); ++it)
												std::cout  << it->first << "    " << it->second << std::endl;

								std::cout << std::endl;
				}
}

void test2(){
				Array<20> times = eq_spaced_times<20>(1);

				auto B = make_rvar(StdBrownian<20>(times),SQRT<20>());

				for(int i=0; i<100; ++i){
								auto val = B();
								for(auto it=val.begin(); it!=val.end(); ++it)
												std::cout  << it->first << "    " << it->second << std::endl;

								std::cout << std::endl;
				}
}

void test3(){
				Array<20> times = eq_spaced_times<20>(1);

				auto BS = make_rvar(Black_Scholes<20>(0.4,1,1,times),Halton<20>());

				for(int i=0; i<100; ++i){
								auto val = BS();
								for(auto it=val.begin(); it!=val.end(); ++it)
												std::cout  << it->first << "    " << it->second << std::endl;

								std::cout << std::endl;
				}
}

void testMC1(){

				Array<20> times = eq_spaced_times<20>(1);

				Last_Value<20> payoff;
				Black_Scholes<20> dist(0.6,1,1,times);


				auto MC = make_mc(dist,payoff);

				int M = 50000;
				std::cout << "Uniform generator and M=" << M << std::endl;

				Uniform_Gen<20> gen;
				MC(gen,M);

				std::cout << MC << std::endl << std::endl;

				std::cout << "Halton generator and M=" << M << std::endl;

				Halton<20> gen2;
				MC(gen2,M);

				std::cout << MC.mean_est() << std::endl << std::endl;

				std::cout << "Faure generator and M=" << M << std::endl;

				Faure<20> gen3;
				MC(gen3,M);

				std::cout << MC.mean_est() << std::endl << std::endl;

				std::cout << "SQRT generator and M=" << M << std::endl;

				SQRT<20> gen4;
				MC(gen4,M);

				std::cout << MC.mean_est() << std::endl << std::endl;

				std::cout << "Real expectation is " << std::exp(0.6) << std::endl;

}

template<unsigned dim>
double func(const Array<dim>& a){
				//double x = 0;

				//for(int i=0; i<dim; ++i)
								//x += a.at(i);

				return a.at(0);
}


double identity(double x){
				return x;
}

int main(){

				//test1();

				//test2();

				//test3();

				testMC1();

    return 0;
}

















