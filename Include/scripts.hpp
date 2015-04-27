#ifndef SCRIPTS_HPP
#define SCRIPTS_HPP


#include "stl_headers.hpp"
#include "distribution.hpp"
#include "generator.hpp"
#include "payoff.hpp"
#include "rand_var.hpp"
#include "uniform_generator.hpp"
#include "tore.hpp"
#include "monte_carlo.hpp"
#include "stratification.hpp"
#include "adaptive_stratification.hpp"



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
	*/
struct  Func_1{

				typedef double result_type;
				static constexpr unsigned dim_alea = 10;

				result_type operator()(const Array<dim_alea>& pt)
				{
								result = 1;
								for(i=0; i<dim_alea; ++i)
												result *= 0.5*M_PI*std::sin(M_PI*pt.at(i));

								return result;
				}

				double result;
				unsigned i;
};

struct  Func_2{

				typedef double result_type;
				static constexpr unsigned dim_alea = 10;

				result_type operator()(const Array<dim_alea>& pt)
				{
								result = 1;
								for(i=0; i<dim_alea; ++i)
												result *= 12*(pt.at(i) - 0.5)*(pt.at(i) - 0.5);

								return result;
				}

				double result;
				unsigned i;
};

struct  Func_3{

				typedef double result_type;
				static constexpr unsigned dim_alea = 10;

				result_type operator()(const Array<dim_alea>& pt)
				{
								result = 0;
								for(i=0; i<dim_alea; ++i)
												result += pt.at(i);

								return (1./5.) * result;
				}

				double result;
				unsigned i;
};

struct  Func_4{

				typedef double result_type;
				static constexpr unsigned dim_alea = 10;

				result_type operator()(const Array<dim_alea>& pt)
				{
								result = 1;
								for(i=0; i<dim_alea; ++i)
												result *= std::exp(-std::abs(pt.at(i) - 0.5)) / (2 - 2*std::exp(0.5));

								return result;
				}

				double result;
				unsigned i;
};

struct  Func_5{

				typedef double result_type;
				static constexpr unsigned dim_alea = 10;

				result_type operator()(const Array<dim_alea>& pt)
				{
								result = 1;
								for(i=0; i<dim_alea; ++i)
												result *= (pt.at(i) > 0.5) ? 2: 0;

								return result;
				}

				double result;
				unsigned i;
};

template<unsigned M,unsigned N_max,typename Func>
void compare_moments(const char* name)
{
				// To store second moment (variance)
				Array<N_max> mc_2, sqmc_2, rsHalt_2;
				// To store third normalized moment (skewness)
				Array<N_max> mc_3, sqmc_3, rsHalt_3;

				// Payoff, QMC generator for randomized QMC and pseudo-random number generator
				Func func;
				SQRT<Func::dim_alea> sqrt_gen;

				std::cout << "Starting simulations..." << std::endl;


				for(int n=1; n<=N_max; n=n+50)
				{
								std::cout << n << std::endl;

								// Same sequence of random variables for each n
								Uniform_Gen_Fixed<Func::dim_alea> u_gen;

								// Monte Carlo object
								auto MC = make_mc(func);
								// Shifted QMC
								auto SQMC = make_mc(make_shifted_qmc(n, func, sqrt_gen));
								// Random-start Halton
								auto RSHALT = make_mc(make_randStart_halton(n, func));

								// Launch Monte Carlo simulations
								MC(u_gen,n*M);
								SQMC(u_gen,M);
								RSHALT(u_gen,M);

								// Record moments
								mc_2.at(n-1) = MC.var_est();
								sqmc_2.at(n-1) = SQMC.var_est();
								rsHalt_2.at(n-1) = RSHALT.var_est();

								mc_3.at(n-1) = MC.moment_3_est() / std::pow(MC.st_dev_est(), 3);
								sqmc_3.at(n-1) = SQMC.moment_3_est() / std::pow(SQMC.st_dev_est(), 3);
								rsHalt_3.at(n-1) = RSHALT.moment_3_est() / std::pow(RSHALT.st_dev_est(), 3);
				}

				std::cout << "Writing to stream..." << std::endl;

				std::string path = std::string("../Data/compare_moments_") +
																							std::string(name) +
																							std::string(".dat");

				std::ofstream stream(path);
				for(int n=1; n<=N_max; n=n+50)
				{
							stream << n*M << " " << mc_2.at(n-1)  << " " << sqmc_2.at(n-1)   << " " << rsHalt_2.at(n-1)
																					<< " " << mc_3.at(n-1)  << " " << sqmc_3.at(n-1)   << " " << rsHalt_3.at(n-1)
																					<< std::endl;
				}
				stream.close();
}

void compare_moments();

/*
	* First test of EJ07
	*
	* Compute average of standard normal gaussian and check convergence of proportion to optimality (tractable).
	*
	* We aim here at reproducing the numerical experiments in the article to be sure our implementation
	* is viable.
	*
	*/


void stratification_reproduce_convergence();

void stratification_reproduce_variance();

void stratification_reproduce_time();

/*
	* Second test of EJ07
	*
	* Application to asian european option
	*
	*/

/* Compute number of iterations for bisection method.
	*/
inline int bisec_nb_iteration(double a, double b, double tolerance)
{
				return std::round( std::log2(b-a) - std::log2(tolerance) - 1.);
}

/* Bisection method to find zero of function
	*/
inline double bisec(std::function<double(double)> func, double a, double b, double tolerance = 1e-12)
{
				// Assume func(a), func(b) are of opposite sign

				double c = 0, y = 0;

				// Find the num of iterations
				int N = bisec_nb_iteration(a,b,tolerance);

				for (int i=0; i<=N; i++) {
								c = (a + b) / 2.;

								y = func(c);

								if(y < 0)
												if(func(a)>0)
																b = c;
												else
																a = c;
								else if(y > 0)
												if(func(a)>0)
																a = c;
												else
																b = c;
								else
												break;

				}

				return c;
}

enum class Type { CALL, PUT };

/* Solving equn_func(y) = 0 allows to compute the stratification direction using method from
	* [GHS99]
	*/
template<unsigned dim>
double eqn_func_call(double y, double r, double sigma, double dt, double K, double S0, Array<dim>& mu)
{
				mu.at(0) = sigma * std::sqrt(dt) * (y + K) / y;
				double S = S0, sum_S = 0.;

				for(int i=1; i<dim; ++i)
				{
								mu.at(i) = mu.at(i-1) - sigma * std::sqrt(dt) * S / (dim * y);
								S *= std::exp( (r - 0.5*sigma*sigma)*dt + sigma*std::sqrt(dt)*mu.at(i) );
								sum_S += S;
				}


				return sum_S/dim - K - y;
}

template<unsigned dim>
double eqn_func_put(double y, double r, double sigma, double dt, double K, double S0, Array<dim>& mu)
{
				mu.at(0) = sigma * std::sqrt(dt) * (y - K) / y;
				double S = S0, sum_S = 0.;

				for(int i=1; i<dim; ++i)
				{
								mu.at(i) = mu.at(i-1) + sigma * std::sqrt(dt) * S / (dim * y);
								S *= std::exp( (r - 0.5*sigma*sigma)*dt + sigma*std::sqrt(dt)*mu.at(i) );
								sum_S += S;
				}


				return K - sum_S/dim - y;
}


/* Payoff
	*
	* Notations from [EJ07]
	*/
template<unsigned dim>
double f_mu(const Array<dim>& x, double r, double sigma, double dt, double K, double S0, const Array<dim>& mu, Type type)
{
				double S = S0, sum_S = 0.;
				Array<dim> x_pl_mu = x + mu;

				for(int i=0; i<dim; ++i)
				{
								S *= std::exp( (r - 0.5*sigma*sigma)*dt + sigma*std::sqrt(dt)*x_pl_mu.at(i) );
								sum_S += S;
				}

				double g = 0.;
				if(type==Type::CALL)
								g = std::exp(-r*dim*dt) * ( (1./dim)*sum_S - K );
				else
								g = std::exp(-r*dim*dt) * ( K - (1./dim)*sum_S );

				g = (g > 0.)? g : 0.;

				return g * std::exp( -scalar_prod(mu,x) - 0.5*scalar_prod(mu,mu) );
}


template<unsigned dim, unsigned K>
void stratification_asian(Type type)
{
				using namespace std::placeholders;

				double S0 = 50.;
				double r = 0.05;
				double sigma = 0.1;
				double T = 1.;
				constexpr unsigned nb_strats = 100;

				double dt = T / dim;

				// Compute strata and probabilities
				Array<nb_strats> probs = makeFill<nb_strats>(1./nb_strats);

				Array<nb_strats-1> y = makeFill<nb_strats-1>(0.);
				for(unsigned i=0; i<nb_strats-1; ++i)
								y.at(i) = inv_normal_cdf((i+1.) / nb_strats);

				// Compute direction of stratification using bisection method from [GHS99]
				Array<dim> mu = makeFill<dim>(0.);

				if(type==Type::CALL)
								bisec(std::bind(eqn_func_call<dim>, _1, r, sigma, dt, K, S0, std::ref(mu)), -1e2, 2e2);
				else
								bisec(std::bind(eqn_func_put<dim>, _1, r, sigma, dt, K, S0, std::ref(mu)), -1e2, 2e2);

				Array<dim> u = mu / std::sqrt(scalar_prod(mu,mu));
				u = ( std::signbit(u.at(0)) )? ((-1.) * u) : u;

				// Payoff function (from importance sampling)
				auto payoff_mu = std::bind(f_mu<dim>, _1, r, sigma, dt, K, S0, std::ref(mu), type);

				// Generator
				Uniform_Gen<Normal_Interval<dim>::dim_alea> gen;

				// Standard sratification
				Stratification<Normal_Interval,dim,nb_strats> strat(payoff_mu, y, u, probs);
				unsigned N_total = 500000;

				strat(gen,N_total);

				// Adaptive stratification
				Adaptive_Stratification<Normal_Interval,dim,nb_strats> adap_stratA(payoff_mu, y, u, probs, Method::A);
				Adaptive_Stratification<Normal_Interval,dim,nb_strats> adap_stratB(payoff_mu, y, u, probs, Method::B);
				std::array<unsigned,4> N;
				N.at(0) = 100000;
				N.at(1) = 400000;
				N.at(2) = 500000;
				N.at(3) = 1000000;

				adap_stratA(gen,N);
				adap_stratB(gen,N);

				// Record results

				if(dim==64 && K==45.)
				{
								// Record u
								std::string path_u = (type==Type::CALL)? "../Data/stratification_asian_call_u.dat" : "../Data/stratification_asian_put_u.dat";
								std::ofstream stream_u(path_u);
								stream_u << u;
								stream_u.close();

								// Record stds estimations
								std::string path_stds = (type==Type::CALL)? "../Data/stratification_asian_call_stds.dat" : "../Data/stratification_asian_put_stds.dat";
								std::ofstream stream_stds(path_stds);
								for(unsigned i=0; i<nb_strats; ++i)
												stream_stds << adap_stratA.stds().at(i)
																								<< " " << adap_stratB.stds().at(i)
																								<< " " << adap_stratA.sum_strata().at(i) / adap_stratA.N_strata().at(i)
																								<< " " << adap_stratB.sum_strata().at(i) / adap_stratB.N_strata().at(i)
																								<< std::endl;
								stream_stds.close();
				}

				// Estimations
				std::string path = (type==Type::CALL)? "../Data/stratification_asian_call.dat" : "../Data/stratification_asian_put.dat";
				std::ofstream stream(path,std::ios_base::app);
				stream << dim << " " << K << " " << adap_stratA.mean_est()
																														<< " " << adap_stratB.mean_est()
																														<< " " << strat.mean_est()
																														<< " " << adap_stratA.var_est()
																														<< " " << adap_stratB.var_est()
																														<< " " << strat.var_est()
																														<< " " << strat.var_est() / adap_stratA.var_est()
																														<< " " << strat.var_est() / adap_stratB.var_est()
																														<< std::endl;
				stream.close();
}

void stratification_asian();



template<unsigned dim, unsigned K, template <unsigned> class Asian_Type>
void rqmc_asian()
{
				double S0 = 50.;
				double r = 0.05;
				double sigma = 0.1;
				double T = 1.;

				constexpr unsigned M = 250;
				constexpr unsigned N = 500;

				Uniform_Gen<Black_Scholes<dim>::dim_alea> u_gen;

				Array<dim> times;
				for(int i=1; i<=dim; ++i)
								times.at(i-1) = i * T / dim;

					auto dist = compose_dist(Black_Scholes<dim>(r,sigma,S0,times), Asian_Type<dim>(K, r, T));

				//	MC
				auto MC = make_mc(dist);
				MC(u_gen, 500000);

				std::cout << MC << std::endl << std::endl;


				//	Shifted QMC
				SQRT<Black_Scholes<dim>::dim_alea> sqrt_gen;
				auto sqmc = make_shifted_qmc(N, dist, sqrt_gen);
				auto MC_sqmc = make_mc(sqmc);

				MC_sqmc(u_gen, M);
				std::cout << MC_sqmc << std::endl << std::endl;

				// Random start Halton RQMC
				auto rdHalt = make_randStart_halton(N, dist);
				auto MC_rdStartHalt = make_mc(rdHalt);

				MC_rdStartHalt(u_gen, M);
				std::cout << MC_rdStartHalt << std::endl << std::endl;
}


void rqmc_asian();


void compare_on_gaussian();



template<unsigned dim, unsigned K, template <unsigned> class Asian_Type>
void compare_on_asian()
{
				using namespace std::placeholders;

				Type type = Type::CALL;
				if(std::is_same<Asian_Type<dim>,Asian_Put<dim>>::value)
								type = Type::PUT;


				double S0 = 50.;
				double r = 0.05;
				double sigma = 0.1;
				double T = 1.;

				double dt = T/dim;

				// Arrays to store computation times and confidence interval
				// (using Berry-Essen bounds for rqmc)
				std::vector<double> I_rqmc;
				std::vector<double> N_rqmc;
				std::vector<double> times_sqmc;
				std::vector<double> times_rdStart;
				std::vector<double> ci_be_sqmc;
				std::vector<double> ci_be_rdStart;
				std::vector<double> mean_sqmc;
				std::vector<double> mean_rdStart;


				std::vector<double> N_stratification;
				std::vector<double> times_stratA;
				std::vector<double> ci_stratA;
				std::vector<double> mean_stratA;

				/*
					* Stratification
					*/
				std::cout << "Stratification..." << std::endl;
				constexpr unsigned nb_strats = 100;

				// Compute strata and probabilities
				Array<nb_strats> probs = makeFill<nb_strats>(1./nb_strats);

				Array<nb_strats-1> y = makeFill<nb_strats-1>(0.);
				for(unsigned i=0; i<nb_strats-1; ++i)
								y.at(i) = inv_normal_cdf((i+1.) / nb_strats);

				// Compute direction of stratification using bisection method from [GHS99]
				Array<dim> mu = makeFill<dim>(0.);

				if(type==Type::CALL)
								bisec(std::bind(eqn_func_call<dim>, _1, r, sigma, dt, K, S0, std::ref(mu)), 1e-2, 2e2);
				else
								bisec(std::bind(eqn_func_put<dim>, _1, r, sigma, dt, K, S0, std::ref(mu)), -1e2, 2e2);

				Array<dim> u = mu / std::sqrt(scalar_prod(mu,mu));
				u = ( std::signbit(u.at(0)) )? ((-1.) * u) : u;

				// Payoff function (from importance sampling)
				auto payoff_mu = std::bind(f_mu<dim>, _1, r, sigma, dt, K, S0, std::ref(mu), type);

				// Generator
				Uniform_Gen<Normal_Interval<dim>::dim_alea> u_gen_strat;

				// Adaptive stratification
				Adaptive_Stratification<Normal_Interval,dim,nb_strats> adap_stratA(payoff_mu, y, u, probs, Method::A);

				// Number of successive drawings

				std::array<unsigned,17> N_strat;
				N_strat.at(0) = 500;
				N_strat.at(1) = 1000;
				N_strat.at(2) = 5000;
				N_strat.at(3) = 10000;
				N_strat.at(4) = 20000;
				N_strat.at(5) = 50000;
				N_strat.at(6) = 100000;
				N_strat.at(7) = 150000;
				N_strat.at(8) = 200000;
				N_strat.at(9) = 300000;
				N_strat.at(10) = 400000;
				N_strat.at(11) = 500000;
				N_strat.at(12) = 600000;
				N_strat.at(13) = 700000;
				N_strat.at(14) = 800000;
				N_strat.at(15) = 900000;
				N_strat.at(16) = 1000000;

				double N_last = 0.;
				for(int k=0; k<N_strat.size(); ++k)
				{
								std::cout << N_strat.at(k) << std::endl;

								adap_stratA(u_gen_strat, N_strat.at(k) - N_last);
								N_last = N_strat.at(k);

								N_stratification.push_back(N_last);
								times_stratA.push_back(adap_stratA.time_span());
								ci_stratA.push_back(adap_stratA.ci());
								mean_stratA.push_back(adap_stratA.mean_est());
				}

				std::cout << std::endl;

				/*
					* Randomized QMC
					*/
				std::cout << "RQMC" << std::endl;

				Uniform_Gen<Black_Scholes<dim>::dim_alea> u_gen_rqmc;

				Array<dim> times;
				for(int i=1; i<=dim; ++i)
								times.at(i-1) = i * T / dim;

					auto dist = compose_dist(Black_Scholes<dim>(r,sigma,S0,times), Asian_Type<dim>(K, r, T));

					std::array<unsigned, 3> I;
					I.at(0) = 800;
					I.at(1) = 900;
					I.at(2) = 1000;

					std::array<unsigned, 13> NI;
					NI.at(0) = 1000;
					NI.at(1) = 2500;
					NI.at(2) = 5000;
					NI.at(3) = 10000;
					NI.at(4) = 20000;
					NI.at(5) = 50000;
					NI.at(6) = 100000;
					NI.at(7) = 200000;
					NI.at(8) = 350000;
					NI.at(9) = 500000;
					NI.at(10) = 800000;
					NI.at(11) = 1000000;
					NI.at(12) = 1200000;

					for(int i=0; i<I.size(); ++i)
					{
									for(int ni=0; ni<NI.size(); ++ni)
									{
													unsigned n = NI.at(ni) / I.at(i);
													std::cout << I.at(i) << " " << n << std::endl;

													SQRT<Black_Scholes<dim>::dim_alea> sqrt_gen;
													auto MC_sqmc = make_mc(make_shifted_qmc(n, dist, sqrt_gen));
													auto MC_rdStart = make_mc(make_randStart_halton(n, dist));

													MC_sqmc(u_gen_rqmc, I.at(i));
													MC_rdStart(u_gen_rqmc, I.at(i));

													N_rqmc.push_back(n);
													I_rqmc.push_back(I.at(i));
													times_sqmc.push_back(MC_sqmc.time());
													times_rdStart.push_back(MC_rdStart.time());
													ci_be_sqmc.push_back(MC_sqmc.ci());
													ci_be_rdStart.push_back(MC_rdStart.ci());
													mean_sqmc.push_back(MC_sqmc.mean_est());
													mean_rdStart.push_back(MC_rdStart.mean_est());
									}
					}



					std::ofstream stream_strat("../Data/compare_asian_strat.dat");
					for(int i=0; i<N_strat.size(); ++i)
					{
									stream_strat << N_strat.at(i) << " " << times_stratA.at(i)
																																							<< " " << ci_stratA.at(i)
																																							<< std::endl;
					}
					stream_strat.close();

					std::ofstream stream_rqmc("../Data/compare_asian_rqmc.dat");
					for(int i=0; i<N_rqmc.size(); ++i)
					{
									stream_rqmc << N_rqmc.at(i) << " " << I_rqmc.at(i) << " " << times_sqmc.at(i)
																																																													<< " " << times_rdStart.at(i)
																																																													<< " " << ci_be_sqmc.at(i)
																																																													<< " " << ci_be_rdStart.at(i)
																																																													<< std::endl;
					}
					stream_rqmc.close();
}



void compare_on_asian();













#endif // SCRIPTS_HPP


