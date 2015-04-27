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
				constexpr unsigned N_max = 1500;

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

void stratification_reproduce_convergence()
{

				// Probabilities for each stratum
				Array<10> probs = makeFill<10>(0.1);
				// Real-valued limits of each stratum
				Array<9> y;
				y.at(0) = -1.282;
				y.at(1) = -0.8416;
				y.at(2) = -0.5244;
				y.at(3) = -0.2533;
				y.at(4) = 0;
				y.at(5) = 0.2533;
				y.at(6) = 0.5244;
				y.at(7) = 0.8416;
				y.at(8) = 1.282;
				// Projection direction
				Array<1> u = makeFill<1>(1.);

				// Generator of pseudo-random variables
				Uniform_Gen<Normal_Interval<1>::dim_alea> gen;

				// Drawing budget for each step (first four values from the article EJ07)
				constexpr int len = 8;
				std::array<unsigned,len+1> N;
				N.at(0) = 0;
				N.at(1) = 300;
				N.at(2) = 700;
				N.at(3) = 1300;
				N.at(4) = 3300;
				N.at(5) = 5300;
				N.at(6) = 11300;
				N.at(7) = 21300;
				N.at(8)= 31300;


				auto identity = [] (const Array<1>& x) { return x.at(0); };

				Adaptive_Stratification<Normal_Interval,1,10>  adap_stratA(identity, y, u, probs, Method::A);
				Adaptive_Stratification<Normal_Interval,1,10>  adap_stratB(identity, y, u, probs, Method::B);

				std::ofstream stream("../Data/stratification_reproduce_convergence.dat");

				for(int i=1; i<=len; ++i)
				{
								std::cout << i << std::endl;
								adap_stratB(gen, N.at(i) - N.at(i-1));
								adap_stratA(gen, N.at(i) - N.at(i-1));
								stream << adap_stratA.sample_size() << " "
															<< (double)adap_stratA.N_strata().at(4) / adap_stratA.sample_size() << " "
															<< (double)adap_stratB.N_strata().at(4) / adap_stratB.sample_size()
															<<  std::endl;
				}

				stream.close();
}

void stratification_reproduce_variance()
{
				// Probabilities for each stratum
				Array<10> probs = makeFill<10>(0.1);
				// Real-valued limits of each stratum
				Array<9> y;
				y.at(0) = -1.282;
				y.at(1) = -0.8416;
				y.at(2) = -0.5244;
				y.at(3) = -0.2533;
				y.at(4) = 0;
				y.at(5) = 0.2533;
				y.at(6) = 0.5244;
				y.at(7) = 0.8416;
				y.at(8) = 1.282;
				// Projection direction
				Array<1> u = makeFill<1>(1.);

				// Generator of pseudo-random variables
				Uniform_Gen<Normal_Interval<1>::dim_alea> gen;

				// Drawing budget for each step (first four values from the article EJ07)
				constexpr int len = 8;
				std::array<unsigned,len+1> N;
				N.at(0) = 0;
				N.at(1) = 300;
				N.at(2) = 700;
				N.at(3) = 1000;
				N.at(4) = 3000;
				N.at(5) = 5000;
				N.at(6) = 11000;
				N.at(7) = 21000;
				N.at(8)= 31000;

				int L = 10000;
				Array<len> sum_est_A = makeFill<len>(0.);
				Array<len> sum_est_B = makeFill<len>(0.);
				Array<len> sum_squares_est_A = makeFill<len>(0.);
				Array<len> sum_squares_est_B = makeFill<len>(0.);
				Array<len> std_est_from_MC_A = makeFill<len>(0.);
				Array<len> std_est_from_algo_A = makeFill<len>(0.);
				Array<len> std_est_from_MC_B = makeFill<len>(0.);
				Array<len> std_est_from_algo_B = makeFill<len>(0.);

				auto identity = [] (const Array<1>& x) { return x.at(0); };

				Adaptive_Stratification<Normal_Interval,1,10>  adap_stratA(identity, y, u, probs, Method::A);
				Adaptive_Stratification<Normal_Interval,1,10>  adap_stratB(identity, y, u, probs, Method::B);


				for(int n=0; n<L; ++n)
				{
								std::cout << n << std::endl;
								adap_stratA.reinit();
								adap_stratB.reinit();

								for(int k=1; k<=len; ++k)
								{
												adap_stratA(gen, N.at(k) - N.at(k-1));
												sum_est_A.at(k-1) += adap_stratA.mean_est();
												sum_squares_est_A.at(k-1) += adap_stratA.mean_est()*adap_stratA.mean_est();

												if(n==0)
																std_est_from_algo_A.at(k-1) = adap_stratA.asymp_std_est();

												adap_stratB(gen, N.at(k) - N.at(k-1));
												sum_est_B.at(k-1) += adap_stratB.mean_est();
												sum_squares_est_B.at(k-1) += adap_stratB.mean_est()*adap_stratB.mean_est();

												if(n==0)
																std_est_from_algo_B.at(k-1) = adap_stratB.asymp_std_est();
								}
				}

				for(int k=0; k<len; ++k)
				{
								std_est_from_MC_A.at(k) = std::sqrt( (1./L) * sum_squares_est_A.at(k) - std::pow((1./L)*sum_est_A.at(k), 2) );
								std_est_from_MC_B.at(k) = std::sqrt( (1./L) * sum_squares_est_B.at(k) - std::pow((1./L)*sum_est_B.at(k), 2) );
				}


				std::ofstream stream("../Data/stratification_reproduce_variance.dat");

				for(int k=0; k<len; ++k)
				{
								stream << 	N.at(k+1) << " " << std_est_from_MC_A.at(k) * std::sqrt(N.at(k+1))
																													<< " " << std_est_from_MC_B.at(k) * std::sqrt(N.at(k+1))
																													<< " " << std_est_from_algo_A.at(k)
																													<< " " << std_est_from_algo_B.at(k)
																													<< std::endl;
				}


				stream.close();
}

void stratification_reproduce_time()
{
				// Probabilities for each stratum
				Array<10> probs = makeFill<10>(0.1);
				// Real-valued limits of each stratum
				Array<9> y;
				y.at(0) = -1.282;
				y.at(1) = -0.8416;
				y.at(2) = -0.5244;
				y.at(3) = -0.2533;
				y.at(4) = 0;
				y.at(5) = 0.2533;
				y.at(6) = 0.5244;
				y.at(7) = 0.8416;
				y.at(8) = 1.282;
				// Projection direction
				Array<1> u = makeFill<1>(1.);

				// Generator of pseudo-random variables
				Uniform_Gen<Normal_Interval<1>::dim_alea> gen;

				// Drawing budget for each step (first four values from the article EJ07)
				std::array<unsigned,4> N;
				N.at(0) = 300;
				N.at(1) = 1300;
				N.at(2) = 11300;
				N.at(3)= 31300;

				auto identity = [] (const Array<1>& x) { return x.at(0); };

				// Stratification objects
				Adaptive_Stratification<Normal_Interval,1,10>  adap_strat(identity, y, u, probs, Method::A);
				Stratification<Normal_Interval,1,10> strat(identity, y, u, probs);

				// Number of estimations
				int L = 1000;
				double sum_est_adap = 0., sum_est = 0.;
				double sum_squares_est_adap = 0., sum_squares_est = 0.;
				double time_adap = 0., time = 0.;

				for(int n=0; n<L; ++n)
				{
								adap_strat(gen,N);
								sum_est_adap += adap_strat.mean_est();
								sum_squares_est_adap += adap_strat.mean_est()*adap_strat.mean_est();
								time_adap += adap_strat.time_span();

								strat(gen,N.at(3));
								sum_est += strat.mean_est();
								sum_squares_est += strat.mean_est()*strat.mean_est();
								time += strat.time_span();
				}

				double var_adap = (1./L) * sum_squares_est_adap - std::pow((1./L)*sum_est_adap, 2);
				double var = (1./L) * sum_squares_est - std::pow((1./L)*sum_est, 2);

				time_adap /= L;
				time /= L;

				std::ofstream stream("../Data/stratification_reproduce_time.dat");
				stream << "Variance of the adaptive stratification" << std::endl
											<< var_adap << std::endl
											<< "Time of the adaptive stratification" << std::endl
											<< time_adap << std::endl
											<< "Variance of the standard stratification" << std::endl
											<< var << std::endl
											<< "Time of the standard stratification" << std::endl
											<< time << std::endl
											<< "Time*Variance for the adaptive stratification" << std::endl
											<< var_adap * time_adap << std::endl
											<< "Time*Variance for the standard stratification" << std::endl
											<< var * time << std::endl;
				stream.close();
}

void stratification_asian()
{
				{ std::ofstream("../Data/stratification_asian_call.dat"); }
				{ std::ofstream("../Data/stratification_asian_put.dat"); }

				std::cout << "Call" << std::endl;
				std::cout << "Dim: 16 and " << "K: 45" << std::endl;
				stratification_asian<16,45>(Type::CALL);
				std::cout << "Dim: 16 and " << "K: 50" << std::endl;
				stratification_asian<16,50>(Type::CALL);
				std::cout << "Dim: 16 and " << "K: 55" << std::endl;
				stratification_asian<16,55>(Type::CALL);

				std::cout << "Dim: 64 and " << "K: 45" << std::endl;
				stratification_asian<64,45>(Type::CALL);
				std::cout << "Dim: 64 and " << "K: 50" << std::endl;
				stratification_asian<64,50>(Type::CALL);
				std::cout << "Dim: 64 and " << "K: 55" << std::endl;
				stratification_asian<64,55>(Type::CALL);

				std::cout << std::endl << "Put" << std::endl;
				std::cout << "Dim: 16 and " << "K: 45" << std::endl;
				stratification_asian<16,45>(Type::PUT);
				std::cout << "Dim: 16 and " << "K: 50" << std::endl;
				stratification_asian<16,50>(Type::PUT);
				std::cout << "Dim: 16 and " << "K: 55" << std::endl;
				stratification_asian<16,55>(Type::PUT);

				std::cout << "Dim: 64 and " << "K: 45" << std::endl;
				stratification_asian<64,45>(Type::PUT);
				std::cout << "Dim: 64 and " << "K: 50" << std::endl;
				stratification_asian<64,50>(Type::PUT);
				std::cout << "Dim: 64 and " << "K: 55" << std::endl;
				stratification_asian<64,55>(Type::PUT);
}

void rqmc_asian()
{
				std::cout << "Call" << std::endl;
				std::cout << "Dim: 16 and " << "K: 50" << std::endl;
				rqmc_asian<16,50,Asian_Call>();

				std::cout << std::endl << "Put" << std::endl;
				std::cout << "Dim: 64 and " << "K: 45" << std::endl;
				rqmc_asian<64,45,Asian_Put>();
}


void compare_on_gaussian()
{
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

				// Direction of stratification using bisection method from [GHS99]
				Array<1> u = makeFill<1>(1.);

				// Generator
				Uniform_Gen<Normal_Interval<1>::dim_alea> u_gen_strat;


				// Adaptive stratification
				auto identity = [] (const Array<1>& x) { return x.at(0); };
				Adaptive_Stratification<Normal_Interval,1,nb_strats> adap_stratA(identity, y, u, probs, Method::A);

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
				for(unsigned k=0; k<N_strat.size(); ++k)
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

				Uniform_Gen<Gaussian_Ind<1>::dim_alea> u_gen_rqmc;

				auto dist = compose_dist(Gaussian_Ind<1>(), identity);

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

					for(unsigned i=0; i<I.size(); ++i)
					{
									for(unsigned ni=0; ni<NI.size(); ++ni)
									{
													unsigned n = NI.at(ni) / I.at(i);
													std::cout << I.at(i) << " " << n << std::endl;

													SQRT<Gaussian_Ind<1>::dim_alea> sqrt_gen;
													auto MC_sqmc = make_mc(make_shifted_qmc(n, dist, sqrt_gen));
													auto MC_rdStart = make_mc(make_randStart_halton(n, dist));

													MC_sqmc(u_gen_rqmc, I.at(i));
													MC_rdStart(u_gen_rqmc, I.at(i));

													N_rqmc.push_back(n);
													I_rqmc.push_back(I.at(i));
													times_sqmc.push_back(MC_sqmc.time());
													times_rdStart.push_back(MC_rdStart.time());
													ci_be_sqmc.push_back(MC_sqmc.ci_berry_essen());
													ci_be_rdStart.push_back(MC_rdStart.ci_berry_essen());
													mean_sqmc.push_back(MC_sqmc.mean_est());
													mean_rdStart.push_back(MC_rdStart.mean_est());
									}
					}



					std::ofstream stream_strat("../Data/compare_gaussian_strat.dat");
					for(unsigned i=0; i<N_strat.size(); ++i)
					{
									stream_strat << N_strat.at(i) << " " << times_stratA.at(i)
																																							<< " " << ci_stratA.at(i)
																																							<< std::endl;
					}
					stream_strat.close();

					std::ofstream stream_rqmc("../Data/compare_gaussian_rqmc.dat");
					for(unsigned i=0; i<N_rqmc.size(); ++i)
					{
									stream_rqmc << N_rqmc.at(i) << " " << I_rqmc.at(i) << " " << times_sqmc.at(i)
																																																													<< " " << times_rdStart.at(i)
																																																													<< " " << ci_be_sqmc.at(i)
																																																													<< " " << ci_be_rdStart.at(i)
																																																													<< std::endl;
					}
					stream_rqmc.close();
}

void compare_on_asian()
{
				compare_on_asian<64,55,Asian_Call>();
}
