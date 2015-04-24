#ifndef STRATIFICATION_HPP
#define STRATIFICATION_HPP

#include <stdio.h>
#include <chrono>
#include <functional>
#include <iostream>
#include <numeric>
#include <cmath>

#include <iostream>

#include "my_array.hpp"
#include "monte_carlo.hpp"
#include "rand_var.hpp"


/********************************************************************************************
	* STRUCT Stratification<Dist,dim,nb_strats>
	*
	* Class to represent Monte Carlo estimation.
	*
	* 'Dist' must be a distribution type
	*
	* TODO: add a way to pass distributions for each stratum
	*
	*
	*******************************************************************************************/
template <typename Dist,unsigned dim,unsigned nb_strats>
struct Stratification {

				typedef std::function<double(const typename Dist::result_type&)> Func_Type;
    
				Stratification(const Dist& dist,
																			const Func_Type& func,
																			const Array<nb_strats-1>& y,
																			const Array<dim>& u,
																			const Array<nb_strats>& probs):
				dist_(dist), func_(func),
				y_(y), u_(u), probs_(probs) {}

				/*
					* Standard non-adaptive stratification with given weights (whose sum
					* must be 1).
					*
					* Does not provide a confidence interval.
					*/
				template<typename Generator>
				double
				operator()(Generator& gen, unsigned N, const Array<nb_strats>& weights)
				{
								reinit();

								auto time_start = std::chrono::steady_clock::now();

								// Compute allocation to each stratum
								std::array<unsigned,nb_strats> Ni;

								double sum = 0.;
								for(int i=0; i<nb_strats; i++)
								{
												Ni.at(i) = std::floor(sum + weights.at(i)*N) - std::floor(sum);
												sum += weights.at(i)*N;
								}

								// Draw realizations and compute estimator
								double partial_sum = 0.;
								for(int i=0; i<nb_strats; ++i)
								{
												partial_sum = 0.;
												for(int n=0; n<Ni.at(i); ++n)
																partial_sum += func_(generate(i,gen));

												mean_est_ += probs_.at(i) * (1/weights.at(i)) * partial_sum;
								}

								sample_size_ += N;
								mean_est_ /= sample_size();
								time_span_ += std::chrono::steady_clock::now() - time_start;

								return mean_est_;
				}
    

				/*
					* 	Standard non-adaptive stratification with proportional allocation.
					*/
				template<typename Generator>
				double
				operator()(Generator& gen, unsigned N)
				{
								return mean_est_ = operator()(gen,N,probs_);
				}


				/*
					* Adaptive optimal stratification.
					*
					* Using method a)
					*/
				template<typename Generator,size_t nb_iter>
				double
				operator()(Generator& gen, Array<nb_iter> N)
				{
								reinit();

								Array<nb_strats> sum_strata = makeFill<nb_strats>(0.);
								Array<nb_strats> sum_squares_strata = makeFill<nb_strats>(0.);
								Array<nb_strats> stds = makeFill<nb_strats>(0.);
								Array<nb_strats> m = makeFill<nb_strats>(0.);

								std::array<unsigned,nb_strats> N_strata;
								for(int i=0; i<nb_strats; ++i)
												N_strata.at(i) = 0;

								auto time_start = std::chrono::steady_clock::now();
								for(int step=0; step<nb_iter; ++step)
								{
												if(step>0)
																compute_empirical_stds(stds,sum_strata,sum_squares_strata,N_strata);

												unsigned N_last = (step>0) ? N.at(step-1) : 0;
												unsigned N_next = N.at(step);

												// Check wether all empirical stds are equal to 0
												if(std::accumulate(stds.begin(),stds.end(),0.)>0.)
																compute_m_A(m, stds, N_last, N_next);
												else
												{
																for(int i=0; i<nb_strats; ++i)
																				m.at(i) = probs_.at(i) * (N_next - N_last - nb_strats);
												}

												draw_realizations(gen, compute_m_tilde(m), N_strata, sum_strata, sum_squares_strata);

								}

								for(int i=0; i<nb_strats; ++i)
								{
												mean_est_ += probs_.at(i) * (1./N_strata.at(i)) * sum_strata.at(i);
												std::cout << mean_est_ << std::endl << std::endl;
												asymp_std_est_ += probs_.at(i) * stds.at(i);
												sample_size_ += N_strata.at(i);
								}

								time_span_ += std::chrono::steady_clock::now() - time_start;


								return mean_est_;
				}

				/*
					* Compute the standard deviation on each stratum
					*/
				void
				compute_empirical_stds(Array<nb_strats>& stds,
																											const Array<nb_strats>& sum_strata,
																											const Array<nb_strats>& sum_squares_strata,
																											const std::array<unsigned,nb_strats>& N_strata)
				{
								for(int i=0; i<nb_strats; ++i)
												stds.at(i) = std::sqrt( sum_squares_strata.at(i)/N_strata.at(i) - std::pow(sum_strata.at(i)/N_strata.at(i), 2) );
				}

				/*
					* Compute the m_i using method a), given that all stds are not equal to 0
					*/
				void
				compute_m_A(Array<nb_strats>& m,
																const Array<nb_strats>& stds,
																unsigned N_last, unsigned N_next)
				{
								for(int i=0; i<nb_strats ;i++)
												m.at(i) = ( probs_.at(i) * stds.at(i) / scalar_prod(probs_,stds) )
																						* (N_next - N_last - nb_strats);
				}

				/*
					* Compute the m_tilde_i. Same for both methods
					*/
				std::array<unsigned,nb_strats>
				compute_m_tilde(const Array<nb_strats>& m)
				{
								std::array<unsigned,nb_strats> m_tilde;
								double m_sum = 0;
								for(int i=0; i<nb_strats; i++)
								{
												m_tilde.at(i) = std::floor(m_sum + m.at(i)) - std::floor(m_sum);
												m_sum += m.at(i);
								}

								return m_tilde;
				}

				/*
					* Draw realizations in each stratum and update the N_i (total number
					* of drawings in each stratum)
					*/
				template<typename Generator>
				void
				draw_realizations(Generator& gen,
																						const std::array<unsigned,nb_strats>& m_tilde,
																						std::array<unsigned,nb_strats>& N_strata,
																						Array<nb_strats>& sum_strata,
																						Array<nb_strats>& sum_squares_strata)
				{
								double x = 0;
								for(int i=0; i<nb_strats; i++)
								{
												for(int n=0; n<1+m_tilde.at(i); n++)
												{
																x = func_(generate(i,gen));
																sum_strata.at(i) += x;
																sum_squares_strata.at(i) += x*x;
												}

												N_strata.at(i) += 1 + m_tilde.at(i);
								}
				}

				/*
					* Check wether vector x is in stratum i
					*/
				bool instrat(Array<dim> x, int i)
				{
								double product = scalar_prod(x,u_);

								if(i==0)
												return product <= y_.at(0);
								else if(i==nb_strats-1)
												return product > y_.at(nb_strats-2);
								else
												return ((product > y_.at(i-1)) && (product <= y_.at(i)));
    }
    
    
				/*
					* Generate realization of distribution 'dist_' given that it is in
					* strata i
					*/
				template<typename Generator>
				Array<dim> generate(int i, Generator& gen)
				{
								Array<dim> x = dist_(gen());

								while ( !(instrat(x, i)) )
												x = dist_(gen());

        return x;
    }

				void reinit()
				{
								mean_est_ = 0.;
								asymp_std_est_ = 0.;
								sample_size_ = 0;
								time_span_ = std::chrono::duration<double>(0);
				}

				double mean_est() const { return mean_est_; }
				double ci() const { return 1.96*asymp_std_est_/std::sqrt(sample_size_); }
				unsigned sample_size() const { return sample_size_; }
				double time() const { return time_span_.count(); }

				template<typename D,unsigned d,unsigned nb>
				friend std::ostream& operator<<(std::ostream& stream, const Stratification<D,d,nb>& MC);
    
protected:
				Dist dist_;
				Func_Type func_;
				const Array<nb_strats-1> y_;						 // Definit les strats
				const Array<dim> u_;												   // Definit les strats
				const Array<nb_strats> probs_;					// Probas d'etre dans les strats

				double mean_est_;
				double asymp_std_est_;
				unsigned sample_size_;
				std::chrono::duration<double> time_span_;
    
};

/*-------------------------------------------------------------------------------------
	* Stream operator
	*------------------------------------------------------------------------------------*/
template<typename Dist,unsigned dim,unsigned nb_strats>
std::ostream&
operator<<(std::ostream& stream, const Stratification<Dist,dim,nb_strats>& strat)
{
				stream << "Stratification mean: " << strat.mean_est() << std::endl;
				stream << "Stratification asymp_std_est: " << strat.asymp_std_est_ << std::endl;
				stream << "Stratification ci: " << strat.ci() << std::endl;
				stream << "Stratification sample size: " << strat.sample_size() << std::endl;
				stream << "Time span: " << strat.time() << std::endl;
				return stream;
}



#endif // STRATIFICATION_HPP
