#ifndef ADAPTIVE_STRATIFICATION_HPP
#define ADAPTIVE_STRATIFICATION_HPP



#include <iostream>

#include "my_array.hpp"
#include "normal_interval.hpp"
#include "monte_carlo.hpp"
#include "rand_var.hpp"
#include "stratification.hpp"



/********************************************************************************************
	* STRUCT Adaptive_Stratification<Dist,dim,nb_strats>
	*
	* Adaptive optimal stratification from EJ07
	*
	*******************************************************************************************/
enum class Method {A, B};

template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
struct Adaptive_Stratification {

				typedef std::function<double(const typename Dist_Strat<dim>::result_type&)> Func_Type;
				typedef std::function<void(Array<nb_strats>& m_,
																															const Array<nb_strats>& stds_,
																															unsigned N_last_, unsigned N_next)> Compute_M_Type;

				Adaptive_Stratification(const Func_Type& func,
																												const Array<nb_strats-1>& y,
																												const Array<dim>& u,
																												const Array<nb_strats>& probs,
																												Method method):
				func_(func), y_(y), u_(u),
				dist_strats_(Stratification<Dist_Strat,dim,nb_strats>::init_dist_strats(u_,y_)),
				probs_(probs), method_(method)
				{ reinit(); }


				/* Adaptive optimal stratification with given drawings budget
					*
					* Using method a)
					*/
				template<typename Generator,size_t nb_iter>
				double
				operator()(Generator& gen, std::array<unsigned,nb_iter> N);

				/* One step of adaptive optimal stratification, adding N to the total
					* number of drawings.
					*
					* May be used to refine previous run.
					*
					* Using method a)
					*/
				template<typename Generator>
				double
				operator()(Generator& gen, unsigned N);

				void reinit()
				{
								sum_strata_ = makeFill<nb_strats>(0.);
								sum_squares_strata_ = makeFill<nb_strats>(0.);
								stds_ = makeFill<nb_strats>(0.);
								m_ = makeFill<nb_strats>(0.);

								for(int i=0; i<nb_strats; ++i)
												N_strata_.at(i) = 0;

								N_last_ = 0;

								mean_est_ = 0.;
								asymp_std_est_ = 0.;
								sample_size_ = 0;
								var_est_ = 0.;
								time_span_ = std::chrono::duration<double>(0.);
				}


				double mean_est()						const { return mean_est_; }
				double asymp_std_est() const { return asymp_std_est_; }
				unsigned sample_size() const { return sample_size_; }
				double ci(double alpha = 0.05)									   const
				{
								return inv_normal_cdf(1-alpha/2.)*asymp_std_est_/std::sqrt(sample_size_);
				}

				double var_est()									   const { return var_est_; }
				double time_span()					const { return time_span_.count(); }

				const std::array<unsigned,nb_strats>& N_strata() { return N_strata_; }
				const Array<nb_strats>& stds() { return stds_; }
				const Array<nb_strats>& sum_strata() { return sum_strata_; }

				template<template <unsigned> class D,unsigned d,unsigned nb>
				friend std::ostream& operator<<(std::ostream& stream, const Adaptive_Stratification<D,d,nb>& MC);

protected:
				Func_Type func_;
				const Array<nb_strats-1> y_;																				// Definit les strats
				const Array<dim> u_;																												// Definit les strats
				std::array<Dist_Strat<dim>,nb_strats> dist_strats_;  // Distribution on each stratum
				const Array<nb_strats> probs_;																		// Probas d'etre dans les strats
				Method method_;																											// Method for computing the mi

				std::function<void(unsigned,unsigned)> compute_m_;

				Array<nb_strats> sum_strata_;
				Array<nb_strats> sum_squares_strata_;
				Array<nb_strats> stds_;
				Array<nb_strats> m_;
				std::array<unsigned,nb_strats> N_strata_;
				unsigned N_last_;

				double mean_est_;
				double asymp_std_est_;
				unsigned sample_size_;
				double var_est_;
				std::chrono::duration<double> time_span_;


				/* Iniitialize distributions on each stratum
					*/
				static std::array<Dist_Strat<dim>,nb_strats>
				init_dist_strats(Array<dim> u, Array<nb_strats-1> y)
				{
								std::array<Dist_Strat<dim>,nb_strats> dist_strats;
								for(int i=0; i<nb_strats; ++i)
								{
												if(i==0)
																dist_strats.at(i) = Dist_Strat<dim>(u,
																																																			-std::numeric_limits<double>::infinity(),
																																																			y.at(0));
												else if(i==nb_strats-1)
																dist_strats.at(i) = Dist_Strat<dim>(u,
																																																			y.at(nb_strats-2),
																																																			std::numeric_limits<double>::infinity());
												else
																dist_strats.at(i) = Dist_Strat<dim>(u,
																																																				y.at(i-1),
																																																				y.at(i));
								}

								return dist_strats;
				}

				/* Compute the standard deviation on each stratum
					*/
				void
				compute_empirical_stds();

				/*
					* Compute the m_i using method a), given that all stds are not equal to 0
					*/
				void
				compute_m_A(unsigned N_last_, unsigned N_next);

				/*
					* Compute the m_i using method b), given that all stds are not equal to 0
					*/
				void
				compute_m_B(unsigned N_last, unsigned N_next);

				/*
					* Compute the m_tilde_i. Same for both methods
					*/
				std::array<unsigned,nb_strats>
				compute_m_tilde(unsigned N_last_, unsigned N_next);

				/*
					* Draw realizations in each stratum and update the N_i (total number
					* of drawings in each stratum)
					*/
				template<typename Generator>
				void
				draw_realizations(Generator& gen,
																						const std::array<unsigned,nb_strats>& m_tilde);

				// One step of the iterative algorithm
				template<typename Generator>
				void
				step(Generator& gen, unsigned N_next);

				// Update estimation information
				void update_info(std::chrono::duration<double> time_span);

};



/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
void
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::update_info(std::chrono::duration<double> time_span)
{
				mean_est_ = 0.;
				asymp_std_est_ = 0.;
				sample_size_ = 0;

				for(int i=0; i<nb_strats; ++i)
				{
								mean_est_ += probs_.at(i) * (1./N_strata_.at(i)) * sum_strata_.at(i);
								asymp_std_est_ += probs_.at(i) * stds_.at(i);
								sample_size_ += N_strata_.at(i);
				}

				var_est_ = asymp_std_est_ * asymp_std_est_ / sample_size_;

				time_span_ += time_span;
}

/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
template<typename Generator>
void
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::step(Generator& gen, unsigned N_next)
{
				// Check wether all empirical stds are equal to 0
				if(std::accumulate(stds_.begin(),stds_.end(),0.)>0.)
				{
								if(method_==Method::A)
												compute_m_A(N_last_,N_next);
								else if(method_==Method::B)
												compute_m_B(N_last_,N_next);
				}
				else
				{
								for(int i=0; i<nb_strats; ++i)
												m_.at(i) = probs_.at(i) * (N_next - N_last_ - nb_strats);
				}

				draw_realizations(gen, compute_m_tilde(N_last_, N_next));
				compute_empirical_stds();

				N_last_ = N_next;
}


/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
template<typename Generator>
double
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::operator()(
								Generator& gen, unsigned N)
{
				auto time_start = std::chrono::steady_clock::now();
				step(gen, N_last_ + N);
				update_info(std::chrono::steady_clock::now() - time_start);

				return mean_est_;
}

/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
template<typename Generator,size_t nb_iter>
double
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::operator()(
								Generator& gen, std::array<unsigned,nb_iter> N)
{
				reinit();

				auto time_start = std::chrono::steady_clock::now();
				for(int s=0; s<nb_iter; ++s)
								step(gen, N.at(s));
				update_info(std::chrono::steady_clock::now() - time_start);

				return mean_est_;
}

/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
void
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::compute_empirical_stds()
{
				for(int i=0; i<nb_strats; ++i)
								stds_.at(i) = std::sqrt( sum_squares_strata_.at(i)/N_strata_.at(i) - std::pow(sum_strata_.at(i)/N_strata_.at(i), 2) );
}

/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
void
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::compute_m_A(unsigned N_last, unsigned N_next)
{
				for(int i=0; i<nb_strats ;i++)
								m_.at(i) = ( probs_.at(i) * stds_.at(i) / scalar_prod(probs_,stds_) )
																		* (N_next - N_last - nb_strats);
}


/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
void
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::compute_m_B(unsigned N_last, unsigned N_next)
{
				// Compute first array (point i) )
				Array<nb_strats> array_1 = makeFill<nb_strats>(0.);
				for(int i=0; i<nb_strats; ++i)
								array_1.at(i) = ( N_strata_.at(i) + 1. ) /
																								( probs_.at(i) * stds_.at(i) );

				// Get index that sort 'array_unordered' in descending order
				std::array<unsigned,nb_strats> idx = sort_idx_desc(array_1);

				// Compute second array ( point ii) )
				Array<nb_strats> array_2 = makeFill<nb_strats>(0.);

//				double sum_nb_drawings = 0., sum_probs_stds = 0.;
//				for(int i=nb_strats-1; i>=0; --i)
//				{
//								array_2.at(i) = ((N_next - N_last - nb_strats) + sum_nb_drawings) / sum_probs_stds;
//								sum_nb_drawings += ( N_strata_.at(idx.at(i)) + 1. );
//								sum_probs_stds += probs_.at(idx.at(i)) * stds_.at(idx.at(i));
//				}

				for(int i=0; i<nb_strats; ++i)
				{
								double aux1 = 0., aux2 = 0.;
								for(int k=i+1; k<nb_strats; ++k)
								{
												aux1 += N_strata_.at(idx.at(k)) + 1.;
												aux2 += probs_.at(idx.at(k)) * stds_.at(idx.at(k));
								}

								array_2.at(i) = ( (N_next - N_last - nb_strats) + aux1 ) / aux2;
				}

				// Find i_star
				int i_star = 0, j = 0;
				while( (j < nb_strats) && (array_1.at(idx.at(j)) >= array_2.at(j)) ) ++j;
				i_star = (j==0)? 0 : (j-1);

				// Compute the mi
				for(int i=0; i<nb_strats; ++i)
				{
								if(i <= i_star)
												m_.at(idx.at(i)) = 0.;
								else
												m_.at(idx.at(i)) = probs_.at(idx.at(i)) * stds_.at(idx.at(i)) * array_2.at(i_star) - N_strata_.at(idx.at(i)) - 1.;

				}
}

/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
std::array<unsigned,nb_strats>
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::compute_m_tilde(unsigned N_last, unsigned N_next)
{
				std::array<unsigned,nb_strats> m_tilde;
				double m_sum = 0.;
				for(int i=0; i<nb_strats; i++)
				{
								m_tilde.at(i) = std::floor(m_sum + m_.at(i)) - std::floor(m_sum);
								m_sum += m_.at(i);
				}

				// Accounts for rounding error with std::floor
				if(std::accumulate(m_tilde.begin(),m_tilde.end(),0) != (N_next - N_last - nb_strats))
								m_tilde.at(nb_strats-1)++;

				return m_tilde;
}

/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
template<typename Generator>
void
Adaptive_Stratification<Dist_Strat,dim,nb_strats>::draw_realizations(
								Generator& gen,
								const std::array<unsigned,nb_strats>& m_tilde)
{
				double x = 0;
				for(int i=0; i<nb_strats; i++)
				{
								for(int n=0; n<1+m_tilde.at(i); n++)
								{
												x = func_(dist_strats_.at(i)(gen()));
												sum_strata_.at(i) += x;
												sum_squares_strata_.at(i) += x*x;
								}

								N_strata_.at(i) += 1 + m_tilde.at(i);
				}
}



/*-------------------------------------------------------------------------------------
	* Stream operator
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
std::ostream&
operator<<(std::ostream& stream, const Adaptive_Stratification<Dist_Strat,dim,nb_strats>& strat)
{
				stream << "Adaptive Stratification mean: " << strat.mean_est() << std::endl;
				stream << "Adaptive Stratification asymp_std_est: " << strat.asymp_std_est() << std::endl;
				stream << "Adaptive Stratification sample size: " << strat.sample_size() << std::endl;
				stream << "Adaptive Stratification ci: " << strat.ci() << std::endl;
				stream << "Adaptive Stratification var: " << strat.var_est() << std::endl;
				stream << "Time span: " << strat.time_span() << std::endl;
				return stream;
}



#endif // ADAPTIVE_STRATIFICATION_HPP
