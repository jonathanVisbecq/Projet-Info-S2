#ifndef STRATIFICATION_HPP
#define STRATIFICATION_HPP

#include "stl_headers.hpp"
#include "my_array.hpp"
#include "normal_interval.hpp"
#include "monte_carlo.hpp"
#include "rand_var.hpp"



/********************************************************************************************
	* STRUCT Stratification<Dist,dim,nb_strats>
	*
	* Standard stratification
	*
	* Does not provide a confidence interval.
	*
	*******************************************************************************************/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
struct Stratification {

				typedef std::function<double(const typename Dist_Strat<dim>::result_type&)> Func_Type;

				Stratification(const Func_Type& func,
																			const Array<nb_strats-1>& y,
																			const Array<dim>& u,
																			const Array<nb_strats>& probs):
				func_(func), y_(y), u_(u), dist_strats_(init_dist_strats(u_,y_)), probs_(probs){ reinit(); }

				/* Standard non-adaptive stratification with given weights (whose sum
					* must be 1).
					*
					*/
				template<typename Generator>
				double
				operator()(Generator& gen, unsigned N, const Array<nb_strats>& weights);


				/* Standard non-adaptive stratification with proportional allocation.
					*
					*/
				template<typename Generator>
				double
				operator()(Generator& gen, unsigned N)
				{
								return operator()(gen,N,probs_);
				}

				void reinit()
				{
								mean_est_ = 0.;
								sample_size_ = 0.;
								time_span_ = std::chrono::duration<double>(0.);
				}

				double mean_est()						const { return mean_est_; }
				double var_est()						const { return var_est_; }
				unsigned sample_size() const { return sample_size_; }
				double time_span()					const { return time_span_.count(); }

				template<template <unsigned> class D,unsigned d,unsigned nb>
				friend std::ostream& operator<<(std::ostream& stream, const Stratification<D,d,nb>& MC);

				template<template <unsigned> class D,unsigned d,unsigned nb>
				friend class Adaptive_Stratification;

protected:
				Func_Type func_;
				Array<nb_strats-1> y_;																																 // Definit les strats
				Array<dim> u_;																																								 // Definit les strats
				std::array<Dist_Strat<dim>,nb_strats> dist_strats_;    // Conditional distribution on each stratum
				Array<nb_strats> probs_;																														 // Probas d'etre dans les strats

				double mean_est_;
				unsigned sample_size_;
				double var_est_;
				std::chrono::duration<double> time_span_;

				/* Initialize distributions on each stratum
					*/
				static std::array<Dist_Strat<dim>,nb_strats>
				init_dist_strats(Array<dim> u, Array<nb_strats-1> y);
};



/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
std::array<Dist_Strat<dim>, nb_strats>
Stratification<Dist_Strat, dim, nb_strats>::init_dist_strats(
								Array<dim> u, Array<nb_strats-1> y)
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

/*-------------------------------------------------------------------------------------
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
template<typename Generator>
double
Stratification<Dist_Strat,dim,nb_strats>::operator()(
								Generator& gen, unsigned N, const Array<nb_strats>& weights)
{
				reinit();

				std::array<unsigned,nb_strats> Ni;
				double s = 0., temp = 0.;
				Array<nb_strats> sum = makeFill<nb_strats>(0.);
				Array<nb_strats> sum_squares = makeFill<nb_strats>(0.);

				auto time_start = std::chrono::steady_clock::now();
				// Compute allocation to each stratum
				for(int i=0; i<nb_strats; i++)
				{
								Ni.at(i) = std::floor(s + weights.at(i)*N) - std::floor(s);
								s += weights.at(i)*N;
				}

				// Draw realizations and compute estimators
				for(int i=0; i<nb_strats; ++i)
				{
								for(int n=0; n<Ni.at(i); ++n)
								{
												temp = func_(dist_strats_.at(i)(gen()));
												sum.at(i) += temp;
												sum_squares.at(i) += temp*temp;
								}
				}
				time_span_ = std::chrono::steady_clock::now() - time_start;

				for(int i=0; i<nb_strats; ++i)
				{
								mean_est_ += ( probs_.at(i)/weights.at(i) ) * sum.at(i);
								var_est_ += probs_.at(i) * ( sum_squares.at(i) / Ni.at(i) - std::pow(sum.at(i)/Ni.at(i), 2) );
				}

				sample_size_ = N;
				var_est_ /= sample_size_;
				mean_est_ /= sample_size_;

				return mean_est_;
}

/*-------------------------------------------------------------------------------------
	* Stream operator
	*------------------------------------------------------------------------------------*/
template <template <unsigned> class Dist_Strat,unsigned dim,unsigned nb_strats>
std::ostream&
operator<<(std::ostream& stream, const Stratification<Dist_Strat,dim,nb_strats>& strat)
{
				stream << "Standard Stratification mean: " << strat.mean_est() << std::endl;
				stream << "Standars Stratification sample size: " << strat.sample_size() << std::endl;
				stream << "Standard Stratification var: " << strat.var_est() << std::endl;
				stream << "Time span: " << strat.time_span() << std::endl;
				return stream;
}


#endif // STRATIFICATION_HPP
