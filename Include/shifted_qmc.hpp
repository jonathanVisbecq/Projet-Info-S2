#ifndef SHIFTED_QMC_HPP
#define SHIFTED_QMC_HPP

#include <type_traits>
#include <vector>

#include "my_array.hpp"

/*************************************************************************************
	* STRUCT Shifted_QMC<Dist,QMC_Generator>
	* Distribution
	*
	* Shifted quasi-monte carlo implemented as a distribution to pass to a Monte Carlo
	* estimator.
	*
	*************************************************************************************/
template<typename Dist,typename QMC_Generator>
struct Shifted_QMC{

				static_assert(std::is_same<double,typename Dist::result_type>::value,
																		"Distribution's result_type typedef should be double for shifted QMC");

				typedef double result_type;
				static constexpr unsigned dim_alea = Dist::dim_alea;

				Shifted_QMC() = delete;
				Shifted_QMC(unsigned N,const Dist& dist,QMC_Generator& gen):
								N_(N), dist_(dist), QMCk_(N)
				{
								for(int k=0; k<N_; ++k)
												QMCk_.at(k) = gen();
				}

				double operator()(const Array<dim_alea>& pt)
				{
								x = 0;
								for(k=0; k<N_; ++k)
								{
												// Add both points and take the fractional part coordinate-wise
												for(l=0; l<dim_alea; ++l)
												{
																a.at(l) = pt.at(l) + QMCk_.at(k).at(l);
																a.at(l) -= (a.at(l)>=1) ? 1 : 0;
												}

												x += dist_(a);
								}

								return x / N_;
				}


protected:
				unsigned N_;
				Dist dist_;
				std::vector<Array<dim_alea> > QMCk_;

				Array<dim_alea> a;
				double x;
				int k,l;
};


/*------------------------------------------------------------------------------------
	* FUNC make_shifted_qmc<Dist,Generator>
	*
	* Make a shifted-QMC distribution out of a distribution, a QMC generator and an
	* integer (number of QMC drawings).
	*
	*-----------------------------------------------------------------------------------*/
template<typename Dist,typename QMC_Generator>
Shifted_QMC<Dist,QMC_Generator>
make_shifted_qmc(unsigned N,
																	const Dist& dist,
																	QMC_Generator& gen)
{
				return Shifted_QMC<Dist,QMC_Generator>(N,dist,gen);
}

#endif // SHIFTED_QMC_HPP
