#ifndef RANDOM_START_HALTON_HPP
#define RANDOM_START_HALTON_HPP

#include <type_traits>

#include "my_array.hpp"
#include "kakutani.hpp"

/*************************************************************************************
	* STRUCT RandStart_Halton<Dist>
	* Distribution
	*
	* Random-start Halton sequence implemented as a distribution to pass to a Monte Carlo
	* estimator.
	*
	*************************************************************************************/
template<typename Dist>
struct RandStart_Halton{

				static_assert(std::is_same<double,typename Dist::result_type>::value,
																		"Distribution's result_type typedef should be double for random-start Halton");

				typedef double result_type;
				static constexpr unsigned dim_alea = Dist::dim_alea;

				RandStart_Halton() = delete;
				RandStart_Halton(unsigned N,const Dist& dist):
								N_(N),dist_(dist),haltF_() {}

				result_type operator()(const Array<dim_alea>& pt)
				{
								haltF_ = Halton_Fast<dim_alea>(pt);

								x = 0;
								for(i=0; i<N_; ++i)
												x += dist_(haltF_());

								return x / N_;
				}



protected:
				unsigned N_;
				Dist dist_;
				Halton_Fast<dim_alea> haltF_;

				double x;
				int i;
};

/*------------------------------------------------------------------------------------
	* FUNC make_randStart_halton<Dist>
	*
	* Make a random-start Halton distribution out of a distribution and ann integer
	* (number of QMC drawings)
	*
	*-----------------------------------------------------------------------------------*/
template<typename Dist>
RandStart_Halton<Dist>
make_randStart_halton(unsigned N,
																						const Dist& dist)
{
				return RandStart_Halton<Dist>(N,dist);
}


#endif // RANDOM_START_HALTON_HPP
