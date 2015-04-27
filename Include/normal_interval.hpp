#ifndef NORMAL_INTERVAL_HPP
#define NORMAL_INTERVAL_HPP

#include "stl_headers.hpp"
#include "my_array.hpp"
#include "gaussian_ind.hpp"
#include "uniform_generator.hpp"

/* Standard normal CDF
	*/
inline double normal_cdf(double x)
{
				return 0.5 + 0.5 * std::erf(x / std::sqrt(2));
}

/* Standard normal inverse CDF
	*
	* Beasley-Springer-Moro algorithm which can be found in [Glasserman 2004].
	*/
inline double inv_normal_cdf(double quantile)
{
				static double a[4] = {2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};

				static double b[4] = {-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};

				static double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209, 0.0276438810333863, 0.0038405729373609,
																										0.0003951896511919, 0.0000321767881768, 0.0000002888167364, 0.0000003960315187};

				double y = quantile - 0.5;
				double r,s,t,result;

				if (fabs(y) < 0.42)
				{
								r = y * y;
								result = y * (a[0] + r * (a[1] + r * (a[2] + r * a[3]))) / (1.0 + r * (b[0] + r * (b[1] + r * (b[2] + r * b[3]))));
				}
				else
				{
								if(y <= 0)
												r = quantile;
								else
												r = 1 - quantile;

								s = log(-log(r));
								t = c[0] + s * (c[1] + s * (c[2] + s * (c[3] + s * (c[4] + s * (c[5] + s * (c[6] + s * (c[7] + s * c[8])))))));

								if(quantile > 0.5)
												result =  t;
								else
												result = -t;
					}

				return result;
}

/*************************************************************************************
	* STRUCT Normal_Interval<dim>
	* Distribution
	*
	* Distribution of a standard multivariate normal random variable x given that u.x
	* lies within an interval. Useful for stratification methods.
	*
	*************************************************************************************/
template<unsigned dim>
struct Normal_Interval{

				typedef Array<dim> result_type;
				static constexpr unsigned dim_alea = Gaussian_Ind<dim>::dim_alea + 1;

				Normal_Interval(Array<dim> u = makeFill<dim>(1.),
																				double min = -std::numeric_limits<double>::infinity(),
																				double max = std::numeric_limits<double>::infinity()):
								u_(u), min_(min), max_(max),
								phi_max_(normal_cdf(max_)), phi_min_(normal_cdf(min_)),
								G_() {}

				result_type operator()(const Array<dim_alea>& pt)
				{
								z = inv_normal_cdf( (phi_max_ - phi_min_)*pt.at(0) + phi_min_ );

								auto it_begin = pt.begin();
								std::copy(++it_begin, pt.end(), part_pt.begin());

								Y = G_(part_pt);

								return z * u_ + Y - scalar_prod(u_,Y) * u_;
				}

protected:
				Array<dim> u_;
				double min_, max_;
				double phi_max_, phi_min_;
				Gaussian_Ind<dim> G_;

				double z;
				Array<dim> Y;
				Array<dim_alea - 1> part_pt;
};


#endif // NORMAL_INTERVAL_HPP
