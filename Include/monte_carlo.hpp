#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include "stl_headers.hpp"
#include "normal_interval.hpp"



/********************************************************************************************
	* STRUCT Monte_Carlo<Generator>
	*
	* Class to represent Monte Carlo estimation.
	*
	* 'Dist' must be a distribution type
	*
	*******************************************************************************************/
template <typename Dist>
struct Monte_Carlo{

				static_assert(std::is_same<double,typename Dist::result_type>::value,
																		"Distribution's result_type typedef should be double for Monte Carlo estimation");

				Monte_Carlo(const Dist& dist): dist_(dist) { reinit(); }


				template<typename G>
				double operator()(G& gen, unsigned M)
				{
								// Needed to estimate moment of order three
								std::vector<double> x(M);

								auto time_start = std::chrono::steady_clock::now();
								for(int m=0; m<M; ++m) {
												auto a = gen();
												x.at(m) = dist_(a);
												sum_ += x.at(m);
												sum_of_squares_ += x.at(m)*x.at(m);
								}

								time_span_ = std::chrono::steady_clock::now() - time_start;

								sample_size_ = M;
								for(int m=0; m<M; ++m)
												sum_moment_3_ += std::pow( std::abs(x.at(m) - mean_est()), 3);

								return mean_est();
				}


				void reinit()
				{
								sum_ = 0;
								sum_of_squares_ = 0;
								sum_moment_3_ = 0;
								sample_size_ = 0;
								time_span_ = std::chrono::duration<double>(0);
				}

				Monte_Carlo& operator+=(const Monte_Carlo& other)
				{
								sum_ += other.sum_;
								sum_of_squares_ += other.sum_of_squares_;
								sum_moment_3_ += other.sum_moment_3_;
								sample_size_ += other.sample_size_;
								time_span_ += other.time_span_;
								return *this;
				}

				double mean_est() const { return sum_ / sample_size_; }
				double mean_of_squares_est() const { return sum_of_squares_ / sample_size_; }
				double var_est() const { return mean_of_squares_est() - mean_est()*mean_est(); }
				double moment_3_est() const { return sum_moment_3_ / sample_size_; }
				double st_dev_est() const { return std::sqrt(var_est()); }
				unsigned sample_size() const { return sample_size_; }
				double time() const { return time_span_.count(); }
				double ci(double alpha = 0.05) const
				{
								return inv_normal_cdf(1.-alpha/2.) * st_dev_est() / std::sqrt(sample_size_);
				}
				double ci_berry_essen(double alpha = 0.05) const
				{
								// Upper bound for Berry-Essen's constant from Shevtsova (2011)
								return ci(alpha) + 0.4748 * moment_3_est() / (var_est() * sample_size_);
				}

				template<typename D>
				friend std::ostream& operator<<(std::ostream& stream, const Monte_Carlo<D>& MC);

protected:
				Dist dist_;

				double sum_;
				double sum_of_squares_;
				double sum_moment_3_;
				unsigned sample_size_;
				std::chrono::duration<double> time_span_;
};

/*-------------------------------------------------------------------------------------
	* Stream operator
	*------------------------------------------------------------------------------------*/
template<typename Dist>
std::ostream&
operator<<(std::ostream& stream, const Monte_Carlo<Dist>& MC)
{
				stream << "MC mean: " << MC.mean_est() << std::endl;
				stream << "MC var: " << MC.var_est() << std::endl;
				stream << "MC sample size: " << MC.sample_size() << std::endl;
				stream << "MC ci: " << MC.ci() << std::endl;
				stream << "Time span: " << MC.time() << std::endl;
				return stream;
}

/*-------------------------------------------------------------------------------------
	* FUNC make_mc<Dist>
	*
	* Make a Monte_Carlo<Dist> object from a distribution and a function
	*
	*------------------------------------------------------------------------------------*/
template<typename Dist>
Monte_Carlo<Dist>
make_mc(const Dist& dist)
{
				return Monte_Carlo<Dist>(dist);
}





















#endif // MONTE_CARLO_HPP
