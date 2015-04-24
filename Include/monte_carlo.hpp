#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include <chrono>
#include <functional>
#include <iostream>
#include <type_traits>
#include <iostream>



/********************************************************************************************
	* STRUCT Monte_Carlo<Generator>
	*
	* Class to represent Monte Carlo estimation.
	*
	* Constraints on type 'Dist':
	*						- typedef 'result_type'
	*						- operator(Generator&) yielding 'result_type'
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
								reinit();
								double x = 0;

								auto time_start = std::chrono::steady_clock::now();
								for(int m=0; m<M; ++m) {
												auto a = gen();
												x = dist_(a);
												sum_ += x;
												sum_of_squares_ += x*x;
												sum_of_cubes_ += x*x*x;
								}

								sample_size_ += M;
								time_span_ = std::chrono::steady_clock::now() - time_start;

								return mean_est();
				}

				void reinit()
				{
								sum_ = 0;
								sum_of_squares_ = 0;
								sum_of_cubes_ = 0;
								sample_size_ = 0;
								time_span_ = std::chrono::duration<double>(0);
				}

				Monte_Carlo& operator+=(const Monte_Carlo& other)
				{
								sum_ += other.sum_;
								sum_of_squares_ += other.sum_of_squares_;
								sum_of_cubes_ += other.sum_of_cubes_;
								sample_size_ += other.sample_size_;
								time_span_ += other.time_span_;
								return *this;
				}

				double mean_est() const { return sum_ / sample_size_; }
				double mean_of_squares_est() const { return sum_of_squares_ / sample_size_; }
				double mean_of_cubes_est() const { return sum_of_cubes_ / sample_size_; }
				double var_est() const { return mean_of_squares_est() - mean_est()*mean_est(); }
				double skew_est() const { return 2*mean_of_cubes_est() - 2*mean_est()*mean_of_squares_est(); }
				double st_dev_est() const { return std::sqrt(var_est()); }
				unsigned sample_size() const { return sample_size_; }
				double ci() const { return 1.96*st_dev_est()/std::sqrt(sample_size_); }
				double time() const { return time_span_.count(); }

				template<typename D>
				friend std::ostream& operator<<(std::ostream& stream, const Monte_Carlo<D>& MC);

protected:
				Dist dist_;

				double sum_;
				double sum_of_squares_;
				double sum_of_cubes_;
				unsigned sample_size_;
				std::chrono::duration<double> time_span_;
};

/*-------------------------------------------------------------------------------------
	* Stream operator
	*------------------------------------------------------------------------------------*/
template<typename Dist>
std::ostream& operator<<(std::ostream& stream, const Monte_Carlo<Dist>& MC)
{
				stream << "MC mean: " << MC.mean_est() << std::endl;
				stream << "MC var: " << MC.var_est() << std::endl;
				stream << "MC skew: " << MC.skew_est() << std::endl;
				stream << "MC std: " << MC.st_dev_est() << std::endl;
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
