#ifndef LINEAR_ESTIMATOR_HPP
#define LINEAR_ESTIMATOR_HPP

#include <chrono>
#include <cmath>


/*************************************************************************************
	*
	* STRUCT Linear_Estimator
	*
	* Base class for estimators.
	*
	* Rk: should rewrite something for QMC estimator because confidence intervals and
	* variances do not make sense
	*
	*************************************************************************************/
struct Linear_Estimator{

				Linear_Estimator() { reinit(); }

				double mean_est() const
				{
								return this->sum_ / this->sample_size_;
				}

				double mean_of_squares_est() const
				{
								return this->sum_of_squares_ / this->sample_size_;
				}

				double var_est() const
				{
								return mean_of_squares_est() - std::pow(mean_est(),2);
				}

				double st_dev_est() const
				{
								return std::sqrt(var_est());
				}

				double ci() const
				{
								return 1.96*st_dev_est()/std::sqrt(sample_size_);
				}

				void reinit()
				{
								sum_ = 0;
								sum_of_squares_ = 0;
								sample_size_ = 0;
								time_span_ = std::chrono::duration<double>(0);
				}

				double time() const
				{
								return time_span_.count();
				}

				Linear_Estimator& operator+=(const Linear_Estimator& other)
				{
								sum_ += other.sum_;
								sum_of_squares_ += other.sum_of_squares_;
								sample_size_ += other.sample_size_;
								time_span_ += other.time_span_;

								return *this;
				}


				friend std::ostream& operator<<(std::ostream& stream, const Linear_Estimator& MC);

protected:
				double sum_;
				double sum_of_squares_;
				unsigned sample_size_;
				std::chrono::duration<double> time_span_;
};


std::ostream& operator<<(std::ostream& stream, const Linear_Estimator& MC)
{
				stream << "MC mean: " << MC.mean_est() << std::endl;
				stream << "MC var: " << MC.var_est() << std::endl;
				stream << "MC std: " << MC.st_dev_est() << std::endl;
				stream << "MC ci: " << MC.ci() << std::endl;
				stream << "Time span: " << MC.time() << std::endl;

				return stream;
}


#endif // LINEAR_ESTIMATOR_HPP
