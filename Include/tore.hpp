#ifndef TORE_HPP
#define TORE_HPP

#include <cmath>

#include "p_adic.hpp"
#include "my_array.hpp"


/*************************************************************************************
	*
	* STRUCT Tore<dim>
	*
	* Generator for Tore sequences
	*
	*************************************************************************************/
template<unsigned dim>
struct Tore{

				typedef Array<dim> result_type;

				Tore(const Array<dim>& alphak): alphak_(alphak), result_(), n_(0)
				{
								operator()();
				}

				result_type operator()()
				{
								n_++;
								for(i=0; i<dim; ++i)
												result_.at(i) = std::modf(n_*alphak_.at(i),&intpart);

								return result_;
				}

				result_type current() { return result_; }


				// Output current vector
				template<unsigned d>
				friend std::ostream& operator<<(std::ostream &stream, const Tore<d>& k);

protected:
				result_type alphak_;
				result_type result_;
				unsigned n_;

				int i;
				double intpart;
};

/*-------------------------------------------------------------------------------------
	* Related non-member functions
	*------------------------------------------------------------------------------------*/
template<unsigned d>
std::ostream& operator<<(std::ostream &stream, const Tore<d>& k){

				for(auto it=k.result_.begin(); it!=k.result_.end(); ++it)
								stream << *it << "    ";

				return stream;
}



/*************************************************************************************
	*
	* STRUCT SQRT<dim>
	*
	* Special case of Tore sequence, initialized with the first prime numbers
	*
	*************************************************************************************/
template<unsigned dim>
struct SQRT: public Tore<dim>{

				SQRT(): Tore<dim>(init_alphaks()) {}


protected:
								static Array<dim> init_alphaks();
};


/*-------------------------------------------------------------------------------------
	* Methods
	*------------------------------------------------------------------------------------*/
template<unsigned dim>
Array<dim>
SQRT<dim>::init_alphaks()
{
				Array<dim> a;
				for(int i=0; i<dim; ++i)
								a.at(i) = std::sqrt(primes[i]);

				return a;
}















#endif // TORE_HPP
