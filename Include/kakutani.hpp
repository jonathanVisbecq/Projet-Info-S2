#ifndef KAKUTANI_HPP
#define KAKUTANI_HPP

#include <iostream>
#include <functional>
#include <algorithm>
#include <cmath>

#include "my_array.hpp"
#include "p_adic.hpp"

/*************************************************************************************
	*
	* STRUCT Kakutani<dim>
	*
	* Generator for Kakutani's sequences
	*
	*************************************************************************************/
template<unsigned dim>
struct Kakutani{

				typedef Array<dim> result_type;

				Kakutani() = delete;
				// Each coordinate of x must have same base 'p_' as the corresponding coordinate of y
				Kakutani(const Vect_P_Adic& x, const Vect_P_Adic& y) :
								xk_(x), yk_(y), result_()
				{
								operator()();
				}

				result_type operator()()
				{
								std::transform(xk_.begin(), xk_.end(), yk_.begin(), result_.begin(), op);
								return result_;
				}

				result_type current() { return result_; }

				// Output current vector
				template<unsigned d>
				friend std::ostream& operator<<(std::ostream &stream, const Kakutani<d>& k);

protected:
				Vect_P_Adic xk_;
				Vect_P_Adic yk_;
				result_type result_;

				static double op(P_Adic& x, const P_Adic& y)
				{
												x = x + y;
												return (double) x;
				}


};


/*-------------------------------------------------------------------------------------
	* Related non-member functions
	*------------------------------------------------------------------------------------*/
template<unsigned d>
std::ostream& operator<<(std::ostream &stream, const Kakutani<d>& k){

				for(auto it=k.result_.begin(); it!=k.result_.end(); ++it)
								stream << *it << "    ";

				return stream;
}



/********************************************************************************************
	*
	* STRUCT Halton<dim>
	*
	* Halton sequences are special cases of Kakutani sequences, but we can make the generation
	* way faster: addition of 'P_Adic' objects is slow.
	*
	* Moreover, with implementation of random-start halton sequences in view, obtaining the
	* p-adic decomposition of an arbitrary number in [0,1) may take a while. That's why we use
	* another implementation for randomized QMC.
	*
	*******************************************************************************************/
template<unsigned dim>
struct Halton: public Kakutani<dim>{

				typedef Array<dim> result_type;

				// Use the first 'dim' prime numbers as increment and an orbit vector
				Halton(const Array<dim> &orbit_pt = makeFill<dim>(0.)):
								Kakutani<dim>(make_init_pt(orbit_pt), make_init_angle()) {}


				result_type operator()()
				{
								for(int i=0; i<dim; ++i)
											this->result_.at(i) = (double) ++(this->xk_).at(i);

								std::cout << this->result_ << std::endl << std::endl;
								return this->result_;
				}


protected:

				// Build the p-adic representation of the starting point (multi-dimensionnal)
				static Vect_P_Adic make_init_pt(const Array<dim> &starting_pt);
				// Build the p-adic representation of the angle (multi-dimensionnal)
				static Vect_P_Adic make_init_angle();

};

/*-------------------------------------------------------------------------------------
	* Methods
	*------------------------------------------------------------------------------------*/

template<unsigned dim>
Vect_P_Adic
Halton<dim>::make_init_pt(const Array<dim> &starting_pt)
{
				Vect_P_Adic v;
				for(int i=0; i<dim; ++i){
								v.push_back(P_Adic(starting_pt.at(i),primes[i]));
								std::cout << v.back() << std::endl << std::endl;
				}

				return v;
}

template<unsigned dim>
Vect_P_Adic
Halton<dim>::make_init_angle()
{
				Vect_P_Adic v;
				for(int i=0; i<dim; ++i)
								v.push_back(P_Adic(1,primes[i]));

				return v;
}

/********************************************************************************************
	*
	* STRUCT Halton_Fast<dim>
	*
	* Faster than the previous version. Better for randomized QMC.
	*
	*******************************************************************************************/
template<unsigned dim>
struct Halton_Fast{

				typedef Array<dim> result_type;

				Halton_Fast(const Array<dim>& starting_pt = makeFill<dim>(0.25)):
								x_(starting_pt)
				{
								for(int i=0; i<dim; ++i)
												bases_.at(i) = primes[i];
				}

				result_type operator()()
				{
								for(i=0; i<dim; ++i){
												base = bases_.at(i);
												k = std::floor(- std::log(1. - x_.at(i)) / std::log(base) );
												x_.at(i) += 1./std::pow((double)base,k) + 1./std::pow((double)base,k+1) - 1.;
								}

								return x_;
				}


protected:
				Array<dim> x_;
				std::array<int,dim> bases_;

				int i,k,base;
};

#endif // KAKUTANI_HPP
