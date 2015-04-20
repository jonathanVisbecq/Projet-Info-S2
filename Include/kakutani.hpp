#ifndef KAKUTANI_HPP
#define KAKUTANI_HPP

#include <iostream>
#include <functional>
#include <algorithm>

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
				Kakutani(const Vect_P_Adic& x, const Vect_P_Adic& y) : xk_(x), yk_(y), result_()
				{
								operator()();
				}

				result_type operator()() {
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

				static double op(P_Adic& x, const P_Adic& y){
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
	* Halton sequences are special cases of Kakutani sequences
	*
	*******************************************************************************************/
template<unsigned dim>
struct Halton: public Kakutani<dim>{


				// Use the first 'dim' prime numbers
				Halton(): Kakutani<dim>(make_init_vect(),make_init_vect()) {}

				Halton(const std::array<unsigned,dim>& primes): Kakutani<dim>(make_init_vect(primes), make_init_vect(primes)) {}


protected:

				// Functions building 'xk_' 'yk_' array to initilize a 'Kakutatni' instance
				static Vect_P_Adic make_init_vect();
				static Vect_P_Adic make_init_vect(const std::array<unsigned, dim> &primes) ;

};

/*-------------------------------------------------------------------------------------
	* Methods
	*------------------------------------------------------------------------------------*/
template<unsigned dim>
Vect_P_Adic
Halton<dim>::make_init_vect(){
				Vect_P_Adic v;

				for(int i=0; i<dim; ++i)
								v.push_back(P_Adic(1,primes[i]));

				return v;
}

template<unsigned dim>
Vect_P_Adic
Halton<dim>::make_init_vect(const std::array<unsigned,dim>& primes){
				Vect_P_Adic v;

				for(auto it=primes.begin(); it!=primes.end(); ++it){
								v.push_back(P_Adic(1,*it));
				}

				return v;
}





#endif // KAKUTANI_HPP
