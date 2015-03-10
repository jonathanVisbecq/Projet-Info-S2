#ifndef FAURE_HPP
#define FAURE_HPP

#include <vector>
#include <iostream>

#include "my_array.hpp"
#include "p_adic.hpp"


/*************************************************************************************
	*
	* STRUCT Faure
	*
	* Generator for Faure's sequences
	*
	*************************************************************************************/
template<unsigned dim>
struct Faure {

				typedef Array<dim> result_type;
				typedef std::vector<std::vector<int>> Coeff_Binom;


				/*
					* 'x.p_' should be the smallest prime number greater than 'dimension'
					*/
				Faure(const P_Adic& x):
								x_(x), result_(), comb_(make_binom_array(32)) {}
				/*
					* Use the first prime number greater than 'dimension' (the dimension should be <= than 1613)
					* and start from the representation of 1
					*/
				Faure(): Faure(P_Adic(1,smallest_greater_prime(dim))) {}


				result_type current() { return result_; }
				result_type operator()();
				P_Adic transform(const P_Adic& y);

				// Output current vector
				template<unsigned d>
				friend std::ostream& operator<<(std::ostream &stream, const Faure<d>& k);

public:
				P_Adic x_;
				result_type result_;
				Coeff_Binom comb_;

				static int smallest_greater_prime(int d);
				static Coeff_Binom make_binom_array(int p);
};


/*-------------------------------------------------------------------------------------
	* Methods
	*------------------------------------------------------------------------------------*/
template<unsigned dim>
auto
Faure<dim>::operator()() -> result_type{

				P_Adic xi = x_++;
				auto it = result_.begin();

				while(it != result_.end()){
								*it++ = (double) xi;
								xi = transform(xi);
				}

				return result_;
}


template<unsigned dim>
P_Adic
Faure<dim>::transform(const P_Adic& y){

				P_Adic::Coeff bk;
				auto it_Comb = comb_.begin();
				auto it = y.ak_.begin();

				while(it != y.ak_.end()){
								bk.push_back(std::inner_product(it, y.ak_.end(), (*it_Comb).begin(), 0) % y.base());
								it++;
								it_Comb++;
				}

				return P_Adic(bk, y.pk_);
}


template<unsigned dim>
int
Faure<dim>::smallest_greater_prime(int d){

				int idx = 0;
				while(primes[idx]<d && idx<255)
								++idx;

				return primes[idx];
}

template<unsigned dim>
auto
Faure<dim>::make_binom_array(int p) -> Coeff_Binom{

				Coeff_Binom comb(p,std::vector<int>(p,0));

				for(int n = 0; n < p; n++)
								comb[0][n] = 1;

				for(int k = 1; k < p; k++) {
								comb[k][0] = 1;

								for (int n = 1; n <= p-k; n++)
												comb[k][n] = comb[k][n-1] + comb[k-1][n];
				}

				return comb;
}



/*-------------------------------------------------------------------------------------
	* Related non-member functions
	*------------------------------------------------------------------------------------*/
template<unsigned d>
std::ostream& operator<<(std::ostream &stream, const Faure<d>& f){

				for(auto it=f.result_.begin(); it!=f.result_.end(); ++it)
								stream << *it << "    ";

				return stream;
}












#endif // FAURE_HPP
