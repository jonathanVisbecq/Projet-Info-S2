#include "p_adic.hpp"


/*************************************************************************************
	*
	* CLASS P_Adic
	*
	*************************************************************************************/

/*-------------------------------------------------------------------------------------
	* Constructors
	*------------------------------------------------------------------------------------*/

P_Adic::P_Adic(double x, int p):
				p_(p)
{

				int puiss = 1;
				double integral = 0;

				while(x>0)
				{
								x = std::modf(x*p_,&integral);
								ak_.push_back((int)integral);
								pk_.push_back(puiss);
								puiss *= p_;
				}

				// One more element because we stored in 'pk_' the power of 'p_' used in
				// the integer representation
				pk_.push_back(puiss);
}


P_Adic::P_Adic(int n, int p):
				p_(p)
{

				int puiss = 1;
				while(n > 0)
				{
								ak_.push_back(n % p_);
								pk_.push_back(puiss);
								puiss *= p_;
								n -= ak_.back();
								n /= p_;
				}

				// One more element to allow later conversion to 'double' from the p-adic representation
				pk_.push_back(puiss);
}

/*-------------------------------------------------------------------------------------
	* Methods
	*------------------------------------------------------------------------------------*/

P_Adic P_Adic::operator++(int)
{
				P_Adic copie = * this;
				increment();
				return copie;
}

P_Adic& P_Adic::operator++()
{
				increment();
				return (* this);
}



// Moves to the p-adic representation of n to the representation
// of n+1
void P_Adic::increment()
{
				Coeff::iterator it = ak_.begin();
				while((it != ak_.end()) && ((*it)+1 == p_))
				{
								(*it) = 0;
								it++;
				}

				if (it == ak_.end())
				{
								ak_.push_back(1);
								pk_.push_back(pk_.back()*p_);
				}
				else
								(*it) += 1;
}



/*-------------------------------------------------------------------------------------
	* Related non-member functions
	*------------------------------------------------------------------------------------*/


P_Adic operator+(const P_Adic& m1, const P_Adic& m2){

				// 'm1' and 'm2' must be decompositions with respect to the same base 'p_'
				int p = m1.p_;
				P_Adic::Coeff ak, pk;

				auto it1 = m1.ak_.begin();
				auto it2 = m2.ak_.begin();

				int r = 0, s = 0;
				int puiss = 1;
				while(it1!=m1.ak_.end() && it2!=m2.ak_.end())
				{
								s = *it1 + *it2 + r;
								ak.push_back(s % p);
								pk.push_back(puiss);
								r = (s >= p) ? 1 : 0;
								puiss *= p;
								++it1;
								++it2;
				}

				// Only one of the following two loops is executed
				while(it1!=m1.ak_.end())
				{
								s = *it1 + r;
								ak.push_back(s % p);
								pk.push_back(puiss);
								r = (s >= p) ? 1 : 0;
								puiss *= p;
								++it1;
				}
				while(it2!=m2.ak_.end())
				{
								s = *it2 + r;
								ak.push_back(s % p);
								pk.push_back(puiss);
								r = (s >= p) ? 1 : 0;
								puiss *= p;
								++it2;
				}

				if(r == 1)
				{
								ak.push_back(1);
								pk.push_back(puiss);
								puiss *= p;
				}

				pk.push_back(puiss);

				return P_Adic(ak,pk);
}



std::ostream& operator<<(std::ostream &stream, const P_Adic& padic)
{
				stream << padic.p_ << "|.";
				for(auto it=padic.ak_.begin(); it!= padic.ak_.end(); ++it)
								stream << *it << " ";

				return stream;
}
