#ifndef P_ADIC_HPP
#define P_ADIC_HPP

#include "stl_headers.hpp"

/*
	* Array of the first 255 prime numbers
	*/
static const int primes[255] = {
													2,    3,    5,    7,   11,   13,   17,   19,   23,
						29,   31,   37,   41,   43,   47,   53,   59,   61,   67,
						71,   73,   79,   83,   89,   97,  101,  103,  107,  109,
					113,  127,  131,  137,  139,  149,  151,  157,  163,  167,
					173,  179,  181,  191,  193,  197,  199,  211,  223,  227,
					229,  233,  239,  241,  251,  257,  263,  269,  271,  277,
					281,  283,  293,  307,  311,  313,  317,  331,  337,  347,
					349,  353,  359,  367,  373,  379,  383,  389,  397,  401,
					409,  419,  421,  431,  433,  439,  443,  449,  457,  461,
					463,  467,  479,  487,  491,  499,  503,  509,  521,  523,
					541,  547,  557,  563,  569,  571,  577,  587,  593,  599,
					601,  607,  613,  617,  619,  631,  641,  643,  647,  653,
					659,  661,  673,  677,  683,  691,  701,  709,  719,  727,
					733,  739,  743,  751,  757,  761,  769,  773,  787,  797,
					809,  811,  821,  823,  827,  829,  839,  853,  857,  859,
					863,  877,  881,  883,  887,  907,  911,  919,  929,  937,
					941,  947,  953,  967,  971,  977,  983,  991,  997, 1009,
				1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063,
				1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129,
				1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217,
				1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289,
				1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367,
				1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447,
				1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
				1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579,
				1583, 1597, 1601, 1607, 1609, 1613
};


/*************************************************************************************
	*
	* STRUCT P_Adic
	*
	* Store and allow manipulation of p-adic representations
	*
	*************************************************************************************/
struct P_Adic {

				typedef std::vector<int> Coeff;

				P_Adic(Coeff ak, Coeff pk) : ak_(ak), pk_(pk), p_(*(++pk.begin())) {}
				P_Adic(int n, int p);
				P_Adic(double x, int p);

				// Cast to a real number in [0,1]^'ak_.size()'
				operator double()
				{
								return std::inner_product(ak_.begin(), ak_.end(), ++pk_.begin(), 0.0, std::plus<double>(), std::divides<double>());
				}

				operator int()
				{
								return std::inner_product(ak_.begin(), ak_.end(), pk_.begin(), 0);
				}


				P_Adic operator++(int);
				P_Adic& operator++();

				int base() const { return p_; }

				// p-adic addition on [0,1], ie Kakutani addition
				friend P_Adic operator+(const P_Adic& m1, const P_Adic& m2);
				// Output double representation under the format 'p_'|.'ak_' (aks separated by single whitespaces)
				friend std::ostream& operator<<(std::ostream &stream, const P_Adic& padic);

public:
				// /!\ 'pk_' should have one more element than 'pk_' to allow conversion to 'double'
				Coeff ak_;
				Coeff pk_;
				int p_;

				void increment();
};

P_Adic operator+(const P_Adic& m1, const P_Adic& m2);
std::ostream& operator<<(std::ostream &stream, const P_Adic& padic);


typedef std::vector<P_Adic> Vect_P_Adic;
std::ostream& operator<<(std::ostream &stream, const Vect_P_Adic& v);

#endif // P_ADIC_HPP
