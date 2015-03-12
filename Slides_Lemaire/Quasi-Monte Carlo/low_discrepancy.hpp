#ifndef LOW_DISCREPANCY_HPP
#define LOW_DISCREPANCY_HPP
#include <iostream>
#include <cmath>
#include <climits>
#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_qrng.h>

class p_adic {
	public: 
		typedef std::list<int> coeff;
		p_adic(coeff ak, coeff pk) : ak(ak), pk(pk), p(*(++pk.begin())) {};
		p_adic(int n, int p = 2);
		p_adic(double x, int p = 2);
		operator double() {
			return std::inner_product(ak.begin(), ak.end(), ++pk.begin(), 
					0.0, std::plus<double>(), std::divides<double>());
		}
		operator int() {
			return std::inner_product(ak.begin(), ak.end(), pk.begin(), 0); 
		}
		friend p_adic operator+(const p_adic &m, const p_adic &o);
		p_adic operator++(int) {
			p_adic copie = *this;
			increment();
			return copie;
		};
		p_adic& operator++() {
			increment();
			return (*this);
		};
		friend struct halton;
		friend struct kakutani;
		friend struct faure;
	protected:
		void increment();
		p_adic add(const p_adic &x) const;
		int p;
		coeff ak, pk;
};

static int primes[255] = {
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

struct halton {
	typedef std::vector<double> result_type;
	typedef std::vector<p_adic> list_p_adic;
	halton(const list_p_adic &x) : nk(x), result(x.size()) {};
	halton(int dimension) : result(dimension) {
		for (int k = 0; k < dimension; k++)
			nk.push_back(p_adic((int) 1, primes[k]));
	};
	halton(int dimension, double x_init[]) : result(dimension) {
		for (int k = 0; k < dimension; k++)
			nk.push_back(p_adic(x_init[k], primes[k]));
	};
	result_type operator()() {
		result_type::iterator ir = result.begin();
		list_p_adic::iterator ink = nk.begin();
		while (ink != nk.end()) {
			*ir++ = (double) (*ink++)++;
		}
		return result;
	}
	protected:
		list_p_adic nk;
		result_type result;
};

struct kakutani {
	typedef std::vector<double> result_type;
	typedef std::vector<p_adic> list_p_adic;
	kakutani(const list_p_adic &x, const list_p_adic &y) : xk(x), yk(y), result(x.size()) {};
	result_type operator()() {
		std::transform(xk.begin(), xk.end(), yk.begin(), result.begin(), op);
		return result;
	}
	protected:
		list_p_adic xk;
		list_p_adic yk;
		result_type result;
		struct t_op {
			double operator()(p_adic &x, const p_adic &y) { 
				x = x + y;
				return (double) x;
			}; 
		} op;
};

struct faure {
	typedef std::vector<double> result_type;
	typedef std::vector< std::vector<int> > coeff_binom;
	faure(int dimension, const p_adic &x) 
		: x(x), result(dimension), Comb(32, std::vector<int>(32,0))  {
		for (int n = 0; n < 32; n++) {
			Comb[0][n] = 1;
		}
		for (int k = 1; k < 32; k++) {
			Comb[k][0] = 1;
			for (int n = 1; n <= 32-k; n++){
				Comb[k][n] = Comb[k][n-1] + Comb[k-1][n];
			}
		}
		std::cout << x.p << std::endl;
	};
	result_type operator()() {
		p_adic xi = x++;
		result_type::iterator it = result.begin();
		while (it != result.end()) {
			*it++ = (double) xi;
			xi = T(xi);
		}
		return result;
	}
	p_adic T(const p_adic &y) {
		p_adic::coeff bk;
		coeff_binom::const_iterator it_Comb = Comb.begin();
		p_adic::coeff::const_iterator it = y.ak.begin();
		while (it != y.ak.end()) {
			bk.push_back(std::inner_product(it, y.ak.end(), (*it_Comb).begin(), 0) % y.p);
			it++; it_Comb++;
		}
		return p_adic(bk, y.pk);
	}
	protected:
		p_adic x;
		result_type result;
		coeff_binom Comb; 
};


struct sobol {
	typedef std::vector<double> result_type;
	sobol(int dimension) : dimension(dimension), result(dimension) {
		q = gsl_qrng_alloc(gsl_qrng_sobol, dimension);
	}
	sobol(sobol const & o) {
		dimension = o.dimension;
		q = gsl_qrng_clone(o.q);
		result = o.result;
	}
	sobol & operator=(sobol const & o) {
		if (this != &o) {
			dimension = o.dimension;
			gsl_qrng_memcpy(q, o.q);
			result = o.result;
		}
		return *this;
	}
	~sobol() { gsl_qrng_free(q); }
	result_type operator()() {
		gsl_qrng_get(q, &(*result.begin()));
		return result;
	}
	protected: 
		int dimension;
		gsl_qrng * q;
		result_type result;
};

#endif
