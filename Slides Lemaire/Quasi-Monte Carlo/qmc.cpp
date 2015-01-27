#include <iostream>
#include <fstream>
#include <vector>
#include "monte-carlo.hpp"
#include "low_discrepancy.hpp"
#include "fcts.hpp"
using namespace std;

struct carre : public unary_function< vector<double> const &, double> {
	double operator()(vector<double> const & x) const { return x[0]*x[0]; }
};

int main() {
	sobol s(2);
	carre f;
	
	vector<double> result = monte_carlo(1e5, f, s);
	cout << result[0] << "\t" << result[1] << endl;

	return 0;
};

