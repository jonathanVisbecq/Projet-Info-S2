#include <iostream>
#include "var_alea.hpp"

int main() {
	init_alea();

	uniform U(0,1);
	std::cout << U() << std::endl;
	
	expo E(1);
	std::cout << E() << std::endl;
	
	gaussian G(0,1);
	std::cout << G() << std::endl;
	
	chi_deux X(1);
	std::cout << X() << std::endl;
	
	inverse_gaussian Y(0.5,1);
	std::cout << Y() << std::endl;
	
	normal_inverse_gaussian Z(0.8,0.1,0,2);
	std::cout << Z() << std::endl;

	return 0;
};
