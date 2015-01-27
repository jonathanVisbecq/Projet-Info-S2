#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP
#include <functional>
#include <vector>
#include <cmath>
#include <list>

/******* MONTE CARLO **********/
template <typename Gen>
std::vector<double> monte_carlo(int n, Gen G)
{
	std::vector<double> result(3,0);
	double x;
	for (int j = 0; j < n; j++) {
		x = G();
		result[0] += x;
		result[1] += x*x;
	}
	result[0] /= (double) n;
	result[1] = (result[1] - n*result[0]*result[0])/(double)(n-1);
	result[2] = 1.96*sqrt(result[1]/(double) n);
	return result;
}


template <typename Gen, typename Fct>
std::vector<double> monte_carlo(int n, Fct f, Gen G)
{
	return monte_carlo(n, compose(f, G)); 
};

template <typename Gen, typename Fct1, typename Fct2>
std::vector<double> monte_carlo(int n, Fct1 f1, Fct2 f2, Gen G)
{
	return monte_carlo(n, compose(f1, compose(f2, G))); 
};


/****** Composition d'une fonction et d'une VA ************/
template <class _Result>
struct generator 
{
	typedef _Result result_type;
};

template <typename Fct, typename VA>
struct compose_t : public generator< typename Fct::result_type >
{
	compose_t(Fct f, VA X) : f(f), X(X) {};
	typename Fct::result_type operator()() { 
		return f(X()); 
	};
	private:
		Fct f; VA X;
};

template <typename Fct, typename VA>
inline compose_t<Fct, VA>
compose(Fct f, VA X) {
    return compose_t<Fct, VA>(f, X);
};

/******** COMPOSITION DE FONCTIONS ***************/
template <typename Fct1, typename Fct2>
struct compo_f_t
	: public std::unary_function<typename Fct2::argument_type, 
								 typename Fct1::result_type>
{
	compo_f_t(Fct1 f, Fct2 g) : f(f), g(g) {};
	typename Fct1::result_type operator()(typename Fct2::argument_type x) { 
		typename Fct1::result_type result(f(g(x)));
		return result;
	};
	private:
		Fct1 f; Fct2 g;
};

template <typename Fct1, typename Fct2>
inline compo_f_t<Fct1, Fct2>
compo_f(Fct1 f, Fct2 g) {
    return compo_f_t<Fct1, Fct2>(f, g);
};

/****** ANTITHETIQUE ********/
template <typename Fct, typename Trans>
struct antithetic_t 
	: public std::unary_function<typename Fct::argument_type, 
								 typename Fct::result_type>
{
    antithetic_t(Fct f, Trans T) 
    	: f(f), T(T) {};
    typename Fct::result_type 
	operator()(const typename Fct::argument_type &x) {
		return 0.5*(f(x) + f(T(x)));
    }
	private:
		Fct f; Trans T;
};

template <typename Fct, typename Trans>
inline antithetic_t<Fct, Trans>
antithet(Fct f, Trans T) { 
    return antithetic_t<Fct, Trans>(f, T);
};

/****** VARIABLE DE CONTROLE ********/
template <typename Fct1, typename Fct2>
struct var_control_t 
	: public std::unary_function<typename Fct1::argument_type, 
								 typename Fct1::result_type>
{
    var_control_t(Fct1 f, Fct2 g) 
    	: f(f), g(g) {};
	typename Fct1::result_type
    operator()(const typename Fct1::argument_type &x) const {
		return (f(x) - g(x));
    }
	private:
		Fct1 f; Fct2 g; 
};

template <typename Fct1, typename Fct2>
inline var_control_t<Fct1, Fct2>
var_control(Fct1 f, Fct2 g) {
    return var_control_t<Fct1, Fct2>(f, g);
}

/****** VARIABLE DE CONTROLE ADAPTATIVE ********/
template <typename Fct1, typename Fct2>
struct adapt_var_control_t 
    : public std::unary_function<typename Fct1::argument_type, 
                                 typename Fct1::result_type>
{
    adapt_var_control_t(Fct1 f, Fct2 g, double &lambda) 
        : f(f), g(g), lambda(lambda), cov(0), var(0) {};
    typename Fct1::result_type
    operator()(const typename Fct1::argument_type &x) {
        double f_x = f(x), g_x = g(x), result = f_x - lambda*g_x;
        cov += f_x*g_x;
        var += g_x*g_x;
        lambda = var > 0 ? cov / var : lambda;
//		std::cout << cov  << "\t" << var << "\t" << lambda << std::endl;
        return result;
    }
    private:
        Fct1 f; Fct2 g;
		double cov, var;
		double &lambda;
};

template <typename Fct1, typename Fct2>
inline adapt_var_control_t<Fct1, Fct2>
var_control(Fct1 f, Fct2 g, double &lambda) {
    return adapt_var_control_t<Fct1, Fct2>(f, g, lambda);
}

/******* STRATIFICATION **********/
template <typename Fct, typename Gen>
std::vector<double> stratification(int n, Fct f, 
								   std::list<double> pk,
								   std::list<Gen> Gk) 
{
	std::vector<double> result(3, 0), res_tmp; 
	std::list<double>::iterator it_pk = pk.begin();
	typename std::list<Gen>::iterator it_Gk = Gk.begin();
	double var_strate_k = 0;
	while (it_pk != pk.end()) {
		unsigned nk = floor(n*(*it_pk));
		res_tmp = monte_carlo(nk, f, *it_Gk);
		result[0] += (*it_pk) * res_tmp[0];
		result[1] += (*it_pk) * res_tmp[1];
		it_pk++; it_Gk++;
	}
	result[2] = 1.96*sqrt(result[1]/(double) n);
	return result;	
}


/****** EXPONENTIAL TILTING **********/
template <typename Fct>
struct expo_tilting_gauss : public std::unary_function<double, double> {
	expo_tilting_gauss(Fct f, double theta) : f(f), theta(theta) {};
	double operator()(double x) {
		return f(x + theta) * exp(-theta*x - 0.5*theta*theta); 
	};
	private:
		Fct f;
		double theta;
};

template <typename Fct>
expo_tilting_gauss<Fct> expo_tilting(Fct f, double theta) {
	return expo_tilting_gauss<Fct>(f, theta);
};

/****** EXPONENTIAL TILTING **********/
/*
template <typename Fct>
struct gauss2expo : public std::unary_function<double, double> {
	gauss2expo(Fct f) : f(f) {};
	double operator()(double x) {
		return f(x); 
	};
	private:
		Fct f;
		double theta;
};

template <typename Fct>
expo_tilting_gauss<Fct> expo_tilting(Fct f, double theta) {
	return expo_tilting_gauss<Fct>(f, theta);
};*/
#endif
