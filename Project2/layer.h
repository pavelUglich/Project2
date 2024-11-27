#pragma once
#include <complex>
#include <functional>
#include <utility>
#include <vector>

#include "boundary_value_problem.h"

template<class Fun>
double secant_method(double a, double b, Fun equation, double eps = 0.1e-6) 
{
	double va = equation(a);
	double vb = equation(b);
	double cs = a;
	double cn = b;
	while (abs(cn - cs) > eps)
		if (va * vb < 0) {
			cs = cn;
			cn = b - (b - a) * vb / (vb - va);
			const double vc = equation(cn);
			if (va * vc < 0) {
				b = cn;
				vb = vc;
			}
			else if (va * vc > 0) {
				a = cn;
				va = vc;
			}
			else
				return cn;
		}
	return cn;
}

class layer
{
	// законы неоднородности
	std::function<double(double)> lambda, mu, rho;

	
	std::vector<
		std::function<
		std::complex<double>(double, const std::vector<std::complex<double>>&)
	>
	> equations(std::complex<double> alpha, double kappa) const;

	// вектор производных правых частей по параметру преобразования Фурье
	std::vector<
	std::function<
		std::complex<double>(double, const std::vector<std::complex<double>>&)
	>
	> extended(std::complex<double> alpha, double kappa) const;
		
	boundary_value_problem<std::complex<double>> bvp(std::complex<double> alpha, double kappa) const;

public:
	layer(std::function<double(double)> lambda, std::function<double(double)> mu,
	      std::function<double(double)> rho)
		: lambda(std::move(lambda)),
		mu(std::move(mu)),
		rho(std::move(rho))
	{

	}

	std::vector<std::complex<double>> evaluate_roots(double kappa) const;


	std::vector<std::complex<double>> transformant(std::complex<double> alpha, double kappa) const;
	std::complex<double> dispersion_equation(std::complex<double> alpha, double kappa) const;
	std::complex<double> dispersion_equation_derivative(std::complex<double> alpha, double kappa) const;

	std::vector<double> real_roots(double kappa, size_t num_roots = 80) const;
	std::vector<double> imaginary_roots(double kappa, size_t num_roots = 80) const;
	std::vector<std::complex<double>> roots(const std::vector<std::complex<double>>& initial_values, double kappa) const;
	std::complex<double> roots(std::complex<double> initial_values, double kappa) const;

	std::vector<double> eigen_frequencies(double max_kappa, double alpha = 0, size_t num_roots = 40) const;
	std::map<double, std::vector<double>> dispersion_set(double max_kappa, size_t freqs) const;
	std::map<double, std::vector<double>> imaginary_dispersion_set(double max_kappa, size_t freqs) const;

	std::map<double, std::vector<double>> dispersion_set_alpha(double max_alpha, double max_kappa, size_t freqs) const;
	std::complex<double> newton_method(std::complex<double> complex, double kappa, double eps = 0.1e-6) const;
	std::map<double, std::complex<double>> dispersional_curve(double max_kappa, std::complex<double> initial_value, size_t freqs) const;
};




