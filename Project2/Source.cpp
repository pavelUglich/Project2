#include<functional>
#include <complex>
#include <fstream>
#include <vector>
#include "boundary_value_problem.h"
#include "layer.h"
#include"plots.h"



const std::function<double(double)> density = [](double t) {return 1.0; };
const std::function<double(double)> poisson_ratio = [](double t) {return 0.3; };
const std::function<double(double)> shear_modulus = [=](double t) {return 1.0/*young_modulus(t) / 2 / (1 + poisson_ratio(t))*/; };
const std::function<double(double)> lambda = [=](double t) {return 2 * poisson_ratio(t) * shear_modulus(t) / (1 - 2 * poisson_ratio(t)); };
const std::function<double(double)> bulk = [=](double t) {return lambda(t) + 2 * shear_modulus(t); };
const double pi = 2.0 * acos(0);






int main() {
	const layer l = { lambda, shear_modulus, density, 4 };
	/*auto v = l.evaluate_roots();
	for (auto x : v)
	{
		std::cout << x << ", ";
	}*/
	system("pause");
}