#include<functional>
#include <complex>
#include <fstream>
#include <vector>
#include "boundary_value_problem.h"
#include "layer.h"
#include"plots.h"

/**
 * \brief Значения комплексных корней дисперсионного уравнения для дальнейшего
 * поиска методом Ньютона
 */
const std::vector<std::complex<double>> initial_root_values = {
	{ 0.0, 0.912784940892 },
	{ 1.45379095058, 2.57038756409 },
	{ 2.23310278473, 5.91307674312 },
	{ 2.64719022914, 9.13968435484 },
	{ 2.93599785383, 12.3311277936 },
	{ 3.15887074448, 15.5061637747 },
	{ 3.34064922434, 18.6720041867 },
	{ 3.49425096472, 21.8321042181 },
	{ 3.62729203504, 24.9883519086 },
	{ 3.74465280440, 28.1418753540 },
	{ 3.84965331757, 31.2933940416 },
	{ 3.94465662336, 34.4433903917 },
	{ 4.03140538392, 37.5922009659 },
	{ 4.11122209138, 40.7400682812 },
	{ 4.18513437513, 43.8871718475 },
	{ 4.25395675138, 47.0336475894 },
	{ 4.31834583247, 50.1796004514 },
	{ 4.37883872007, 53.3251128379 },
	{ 4.43588038346, 56.4702504167 },
	{ 4.48984361388, 59.6150662064 },
	{ 4.54104384863, 62.7596035120 }
};

const std::function<double(double)> density = [](double t) {return 1.0; };
const std::function<double(double)> poisson_ratio = [](double t) {return 0.3; };
const std::function<double(double)> shear_modulus = [=](double t) {return 1.0/*young_modulus(t) / 2 / (1 + poisson_ratio(t))*/; };
const std::function<double(double)> lambda = [=](double t) {return 2 * poisson_ratio(t) * shear_modulus(t) / (1 - 2 * poisson_ratio(t)); };
const std::function<double(double)> bulk = [=](double t) {return lambda(t) + 2 * shear_modulus(t); };
const double pi = 2.0 * acos(0);






int main() {
	const layer l = { lambda, shear_modulus, density };
	double kappa = 0.01;
	//std::map<double, std::vector<std::complex<double>>> ds;	
	std::map<double, std::complex<double>> ds;
	auto r = l.dispersion_set_alpha(9, 9, 50);
	auto r = l.imaginary_dispersion_set(9, 100);
	/*while (kappa < 10)
	{
		ds[kappa] = r;
		kappa += 0.1;
		try
		{
			r = l.roots(r, kappa);
		}
		catch (const std::exception&)
		{
			break;
		}
	}*/
	plot_the_imaginary_roots(r, "ds.txt");
	system("pause");
}