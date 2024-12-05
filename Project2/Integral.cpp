#include "Integral.h"
#include<functional>
using namespace std;
const double pi = 3.1415926538;
const std::complex<double> I = { 0,1 };


complex<double> Integral::gauss(double a, double b, double x,
	const function<complex<double>(double)>& contour,
	const function<complex<double>(double)>& derivative) const
{
	const double m = 0.5 * (a + b);
	const double n = 0.5 * (b - a);
	complex<double> result;
	for (size_t i = 0; i < 16; i++)
	{
		const double y = m + n * xi[i];
		const complex<double> alpha = contour(y);
		const complex<double> da = derivative(y);
		const auto symbol = _integrand(-alpha) * exp(I * alpha * x) * da;
		const auto symbol1 = _integrand(alpha) * exp(-I * alpha * x) * da;
		result += eta[i] * (symbol + symbol1);
	}
	return 0.5 * result;
}

complex<double> Integral::gauss(double a, double b, double x) const
{
	const double m = 0.5 * (a + b);
	const double n = 0.5 * (b - a);
	complex<double> result;
	for (size_t i = 0; i < 16; i++)
	{
		const double y = m + n * xi[i];
		const auto symbol = _integrand(-y) * exp(I * y * x) + _integrand(y)
			* exp(-I * y * x);
		result += 0.5 * eta[i] * symbol;
	}
	return result;
}


complex<double> Integral::evaluate(double x) const {
	complex<double> result(0, 0);
	double delta = _kappa / 10;
	for (size_t j = 0; j < 8; j++) {
		result += gauss(0.25 * _kappa * j, 0.25 * _kappa * (j + 1), x,
			[delta, this](double y) {return complex<double>(y, -1 * delta * (1 - (y - _kappa) * (y - _kappa) / (_kappa * _kappa))); },
			[delta, this](double y) {return complex<double>(1, 2 * delta * (y - _kappa) / (_kappa * _kappa)); });
	}
	result *= 0.125 * _kappa;
	double segment = x < 0.1 ? 1.0 : abs(2 * pi / x);
	double x0 = 2 * _kappa;
	complex<double> term;
	do
	{
		term = gauss(x0, x0 + segment, x);
		result += term * 0.5 * segment;
		x0 += segment;
	} while (abs(term) > 0.3e-4);//*/
	return -result / pi;
}

Integral::~Integral() = default;
