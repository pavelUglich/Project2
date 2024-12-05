#pragma once
#include <utility>
//#include"BoundaryValueProblem.h"
#include <functional>
#include <vector>
#include <complex>
using namespace std;

const vector<double> xi{ -0.989400934991650 , -0.944575023073233, -0.865631202387832 ,
-0.755404408355003 , -0.617876244402644 , -0.458016777657227 , -0.281603550779259, -0.950125098376374e-1 ,
0.950125098376374e-1 ,0.281603550779259, 0.458016777657227, 0.617876244402644, 0.755404408355003,
0.865631202387832, 0.944575023073233, 0.989400934991650 };


const vector<double> eta{ 0.0271524594117540, 0.0622535239386475, 0.0951585116824939 ,
0.124628971255534 , 0.149595988816577 , 0.169156519395002 , 0.182603415044922, 0.189450610455068,
0.189450610455068 ,0.18260341504492, 0.16915651939500, 0.14959598881657, 0.124628971255534,
0.0951585116824939, 0.0622535239386475, 0.0271524594117540 };

template<class Fun>
std::complex<double> integral_gauss(double a, double b, Fun f) {
	const double m = 0.5 * (a + b);
	const double n = 0.5 * (b - a);
	complex<double> result;
	for (size_t i = 0; i < 16; i++)
	{
		double y = m + n * xi[i];
		result += eta[i] * f(y);
	}
	return 0.5 * (b - a) * result;
}


class Integral
{
	function<complex<double>(complex<double>)> _integrand;
	double _kappa;

	complex<double> gauss(double a, double b, double x) const;

	complex<double> gauss(double a, double b, double x,
		const function<complex<double>(double)>& contour,
		const function<complex<double>(double)>& derivative) const;

public:
	complex<double> evaluate(double x) const;
	explicit Integral(const function<complex<double>(complex<double>)>& integrand,
		double kappa) :_integrand(integrand), _kappa(kappa) {}
	~Integral();
};

