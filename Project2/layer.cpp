#include "layer.h"

const std::complex<double> I = { 0,1 };

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


/**
 * \brief вектор правых частей однородной системы уравнений
 * \param alpha параметр преобразования Фурье
 * \param kappa частота
 * \return значения правых частей системы, найденные при alpha, kappa
 */
std::vector<std::function<std::complex<double>(double,
	const std::vector<std::complex<double>>&)>> layer::equations(
		std::complex<double> alpha, double kappa) const
{
	const std::function<double(double)> bulk = [=](double t)
	{
		return _lambda(t) + 2 * _mu(t);
	};
	return {
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return I * alpha * v[1] + v[2] / _mu(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return I * alpha * _lambda(x) * v[0] / bulk(x) + v[3] / bulk(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return (4.0 * alpha * alpha * _mu(x) * (_lambda(x) + _mu(x)) / bulk(x) -
				_rho(x) * kappa * kappa) * v[0] + I * alpha * _lambda(x) * v[3] / bulk(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return -_rho(x) * kappa * kappa * v[1] + I * alpha * v[2];
		},
	};
}

/**
 * \brief вектор правых частей, дополненный производными по alpha
 * \param alpha параметр преобразования Фурье
 * \param kappa частота
 * \return значения правых частей системы, найденные при alpha, kappa
 */
std::vector<std::function<std::complex<double>(double,
	const std::vector<std::complex<double>>&)>> layer::extended(
		std::complex<double> alpha, double kappa) const
{
	auto equations = this->equations(alpha, kappa);
	const std::function<double(double)> bulk = [=](double t)
	{
		return _lambda(t) + 2 * _mu(t);
	};
	std::vector<std::function<std::complex<double>(double,
		const std::vector<std::complex<double>>&)>> added = {
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return I * v[1] + I * alpha * v[5] + v[6] / _mu(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return I * _lambda(x) * v[0] / bulk(x) + I * alpha * _lambda(x) * v[4] / bulk(x) + v[7] / bulk(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return 8.0 * alpha * _mu(x) * (_lambda(x) + _mu(x)) / bulk(x) * v[0] + I * _lambda(x) * v[3] / bulk(x)
			+ (4.0 * alpha * alpha * _mu(x) * (_lambda(x) + _mu(x)) / bulk(x) -
				_rho(x) * kappa * kappa) * v[4] + I * alpha * _lambda(x) * v[7] / bulk(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return  I * v[2] - _rho(x) * kappa * kappa * v[5] + I * alpha * v[6];
		},
	};
	equations.insert(equations.end(), added.begin(), added.end());
	return equations;
}

boundary_value_problem<std::complex<double>> layer::bvp(std::complex<double> alpha, double kappa) const
{
	return
	{
		equations(alpha, kappa), {{0,{0,0}}, {1,{0,0}}},{{2,{0,0}}, {3,{1,0}}}
	};
}

std::vector<std::complex<double>> layer::evaluate_roots() const
{
	auto roots = real_roots(_kappa);
	std::vector<std::complex<double>> result(roots.size());
	std::transform(roots.begin(), roots.end(), result.begin(), 
		[](std::complex<double> x)->double {return x.real(); });
	roots = imaginary_roots(_kappa);
	result.insert(result.end(), roots.begin(), roots.end());
	auto complex_roots = this->roots(initial_root_values, _kappa);
	result.insert(result.end(), complex_roots.begin(), complex_roots.end());
	return result;
}

std::vector<std::complex<double>> layer::transformant(
	std::complex<double> alpha, double kappa) const
{
	return bvp(alpha, kappa).solve();
}

std::complex<double> layer::dispersion_equation(std::complex<double> alpha,
	double kappa) const
{
	return bvp(alpha, kappa).determinant();
}

std::complex<double> layer::dispersion_equation_derivative(
	std::complex<double> alpha, double kappa) const
{
	const auto bvp = this->bvp(alpha, kappa);
	auto initials = bvp.initial_conditions();
	const auto extended = this->extended(alpha, kappa);
	// транспозиция
	matrix<std::complex<double>> transposed = { initials.size(), initials.front().size(), initials };
	transposed = transposed.transpose();
	square_matrix<std::complex<double>> matrix(transposed.get_rows());
	square_matrix<std::complex<double>> derivative(transposed.get_rows());
	for (size_t i = 0; i < transposed.get_rows(); i++)
	{
		auto& initial_conditions = transposed[i];
		initial_conditions.resize(initial_conditions.size() * 2, 0);
		OdeSolver<std::complex<double>> ode_solver = { extended };
		auto solution = ode_solver.solve(0, 1, initial_conditions);
		size_t counter = 0;
		for (const auto& item : bvp.right_conditions())
		{
			matrix[i][counter] = solution[item.first];
			derivative[i][counter++] = solution[item.first + 4];
		}
	}
	square_matrix<std::complex<double>> common(matrix);
	std::complex<double> result = { 0,0 };
	for (size_t i = 0; i < common.get_rows(); i++)
	{
		common[i] = derivative[i];
		result += common.det();
		common[i] = matrix[i];
	}
	return result;
}

std::vector<double> layer::real_roots(double kappa, size_t num_roots) const
{
	std::vector<double> result;
	const double top = 2.0 * kappa;
	const auto h = top / num_roots;
	double fa = this->dispersion_equation({ 0,0 }, kappa).real();
	for (size_t i = 0; i < num_roots; i++)
	{
		const double fb = this->dispersion_equation({ (i + 1) * h,0 }, kappa).real();
		if (fa * fb < 0)
		{
			result.push_back(secant_method(i * h, (i + 1) * h, [=](double x) {
				return dispersion_equation({ x,0 }, kappa).real(); }));
		}
		fa = fb;
	}
	return result;
}

std::vector<double> layer::imaginary_roots(double kappa, size_t num_roots) const
{
	std::vector<double> result;
	const double top = 5.0;
	const auto h = top / num_roots;
	double fa = dispersion_equation({ 0,0 }, kappa).real();
	for (size_t i = 0; i < num_roots; i++)
	{
		const auto fb = dispersion_equation({ 0, (i + 1) * h }, kappa).real();
		if (fa * fb < 0)
		{
			const auto root = secant_method(i * h, (i + 1) * h,
				[=](double x)
				{
					return dispersion_equation({ 0,x }, kappa).real();
				});
			std::cout << kappa << ", " << root << std::endl;
			result.push_back(root);
		}
		fa = fb;
	}
	return result;
}

std::vector<std::complex<double>> layer::roots(
	const std::vector<std::complex<double>>& initial_values, double kappa) const
{
	std::vector<std::complex<double>> result;
	size_t j = 0;
	for (size_t i = 0; i < initial_root_values.size(); i++)
	{
		try
		{
			const auto new_root = newton_method(initial_values[i], kappa);
			result.push_back(new_root);
		}
		catch (const std::exception&)
		{
			continue;
		}
	}
	return result;
}

std::complex<double> layer::roots(std::complex<double> initial_value, double kappa) const
{
	return newton_method(initial_value, kappa);
}


std::complex<double> layer::newton_method(std::complex<double> complex,
	double kappa, double eps) const
{
	size_t iter = 0;
	auto numerator = [=](std::complex<double> x) {return this->dispersion_equation(x, kappa); };
	auto denumerator = [=](std::complex<double> x) {return this->dispersion_equation_derivative(x, kappa); };
	auto a = numerator(complex) / denumerator(complex);
	while (abs(a) > eps)
	{
		complex -= a;
		a = numerator(complex) / denumerator(complex);
		if (++iter > 100)
		{
			throw std::exception();
		}
	}
	return complex;
}
