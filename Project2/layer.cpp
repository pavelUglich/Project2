#include "layer.h"

const std::complex<double> I = { 0,1 };


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
			return lambda(t) + 2 * mu(t);
		};
	return {
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return I * alpha * v[1] + v[2] / mu(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return I * alpha * lambda(x) * v[0] / bulk(x) + v[3] / bulk(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return (4.0 * alpha * alpha * mu(x) * (lambda(x) + mu(x)) / bulk(x) -
				rho(x) * kappa * kappa) * v[0] + I * alpha * lambda(x) * v[3] / bulk(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return -rho(x) * kappa * kappa * v[1] + I * alpha * v[2];
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
			return lambda(t) + 2 * mu(t);
		};
	std::vector<std::function<std::complex<double>(double,
		const std::vector<std::complex<double>>&)>> added = {
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return I * v[1] + I * alpha * v[5] + v[6] / mu(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return I * lambda(x) * v[0] / bulk(x) + I * alpha * lambda(x) * v[4] / bulk(x) + v[7] / bulk(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return 8.0 * alpha * mu(x) * (lambda(x) + mu(x)) / bulk(x) * v[0] + I * lambda(x) * v[3] / bulk(x)
			+ (4.0 * alpha * alpha * mu(x) * (lambda(x) + mu(x)) / bulk(x) -
				rho(x) * kappa * kappa) * v[4] + I * alpha * lambda(x) * v[7] / bulk(x);
		},
		[=](double x, const std::vector<std::complex<double>>& v)
		{
			return  I * v[2] - rho(x) * kappa * kappa * v[5] + I * alpha * v[6];
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

std::vector<std::complex<double>> layer::evaluate_roots(double kappa) const
{
	auto roots = real_roots(kappa);
	std::vector<std::complex<double>> result(roots.size());
	std::transform(roots.begin(), roots.end(), result.begin(), [](auto x) {return x.real(); });
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

std::complex<double> layer::dispersion_equation_derivative(std::complex<double> alpha, double kappa) const
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
	std::vector<std::complex<double>> result(initial_values.size());
	size_t j = 0;
	for (size_t i = 0; i < result.size(); i++)
	{
		try
		{
			const auto new_root = newton_method(initial_values[i], kappa);
			result[j++] = new_root;
		}
		catch (const std::exception&)
		{
			result.erase(result.begin() + i);
		}
	}
	return result;
}

std::complex<double> layer::roots(std::complex<double> initial_value, double kappa) const
{
	return newton_method(initial_value, kappa);
}

std::vector<double> layer::eigen_frequencies(double max_kappa, double alpha,
	size_t num_roots) const
{
	std::vector<double> result;
	const auto h = max_kappa / num_roots;
	double fa = dispersion_equation({ alpha,0 }, 0).real();
	for (size_t i = 0; i < num_roots; i++)
	{
		const double fb = dispersion_equation({ alpha,0 }, (i + 1) * h).real();
		if (fa * fb < 0)
		{
			result.push_back(secant_method(i * h, (i + 1) * h,
				[=](double x) {
					return dispersion_equation({ alpha,0 }, x).real();
				}));
		}
		fa = fb;
	}
	return result;
}


std::map<double, std::vector<double>> layer::dispersion_set(double max_kappa,
	size_t freqs) const
{
	std::map<double, std::vector<double>> result;
	const auto h = max_kappa / freqs;
	for (auto kappa = h; kappa < max_kappa; kappa += h)
	{
		std::cout << kappa << std::endl;
		result[kappa] = real_roots(kappa);
	}
	return result;
}

std::map<double, std::vector<double>> layer::imaginary_dispersion_set(
	double max_kappa, size_t freqs) const
{
	std::map<double, std::vector<double>> result;
	const auto h = max_kappa / freqs;
	for (auto kappa = h; kappa < max_kappa; kappa += h)
	{
		std::cout << kappa << std::endl;
		result[kappa] = imaginary_roots(kappa);
	}
	return result;
}


std::map<double, std::vector<double>> layer::dispersion_set_alpha(
	double max_alpha, double max_kappa, size_t freqs) const
{
	std::map<double, std::vector<double>> result;
	const auto h = max_alpha / freqs;
	for (auto kappa = 0.0; kappa < max_alpha; kappa += h)
	{
		result[kappa] = eigen_frequencies(max_kappa, kappa, 40);
	}
	return result;
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

std::map<double, std::complex<double>> layer::dispersional_curve(double max_kappa,
	std::complex<double> initial_value, size_t freqs) const
{
	std::map<double, std::complex<double>> result = { {0.0, initial_value } };
	const auto h = max_kappa / freqs;
	for (double kappa = h; kappa < max_kappa; kappa += h)
	{
		try
		{
			initial_value = newton_method(initial_value, kappa);
		}
		catch (...)
		{
			break;
		}
		result[kappa] = initial_value;
	}
	return result;
}
