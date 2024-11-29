#pragma once
#include <vector>
#include <functional>
#include <map>
#include "OdeSolver.h"


/**
 * \brief краевая задача
 * \tparam T тип результата
 */
template<class T>
class boundary_value_problem
{
	std::vector<std::function<T(double, const std::vector<T>&)>>
		_equations;                        // правые части ОДУ
	std::map<size_t, T> _left_conditions;  // граничные условия слева
	std::map<size_t, T> _right_conditions; // граничные условия справа
	double _epsilon; // погрешность

	std::vector<std::vector<T>> cauchy_problem_solutions() const;
	std::vector<T> get_the_initial_conditions(
		const std::vector<std::vector<T>>& solutions) const;
	square_matrix<T> slau_matrix(
		const std::vector<std::vector<T>>& solutions) const;
	std::vector<T> right_part() const;

public:
	boundary_value_problem(
		std::vector<std::function<T(double,
			const std::vector<T>&)>> functions,
		std::map<size_t, T> left_conditions,
		std::map<size_t, T> right_conditions, double epsilon = 0.1e-6);
	std::vector<T> solve() const;
	std::vector<std::vector<T>> solve(
		const std::vector<double>& points) const;
	T determinant() const;
	std::vector<std::vector<T>> initial_conditions() const;
	std::map<size_t, T> left_conditions() const;
	std::map<size_t, T> right_conditions() const;
	std::vector<T> numerator() const;
};

/**
 * \brief находит набор начальных условия для вспомогательных задач Коши
 * \return вектор, состоящий из набора начальных условий
 */
template <class T>
std::vector<std::vector<T>>
boundary_value_problem<T>::initial_conditions() const
{
	const auto size = _equations.size();
	const auto cols = size - _left_conditions.size();
	std::vector<std::vector<T>> result(size);
	size_t counter = 0;
	for (size_t i = 0; i < result.size(); i++)
	{
		if (_left_conditions.find(i) == _left_conditions.end())
		{
			result[i] = std::vector<T>(cols, 0);
			result[i][counter++] = 1;
		}
		else
		{
			result[i] = std::vector<T>(cols);
			result[i][0] = _left_conditions.at(i);
		}
	}
	return result;
}

template <class T>
std::map<size_t, T> boundary_value_problem<T>::left_conditions() const
{
	return _left_conditions;
}

template <class T>
std::map<size_t, T> boundary_value_problem<T>::right_conditions() const
{
	return _right_conditions;
}

template<class T>
std::vector<T> boundary_value_problem<T>::numerator() const
{
	const auto solutions = cauchy_problem_solutions();
	const auto mat = slau_matrix(solutions);
	const auto rp = right_part();
	std::vector<T> minors(rp.size());
	for (size_t i = 0; i < rp.size(); i++)
	{
		auto copy = mat;
		for (size_t ii = 0; ii < rp.size(); ii++)
		{
			copy[ii][i] = rp[ii];
		}
		minors[i] = copy.det();
	}
	return  matrix<T>(solutions.size(), solutions.front().size(), solutions).transpose() * minors;
}

/**
 * \brief находит набор решений вспомогательных задач Коши
 * \return значения решений в точке x = 1 в виде вектора
 */
template <class T>
std::vector<std::vector<T>>
boundary_value_problem<T>::cauchy_problem_solutions() const
{
	auto conditions = this->initial_conditions();
	const auto size = conditions[0].size();
	std::vector<std::vector<T>> solutions;
	OdeSolver<T> ode_solver(_equations, this->_epsilon, RKF_78);
	for (size_t i = 0; i < size; i++)
	{
		std::vector<T> initials(_equations.size());
		for (size_t ii = 0; ii < initials.size(); ii++)
		{
			initials[ii] = conditions[ii][i];
		}
		auto solution = ode_solver.solve(0, 1, initials);
		solutions.push_back(solution);
	}
	return solutions;
}

/**
 * \brief реализация метода пристрелки, решение СЛАУ
 * \param solutions набор решений вспомогательных задач Коши
 * \return вектор начальных условий для решения краевой задачи
 */
template<class T>
std::vector<T> boundary_value_problem<T>::get_the_initial_conditions(
	const std::vector<std::vector<T>>& solutions) const
{
	return slau_matrix(solutions).linSolve(right_part());
}

template<class T>
square_matrix<T> boundary_value_problem<T>::slau_matrix(
	const std::vector<std::vector<T>>& solutions) const
{
	std::vector<std::vector<T>> matrix;
	for (const auto& x : _right_conditions)
	{
		std::vector<T> row(_right_conditions.size());
		for (size_t i = 0; i < row.size(); i++)
		{
			row[i] = solutions[i][x.first];
		}
		matrix.push_back(row);
	}
	return square_matrix<T>(_right_conditions.size(), matrix);
}

template<class T>
std::vector<T> boundary_value_problem<T>::right_part() const
{
	std::vector<T> result;
	for (const auto& x : _right_conditions)
		result.push_back(x.second);
	return result;
}

/**
 * \brief определитель системы метода пристрелки
 * \return определитель системы метода пристрелки
 */
template <class T>
T boundary_value_problem<T>::determinant() const
{
	const auto solutions = cauchy_problem_solutions();
	return slau_matrix(solutions).det();
}

/**
 * \brief конструктор
 * \param functions правые части однородных уравнений системы
 * \param left_conditions условия на правом конце отрезка (условия линейные и разделённые)
 * \param right_conditions условия на правом конце отрезка
 * \param epsilon погрешность вычислений
 */
template<class T>
boundary_value_problem<T>::boundary_value_problem(
	std::vector<std::function<T(double, const std::vector<T>&)>> functions,
	std::map<size_t, T> left_conditions, std::map<size_t, T> right_conditions,
	double epsilon) :
	_equations(std::move(functions)),
	_left_conditions(std::move(left_conditions)),
	_right_conditions(std::move(right_conditions)),
	_epsilon(epsilon) {
}

/**
 * \brief решение краевой задачи на правом конце отрезка
 * \return решение краевой задачи на правом конце отрезка
 */
template<class T>
std::vector<T> boundary_value_problem<T>::solve() const
{
	std::vector<std::vector<T>> solutions = cauchy_problem_solutions();
	std::vector<T> initial_conditions = get_the_initial_conditions(solutions);
	std::vector<T> result(_equations.size(), 0);
	for (size_t i = 0; i < solutions.size(); i++)
	{
		result = result + initial_conditions[i] * solutions[i];
	}
	return result;
}

/**
 * \brief решение краевой задачи в наборе точек
 * \param points набор точек
 * \return решение краевой задачи в наборе точек
 */
template<class T>
std::vector<std::vector<T>> boundary_value_problem<T>::solve(
	const std::vector<double>& points) const
{
	const auto solutions = cauchy_problem_solutions();
	const auto initial_conditions = get_the_initial_conditions(solutions);
	std::vector<T> initials;
	size_t counter = 0;
	for (size_t i = 0; i < _equations.size(); i++)
	{
		if (_left_conditions.find(i) == _left_conditions.end())
		{
			initials.push_back(initial_conditions[counter++]);
		}
		else
		{
			initials.push_back(_left_conditions.at(i));
		}
	}
	OdeSolver<T> ode_solver(_equations, this->_epsilon, RKF_78);
	return ode_solver.solve(points, initials);
}