#pragma once
#include "Matrix.h"

template<class T>
class square_matrix : public matrix<T>
{
	std::vector<T> transform_the_right_part(
		const std::vector<T>& rightPart, const std::vector<size_t>& p) const;
	std::vector<T> there_and_back_again(const std::vector<std::vector<T>>& lu,
		const std::vector<T>& rp) const;


public:
	std::vector<std::vector<T>> lu_decompose(std::vector<size_t>& p) const;

	/**
	 * \brief конструктор задаёт размеры квадратной матрицы (все конструкторы
	 * должны обращаться к конструкторам базового класса)
	 * \param rows количество строк
	 */
	square_matrix(size_t rows);

	/**
	 * \brief конструктор задаёт размеры квадратной матрицы и заполняет её
	 * значеним value
	 * \param rows количество строк
	 * \param value значение
	 */
	square_matrix(size_t rows, const T& value);

	/**
	 * \brief конструктор задаёт размеры квадратной матрицы и заполняет её
	 * значениями из вектора
	 * \param rows размер матрицы
	 * \param elements вектор значений
	 */
	square_matrix(size_t rows, const std::vector<std::vector<T>>& elements);

	square_matrix(const matrix<T>& matrix);


	/**
	 * \brief определитель матрицы
	 * \return определитель матрицы
	 */
	T det() const;

	/**
	 * \brief решение системы линейных алгебраических уравнений
	 * \param rightPart вектор правой части
	 * \return решение СЛАУ
	 */
	std::vector<T> linSolve(const std::vector<T>& rightPart) const;

	/**
	 * \brief обращение матрицы
	 * \return обратная матрица
	 */
	square_matrix operator~() const;

	/**
	 * \brief деление на другую матрицу (умножение справа на матрицу, обратную
	 * к другой)
	 * \param squareMatrix
	 * \return
	 */
	square_matrix& operator/=(const square_matrix<T>& squareMatrix);
	square_matrix operator/(const square_matrix<T>& squareMatrix) const;
};

template <class T>
square_matrix<T>::square_matrix(size_t rows) : matrix<T>(rows, rows)
{
}

template <class T>
square_matrix<T>::square_matrix(size_t rows, const T& value) :
	matrix<T>(rows, rows, value)
{
}

template <class T>
square_matrix<T>::square_matrix(size_t rows,
	const std::vector<std::vector<T>>& elements) :
	matrix<T>(rows, rows, elements)
{
}

template<class T>
square_matrix<T>::square_matrix(const matrix<T>& matrix) :matrix<T>(
	matrix.get_rows(), matrix.getColumns(), matrix.get_elements())
{
	if (matrix.get_rows() != matrix.getColumns())
	{
		throw std::invalid_argument("");
	}
}

template<class T>
T square_matrix<T>::det() const
{
	std::vector<size_t> p(this->get_rows());
	auto lu = lu_decompose(p);
	T result = 1;
	for (size_t i = 0; i < p.size(); i++)
	{
		result *= lu[i][i];
	}
	int parity = 0;
	for (size_t i = 0; i < p.size(); i++)
	{
		for (size_t ii = i + 1; ii < p.size(); ii++)
		{
			if (p[i] > p[ii])
			{
				++parity;
			}
		}
	}
	return result * (parity % 2 == 0 ? 1.0 : -1.0);
}

template <class T>
std::vector<std::vector<T>> square_matrix<T>::lu_decompose(
	std::vector<size_t>& p) const
{
	auto elems = this->get_elements();
	p.resize(this->get_rows());
	for (size_t i = 0; i < this->get_rows(); i++)
		p[i] = i;
	for (size_t i = 0; i < this->get_rows(); i++) {
		auto max = abs(elems[i][i]);
		size_t imax = i;
		for (size_t k = i; k < this->get_rows(); k++) {
			if (abs(elems[k][i]) > max) {
				max = abs(elems[k][i]);
				imax = k;
			}
		}
		if (max < DBL_EPSILON) throw std::invalid_argument("");
		if (imax != i) {
			std::swap(p[i], p[imax]);
			std::swap(elems[i], elems[imax]);
		}
		for (size_t j = i + 1; j < this->get_rows(); j++) {
			elems[j][i] /= elems[i][i];
			for (size_t k = i + 1; k < this->get_rows(); k++)
				elems[j][k] -= elems[j][i] * elems[i][k];
		}
	}
	return elems;
}

template <class T>
std::vector<T> square_matrix<T>::transform_the_right_part(
	const std::vector<T>& rightPart, const std::vector<size_t>& p) const
{
	std::vector<T> rp(this->get_rows());
	for (size_t i = 0; i < this->get_rows(); i++)
	{
		rp[i] = rightPart[p[i]];
	}
	return rp;
}

template <class T>
std::vector<T> square_matrix<T>::there_and_back_again(
	const std::vector<std::vector<T>>& lu,
	const std::vector<T>& rp) const
{
	std::vector<T> y(this->get_rows());
	for (size_t i = 0; i < this->get_rows(); i++)
	{
		T sum = 0;
		for (size_t j = 0; j < i; j++)
		{
			sum += y[j] * lu[i][j];
		}
		y[i] = rp[i] - sum;
	}
	for (int i = static_cast<int>(this->get_rows() - 1); i >= 0; --i)
	{
		T sum = 0;
		for (size_t j = i + 1; j < this->get_rows(); ++j)
		{
			sum += y[j] * lu[i][j];
		}
		y[i] = (y[i] - sum) / lu[i][i];
	}
	return y;
}

template<class T>
std::vector<T>
square_matrix<T>::linSolve(const std::vector<T>& rightPart) const
{
	std::vector<size_t> p(this->get_rows());
	const auto lu = lu_decompose(p);
	const auto rp = transform_the_right_part(rightPart, p);
	return there_and_back_again(lu, rp);
}

template<class T>
square_matrix<T> square_matrix<T>::operator~() const
{
	square_matrix<T> result(this->get_rows());
	std::vector<size_t> p(this->get_rows());
	auto lu = lu_decompose(p);
	for (size_t i = 0; i < this->get_rows(); i++)
	{
		std::vector<T> right_part(this->get_rows(), 0);
		right_part[i] = 1;
		auto rp = transform_the_right_part(right_part, p);
		auto column = there_and_back_again(lu, rp);
		for (size_t ii = 0; ii < this->get_rows(); ii++)
		{
			result[ii][i] = column[ii];
		}
	}
	return result;
}

template<class T>
square_matrix<T>& square_matrix<T>::operator/=(
	const square_matrix<T>& squareMatrix)
{
	*this *= ~squareMatrix;
	return *this;
}

template<class T>
square_matrix<T> square_matrix<T>::operator/(
	const square_matrix<T>& squareMatrix) const
{
	auto temp(*this);
	return temp /= squareMatrix;
}