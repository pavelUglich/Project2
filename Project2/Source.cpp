#include<functional>
#include <complex>
#include <fstream>
#include <vector>
#include "boundary_value_problem.h"
#include "layer.h"
#include"plots.h"
#include "Integral.h"



const std::function<double(double)> density = [](double t) {return 1.0; };
const std::function<double(double)> poisson_ratio = [](double t) {return 0.3; };
const std::function<double(double)> shear_modulus = [=](double t) {return 1.0/*young_modulus(t) / 2 / (1 + poisson_ratio(t))*/; };
const std::function<double(double)> lambda = [=](double t) {return 2 * poisson_ratio(t) * shear_modulus(t) / (1 - 2 * poisson_ratio(t)); };
const std::function<double(double)> bulk = [=](double t) {return lambda(t) + 2 * shear_modulus(t); };
const double pi = 2.0 * acos(0);


void showVector(const std::vector<double>& nodes,
	const std::vector<std::complex<double>>& complexValuedVector,
	const std::function<double(std::complex<double>)>& mapping,
	std::ostream& stream)
{
	for (size_t i = 0; i < complexValuedVector.size(); i++)
	{
		stream << "(" << nodes[i] << ", " << mapping(complexValuedVector[i]) << ") ";
	}
}

void addTheCurve(const std::vector<double>& nodes,
	const std::vector<std::complex<double>>& fieldHomogeneous,
	std::ofstream& stream, const std::string& color)
{
	stream << " \\addplot[line width = 0.25mm, smooth, ";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(nodes, fieldHomogeneous, [](auto x) { return x.real(); }, stream);
	stream << " };\n";
	stream << " \\addplot[smooth, dashed,";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(nodes, fieldHomogeneous, [](auto x) { return x.imag(); }, stream);
	stream << " };\n";
}


void plot_the_wave_field(const std::vector<double>& nodes,
	const std::map<std::string, std::vector<std::complex<double>>>& vectors,
	const std::string& fileName)
{
	std::ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.1]\\tiny\n";
	stream << "\\begin{axis}[grid]\n";
	for (const auto& item : vectors)
	{
		addTheCurve(nodes, item.second, stream, item.first);
	}
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}


int main() {
	const double kappa = 4;
	const layer l = { lambda, shear_modulus, density, kappa };
	Integral integral([=](auto alpha) {return l.transformant(alpha, kappa)[1]; }, kappa);
	size_t num_points = 50;
	double max = 5.0;
	double h = max / num_points;
	std::vector<std::complex<double>> w;
	std::vector<double> nodes;
	for (size_t i = 0; i < num_points; i++)
	{
		const auto x = (i + 0.5) * h;
		nodes.push_back(x);
		w.push_back(integral.evaluate(x));
		cout << w.back() << endl;
		// w.push_back(l.displacement(x)[1]);
	}
	plot_the_wave_field(nodes, { {"red", w} }, "horiz.txt");
	system("pause");
}