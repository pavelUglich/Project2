#include "plots.h"

std::ofstream initialize(const std::string& fileName) 
{
	std::ofstream os(fileName);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	return os;
}

void finalize(std::ofstream& os) 
{
	os << "\\end{axis}\n";
	os << "\\end{tikzpicture}\n";
	os.close();
}

void plot_the_dispersional_curves(
	std::map<double, std::vector<double>>& dispersionSet, 
	const std::string& fileName)
{
	std::ofstream os = initialize(fileName);
	const auto numberOfCurves = dispersionSet.rbegin()->second.size();
	for (auto& item : dispersionSet) {
		sort(item.second.begin(), item.second.end(),
			[](auto x, auto y) { return x > y; });
	}
	const auto maxFreq = dispersionSet.rbegin()->first;
	for (size_t i = 0; i < numberOfCurves; i++) {
		os << "\\addplot[smooth, black] plot coordinates{\n";
		for (auto x : dispersionSet) {
			if (i < x.second.size())
			{
				os << "(" << x.first << ", " << x.second[i] << ") ";
			}
		}
		os << "};\n";
	}
	finalize(os);
}


void plot_the_dispersional_curves(
	std::map<double, std::vector<std::complex<double>>>& dispersionSet,
	const std::string& fileName)
{
	std::ofstream os = initialize(fileName);
	const auto numberOfCurves = dispersionSet.begin()->second.size();
	for (auto& item : dispersionSet) {
		sort(item.second.begin(), item.second.end(),
			[](auto x, auto y)
			{
				if (x.imag() < 0.01)
				{
					if (y.imag() > 0.01)
					{
						return true;
					}
					return x.real() > y.real();
				}
				return x.imag() > y.imag();
			});
	}
	const auto maxFreq = dispersionSet.rbegin()->first;
	for (size_t i = 0; i < numberOfCurves; i++) {
		os << "\\addplot[smooth, red] plot coordinates{\n";
		for (auto x : dispersionSet) {
			if (i < x.second.size())
			{
				os << "(" << x.first << ", " << x.second[i].real() << ") ";
			}
		}
		os << "};\n";
		os << "\\addplot[smooth, blue] plot coordinates{\n";
		for (auto& x : dispersionSet) {
			if (i < x.second.size())
			{
				os << "(" << x.first << ", " << x.second[i].imag() << ") ";
			}
		}
		os << "};\n";
	}
	finalize(os);
}

void plot_the_dispersional_curves_alpha(
	std::map<double, std::vector<double>>& dispersionSet, 
	const std::string& fileName)
{
	std::ofstream os = initialize(fileName);
	const auto numberOfCurves = dispersionSet.begin()->second.size();
	for (auto& item : dispersionSet) {
		sort(item.second.begin(), item.second.end(),
			[](auto x, auto y) { return x < y; });
	}
	//const auto maxFreq = dispersionSet.rbegin()->first;
	for (size_t i = 0; i < numberOfCurves; i++) {
		os << "\\addplot[smooth, black] plot coordinates{\n";
		for (auto x : dispersionSet) {
			if (i < x.second.size())
			{
				os << "(" << x.first << ", " << x.second[i] << ") ";
			}
		}
		os << "};\n";
	}
	finalize(os);
}


void plot_the_dispersional_curves(
	std::map<double, std::complex<double>>& dispersionSet,
	const std::string& fileName)
{
	std::ofstream os = initialize(fileName);
	os << "\\addplot[smooth, blue] plot coordinates{\n";
	for (const auto& x : dispersionSet) {
		os << "(" << -x.second.imag() << ", " << x.first << ") ";
	}
	os << "};\n";
	os << "\\addplot[smooth, red] plot coordinates{\n";
	for (const auto& x : dispersionSet) {
		os << "(" << x.second.real() << ", " << x.first << ") ";
	}
	os << "};\n";
	finalize(os);
}

void plot_the_imaginary_roots(
	std::map<double, std::vector<double>>& dispersionSet, 
	const std::string& fileName)
{
	bool flag = false;
	std::ofstream os = initialize(fileName);
	for (const auto& x : dispersionSet) {
		if (x.second.size())
		{
			if (!flag)
			{
				os << "\\addplot[smooth, blue] plot coordinates{\n";
				flag = true;
			}
			os << "(" << -x.second.front() << ", " << x.first << ") ";
		}
		else
		{
			if (flag)
			{
				os << "};\n";
				flag = false;
			}
		}
	}
	finalize(os);
}
