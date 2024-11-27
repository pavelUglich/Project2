#pragma once
#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<map>
#include <algorithm>
#include <complex>

void plot_the_dispersional_curves(
	std::map<double, std::vector<double>>& dispersionSet,
	const std::string& fileName);

void plot_the_dispersional_curves(
	std::map<double, std::vector<std::complex<double>>>& dispersionSet,
	const std::string& fileName);

void plot_the_dispersional_curves_alpha(
	std::map<double, std::vector<double>>& dispersionSet,
	const std::string& fileName);

void plot_the_dispersional_curves(
	std::map<double, std::complex<double>>& dispersionSet,
	const std::string& fileName);

void plot_the_imaginary_roots(
	std::map<double, std::vector<double>>& dispersionSet,
	const std::string& fileName);
