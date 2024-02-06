#pragma once
#include <vector>

class Data {
public:
	std::vector<double> x;
	std::vector<double> t;
	std::vector<std::vector<double>> rho;
	std::vector<std::vector<double>> V;
	std::vector<std::vector<double>> T;
	std::vector<std::vector<double>> p;
	std::vector<std::vector<double>> Ma;
	std::vector<std::vector<double>> m;
	std::vector<std::vector<double>> U1;
	std::vector<std::vector<double>> U2;
	std::vector<std::vector<double>> U3;
};
