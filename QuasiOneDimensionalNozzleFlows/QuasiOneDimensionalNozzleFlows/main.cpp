#include <fstream>
#include "QuasiOneDimensionalNozzleFlows.h"

void save(std::vector<double>& data, std::string path, std::string filename) {
	std::ofstream file("./" + path + "/" + filename + ".txt");
	for (int i = 0; i < data.size(); i++) {
		file << data[i];
		if (i < data.size() - 1) {
			file << " ";
		}
	}
	file.close();
}

void save(std::vector<std::vector<double>>& data, std::string path, std::string filename) {
	std::ofstream file("./" + path + "/" + filename + ".txt");
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data[i].size(); j++) {
			file << data[i][j];
			if (j < data[i].size() - 1) {
				file << " ";
			}
		}
		if (i < data.size() - 1) {
			file << "\n";
		}
	}
	file.close();
}

void save(Data& data, std::string path) {
	save(data.x, path, "x");
	save(data.t, path, "time");
	save(data.rho, path, "rho");
	save(data.V, path, "V");
	save(data.T, path, "T");
	save(data.p, path, "p");
	save(data.Ma, path, "Ma");
	save(data.m, path, "m");
	save(data.U1, path, "U1");
	save(data.U2, path, "U2");
	save(data.U3, path, "U3");
}

double A(double x) {
	return 1 + 2.2 * pow(x - 1.5, 2);
}

double rho1(double x) {
	if (x <= 0.5) {
		return 1;
	}
	else if (0.5 < x && x <= 1.5) {
		return 1 - 0.366 * (x - 0.5);
	}
	else {
		return 0.634 - 0.3879 * (x - 1.5);
	}
}

double V1(double x) {
	return 0.59 / (rho1(x) * A(x));
}

double T1(double x) {
	if (x <= 0.5) {
		return 1;
	}
	else if (0.5 < x && x <= 1.5) {
		return 1 - 0.167 * (x - 0.5);
	}
	else {
		return 0.833 - 0.3507 * (x - 1.5);
	}
}

double rho2(double x) {
	if (x <= 0.5) {
		return 1;
	}
	else if (0.5 < x && x <= 1.5) {
		return 1 - 0.366 * (x - 0.5);
	}
	else if (1.5 < x && x <= 2.1) {
		return 0.634 - 0.702 * (x - 1.5);
	}
	else {
		return 0.5892 - 0.10228 * (x - 2.1);
	}
}

double V2(double x) {
	return 0.59 / (rho2(x) * A(x));
}

double T2(double x) {
	if (x <= 0.5) {
		return 1;
	}
	else if (0.5 < x && x <= 1.5) {
		return 1 - 0.167 * (x - 0.5);
	}
	else if (1.5 < x && x <= 2.1) {
		return 0.833 - 0.4908 * (x - 1.5);
	}
	else {
		return 0.93968 - 0.0622 * (x - 2.1);
	}
}

void subsonicSupersonicIsentropicNozzleFlow() {
	QuasiOneDimensionalNozzleFlows QODNF(1.4, 3, A, rho1, V1, T1, 30, 0.5, 0, 1e-6);
	Data data;
	QODNF.MacCormack(data);
	save(data, "subsonicSupersonicIsentropicNozzleFlow");
}

void purelySubsonicIsentropicNozzleFlow() {
	QuasiOneDimensionalNozzleFlows QODNF(1.4, 3, 0.6784, A, rho2, V2, T2, 60, 0.5, 0.2, 1e-6);
	Data data;
	QODNF.MacCormack(data);
	save(data, "purelySubsonicIsentropicNozzleFlow");
}

int main() {
	subsonicSupersonicIsentropicNozzleFlow();
	purelySubsonicIsentropicNozzleFlow();
	return 0;
}
