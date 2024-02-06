#pragma once
#include <iostream>
#include "Data.h"

class QuasiOneDimensionalNozzleFlows {
	double gamma;
	double pe;
	double (*A)(double);
	double (*rho)(double);
	double (*V)(double);
	double (*T)(double);
	int N;
	double dx;
	double C;
	double Cx;
	double tol;
	bool supersonicOutflow;
	void initialize(Data& data, std::vector<double>& rho, std::vector<double>& V, std::vector<double>& T, std::vector<double>& p, std::vector<double>& Ma, std::vector<double>& m, std::vector<double>& U1, std::vector<double>& U2, std::vector<double>& U3, std::vector<double>& dlogAdx_forward, std::vector<double>& dlogAdx_rearward) {
		data.x.resize(N + 1);
		for (int i = 0; i < N + 1; i++) {
			data.x[i] = i * dx;
		}
		for (int i = 0; i < N + 1; i++) {
			rho[i] = this->rho(data.x[i]);
		}
		for (int i = 0; i < N + 1; i++) {
			V[i] = this->V(data.x[i]);
		}
		for (int i = 0; i < N + 1; i++) {
			T[i] = this->T(data.x[i]);
		}
		for (int i = 0; i < N + 1; i++) {
			p[i] = rho[i] * T[i];
		}
		for (int i = 0; i < N + 1; i++) {
			Ma[i] = V[i] / sqrt(T[i]);
		}
		for (int i = 0; i < N + 1; i++) {
			m[i] = rho[i] * A(data.x[i]) * V[i];
		}
		for (int i = 0; i < N + 1; i++) {
			U1[i] = rho[i] * A(data.x[i]);
		}
		for (int i = 0; i < N + 1; i++) {
			U2[i] = rho[i] * A(data.x[i]) * V[i];
		}
		for (int i = 0; i < N + 1; i++) {
			U3[i] = rho[i] * (T[i] / (gamma - 1) + gamma / 2 * pow(V[i], 2)) * A(data.x[i]);
		}
		for (int i = 0; i < N; i++) {
			dlogAdx_forward[i] = (log(A(data.x[i + 1])) - log(A(data.x[i]))) / dx;
		}
		for (int i = 1; i < N + 1; i++) {
			dlogAdx_rearward[i] = (log(A(data.x[i])) - log(A(data.x[i - 1]))) / dx;
		}
		data.t.push_back(0);
		data.rho.push_back(rho);
		data.V.push_back(V);
		data.T.push_back(T);
		data.p.push_back(p);
		data.Ma.push_back(Ma);
		data.m.push_back(m);
		data.U1.push_back(U1);
		data.U2.push_back(U2);
		data.U3.push_back(U3);
	}
	double timeStep(std::vector<double>& V, std::vector<double>& T) {
		std::vector<double> dt(N + 1);
		for (int i = 1; i < N; i++) {
			dt[i] = C * dx / (V[i] + sqrt(T[i]));
		}
		return *std::min_element(dt.begin() + 1, dt.end() - 1);
	}
	double residual(std::vector<std::vector<double>>& U) {
		std::vector<double> res(N + 1);
		for (int i = 0; i < N + 1; i++) {
			res[i] = abs(U[U.size() - 1][i] - U[U.size() - 2][i]);
		}
		return *std::max_element(res.begin(), res.end());
	}
	double residual(Data& data) {
		double res1 = residual(data.U1);
		double res2 = residual(data.U2);
		double res3 = residual(data.U3);
		return std::max(std::max(res1, res2), res3);
	}
public:
	QuasiOneDimensionalNozzleFlows(double gamma, double L, double (*A)(double), double (*rho)(double), double (*V)(double), double (*T)(double), int N, double C, double Cx, double tol) :gamma(gamma), A(A), rho(rho), V(V), T(T), N(N), dx(L / N), C(C), Cx(Cx), tol(tol), supersonicOutflow(true) {}
	QuasiOneDimensionalNozzleFlows(double gamma, double L, double pe, double (*A)(double), double (*rho)(double), double (*V)(double), double (*T)(double), int N, double C, double Cx, double tol) :gamma(gamma), pe(pe), A(A), rho(rho), V(V), T(T), N(N), dx(L / N), C(C), Cx(Cx), tol(tol), supersonicOutflow(false) {}
	void MacCormack(Data& data) {
		double res;
		std::vector<double> rho(N + 1), V(N + 1), T(N + 1), p(N + 1), Ma(N + 1), m(N + 1), U1(N + 1), U2(N + 1), U3(N + 1), dlogAdx_forward(N + 1), dlogAdx_rearward(N + 1), F1(N + 1), F2(N + 1), F3(N + 1), J2(N + 1), dU1dt(N + 1), dU2dt(N + 1), dU3dt(N + 1), U1_bar(N + 1), U2_bar(N + 1), U3_bar(N + 1), F1_bar(N + 1), F2_bar(N + 1), F3_bar(N + 1), J2_bar(N + 1), dU1dt_bar(N + 1), dU2dt_bar(N + 1), dU3dt_bar(N + 1), dU1dt_av(N + 1), dU2dt_av(N + 1), dU3dt_av(N + 1), S1(N + 1), S2(N + 1), S3(N + 1), S1_bar(N + 1), S2_bar(N + 1), S3_bar(N + 1), rho_bar(N + 1), V_bar(N + 1), T_bar(N + 1), p_bar(N + 1);
		initialize(data, rho, V, T, p, Ma, m, U1, U2, U3, dlogAdx_forward, dlogAdx_rearward);
		do {
			double dt = timeStep(V, T);
			for (int i = 0; i < N + 1; i++) {
				F1[i] = U2[i];
			}
			for (int i = 0; i < N + 1; i++) {
				F2[i] = pow(U2[i], 2) / U1[i] + (gamma - 1) / gamma * (U3[i] - gamma / 2 * pow(U2[i], 2) / U1[i]);
			}
			for (int i = 0; i < N + 1; i++) {
				F3[i] = gamma * U2[i] * U3[i] / U1[i] - gamma * (gamma - 1) / 2 * pow(U2[i], 3) / pow(U1[i], 2);
			}
			for (int i = 0; i < N; i++) {
				J2[i] = (gamma - 1) / gamma * (U3[i] - gamma / 2 * pow(U2[i], 2) / U1[i]) * dlogAdx_forward[i];
			}
			for (int i = 0; i < N; i++) {
				dU1dt[i] = -(F1[i + 1] - F1[i]) / dx;
			}
			for (int i = 0; i < N; i++) {
				dU2dt[i] = -(F2[i + 1] - F2[i]) / dx + J2[i];
			}
			for (int i = 0; i < N; i++) {
				dU3dt[i] = -(F3[i + 1] - F3[i]) / dx;
			}
			for (int i = 1; i < N; i++) {
				S1[i] = Cx * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (U1[i + 1] - 2 * U1[i] + U1[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				S2[i] = Cx * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (U2[i + 1] - 2 * U2[i] + U2[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				S3[i] = Cx * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (U3[i + 1] - 2 * U3[i] + U3[i - 1]);
			}
			for (int i = 0; i < N; i++) {
				U1_bar[i] = U1[i] + dU1dt[i] * dt + S1[i];
			}
			for (int i = 0; i < N; i++) {
				U2_bar[i] = U2[i] + dU2dt[i] * dt + S2[i];
			}
			for (int i = 0; i < N; i++) {
				U3_bar[i] = U3[i] + dU3dt[i] * dt + S3[i];
			}
			for (int i = 0; i < N; i++) {
				rho_bar[i] = U1_bar[i] / A(data.x[i]);
			}
			for (int i = 0; i < N; i++) {
				V_bar[i] = U2_bar[i] / U1_bar[i];
			}
			for (int i = 0; i < N; i++) {
				T_bar[i] = (gamma - 1) * (U3_bar[i] / U1_bar[i] - gamma / 2 * pow(V_bar[i], 2));
			}
			for (int i = 0; i < N; i++) {
				p_bar[i] = rho_bar[i] * T_bar[i];
			}
			if (!supersonicOutflow) {
				p_bar[N] = pe;
			}
			for (int i = 0; i < N; i++) {
				F1_bar[i] = U2_bar[i];
			}
			for (int i = 0; i < N; i++) {
				F2_bar[i] = pow(U2_bar[i], 2) / U1_bar[i] + (gamma - 1) / gamma * (U3_bar[i] - gamma / 2 * pow(U2_bar[i], 2) / U1_bar[i]);
			}
			for (int i = 0; i < N; i++) {
				F3_bar[i] = gamma * U2_bar[i] * U3_bar[i] / U1_bar[i] - gamma * (gamma - 1) / 2 * pow(U2_bar[i], 3) / pow(U1_bar[i], 2);
			}
			for (int i = 1; i < N; i++) {
				J2_bar[i] = (gamma - 1) / gamma * (U3_bar[i] - gamma / 2 * pow(U2_bar[i], 2) / U1_bar[i]) * dlogAdx_rearward[i];
			}
			for (int i = 1; i < N; i++) {
				dU1dt_bar[i] = -(F1_bar[i] - F1_bar[i - 1]) / dx;
			}
			for (int i = 1; i < N; i++) {
				dU2dt_bar[i] = -(F2_bar[i] - F2_bar[i - 1]) / dx + J2_bar[i];
			}
			for (int i = 1; i < N; i++) {
				dU3dt_bar[i] = -(F3_bar[i] - F3_bar[i - 1]) / dx;
			}
			for (int i = 1; i < N; i++) {
				dU1dt_av[i] = (dU1dt[i] + dU1dt_bar[i]) / 2;
			}
			for (int i = 1; i < N; i++) {
				dU2dt_av[i] = (dU2dt[i] + dU2dt_bar[i]) / 2;
			}
			for (int i = 1; i < N; i++) {
				dU3dt_av[i] = (dU3dt[i] + dU3dt_bar[i]) / 2;
			}
			for (int i = 1; i < N; i++) {
				S1_bar[i] = Cx * abs(p_bar[i + 1] - 2 * p_bar[i] + p_bar[i - 1]) / (p_bar[i + 1] + 2 * p_bar[i] + p_bar[i - 1]) * (U1_bar[i + 1] - 2 * U1_bar[i] + U1_bar[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				S2_bar[i] = Cx * abs(p_bar[i + 1] - 2 * p_bar[i] + p_bar[i - 1]) / (p_bar[i + 1] + 2 * p_bar[i] + p_bar[i - 1]) * (U2_bar[i + 1] - 2 * U2_bar[i] + U2_bar[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				S3_bar[i] = Cx * abs(p_bar[i + 1] - 2 * p_bar[i] + p_bar[i - 1]) / (p_bar[i + 1] + 2 * p_bar[i] + p_bar[i - 1]) * (U3_bar[i + 1] - 2 * U3_bar[i] + U3_bar[i - 1]);
			}
			for (int i = 1; i < N; i++) {
				U1[i] += dU1dt_av[i] * dt + S1_bar[i];
			}
			for (int i = 1; i < N; i++) {
				U2[i] += dU2dt_av[i] * dt + S2_bar[i];
			}
			for (int i = 1; i < N; i++) {
				U3[i] += dU3dt_av[i] * dt + S3_bar[i];
			}
			U1[0] = A(data.x[0]);
			U2[0] = 2 * U2[1] - U2[2];
			U3[0] = U1[0] * (1 / (gamma - 1) + gamma / 2 * pow(U2[0] / U1[0], 2));
			if (supersonicOutflow) {
				U1[N] = 2 * U1[N - 1] - U1[N - 2];
				U2[N] = 2 * U2[N - 1] - U2[N - 2];
				U3[N] = 2 * U3[N - 1] - U3[N - 2];
			}
			else {
				U1[N] = 2 * U1[N - 1] - U1[N - 2];
				U2[N] = 2 * U2[N - 1] - U2[N - 2];
				U3[N] = pe * A(data.x[N]) / (gamma - 1) + gamma / 2 * pow(U2[N], 2) / U1[N];
			}
			for (int i = 0; i < N + 1; i++) {
				rho[i] = U1[i] / A(data.x[i]);
			}
			for (int i = 0; i < N + 1; i++) {
				V[i] = U2[i] / U1[i];
			}
			for (int i = 0; i < N + 1; i++) {
				T[i] = (gamma - 1) * (U3[i] / U1[i] - gamma / 2 * pow(V[i], 2));
			}
			for (int i = 0; i < N + 1; i++) {
				p[i] = rho[i] * T[i];
			}
			for (int i = 0; i < N + 1; i++) {
				Ma[i] = V[i] / sqrt(T[i]);
			}
			for (int i = 0; i < N + 1; i++) {
				m[i] = rho[i] * A(data.x[i]) * V[i];
			}
			data.t.push_back(data.t[data.t.size() - 1] + dt);
			data.rho.push_back(rho);
			data.V.push_back(V);
			data.T.push_back(T);
			data.p.push_back(p);
			data.Ma.push_back(Ma);
			data.m.push_back(m);
			data.U1.push_back(U1);
			data.U2.push_back(U2);
			data.U3.push_back(U3);
			res = residual(data);
			printf("time=%.6f res=%.6f\n", data.t[data.t.size() - 1], res);
		} while (res > tol);
	}
};
