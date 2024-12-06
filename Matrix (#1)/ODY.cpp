#include "ODY.h"
#include <cmath>

ODY::ODY(double left_border, double right_border, int number_of_nodes) {
	a = left_border;
	b = right_border;
	n = number_of_nodes;
	h = (b - a) / n;
	for (int i = 0; i < n + 1; ++i) {
		x.push_back(a + i * h);
	}
	double a_shift = a - h / 2;
	for (int i = 0; i < n + 2; ++i) {
		x_shift.push_back(a_shift + i * h);
	}
}


Matrix ODY::getMatrix(int order, double alpha_1, double alpha_2, double beta_1, double beta_2) {
	int m = (order == 1) ? n : n + 1;
	vector<double> v = (order == 1) ? x : x_shift;
	Matrix a(m+1, m+1);
	for (int i = 1; i < m; ++i) {
		a(i, i) = B(v[i]);
		a(i, i - 1) = A(v[i]);
		a(i, i + 1) = C(v[i]);
	}
	if (order == 1) {
		a(0, 0) = alpha_1 + alpha_2 / h;
		a(0, 1) = -alpha_2 / h;
		a(m, m - 1) = -beta_2 / h;
		a(m, m) = beta_1 + beta_2 / h;
	}
	else {
		a(0, 0) = alpha_1 / 2 + alpha_2 / h;
		a(0, 1) = alpha_1 / 2 - alpha_2 / h;
		a(m, m - 1) = beta_1 / 2 - beta_2 / h;
		a(m, m) = beta_1 / 2 + beta_2 / h;
	}
	return a;
}

vctr ODY::getVctr(int order, double alpha, double beta) {
	int m = (order == 1) ? n : n + 1;
	vector<double> v = (order == 1) ? x : x_shift;
	vctr b(0);
	b.vec.push_back(alpha);
	for (int i = 1; i < m; ++i) {
		b.vec.push_back(G(v[i]));
	}
	b.vec.push_back(beta);
	return b;
}