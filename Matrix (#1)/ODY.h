#pragma once
#ifndef ODY_H
#define ODY_H

#include"Matrix.h"
#include<iomanip>
/*			Разностный метод решения краевой задачи для ОДУ второго порядка		
Решение возможно с порядком O(h) и O(h^2) (order = 1 и order = 2 соответственно)
*/

class ODY {
	public:
		ODY(double left_border, double right_border, int number_of_nodes);
		vector<double> x;
		vector<double> x_shift;

		//double p(double x) { return 2/(2+x); }
		//double q(double x) { return (1+x)/2 + 2/pow(x+2,2); }
		//double r(double x) { return cos(x/2); }
		//double f(double x) { return 1 + x/2; }

		//double A(double x) { return -p(x)/pow(h, 2) - q(x)/(2*h); }
		//double B(double x) { return (2*p(x))/pow(h, 2) + r(x); }
		//double C(double x) { return -p(x)/pow(h, 2) + q(x)/(2*h); }
		//double G(double x) { return f(x); }
		// -1, -1, 0, 1

		//double p(double x) { return (2 + x) / (3 + x); }
		//double q(double x) { return (2+x)/pow(3+x,2) - 1/(3+x); }
		//double r(double x) { return 1 + sin(x); }
		//double f(double x) { return 1 - x; }
		// 0, -1, 1, 1
		double p(double x) { return 1 / (2 + x); }
		double q(double x) { return 1 / pow(2 + x, 2); }
		double r(double x) { return cos(x); }
		double f(double x) { return 1 + x; }

		double A(double x) { return -p(x) / pow(h, 2) - q(x) / (2 * h); }
		double B(double x) { return (2*p(x))/pow(h, 2) + r(x); }
		double C(double x) { return -p(x)/pow(h, 2) + q(x)/(2*h); }
		//double A(double x) { return -p(x - h/2)/pow(h, 2); }
		//double B(double x) { return p(x - h/2)/pow(h, 2) + p(x + h/2)/pow(h, 2) + r(x); }
		//double C(double x) { return -p(x + h/2)/pow(h,2); }
		double G(double x) { return f(x); }

		Matrix getMatrix(int order, double alpha_1, double alpha_2, double beta_1, double beta_2);
		vctr getVctr(int order, double alpha, double beta);

	private:
		double a;
		double b;
		int n;
		double h;
};

#endif // !ODY_H