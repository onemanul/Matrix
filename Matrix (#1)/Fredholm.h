#pragma once
#ifndef FREDHOLM_H
#define FREDHOLM_H
#include "Matrix.h"

/*
Решение уравнения Фредгольма вида:		u(x) - \lambda \int_a^b K(x,s)u(s)ds = f(x)			с использованием составной формулы 3/8. 
В данном примере: для формулы			u(x) - \lambda \int_a^b sinh(xs^2)u(s)ds = x + 0.1
*/

class Fredholm {
	public:
		Fredholm(double left_border, double right_border, double Lambda) : a(left_border), b(right_border), lambda(Lambda) {}
		double f(double x) { return x + 0.1; }
		double K(double x, double s) { return sinh(x * s * s); }

		Matrix getMatrix_A(int n);
		vctr getVector_b(int n);
		vctr getVector_u(double epsilon);

	private:
		double a;
		double b;
		double lambda;
};

#endif // !FREDHOLM_H