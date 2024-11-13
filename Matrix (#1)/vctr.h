#ifndef VCTR_H
#define VCTR_H

#include<iostream>
#include<vector>
using namespace std;

class vctr {
	public:
		vector <double> vec;
		bool as_row;
		vctr() : as_row(true) {}
		vctr(bool row) : as_row(row) {}
		vctr(vector <double> v) : vec(v), as_row(true) {}
		vctr(vector <double> v, bool row) : vec(v), as_row(row) {}
		void trans() { as_row = !as_row; }

		double norm_first();
		double norm_second();
		double norm_infinity();

		friend vctr transpos(vctr& v);
		friend vctr operator*(double A, vctr& v);
		friend vctr operator*(vctr& v, double A);
		friend vctr operator+(vctr& l, vctr& r);
		friend vctr operator-(vctr& l, vctr& r);
		friend double operator*(vctr& l, vctr& r);
		friend bool operator==(const vctr& left, const vctr& right);
		friend ostream& operator<<(ostream& out, vctr& v);
		friend vctr operator/(vctr& v, double A);
};

#endif // VCTR_H