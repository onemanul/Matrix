#include "vctr.h"

double norm_first(vctr& v) {
	double s = 0;
	for (int i = 0; i < v.vec.size(); ++i)
		s += fabs(v.vec[i]);
	return s;
};
double norm_second(vctr& v) {
	double s = 0;
	for (int i = 0; i < v.vec.size(); ++i)
		s += v.vec[i] * v.vec[i];
	return sqrt(s);
}
double norm_infinity(vctr& v) {
	double m = 0;
	for (int i = 0; i < v.vec.size(); ++i)
		if (fabs(v.vec[i]) > m) m = fabs(v.vec[i]);
	return m;
}

vctr transpos(vctr& v) {
	vctr tmp = v;
	tmp.as_row = !tmp.as_row;
	return tmp;
}
vctr operator*(double A, vctr& v) {
	vctr tmp = v;
	for (int i = 0; i < tmp.vec.size(); ++i)
		tmp.vec[i] *= A;
	return tmp;
}
vctr operator*(vctr& v, double A) {
	return (A * v);
}
vctr operator+(vctr& l, vctr& r) {
	if (l.vec.size() != r.vec.size() || l.as_row != r.as_row) return NULL;
	vctr a = l, b = r, c(a.as_row);
	for (int i = 0; i < a.vec.size(); ++i)
		c.vec.push_back(a.vec[i] + b.vec[i]);
	return c;
}
double operator*(vctr& l, vctr& r) {
	if (l.vec.size() == r.vec.size() && r.vec.size() == 1) {
		return l.vec[0] * r.vec[0];
	}
	if (l.vec.size() != r.vec.size() || !l.as_row || r.as_row) {
		return NULL;
	}
	double s = 0;
	for (int i = 0; i < l.vec.size(); ++i)
		s += l.vec[i] * r.vec[i];
	return s;
}
bool operator==(const vctr& left, const vctr& right) {
	if (left.vec.size() != right.vec.size() || left.as_row != right.as_row) return false;
	for (int i = 0; i < left.vec.size(); ++i) {
		if (!(left.vec[i] == right.vec[i])) return false;
	}
	return true;
}
ostream& operator<<(ostream& out, vctr& v) {
	for (int i = 0; i < v.vec.size(); ++i)
		out << v.vec[i];
	if (v.as_row) out << " row \n";
	else out << " column \n";
	return out;
}