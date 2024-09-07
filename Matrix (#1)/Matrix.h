#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include"vctr.h"

class Matrix {
	public:
		Matrix() {}
		Matrix(size_t n) : data(n, vector<double>(n, 0.0)) {}
		Matrix(size_t rows, size_t cols) : data(rows, vector<double>(cols, 0.0)) {}
		Matrix(vector<vector<double>> dt) : data(dt) {}
		Matrix(vctr v);
		
		size_t getRows() const { return data.size(); }
		size_t getCols() const { return data[0].size(); }
		double& operator()(size_t row, size_t col) { return data[row][col]; }

		friend double norm_first(Matrix& m);
		friend double norm_infinity(Matrix& m);
		friend double norm_evkl(Matrix& m);
		//friend double conditioning_number(Matrix& m);

		friend Matrix inverse(Matrix& m);
		friend Matrix transpos(Matrix& m);
		friend Matrix unit_matr(size_t n);
		friend Matrix rand_matr(size_t rows, size_t cols);

		friend Matrix operator*(double A, Matrix& m);
		friend Matrix operator*(Matrix& m, double A);
		friend Matrix operator*(vctr& v, Matrix& m);
		friend Matrix operator*(Matrix& m, vctr& v);
		friend Matrix operator+(Matrix& l, Matrix& r);
		friend Matrix operator*(Matrix& l, Matrix& r);
		friend bool operator==(const Matrix& left, const Matrix& right);
		friend ostream& operator<<(ostream& out, Matrix& m);

	private:
		vector<vector<double>> data;
};

#endif // !MATRIX_H