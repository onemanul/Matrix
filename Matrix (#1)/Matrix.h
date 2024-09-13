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
		static Matrix unit_matr(size_t n);						// tmp = Matrix::unit_matr(3)
		static Matrix rand_matr(size_t rows, size_t cols);		// tmp = Matrix::rand_matr(3, 5)
		
		size_t getRows() const { return data.size(); }
		size_t getCols() const { return data[0].size(); }
		double& operator()(size_t row, size_t col) { return data[row][col]; }

		friend double norm_first(Matrix& m);
		friend double norm_infinity(Matrix& m);
		friend double norm_evkl(Matrix& m);
		friend double conditioning_first(Matrix& m);
		friend double conditioning_infinity(Matrix& m);
		friend double conditioning_evkl(Matrix& m);

		friend Matrix inverse(Matrix& m);
		friend double determinant(Matrix& m);
		friend vctr SLAE(Matrix& mtr);
		friend Matrix transpos(Matrix& m);
		

		friend Matrix operator*(double A, Matrix& m);
		friend Matrix operator*(Matrix& m, double A);
		friend Matrix operator*(vctr& v, Matrix& m);
		friend Matrix operator*(Matrix& m, vctr& v);
		friend Matrix operator+(Matrix& l, Matrix& r);
		friend Matrix operator-(Matrix& l, Matrix& r);
		friend Matrix operator*(Matrix& l, Matrix& r);
		friend bool operator==(Matrix& left, Matrix& right);
		friend ostream& operator<<(ostream& out, Matrix& m);

	private:
		vector<vector<double>> data;
};

#endif // !MATRIX_H