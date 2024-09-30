#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include"vctr.h"
#include<iomanip>

class Matrix {
	public:
		size_t operation_in_slae = 0;

		Matrix() {}
		Matrix(size_t n) : data(n, vector<double>(n, 0.0)) {}
		Matrix(size_t rows, size_t cols) : data(rows, vector<double>(cols, 0.0)) {}
		Matrix(vector<vector<double>> dt) : data(dt) {}
		Matrix(vctr v);
		static Matrix unit_matr(size_t n);						// tmp = Matrix::unit_matr(3)
		static Matrix rand_matr(size_t rows, size_t cols);		// tmp = Matrix::rand_matr(3, 5)
		
		size_t getRows() const { return data.size(); }
		size_t getCols() const { return data[0].size(); }
		void changeRow(size_t r);								//  ���������� ����� = r. ��������� ������, ���� ������� �����������
		void changeCol(size_t �);								//  ���������� �������� = �. ��������� ������, ���� ������� �����������
		void swap_rows(size_t a, size_t b);
		double& operator()(size_t row, size_t col) { return data[row][col]; }

		friend double norm_first(Matrix& m);
		friend double norm_infinity(Matrix& m);
		friend double norm_evkl(Matrix& m);
		friend double conditioning_first(Matrix& m);
		friend double conditioning_infinity(Matrix& m);
		friend double conditioning_evkl(Matrix& m);

		friend Matrix inverse(Matrix& m);						// �������� �������
		friend Matrix transpos(Matrix& m);						// ������������
		friend double determinant(Matrix& m);					// ������������
		friend double determinant_with_LU(Matrix& m);			// ������������ � ������� LU
		friend vctr SLAE(Matrix& mtr);							// ������� ����, ����� �������� �������� �� �������
		friend vctr SLAE(Matrix& mtr, int metod);				// ������� ���� � ������� ������ ������ �������� ��������: 0 - �� �������, 1 - �� ��������, 2 - �� ������� � ��������
		friend bool find_LU_matrix (Matrix& A, Matrix& L, Matrix& U);
		friend bool find_P_matrix(Matrix& A, Matrix& P);
		friend vctr SLAE_LU(Matrix& mtr, vctr& b);
		friend vctr SLAE_LUP(Matrix& mtr, vctr& b);

		friend Matrix operator*(double A, Matrix& m);
		friend Matrix operator*(Matrix& m, double A);
		friend Matrix operator*(vctr& v, Matrix& m);
		friend Matrix operator*(Matrix& m, vctr& v);
		friend Matrix operator+(Matrix& l, Matrix& r);
		friend Matrix operator-(Matrix& l, Matrix& r);
		friend Matrix operator*(Matrix& l, Matrix& r);
		friend bool operator==(Matrix& left, Matrix& right);
		friend ostream& operator<<(ostream& out, Matrix& m);
		friend Matrix operator+(Matrix& m, vctr& v);			// ��������� � ������� ������ ��� ������� �������

	private:
		vector<vector<double>> data;
};

#endif // !MATRIX_H