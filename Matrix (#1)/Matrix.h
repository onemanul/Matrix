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
		void changeRow(size_t r);								//  количество строк = r. Заполняет нулями, если матрица увеличилась
		void changeCol(size_t с);								//  количество столбцов = с. Заполняет нулями, если матрица увеличилась
		void swap_rows(size_t a, size_t b);
		double& operator()(size_t row, size_t col) { return data[row][col]; }

		friend double norm_first(Matrix& m);
		friend double norm_infinity(Matrix& m);
		friend double norm_evkl(Matrix& m);
		friend double conditioning_first(Matrix& m);
		friend double conditioning_infinity(Matrix& m);
		friend double conditioning_evkl(Matrix& m);

		friend Matrix inverse(Matrix& m);						// обратная матрица
		friend Matrix transpos(Matrix& m);						// транспозиция
		friend double determinant(Matrix& m);					// определитель
		friend double determinant_with_LU(Matrix& m);			// определитель с помощью LU
		friend vctr SLAE(Matrix& mtr);							// решение СЛАУ, поиск главного элемента по строкам
		friend vctr SLAE(Matrix& mtr, int metod);				// решение СЛАУ с выбором метода поиска главного элемента: 0 - по строкам, 1 - по столбцам, 2 - по строкам и столбцам
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
		friend Matrix operator+(Matrix& m, vctr& v);			// добавляет к матрице строку или столбец вектора

	private:
		vector<vector<double>> data;
};

#endif // !MATRIX_H