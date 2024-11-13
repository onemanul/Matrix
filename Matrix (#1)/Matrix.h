#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include"vctr.h"
#include<iomanip>

class Matrix {
	public:
		Matrix() : data(1, vector<double>(1, 0.0)) {}
		Matrix(size_t n) : data(n, vector<double>(n, 0.0)) {}
		Matrix(size_t rows, size_t cols) : data(rows, vector<double>(cols, 0.0)) {}
		Matrix(vector<vector<double>> dt) : data(dt) {}
		Matrix(vctr v);
		static Matrix unit_matr(size_t n);						// tmp = Matrix::unit_matr(3)
		static Matrix rand_matr(size_t rows, size_t cols);		// tmp = Matrix::rand_matr(3, 5)
		
		size_t getRows() const { return data.size(); }
		size_t getCols() const { return data[0].size(); }
		double& operator()(size_t row, size_t col) { return data[row][col]; }
		vctr to_vector();										// перевод матрицы в вектор (строку или столбец), если это возможно

		void changeRow(size_t r);								// количество строк = r. Заполняет нулями, если матрица увеличилась
		void changeCol(size_t с);								// количество столбцов = с. Заполняет нулями, если матрица увеличилась
		void swap_rows(size_t a, size_t b);						// меняет местами строки a и b

		bool is_symmetry();										// проверка матрицы на симметричность
		bool is_positive_definite();							// проверка матрицы на положительную определённость
		bool is_tridiagonal();									// проверка матрицы на тридиагональность

		double norm_first();												// вычисление абсолютной (1-нормы)
		double norm_infinity();												// вычисление  максимальной (inf-нормы)
		double norm_evkl();													// вычисление евклидовой нормы

		double determinant();												// определитель
		double determinant_with_LU();										// определитель с помощью LU-разложения

		bool find_LU_matrix(Matrix& L, Matrix& U);							// нахождение LU разложения матрицы А
		bool find_QR_matrix(Matrix& Q, Matrix& R);							// нахождение QR разложения матрицы А
		bool find_P_matrix(Matrix& P);										// нахождение матрицы перестановок Р для матрицы А
		bool find_L_Lt_matrix(Matrix& L, Matrix& L_t);						// нахождение матриц L и L_t для метода Холецкого


		friend Matrix inverse(Matrix& m);									// обратная матрица
		friend Matrix transpos(Matrix& m);									// транспозиция
		friend vctr SLAE(Matrix& mtr);										// решение СЛАУ, поиск главного элемента по строкам
		friend vctr SLAE(Matrix& mtr, int metod);							// решение СЛАУ с выбором метода поиска главного элемента: 0 - по строкам, 1 - по столбцам, 2 - по строкам и столбцам
		friend vctr SLAE_LU(Matrix& A, vctr& b);							// решение СЛАУ с помощью LU разложения
		friend vctr SLAE_LU(Matrix& L, Matrix& U, vctr& b);					// решение СЛАУ с помощью LU разложения
		friend vctr SLAE_LUP(Matrix& A, vctr& b);							// решение СЛАУ с помощью LUР разложения
		friend vctr SLAE_Cholesky(Matrix& A, vctr& b);						// решение СЛАУ методом Холецкого
		friend vctr SLAE_Cholesky(Matrix& L, Matrix& L_t, vctr& b);			// решение СЛАУ методом Холецкого
		friend vctr SLAE_Tomas(Matrix& A, vctr& b);							// решение СЛАУ для тридиагональных матриц методом Томаса (прогонки)
		friend vctr SLAE_QR(Matrix& A, vctr& b);							// решение СЛАУ с помощью QR разложения
		friend vctr SLAE_QR(Matrix& Q, Matrix& R, vctr& b);					// решение СЛАУ с помощью QR разложения
		friend vctr SLAE_Iteration(Matrix& A, vctr& b, double accuracy, vctr& initial_approximation);		// решение СЛАУ методом простой итерации
		friend vctr SLAE_Seidel(Matrix& A, vctr& b, double accuracy, vctr& initial_approximation);			// решение СЛАУ методом Зейделя

		friend double power_law_method(Matrix& A, double accuracy, vctr& x_next);							// нахождения собственного значения с максимальной абсолютной величиной; условие: |λ| < 1
		friend double power_law_method_with_normalization(Matrix& A, double accuracy, vctr& x_next);		// нахождения собственного значения с максимальной абсолютной величиной; для любых λ


		//------------------------------Переопределение операторов------------------------------
		friend Matrix operator*(double A, Matrix& m);
		friend Matrix operator*(Matrix& m, double A);
		friend Matrix operator*(vctr& v, Matrix& m);
		friend Matrix operator*(Matrix& m, vctr& v);
		friend Matrix operator+(Matrix& l, Matrix& r);
		friend Matrix operator-(Matrix& l, Matrix& r);
		friend Matrix operator*(Matrix& l, Matrix& r);
		friend bool operator==(Matrix& left, Matrix& right);
		friend ostream& operator<<(ostream& out, Matrix& m);
		friend Matrix operator+(Matrix& m, vctr& v);						// добавляет к матрице строку или столбец вектора

	private:
		vector<vector<double>> data;
};

#endif // !MATRIX_H