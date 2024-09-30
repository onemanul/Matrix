#include<iostream>
#include"Matrix.h"
#include <fstream>
using namespace std;

double pre_estimate_slae(Matrix& m) {
    if (m.getCols() != m.getRows() + 1) {
        cout << "    Error: wrong size for SLAE\n\n";
        return NULL;
    }
    else return m.getRows() * m.getRows() * m.getRows() * 2. / 3 + m.getRows() * m.getRows() * 3. / 2 - m.getRows() * 7. / 6;
}

void Zadanie_2() {
    int n;
    cout << "----------------Задание 2---------------\n Введите размер матрицы А: ";
    cin >> n;
    Matrix A(n), B;
    vctr x(0), b(0);
    for (int i = 0; i < n; ++i) {
        x.vec.push_back(1);
        for (int j = 0; j < n; ++j) {
            A(i, j) = 1. / (i + j + 1); // = i+j-1, если i,j начинаются с 1
        }
    }
    cout << "Матрица А:\n" << A << "Вектор b:\n" << (B = A * x);
    for (int i = 0; i < n; ++i) {
        b.vec.push_back(B(i, 0));
        A(i, i) += 1e-5; // помогает сильно повысить точность
    }
    cout << "Числа обусловленности матрицы A: " << conditioning_first(A) << "   " << conditioning_evkl(A) << "    " << conditioning_infinity(A) << "\n";
    A = A + b;
    cout << "Решение х:\n" << (x = SLAE(A)) << "---------------------------\n\n";
}
void SLAE_from_file() {
    ifstream File("variant5.txt");
    if (!File) cout << "Unable to open file variant.txt\n";
    Matrix Var(5, 6);
    vctr Ans;
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            File >> Var(i, j);
        }
    }
    for (size_t i = 0; i < 5; ++i) 
        File >> Var(i, 5);
    File.close();
    cout << "Предварительная оценка:" << pre_estimate_slae(Var) << "\n";
    cout << "Введённая матрица:\n" << Var << "Ответ:\n" << (Ans = SLAE(Var));
    cout << "Количество операций:" << Var.operation_in_slae << "\n";
    Var.changeCol(5);
    cout << "Проверка:\n" << (Var = Var*Ans);
}

int main() {
    setlocale(LC_ALL, "Russian");

    vector<double> V = { 1,2,-3 }, K = { 4,-1,6 }, B = { 1,3,3 };
    vctr v(V, 1), k(K, 1), x, b(B,0);
    vector<vector<double>> z = { {2,-3,1},{4,3,-3},{1,7,-4} }, s = { {9,0,2},{3,8,8},{0,2,1} };
    Matrix a(z), c(s), m, l(3), u(3), p;

    /*
    cout << "---Операции над векторами---\n";
    cout << v << k << "Нормы вектора v: " << norm_first(v) << " " << norm_second(v) << " " << norm_infinity(v) << "\n";
    cout <<  "Нормы вектора k: " << norm_first(k) << " " << norm_second(k) << " " << norm_infinity(k) << "\n";
    cout << "5*v : " << (t = 5 * v) << "k*(-3): " << (t = k * (-3)) << "v+k: " << (t = v + k) << "k-v: " << (t = k - v);
    k.trans();
    cout << "k transpos " << k << "v*k: " << v*k << "\n---------------------------\n\n";

    cout << "---Операции над матрицами---\n";
    cout << "Матрица а:\n" << a << "Нормы матрицы a: " << norm_first(a) << "   " << norm_evkl(a) << "    " << norm_infinity(a) << "\n";
    m = inverse(a);
    cout << "Матрица a(-1):\n" << m << "Нормы матрицы a(-1): " << norm_first(m) << "   " << norm_evkl(m) << "    " << norm_infinity(m) << "\n";
    cout << "Числа обусловленности матрицы а: " << conditioning_first(a) << "   " << conditioning_evkl(a) << "    " << conditioning_infinity(a) << "\n";
    cout << "Матрица c:\n" << c << "Нормы матрицы c: " << norm_first(c) << "   " << norm_evkl(c) << "    " << norm_infinity(c) << "\n";
    cout << "Числа обусловленности матрицы c: " << conditioning_first(c) << "   " << conditioning_evkl(c) << "    " << conditioning_infinity(c) << "\n";
    cout << "5*a:\n" << (m = 5 * a) << "c*(-3):\n" << (m = c * (-3)) << "a+c:\n" << (m = a+c) << "c-a:\n" << (m=c-a);
    cout << "v*a:\n" << (m = v*a) << "c*k:\n" << (m = c*k) << "a*c:\n" << (m = a * c);
    cout << "Единичная матрица рахмера 3:\n" << (m = Matrix::unit_matr(3)) << "Случайная матрица размера 3х5:\n" << (m = Matrix::rand_matr(3, 5));
    cout << "Обратная а:\n" << (m = inverse(a)) << "Транспонирование с:\n" << (m = transpos(c)) << "Определитель а:   " << determinant(a) << "\n";
    cout << "СЛАУ аx=b:\n" << (m = a + b) << "Решение СЛАУ с поиском max элемента по:\nCтроке: " << (t = SLAE(m, 0));
    cout << "Столбцу: " << (t = SLAE(m, 1)) << "Cтроке и столбцу: " << (t = SLAE(m, 2)) << "-------------------------- - \n\n";*/

    find_LU_matrix(a, l, u);
    find_P_matrix(a, p);
    cout << "Матрица A:\n" << a << "Матрица L: \n" << l << "Матрица U : \n" << u << "Матрица LU : \n" << (m = l * u) << "Матрица P:\n" << p;
    m = a + b;
    cout << (x = SLAE_LU(a, b)) << "Oper LU: " << a.operation_in_slae << "\n\n" << (x = SLAE_LUP(a, b)) << "Oper LUP: " << a.operation_in_slae << "\n\n";
    cout << (x = SLAE(m)) << "Oper slae: " << m.operation_in_slae << "\n\ndet and det_LU: " << determinant(a) << " " << determinant_with_LU(a);


    return 0;
}