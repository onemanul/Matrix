#include<iostream>
#include"Matrix.h"
#include <fstream>
#include <random>
using namespace std;
/*
double pre_estimate_slae(Matrix& m) {
    if (m.getCols() != m.getRows() + 1) {
        cout << "    Error: wrong size for SLAE\n\n";
        return NULL;
    }
    else return m.getRows() * m.getRows() * m.getRows() * 2. / 3 + m.getRows() * m.getRows() * 3. / 2 - m.getRows() * 7. / 6;
}

void Zadanie_1() {
    vector<double> V = { 1,2,-3 }, K = { 4,-1,6 }, B = { 1,3,3 };
    vctr v(V, 1), k(K, 1), t, b(B, 0);
    vector<vector<double>> A = { {2,-3,1},{4,3,-3},{1,7,-4} }, C = { {9,0,2},{3,8,8},{0,2,1} };
    Matrix a(A), c(C), m;

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
    cout << "Столбцу: " << (t = SLAE(m, 1)) << "Cтроке и столбцу: " << (t = SLAE(m, 2)) << "-------------------------- - \n\n";
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
    Var.changeCol(5);
    cout << "Проверка:\n" << (Var = Var*Ans);
}

vctr generate_column(size_t size, double a, double b) {
    vector<double> x(size);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(a, b);
    for (int i = 0; i < size; i++) 
        x[i] = dis(gen);
    vctr v(x, 0);
    return v;
}
void Zadanie_3() {
    vector<double> B = { 1,3,3 };
    vctr x, b(B, 0); 
    vector<vector<double>> A = { {2,-3,1},{4,3,-3},{1,7,-4} };
    Matrix a(A), l, u, m;
    find_LU_matrix(a, l, u);
    cout << "Исходная система:\n" << (m = a+b) << "Решение:\n" << (x = SLAE_LU(a, b)) << "Матрица L: \n" << l << "Матрица U : \n" << u;
    for (int i = 0; i < 5; ++i) {
        b = generate_column(a.getRows(), -10, 10);
        m = a + b;
        cout << "Столбец b:\n" << b << "Решение:LU\n" << (x = SLAE_LU(a, b)) << "\n\n";
    }
}

void Zadanie_4() {
    vector<double>  V = { 18,2,14,5,-21 }, B_TOM = { 1,0,0,0 }, B_TOM2 = { 4.5,8.45,10.63,38.123,52.12365,69.4125,26.52,7.002,22.01203,9.1,3.999,-1.156 };
    vctr x, v(V, 0), b_tom(B_TOM, 0), b_tom2(B_TOM2, 0);
    vector<vector<double>> C = { {4,1,2,1,1},{1,5,3,2,2},{2,3,6,3,3},{1,2,3,7,4},{1,2,3,4,8} },
        M = { {2,1,1,1,1},{1,1,0,0,0},{1,0,1,0,0},{1,0,0,1,0},{1,0,0,0,1} }, TOM = { {2,-1,0,0},{-1,2,-1,0},{0,-1,2,-1}, {0,0,-1,2} },
        TOM2 = { {2,1,0,0,0,0,0,0,0,0,0,0},
                {-1,3,1,0,0,0,0,0,0,0,0,0},
                {0,-3,4,1,0,0,0,0,0,0,0,0},
                {0,0,1,5,3,0,0,0,0,0,0,0},
                {0,0,0,4,6,1,0,0,0,0,0,0},
                {0,0,0,0,3,7,2,0,0,0,0,0},
                {0,0,0,0,0,-2,5,2,0,0,0,0},
                {0,0,0,0,0,0,-1,2,1,0,0,0},
                {0,0,0,0,0,0,0,1,3,1,0,0},
                {0,0,0,0,0,0,0,0,-3,6,3,0},
                {0,0,0,0,0,0,0,0,0,-2,4,1},
                {0,0,0,0,0,0,0,0,0,0,-5,9} };
    Matrix c(C), m(M), l, t, tom(TOM), tom2(TOM2), tmp;

    cout << "Positive definite: c = " << c.is_positive_definite() << ", m = " << m.is_positive_definite() << "\nDeterminant: c = " << determinant(c) << ", m = " << determinant(m) << "\n";
    Cholesky_metod(c, l, t);
    cout << "Matrix c: \n" << c << l << t << (tmp = (tmp = l * t) + v) << (x = SLAE_Cholesky(c, v)) << "Операций Холецкий: " << c.operation_in_slae;
    c.operation_in_slae = 0;
    cout << "\n\nРешение LU: " << (x = SLAE_LU(c, v)) << "Операций LU: " << c.operation_in_slae;
    cout << "\n\nРешение Гаусс: " << (x = SLAE(tmp)) << "Операций Гаусс: " << tmp.operation_in_slae;

    cout << "\n\nMatrix m: \n" << m;
    Cholesky_metod(m, l, t);

    cout << "\n\n\n-----МЕТОД ПРОГОНКИ-----\nCЛАУ:\n" << (tmp = tom + b_tom) << "\nРешение методом прогонки: " << (x = SLAE_Tomas(tom, b_tom)) << "Операций: " << tom.operation_in_slae;
    tom.operation_in_slae = 0;
    cout << "\n\nРешение LU: " << (x = SLAE_LU(tom, b_tom)) << "Операций LU: " << tom.operation_in_slae;
    tom.operation_in_slae = 0;
    cout << "\n\nРешение методом Холецкого: " << (x = SLAE_Cholesky(tom, b_tom)) << "Операций Холецкий: " << tom.operation_in_slae;
    cout << "\n\nРешение Гаусс: " << (x = SLAE(tmp)) << "Операций Гаусс: " << tmp.operation_in_slae;

    cout << "\n\n\n-----МЕТОД ПРОГОНКИ ЕЩЁ РАЗ-----\nCЛАУ:\n" << (tmp = tom2 + b_tom2) << "\nРешение методом прогонки: " << (x = SLAE_Tomas(tom2, b_tom2)) << "Операций: " << tom2.operation_in_slae;
    tom2.operation_in_slae = 0;
    cout << "\n\nРешение LU: " << (x = SLAE_LU(tom2, b_tom2)) << "Операций LU: " << tom2.operation_in_slae;
    cout << "\n\nРешение Гаусс: " << (x = SLAE(tmp)) << "Операций Гаусс: " << tmp.operation_in_slae << "\n\n\n";

    tmp = tom2 * x;
    for (int i = 0; i < 12; i++) x.vec[i] = tmp(i, 0);
    cout << "b найденный умножением на найденный х: " << x << "Нормы вектора b найденного: " << norm_first(x) << " " << norm_second(x) << " " << norm_infinity(x) << "\n";
    cout << "\nb известный: " << b_tom2 << "Нормы вектора b известного: " << norm_first(b_tom2) << " " << norm_second(b_tom2) << " " << norm_infinity(b_tom2) << "\n";
}
void Zadanie_5() {
    vector<double> B = { 1.058, 2.155626, -0.15626, -8.100005, 5.102034, 1.99958, -2.159875 };
    vctr x, b(B, 0);
    Matrix a = Matrix::rand_matr(7, 7), q, r, tmp;

    find_QR_matrix(a, q, r);
    cout << "Matrix A: \n" << a << "Matrix Q: \n" << q << "Matrix R: \n" << r << "Q * R =\n" << (tmp = q * r) << "SLAE Ax = b\n" << (tmp = a + b);
    cout << "\n\nОпераций для QR : " << q.operation_in_slae;
    q.operation_in_slae = 0;
    cout << "\n\nРешение QR: " << (x = SLAE_QR(q, r, b)) << "Операций QR: " << q.operation_in_slae;
    cout << "\n\nРешение LU: " << (x = SLAE_LU(a, b)) << "Операций LU: " << a.operation_in_slae;
    cout << "\n\nРешение Гаусс: " << (x = SLAE(tmp = a + b)) << "Операций Гаусс: " << tmp.operation_in_slae;

    tmp = a * x;
    for (int i = 0; i < 7; i++) x.vec[i] = tmp(i, 0);
    cout << "\n\n\n\nb найденный умножением на найденный х: " << x << "Нормы вектора b найденного: " << norm_first(x) << " " << norm_second(x) << " " << norm_infinity(x) << "\n";
    cout << "\n\nb из уравнения: " << b << "Нормы вектора b известного: " << norm_first(b) << " " << norm_second(b) << " " << norm_infinity(b) << "\n";

}
*/
void Zadanie_6() {
    vector<double> B = { 1.058, 2.155626, -0.15626, -8.100005, 5.102034 }, init1 = { 0,0,0,0,0 }, init2 = { 1,1,1,1,1 }, init3 = { -0.19505, 1.15, 0.845, 5.426048, -2.1002 };
    vctr x, b(B, 0), appr1(init1, 0), appr2(init2, 0), appr3(init3, 0);
    vector<vector<double>> A = { {4.82, -0.31, -0.49, 0.13, -0.27},
                                {-0.55, 3.19, -0.11, -0.47, 0.38},
                                {-0.29, -0.43, 2.91, -0.19, 0.51},
                                {0.49, -0.17, -0.35, 4.13, -0.41},
                                {-0.11, 0.53, -0.29, -0.31, 3.67} },
        C = {
    {4.0, 1.0, 0.0, 0.0, 0.0},
    {2.0, 3.0, 1.0, 0.0, 0.0},
    {0.0, 1.0, 3.0, 1.0, 0.0},
    {0.0, 0.0, 1.0, 4.0, 1.0},
    {0.0, 0.0, 0.0, 2.0, 5.0}
    },
        H = {
    {3.2, 0.8, 0.2, 0.0, 0.0},
    {1.5, 4.5, 1.2, 0.3, 0.0},
    {0.7, 1.9, 5.1, 1.5, 0.4},
    {0.2, 0.9, 2.7, 6.3, 1.2},
    {0.1, 0.3, 1.4, 2.9, 7.5}
    };
    Matrix a(A), c(C), h(H), tmp;

    cout << "\n-------------SLAE Ax = b-------------\n-------------Метод решения: простая итерация-------------\n" << (tmp = a + b);
    SLAE_Iteration(a, b, 0.0001, b);
    SLAE_Iteration(a, b, 0.0001, appr1);
    SLAE_Iteration(a, b, 0.0001, appr2);
    SLAE_Iteration(a, b, 0.0001, appr3);
    
    cout << "\n-------------SLAE Аx = b-------------\n-------------Метод решения: Зейдель-------------\n" << (tmp = a + b);
    SLAE_Seidel(a, b, 0.0001, b);
    SLAE_Seidel(a, b, 0.0001, appr1);
    SLAE_Seidel(a, b, 0.0001, appr2);
    SLAE_Seidel(a, b, 0.0001, appr3);

    cout << "\n-------------SLAE Сx = b-------------\n-------------Метод решения: простая итерация-------------\n" << (tmp = c + b);
    SLAE_Iteration(c, b, 0.0001, b);
    SLAE_Iteration(c, b, 0.0001, appr1);
    SLAE_Iteration(c, b, 0.0001, appr2);
    SLAE_Iteration(c, b, 0.0001, appr3);
    
    cout << "\n-------------SLAE Сx = b-------------\n-------------Метод решения: Зейдель-------------\n" << (tmp = c + b);
    SLAE_Seidel(c, b, 0.0001, b);
    SLAE_Seidel(c, b, 0.0001, appr1);
    SLAE_Seidel(c, b, 0.0001, appr2);
    SLAE_Seidel(c, b, 0.0001, appr3);

    cout << "\n-------------SLAE Нx = b-------------\n-------------Метод решения: простая итерация-------------\n" << (tmp = h + b);
    SLAE_Iteration(h, b, 0.0001, b);
    SLAE_Iteration(h, b, 0.0001, appr1);
    SLAE_Iteration(h, b, 0.0001, appr2);
    SLAE_Iteration(h, b, 0.0001, appr3);
    
    cout << "\n-------------SLAE Нx = b-------------\n-------------Метод решения: Зейдель-------------\n" << (tmp = h + b);
    SLAE_Seidel(h, b, 0.0001, b);
    SLAE_Seidel(h, b, 0.0001, appr1);
    SLAE_Seidel(h, b, 0.0001, appr2);
    SLAE_Seidel(h, b, 0.0001, appr3);
}
void Zadanie_7() {
    vector<vector<double>> A = {
        {3.2, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.5, 2.8, 0.7, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.7, 2.1, 0.9, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.9, 3.5, 1.1, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.1, 2.9, 0.8, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.8, 3.2, 0.6},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 4.1}
    },
        C = {
            {0.5, 0.2, 0.1, 0.1, 0.1},
            {0.1, 0.6, 0.1, 0.1, 0.1},
            {0.1, 0.1, 0.5, 0.2, 0.1},
            {0.1, 0.1, 0.1, 0.4, 0.2},
            {0.1, 0.1, 0.1, 0.1, 0.5}
    };
    vctr x;
    Matrix a(A), c(C), tmp;

    cout << "-----СОБСТВЕННОЕ ЧИСЛО: 4.78303228 --------\n";
    double lambda = power_law_method(a, 0.000001, x);
    cout << "\n\nВектор x: " << x << "Вторая норма вектора x: " << norm_second(x) << "\n";
    cout << (tmp = a * x) << (tmp = x * lambda);

    lambda = power_law_method_with_normalization(a, 0.000001, x);
    cout << "\n\nВектор x: " << x << "Вторая норма вектора x: " << norm_second(x) << "\n";
    cout << (tmp = a * x) << (tmp = x * lambda);


    cout << "\n\n-----СОБСТВЕННОЕ ЧИСЛО: 0.96118587 --------\n";
    lambda = power_law_method(c, 0.000001, x);
    cout << "\n\nВектор x: " << x << "Вторая норма вектора x: " << norm_second(x) << "\n";
    cout << (tmp = c * x) << (tmp = x * lambda);

    lambda = power_law_method_with_normalization(c, 0.000001, x);
    cout << "\n\nВектор x: " << x << "Вторая норма вектора x: " << norm_second(x) << "\n";
    cout << (tmp = c * x) << (tmp = x * lambda);
}

int main() {
    setlocale(LC_ALL, "Russian");
    vctr x, b({ 1.058, 2.155626, -0.15626, -8.100005, 5.102034 }, 0);
    vector<vector<double>> A = {
        {3.2, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.5, 2.8, 0.7, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.7, 2.1, 0.9, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.9, 3.5, 1.1, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.1, 2.9, 0.8, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.8, 3.2, 0.6},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 4.1}
    },
        C = {
            {0.5, 0.2, 0.1, 0.1, 0.1},
            {0.1, 0.6, 0.1, 0.1, 0.1},
            {0.1, 0.1, 0.5, 0.2, 0.1},
            {0.1, 0.1, 0.1, 0.4, 0.2},
            {0.1, 0.1, 0.1, 0.1, 0.5}
    };
    Matrix a(A), c(C), l, u, L, U, tmp;
    
    return 0;
}