#include<iostream>
#include"Matrix.h"
#include <fstream>
using namespace std;

void Zadanie_2() {
    int n;
    cout << "    Задание 2\n Введите размер матрицы А: ";
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
    for (int i = 0; i < n; ++i) b.vec.push_back(B(i, 0));
    A = A + b;
    cout << "Решение х:\n" << (x = SLAE(A));
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
    cout << "Введённая матрица:\n" << Var << "Ответ:\n" << (Ans = SLAE(Var));
}

int main() {
    setlocale(LC_ALL, "Russian");
    SLAE_from_file();
    Zadanie_2();


    vector<vector<double>> z = { {3,2,-1,4},{2,-1,5,23},{1,7,-1,5} }, s = { {9,0,2},{3,8,8},{0,2,1} };
    Matrix a(z);
    vector<double> v = { 1,2,-3 }, k = { 4,3,6 };
    vctr V(v, 1), K(k, 1), Z = SLAE(a); 
    cout << a << Z; 
   /* 
    cout << norm_first(V) << "   " << norm_second(V) << "    " << norm_infinity(V) << "\n";
    cout << norm_first(a) << "   " << norm_evkl(a) << "    " << norm_infinity(a) << "\n";
    cout << norm_first(tmp) << "   " << norm_evkl(tmp) << "    " << norm_infinity(tmp) << "\n";
    cout << conditioning_first(a) << "   " << conditioning_evkl(a) << "    " << conditioning_infinity(a) << "\n";
    cout << conditioning_first(tmp) << "   " << conditioning_evkl(tmp) << "    " << conditioning_infinity(tmp) << "\n";*/

return 0;
}