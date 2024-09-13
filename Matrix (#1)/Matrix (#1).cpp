#include<iostream>
#include"Matrix.h"
using namespace std;

int main() {
    vector<vector<double>> z = { {3,2,-1,4},{2,-1,5,23},{1,7,-1,5} }, s = { {9,0,2},{3,8,8},{0,2,1} };
    Matrix a(z), b(s), tmp = inverse(b);
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