#include<iostream>
#include"Matrix.h"
using namespace std;

int main() {
    vector<double> v = { 1,2,-3 };
    vctr V(v, 1);
    
    cout << norm_first(V) << "   " << norm_second(V) << "    " << norm_infinity(V) << "\n";

return 0;
}