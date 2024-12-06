#include "Fredholm.h"

Matrix Fredholm::getMatrix_A(int n) {
    Matrix A(n + 1, n + 1);
    double h = (b - a) / n;
    vctr x(0);
    for (int i = 0; i <= n; i++) {
        x.vec.push_back(a + i * h);
    }
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            double coef = (j == 0 || j == n) ? 1 : (j % 3 == 0) ? 2 : 3;
            coef *= lambda * 3.0 * h * K(x.vec[i], x.vec[j]) / 8.0;
            if (i == j) A(i, j) = 1 - coef;
            else A(i, j) = -coef;
        }
    }
    return A;
}

vctr Fredholm::getVector_b(int n) {
    vctr B(0);
    double h = (b - a) / n;
    for (int i = 0; i <= n; i++) {
        B.vec.push_back(f(a + i * h));
    }
    return B;
}

vctr Fredholm::getVector_u(double epsilon) {
    int n = 3;
    Matrix A = getMatrix_A(n);
    vctr B = getVector_b(n), u, u_new = SLAE_QR(A, B);
    cout << u_new;
    bool done = false;
    while (!done) {
        n *= 2;
        A = getMatrix_A(n);
        B = getVector_b(n);
        u = u_new;
        u_new = SLAE_QR(A, B);
        cout << u_new;
        for (int i = 0; i < u.vec.size(); i++) {
            done = true;
            if (fabs(u_new.vec[i * 2] - u.vec[i]) >= epsilon) {
                done = false;
                break;
            }
        }
        u = u_new;
    }
    return u;
}