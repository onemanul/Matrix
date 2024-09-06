#include "Matrix.h"

Matrix::Matrix(vctr v) {
    if (v.as_row) data.push_back(v.vec);
    else {
        data.resize(v.vec.size(), vector<double>(1, 0.0));
        for (int i = 0; i < v.vec.size(); ++i)
            data[i][0] = v.vec[i];
    }
}

double norm_first(Matrix& m) {
    double mx = 0, s;
    for (int j = 0; j < m.getCols(); ++j) {
        s = 0;
        for (int i = 0; i < m.getRows(); ++i) {
            s += fabs(m(i,j));
        }
        if (mx<s) mx = s;
    }
    return mx;
}
double norm_infinity(Matrix& m) {
    double mx = 0, s;
    for (int i = 0; i < m.getRows(); ++i) {
        s = 0;
        for (int j = 0; j < m.getCols(); ++j) {
            s += fabs(m(i, j));
        }
        if (mx < s) mx = s;
    }
    return mx;
}
double norm_evkl(Matrix& m) {
    double s = 0;
    for (int i = 0; i < m.getRows(); ++i) {
        for (int j = 0; j < m.getCols(); ++j) {
            s += m(i, j)*m(i,j);
        }
    }
    return sqrt(s);
}

friend Matrix inverse(Matrix& m);
Matrix transpos(Matrix& m) {
    Matrix tmp = m;
    for (int i = 0; i < tmp.getRows(); ++i) {
        for (int j = 0; j < tmp.getCols(); ++j) {
            tmp(i,j) = m(j,i);
        }
    }
    return tmp;
}
Matrix unit_matr(size_t n) {
    Matrix m(n);
    for (int i = 0; i < n; ++i)
        m(i, i) = 1;
    return m;
}
Matrix rand_matr(size_t rows, size_t cols) {
    Matrix m(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            m(i, j) = rand() % 10;
        }
    }
    return m;
}

Matrix operator*(double A, Matrix& m) {
    Matrix tmp = m;
    for (int i = 0; i < tmp.getRows(); ++i) {
        for (int j = 0; j < tmp.getCols(); ++j) {
            tmp(i, j) *= A;
        }
    }
    return tmp;
}
Matrix operator*(Matrix& m, double A) {
    return A * m;
}
Matrix operator*(vctr& v, Matrix& m) {
    Matrix tmp(v);
    return tmp * m;
}
Matrix operator*(Matrix& m, vctr& v) {
    Matrix tmp(v);
    return m*tmp;
}
Matrix operator+(Matrix& l, Matrix& r) {
    if (l.getRows() != r.getRows() || l.getCols() != r.getCols()) return NULL;
    Matrix tmp(l.getRows(), l.getCols());
    for (int i = 0; i < tmp.getRows(); ++i) {
        for (int j = 0; j < tmp.getCols(); ++j) {
            tmp(i, j) = l(i, j) + r(i, j);
        }
    }
    return tmp;
}
friend Matrix operator*(Matrix& l, Matrix& r) {
    if (l.getCols() != r.getRows()) return NULL;
    Matrix tmp(l.getRows(), r.getCols());

    return tmp;
}
bool operator==(Matrix& left, Matrix& right) {
    if (left.getRows() != right.getRows() || left.getCols() != right.getCols()) return false;
    for (int i = 0; i < left.getRows(); ++i) {
        for (int j = 0; j < left.getCols(); ++j) {
            if (left(i,j) != right(i,j)) return false;
        }
    }
    return true;
}
ostream& operator<<(ostream& out, Matrix& m) {
    for (int i = 0; i < m.getRows(); ++i) {
        for (int j = 0; j < m.getCols(); ++j) {
            out << m(i, j);
        }
        out << "\n";
    }
    out << "\n";
    return out;
}



void swap(double* a, double *b) {
    double t = *a;
    *a = *b;
    *b = t;
}

void Obratnaya(double** a, int n) {
    int Nmax;
    double del;
    // дозаполнение матрицы Е справа
    for (int i = 0; i < n; i++)
        for (int j = n; j < n + n; j++) {
            if (j == n + i) a[i][j] = 1;
            else a[i][j] = 0;
        }

    for (int K = 0; K < n; K++) {
        // |max| по главной диагонали
        Nmax = K;
        for (int i = K + 1; i < n; i++) {
            if (fabs(a[i][K]) > fabs(a[Nmax][K]))
                Nmax = i;
        }
        if (Nmax != K) {
            for (int j = 0; j < n + n; j++)
                swap(&a[K][j], &a[Nmax][j]);
        }
        // обнуление ниже главной диагонали
        for (int i = K; i < n; i++) {
            if (i == K) del = a[K][K];
            else del = a[i][K];
            for (int j = K; j < n + n; j++) {
                if (i == K)
                    a[K][j] /= del;
                else
                    a[i][j] = a[i][j] - a[K][j] * del;
            }
        }
    }

    // обнуление выше главной диагонали
    for (int K = n - 1; K > 0; K--) {
        for (int i = K - 1; i >= 0; i--) {
            del = a[i][K];
            for (int j = K; j < n + n; j++)
                a[i][j] = a[i][j] - a[K][j] * del;
        }
    }
}