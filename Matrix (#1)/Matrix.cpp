#include "Matrix.h"

Matrix::Matrix(vctr v) {
    if (v.as_row) data.push_back(v.vec);
    else {
        data.resize(v.vec.size(), vector<double>(1, 0.0));
        for (int i = 0; i < v.vec.size(); ++i)
            data[i][0] = v.vec[i];
    }
}
Matrix Matrix::unit_matr(size_t n) {
    Matrix m(n);
    for (int i = 0; i < n; ++i)
        m(i, i) = 1;
    return m;
}
Matrix Matrix::rand_matr(size_t rows, size_t cols) {
    Matrix m(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            m(i, j) = rand() % 10;
        }
    }
    return m;
}
void Matrix::moreRow(size_t r) {
    vector<double> v(data[0].size(), 0);
    for (size_t i = 0; i < r; ++i) data.push_back(v);
}
void Matrix::moreCol(size_t c) {
    size_t n = data.size();
    for (size_t i = 0; i < n; ++i)
        for (size_t k = 0; k < c; ++k)
            data[i].push_back(0);
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
double conditioning_first(Matrix& m) {
    Matrix tmp = inverse(m);
    return norm_first(m) * norm_first(tmp);
}
double conditioning_infinity(Matrix& m) {
    Matrix tmp = inverse(m);
    return norm_infinity(m) * norm_infinity(tmp);
}
double conditioning_evkl(Matrix& m) {
    Matrix tmp = inverse(m);
    return norm_evkl(m) * norm_evkl(tmp);
}

// Нужны для обратной, определителя, решения СЛАУ
bool mod_max_row(Matrix& tmp, size_t K) {                   // |max| по главной диагонали, меняя строки. K - от какой строки вниз. При смене строк det меняет знак
    size_t Nmax = K, n = tmp.getRows(), m = tmp.getCols();
    for (size_t i = K + 1; i < n; i++) {
        if ( fabs(tmp(i,K)) > fabs(tmp(Nmax,K)) )
            Nmax = i;
    }
    if (Nmax != K) {
        for (size_t j = 0; j < m; j++)
            swap(tmp(K,j), tmp(Nmax,j));
        return true;
    }
    return false;
}
void do_zero_under_diag(Matrix& tmp, size_t K) {            // обнуление ниже главной диагонали, на диагонали единица. К - от какой строки вниз
    size_t n = tmp.getRows(), m = tmp.getCols();
    double del;
    for (size_t i = K; i < n; i++) {
        if (i == K) del = tmp(K, K);
        else del = tmp(i, K);
        if (!del || (i == K && del == 1)) continue;
        for (size_t j = K; j < m; j++) {
            if (i == K)
                tmp(K,j) /= del;
            else
                tmp(i,j) = tmp(i,j) - tmp(K,j) * del;
        }
    }
}
void do_zero_upper_diag(Matrix& tmp) {                      // обнуление выше главной диагонали
    size_t n = tmp.getRows(), m = tmp.getCols();
    double del;
    for (int K = n - 1; K > 0; K--) {
        for (int i = K - 1; i >= 0; i--) {
            del = tmp(i,K);
            for (int j = K; j < m; j++)
                tmp(i,j) = tmp(i,j) - tmp(K,j) * del;
        }
    }
}

Matrix inverse(Matrix& m) {
    if (m.getCols() != m.getRows()) {
        cout << "    Error: rows != cols for inverse\n\n";
        return NULL;
    }
    Matrix tmp = m;
    size_t n = tmp.getCols();
    // создание матрицы Е справа
    for (size_t i = 0; i < n; i++) {
        for (size_t j = n; j < n + n; j++) {
            if (j == n + i) tmp.data[i].push_back(1);
            else tmp.data[i].push_back(0);
        }
    }
    for (size_t K = 0; K < n; K++) {
        mod_max_row(tmp, K);                // |max| по главной диагонали, меняя строки
        do_zero_under_diag(tmp, K);         // обнуление ниже главной диагонали, на диагонали единица
    }
    do_zero_upper_diag(tmp);                // обнуление выше главной диагонали
    // удаление исходной матрицы
    for (int i = 0; i < n; ++i)
        tmp.data[i].erase(tmp.data[i].begin(), tmp.data[i].begin() + n);
    return tmp;
}
double determinant(Matrix& m) {
    if (m.getCols() != m.getRows()) return 0;
    Matrix tmp = m;
    size_t n = m.getCols();
    double del, det = 1;
    for (size_t K = 0; K < n; K++) {
        if (mod_max_row(tmp, K)) det = -det;
        det *= tmp.data[K][K];
        do_zero_under_diag(tmp, K);
    }
    return det;
}

vctr SLAE(Matrix& mtr) {
    if (mtr.getCols() != mtr.getRows()+1) {
        cout << "    Error: wrong size for SLAE\n\n";
        return NULL;
    }
    Matrix tmp = mtr;
    size_t n = tmp.getRows();
    vctr ans(0);
    for (size_t K = 0; K < n; K++) {
        mod_max_row(tmp, K);                // |max| по главной диагонали, меняя строки
        do_zero_under_diag(tmp, K);         // обнуление ниже главной диагонали, на диагонали единица
    }
    do_zero_upper_diag(tmp);                // обнуление выше главной диагонали
    for (size_t i = 0; i < n; ++i)
        ans.vec.push_back(tmp(i, n));
    return ans;
}

Matrix transpos(Matrix& m) {
    Matrix tmp = m;
    for (int i = 0; i < tmp.getRows(); ++i) {
        for (int j = 0; j < tmp.getCols(); ++j) {
            tmp(i, j) = m(j, i);
        }
    }
    return tmp;
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
    if (l.getRows() != r.getRows() || l.getCols() != r.getCols()) {
        cout << "    Error: wrong sizes (matrix +(-) matrix)\n\n";
        return NULL;
    }
    Matrix tmp(l.getRows(), l.getCols());
    for (int i = 0; i < tmp.getRows(); ++i) {
        for (int j = 0; j < tmp.getCols(); ++j) {
            tmp(i, j) = l(i, j) + r(i, j);
        }
    }
    return tmp;
}
Matrix operator-(Matrix& l, Matrix& r) {
    Matrix tmp = -1 * r;
    return l + tmp;
}
Matrix operator*(Matrix& l, Matrix& r) {
    if (l.getCols() != r.getRows()) {
        cout << "    Error: wrong sizes (matrix * matrix)\n\n";
        return NULL;
    }
    size_t n = l.getRows(), m = r.getCols(), p = l.getCols();
    Matrix tmp(n,m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < p; ++k) {
                tmp(i, j) += l.data[i][k] * r.data[k][j];
            }
            if (fabs(tmp(i, j)) < 1e-14) tmp(i, j) = 0;
        }
    }
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
            out << setw(12) << setprecision(10) << m(i, j) << " ";
        }
        out << "\n";
    }
    out << "\n";
    return out;
}
Matrix operator+(Matrix& m, vctr& v) {
    if (v.as_row && v.vec.size() == m.getCols()) {
        Matrix tmp = m;
        tmp.moreRow(1);
        size_t n = tmp.getCols(), i = m.getRows();
        for (int j = 0; j < n; ++j)
            tmp(i, j) = v.vec[j];
        return tmp;
    }
    if (!v.as_row && v.vec.size() == m.getRows()) {
        Matrix tmp = m;
        tmp.moreCol(1);
        size_t n = tmp.getRows(), j = m.getCols();
        for (int i = 0; i < n; ++i)
            tmp(i, j) = v.vec[i];
        return tmp;
    }
    cout << "    Error: wrong sizes or position (matrix + vctr)\n\n";
    return NULL;
}