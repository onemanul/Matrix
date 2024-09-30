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
void Matrix::changeRow(size_t r) {
    vector<double> v(data[0].size(), 0);
    data.resize(r, v);
}
void Matrix::changeCol(size_t c) {
    size_t n = data.size();
    for (size_t i = 0; i < n; ++i)
        data[i].resize(c, 0);
}
void Matrix::swap_rows(size_t a, size_t b) {
    if (a >= data.size() || b >= data.size()) return;
    vector<double> tmp = data[a];
    data[a] = data[b];
    data[b] = tmp;
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
int mod_max_row(Matrix& tmp, size_t K) {                   // |max| по главной диагонали, меняя строки. K - от какой строки вниз. При смене строк det меняет знак
    size_t Nmax = K, n = tmp.getRows(), m = tmp.getCols();
    for (size_t i = K + 1; i < n; i++) {
        if ( fabs(tmp(i,K)) > fabs(tmp(Nmax,K)) )
            Nmax = i;
    }
    if (Nmax != K) {
        for (size_t j = 0; j < m; j++)
            swap(tmp(K,j), tmp(Nmax,j));
        return Nmax;
    }
    return 0;
}
int mod_max_col(Matrix& tmp, size_t K) {                    // НУЖЕН ТОЛЬКО ДЛЯ МЕТОДОВ СЛАУ! |max| по главной диагонали, меняя столбцы. K - от какого столбца вправо.
    size_t Nmax = K, n = tmp.getRows();                     // При смене столбца надо учитывать порядок Х в ответах; возращает номер столбца, с которым был обмен
    for (size_t j = K + 1; j < n; ++j) {
        if (fabs(tmp(K, j)) > fabs(tmp(K, Nmax)))
            Nmax = j;
    }
    if (Nmax != K) {
        for (size_t i = 0; i < n; i++)
            swap(tmp(i, K), tmp(i, Nmax));
        return Nmax;
    }
    return 0;
}
int both_metods(Matrix& tmp, size_t K) {                    // НУЖЕН ТОЛЬКО ДЛЯ МЕТОДОВ СЛАУ!
    size_t n = tmp.getRows(), I = K, J = K;
    for (size_t i = K; i < n; i++) {
        for (size_t j = K; j < n; j++) {
            if (fabs(tmp(i, j)) > fabs(tmp(I, J))) {
                I = i;
                J = j;
            }
        }
    }
    if (I != K) {
        for (size_t j = 0; j <= n; j++)
            swap(tmp(K, j), tmp(I, j));
    }
    if (J != K) {
        for (size_t i = 0; i < n; i++)
            swap(tmp(i, K), tmp(i, J));
        return J;
    }
    return 0;
}
void do_zero_under_diag(Matrix& tmp, size_t K) {            // обнуление ниже главной диагонали, на диагонали единица. К - от какой строки вниз
    size_t n = tmp.getRows(), m = tmp.getCols();
    double del;
    for (size_t i = K; i < n; i++) {
        if (i == K) del = tmp(K, K);
        else del = tmp(i, K);
        if (!del || (i == K && del == 1)) continue;
        for (size_t j = K; j < m; j++) {
            if (i == K) {
                tmp(K, j) /= del;
                ++tmp.operation_in_slae;
            }
            else {
                if (tmp(K, j)) {
                    tmp(i, j) = tmp(i, j) - tmp(K, j) * del;
                    tmp.operation_in_slae += 2;
                }
            }
        }
    }
}
bool do_zero_under_diag_for_LU(Matrix& U, size_t K, Matrix& L) {            // обнуление ниже главной диагонали БЕЗ СМЕНЫ СТРОК (где-то можно сломаться). К - от какой строки вниз
    size_t n = U.getRows(), m = U.getCols();
    double mu;
    if (!U(K, K)) return false;
    for (size_t i = K+1; i < n; i++) {
        if (!U(i, K)) continue;
        else mu = U(i, K)/U(K,K);
        for (size_t j = K; j < n; j++) {
            if (U(K, j)) {
                U(i, j) = U(i, j) - U(K, j) * mu;
            }
        }
        L(i, K) = mu;
    }
    return true;
}
void do_zero_upper_diag(Matrix& tmp) {                      // обнуление выше главной диагонали
    size_t n = tmp.getRows(), m = tmp.getCols();
    double del;
    for (int K = n - 1; K > 0; K--) {
        for (int i = K - 1; i >= 0; i--) {
            del = tmp(i,K);
            if (del) {
                for (int j = K; j < m; j++) {
                    if (tmp(K, j)) {
                        tmp(i, j) = tmp(i, j) - tmp(K, j) * del;
                        tmp.operation_in_slae += 2;
                    }
                }
            }
        }
    }
}
void do_zero_upper_diag_without_ones(Matrix& tmp) {                      // сделать единицы на диагонали
    size_t n = tmp.getRows(), m = tmp.getCols();
    double del;
    for (int K = n - 1; K >= 0; K--) {
        if (tmp(K, K) != 1) {
            del = tmp(K,K);
            for (int j = K; j < m; ++j) {
                if (tmp(K, j)) {
                    tmp(K, j) /= del;
                    ++tmp.operation_in_slae;
                }
            }
        }
        for (int i = K - 1; i >= 0; i--) {
            del = tmp(i, K);
            if (del) {
                for (int j = K; j < m; j++) {
                    if (tmp(K, j)) {
                        tmp(i, j) = tmp(i, j) - tmp(K, j) * del;
                        tmp.operation_in_slae += 2;
                    }
                }
            }
        }
    }
}
// -----------------------------------------------

Matrix inverse(Matrix& m) {
    if (m.getCols() != m.getRows()) {
        cout << "    Error: rows != cols for inverse\n\n";
        return NULL;
    }
    Matrix tmp = m;
    size_t n = tmp.getCols();
    // создание матрицы Е справа
    tmp.changeCol(n * 2);                                   // заполнены нулями
    for (size_t i = 0; i < n; ++i) tmp(i, i + n) = 1;       // единицы на диагонали матрицы Е
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
Matrix transpos(Matrix& m) {
    Matrix tmp = m;
    for (int i = 0; i < tmp.getRows(); ++i) {
        for (int j = 0; j < tmp.getCols(); ++j) {
            tmp(i, j) = m(j, i);
        }
    }
    return tmp;
}
double determinant(Matrix& m) {
    if (m.getCols() != m.getRows()) return 0;
    Matrix tmp = m;
    size_t n = m.getCols();
    double det = 1;
    for (size_t K = 0; K < n; K++) {
        if (mod_max_row(tmp, K)) det = -det;
        det *= tmp.data[K][K];
        do_zero_under_diag(tmp, K);
    }
    return det;
}
double determinant_with_LU(Matrix& m) {
    Matrix L, U;
    if (!find_LU_matrix(m, L, U)) return 0;
    else return determinant(U);
}
vctr SLAE(Matrix& mtr) {
    if (mtr.getCols() != mtr.getRows() + 1) {
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
    mtr.operation_in_slae = tmp.operation_in_slae - mtr.operation_in_slae;
    return ans;
}
vctr SLAE(Matrix& mtr, int metod) {             // 0 - строки, 1 - столбцы, 2 - оба метода
    if (!metod) return SLAE(mtr);
    if (mtr.getCols() != mtr.getRows()+1) {
        cout << "    Error: wrong size for SLAE\n\n";
        return NULL;
    }
    Matrix tmp = mtr;
    size_t n = tmp.getRows(), changed;
    vector<int> index(n);
    vector<double> ans(n);
    for (size_t i = 0; i < n; ++i) index[i] = i;
    for (size_t K = 0; K < n; K++) {
        if (metod==1) changed = mod_max_col(tmp, K);                      // |max| по главной диагонали, меняя столбцы
        else changed = both_metods(tmp, K);
        if (changed) swap(index[K], index[changed]);
        do_zero_under_diag(tmp, K);         // обнуление ниже главной диагонали, на диагонали единица
    }
    do_zero_upper_diag(tmp);                // обнуление выше главной диагонали
    for (size_t i = 0; i < n; ++i) ans[index[i]] = tmp(i, n);
    vctr v(ans, 0);
    mtr.operation_in_slae = tmp.operation_in_slae;
    return v;
}
bool find_LU_matrix(Matrix& A, Matrix& L, Matrix& U) {
    if (A.getCols() != A.getRows()) {
        cout << "    Error: wrong size for LU_matrix\n\n";
        return false;
    }
    size_t n = A.getRows();
    Matrix tmp_U = A, tmp_L = Matrix::unit_matr(n);
    for (size_t K = 0; K < n-1; K++) {
        if (!do_zero_under_diag_for_LU(tmp_U, K, tmp_L)) {           // обнуление ниже главной диагонали
            cout << "    Error: LU is imposible\n\n";
            return false;
        }
    }
    L = tmp_L;
    U = tmp_U;
    return true;
}
bool find_P_matrix(Matrix& A, Matrix& P) {
    if (A.getCols() != A.getRows()) {
        cout << "    Error: wrong size for P_matrix\n\n";
        return false;
    }
    size_t n = A.getRows(), Nmax;
    Matrix tmp_A = A, tmp_P = Matrix::unit_matr(n);
    for (size_t K = 0; K < n - 1; K++) {
        Nmax = mod_max_row(tmp_A, K);
        if (Nmax) tmp_P.swap_rows(K, Nmax);
    }
    P = tmp_P;
    return true;
}
vctr SLAE_LU(Matrix& mtr, vctr& b) {
    if (mtr.getCols() != mtr.getRows()) {
        cout << "    Error: wrong size for SLAE_LU without vector b\n\n";
        return NULL;
    }
    Matrix A = mtr, L, U;
    if (!find_LU_matrix(A, L, U)) return NULL;
    Matrix tmp = L + b;
    size_t n = tmp.getRows();
    vctr y(0), x(0);
    for (size_t K = 0; K < n; K++) {
        do_zero_under_diag(tmp, K);          // обнуление ниже главной диагонали, на диагонали УЖЕ единица
        y.vec.push_back(tmp(K, n));
    }
    mtr.operation_in_slae = tmp.operation_in_slae;
    tmp = U + y;
    do_zero_upper_diag_without_ones(tmp);                               // обнуление выше главной диагонали, на диагонали НЕ нули
    for (size_t i = 0; i < n; ++i) x.vec.push_back(tmp(i, n));
    mtr.operation_in_slae += tmp.operation_in_slae;
    return x;
}
vctr SLAE_LUP(Matrix& mtr, vctr& b) {
    if (mtr.getCols() != mtr.getRows()) {
        cout << "    Error: wrong size for SLAE_LUP without vector b\n\n";
        return NULL;
    }
    Matrix A = mtr, L, U, P; 
    if (!find_P_matrix(A, P) || !find_LU_matrix((A = P*A), L, U)) return NULL;
    size_t n = A.getRows();
    Matrix B = P * b;
    vctr pb(0), y(0), x(0);
    for (size_t i = 0; i < n; ++i) pb.vec.push_back(B(i, 0));
    Matrix tmp = L + pb;
    for (size_t K = 0; K < n; K++) {
        do_zero_under_diag(tmp, K);          // обнуление ниже главной диагонали, на диагонали УЖЕ единица
        y.vec.push_back(tmp(K, n));
    }
    mtr.operation_in_slae = tmp.operation_in_slae;
    tmp = U + y;
    do_zero_upper_diag_without_ones(tmp);                               // обнуление выше главной диагонали, на диагонали НЕ нули
    for (size_t i = 0; i < n; ++i) x.vec.push_back(tmp(i, n));
    mtr.operation_in_slae += tmp.operation_in_slae;
    return x;
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
        cout << "    Error: wrong sizes (matrix * matrix)\n\n" << l.getCols() << " " << r.getRows() << "\n";
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
            out << setw(12) << setprecision(8) << m(i, j) << " ";
        }
        out << "\n";
    }
    out << "\n";
    return out;
}
Matrix operator+(Matrix& m, vctr& v) {
    if (v.as_row && v.vec.size() == m.getCols()) {
        Matrix tmp = m;
        tmp.changeRow(tmp.getRows()+1);
        size_t n = tmp.getCols(), i = m.getRows();
        for (int j = 0; j < n; ++j)
            tmp(i, j) = v.vec[j];
        return tmp;
    }
    if (!v.as_row && v.vec.size() == m.getRows()) {
        Matrix tmp = m;
        tmp.changeCol(tmp.getCols()+1);
        size_t n = tmp.getRows(), j = m.getCols();
        for (int i = 0; i < n; ++i)
            tmp(i, j) = v.vec[i];
        return tmp;
    }
    cout << "    Error: wrong sizes or position (matrix + vctr)\n\n";
    return NULL;
}