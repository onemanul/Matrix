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
            m(i, j) = (double)(rand()) / (double)(RAND_MAX);
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

bool Matrix::is_symmetry() {
    if (data.size() != data[0].size()) return false;
    size_t n = data.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (data[i][j] != data[j][i]) return false;
        }
    }
    return true;
}
bool Matrix::is_positive_definite() {
    if (data.size() != data[0].size()) return false;
    Matrix tmp(data);
    for (size_t i = tmp.getCols()-1; i > 0; --i) {
        if (tmp.determinant() <= 0) return false;
        tmp.changeRow(i);
        tmp.changeCol(i);
    }
    if (tmp(0, 0) > 0) return true;
    else return false;
}
bool Matrix::is_tridiagonal() {
    if (data.size() != data[0].size() || data.size() < 3) return false;
    size_t n = data.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            if (j < i + 2) {
                if (!data[i][j] || !data[j][i]) return false;
            } 
            else {
                if (data[i][j] || data[j][i]) return false;
            }
        }
    }
    return true;
}
vctr Matrix::to_vector() {
    if (data.size() == 1) {
        vector<double> v = data[0];
        vctr V(v, 1);
        return V;
    }
    else {
        if (data[0].size() == 1) {
            vector<double> v;
            for (size_t i = 0; i < data.size(); ++i)
                v.push_back(data[i][0]);
            vctr V(v, 0);
            return V;
        }
        else return NULL;
    }
}


// Нужны для обратной, определителя, решений СЛАУ-------------------
int mod_max_row(Matrix& tmp, size_t K) {                   // |max| по главной диагонали, меняя строки. K - от какой строки вниз. При смене строк det меняет знак
    size_t Nmax = K, n = tmp.getRows(), m = tmp.getCols();
    for (size_t i = K + 1; i < n; i++) {
        if (fabs(tmp(i, K)) > fabs(tmp(Nmax, K)))
            Nmax = i;
    }
    if (Nmax != K) {
        for (size_t j = 0; j < m; j++)
            swap(tmp(K, j), tmp(Nmax, j));
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
            }
            else {
                if (tmp(K, j)) {
                    tmp(i, j) = tmp(i, j) - tmp(K, j) * del;
                }
            }
        }
    }
}
bool do_zero_under_diag_for_LU(Matrix& U, size_t K, Matrix& L) {            // обнуление ниже главной диагонали БЕЗ СМЕНЫ СТРОК (где-то можно сломаться). К - от какой строки вниз
    size_t n = U.getRows(), m = U.getCols();
    double mu;
    if (!U(K, K)) return false;
    for (size_t i = K + 1; i < n; i++) {
        if (!U(i, K)) continue;
        else mu = U(i, K) / U(K, K);
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
            del = tmp(i, K);
            if (del) {
                for (int j = K; j < m; j++) {
                    if (tmp(K, j)) {
                        tmp(i, j) = tmp(i, j) - tmp(K, j) * del;
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
            del = tmp(K, K);
            for (int j = K; j < m; ++j) {
                if (tmp(K, j)) {
                    tmp(K, j) /= del;
                }
            }
        }
        for (int i = K - 1; i >= 0; i--) {
            del = tmp(i, K);
            if (del) {
                for (int j = K; j < m; j++) {
                    if (tmp(K, j)) {
                        tmp(i, j) = tmp(i, j) - tmp(K, j) * del;
                    }
                }
            }
        }
    }
}
Matrix B_matrix(Matrix& A) {
    double n = A.getRows(), divisor;
    Matrix B(n, n);
    for (size_t i = 0; i < n; ++i) {
        divisor = A(i, i);
        for (size_t j = 0; j < n; ++j) {
            B(i, j) = -A(i, j) / divisor;
        }
        B(i, i) = 0;
    }
    return B;
}
void B1_B2_matrix(Matrix& A, Matrix& B1, Matrix& B2) {
    Matrix B = B_matrix(A);
    double n = B.getRows();
    B2.changeCol(n);
    B2.changeRow(n);
    for (size_t i = 0; i < n - 1; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            B2(i, j) = B(i, j);
            B(i, j) = 0;
        }
    }
    B1 = B;
}
// ----------------------------------------------------------------


double Matrix::norm_first() {
    if (!data.size()) return 0;
    double mx = 0, s;
    for (int j = 0; j < data[0].size(); ++j) {
        s = 0;
        for (int i = 0; i < data.size(); ++i) {
            s += fabs(data[i][j]);
        }
        if (mx < s) mx = s;
    }
    return mx;
}
double Matrix::norm_infinity() {
    double mx = 0, s;
    for (int i = 0; i < data.size(); ++i) {
        s = 0;
        for (int j = 0; j < data[0].size(); ++j) {
            s += fabs(data[i][j]);
        }
        if (mx < s) mx = s;
    }
    return mx;
}
double Matrix::norm_evkl() {
    double s = 0;
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[0].size(); ++j) {
            s += data[i][j] * data[i][j];
        }
    }
    return sqrt(s);
}

double Matrix::determinant() {
    if (getCols() != getRows()) return 0;
    Matrix tmp(data);
    size_t n = getCols();
    double det = 1;
    for (size_t K = 0; K < n; K++) {
        if (mod_max_row(tmp, K)) det = -det;
        det *= tmp(K,K);
        do_zero_under_diag(tmp, K);
    }
    return det;
}
double Matrix::determinant_with_LU() {
    Matrix L, U;
    if (!find_LU_matrix(L, U)) return 0;
    else {
        double s = 1;
        for (size_t i = 0; i < U.getRows(); ++i)
            s *= U(i, i);
        return s;
    }
}

bool Matrix::find_LU_matrix(Matrix& L, Matrix& U) {
    if (getCols() != getRows()) {
        cout << "    Error: wrong size for LU_matrix\n\n";
        return false;
    }
    size_t n = getRows();
    Matrix tmp_U(data), tmp_L = Matrix::unit_matr(n);
    for (size_t K = 0; K < n - 1; K++) {
        if (!do_zero_under_diag_for_LU(tmp_U, K, tmp_L)) {
            cout << "    Error: LU is impossible\n\n";
            return false;
        }
    }
    L = tmp_L;
    U = tmp_U;
    return true;
}
bool Matrix::find_P_matrix(Matrix& P) {
    if (getCols() != getRows()) {
        cout << "    Error: wrong size for P_matrix\n\n";
        return false;
    }
    size_t n = getRows(), Nmax;
    Matrix tmp_A(data), tmp_P = Matrix::unit_matr(n);
    for (size_t K = 0; K < n - 1; K++) {
        Nmax = mod_max_row(tmp_A, K);
        if (Nmax) tmp_P.swap_rows(K, Nmax);
    }
    P = tmp_P;
    return true;
}
bool Matrix::find_L_Lt_matrix(Matrix& L, Matrix& L_t) {
    if (getRows() == 1 || !is_symmetry() || !is_positive_definite()) {
        cout << "    Error: Cholesky decomposition is impossible\n\n";
        return false;
    }
    size_t n = getCols();
    Matrix tmp(n);
    double sum;
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = j; i < n; ++i) {
            sum = 0;
            for (int k = j - 1; k >= 0; --k) {
                if (i == j) sum += tmp(i, k) * tmp(i, k);
                else sum += tmp(i, k) * tmp(j, k);
            }
            if (i == j) tmp(j, j) = sqrt(data[j][j] - sum);
            else tmp(i, j) = (data[i][j] - sum) / tmp(j, j);
        }
    }
    L = tmp;
    L_t = transpos(tmp);
    return true;
}
bool Matrix::find_QR_matrix(Matrix& Q, Matrix& R) {
    if (getCols() != getRows()) {
        cout << "    Error: wrong size for QR_matrix\n\n";
        return false;
    }
    size_t n = getRows();
    Matrix T = Matrix::unit_matr(n), T_ij, tmp_A(data);
    double c, s, a;
    for (size_t j = 0; j < n - 1; j++) {
        for (size_t i = j + 1; i < n; i++) {
            c = tmp_A(j, j) / sqrt(tmp_A(j, j) * tmp_A(j, j) + tmp_A(i, j) * tmp_A(i, j));
            s = tmp_A(i, j) / sqrt(tmp_A(j, j) * tmp_A(j, j) + tmp_A(i, j) * tmp_A(i, j));
            for (size_t k = j; k < n; ++k) {
                a = tmp_A(j, k);
                tmp_A(j, k) = c * tmp_A(j, k) + s * tmp_A(i, k);
                tmp_A(i, k) = -s * a + c * tmp_A(i, k);
            }
            T_ij = Matrix::unit_matr(n);
            T_ij(j, j) = c;
            T_ij(i, i) = c;
            T_ij(j, i) = s;
            T_ij(i, j) = -s;
            T = T_ij * T;
        }
    }
    tmp_A.data = data;
    R = T * tmp_A;
    Q = transpos(T);
    return true;
}


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
    return ans;
}
vctr SLAE(Matrix& mtr, int metod) {             // 0 - строки, 1 - столбцы, 2 - оба метода
    if (!metod) return SLAE(mtr);
    if (mtr.getCols() != mtr.getRows() + 1) {
        cout << "    Error: wrong size for SLAE\n\n";
        return NULL;
    }
    Matrix tmp = mtr;
    size_t n = tmp.getRows(), changed;
    vector<int> index(n);
    vector<double> ans(n);
    for (size_t i = 0; i < n; ++i) index[i] = i;
    for (size_t K = 0; K < n; K++) {
        if (metod == 1) changed = mod_max_col(tmp, K);                      // |max| по главной диагонали, меняя столбцы
        else changed = both_metods(tmp, K);
        if (changed) swap(index[K], index[changed]);
        do_zero_under_diag(tmp, K);         // обнуление ниже главной диагонали, на диагонали единица
    }
    do_zero_upper_diag(tmp);                // обнуление выше главной диагонали
    for (size_t i = 0; i < n; ++i) ans[index[i]] = tmp(i, n);
    vctr v(ans, 0);
    return v;
}
vctr SLAE_LU(Matrix& mtr, vctr& b) {
    if (mtr.getCols() != mtr.getRows() || b.vec.size() != mtr.getRows()) {
        cout << "    Error: wrong size for SLAE_LU without vector b or wrong vector b\n\n";
        return NULL;
    }
    Matrix A = mtr, L, U;
    if (!A.find_LU_matrix(L, U)) return NULL;
    else return SLAE_LU(L, U, b);
}
vctr SLAE_LU(Matrix& L, Matrix& U, vctr& b) {
    if (L.getCols() != L.getRows() || U.getCols() != U.getRows() || L.getCols() != U.getCols() || b.vec.size() != L.getRows()) {
        cout << "    Error: wrong size of L or U or b for SLAE_LU\n\n";
        return NULL;
    }
    Matrix tmp = L + b;
    size_t n = tmp.getRows();
    vctr y(0), x(0);
    for (size_t K = 0; K < n; K++) {
        do_zero_under_diag(tmp, K);          // обнуление ниже главной диагонали, на диагонали УЖЕ единица
        y.vec.push_back(tmp(K, n));
    }
    tmp = U + y;
    do_zero_upper_diag_without_ones(tmp);                               // обнуление выше главной диагонали, на диагонали НЕ нули
    for (size_t i = 0; i < n; ++i) x.vec.push_back(tmp(i, n));
    return x;
}
vctr SLAE_LUP(Matrix& mtr, vctr& b) {
    if (mtr.getCols() != mtr.getRows() || b.vec.size() != mtr.getRows()) {
        cout << "    Error: wrong size for SLAE_LUP without vector b or wrong vector b\n\n";
        return NULL;
    }
    Matrix A = mtr, L, U, P; 
    if (!A.find_P_matrix(P) || !(P*A).find_LU_matrix(L, U)) return NULL;
    size_t n = A.getRows();
    Matrix B = P * b;
    vctr pb(0), y(0), x(0);
    for (size_t i = 0; i < n; ++i) pb.vec.push_back(B(i, 0));
    Matrix tmp = L + pb;
    for (size_t K = 0; K < n; K++) {
        do_zero_under_diag(tmp, K);          // обнуление ниже главной диагонали, на диагонали УЖЕ единица
        y.vec.push_back(tmp(K, n));
    }
    tmp = U + y;
    do_zero_upper_diag_without_ones(tmp);                               // обнуление выше главной диагонали, на диагонали НЕ нули
    for (size_t i = 0; i < n; ++i) x.vec.push_back(tmp(i, n));
    return x;
}
vctr SLAE_Cholesky(Matrix& mtr, vctr& b) {
    if (mtr.getCols() != mtr.getRows() || b.vec.size() != mtr.getRows()) {
        cout << "    Error: wrong size for SLAE_Cholesky without vector b or wrong vector b\n\n";
        return NULL;
    }
    Matrix A = mtr, L, L_t;
    if (!A.find_L_Lt_matrix(L, L_t)) return NULL;
    else return SLAE_Cholesky(L, L_t, b);
}
vctr SLAE_Cholesky(Matrix& L, Matrix& L_t, vctr& b) {
    Matrix tmp = transpos(L_t);
    if (tmp != L) {
        cout << "    Error: wrong L/L_t for SLAE_Cholesky\n\n";
        return NULL;
    }
    tmp = L + b;
    size_t n = tmp.getRows();
    vctr y(0), x(0);
    for (size_t K = 0; K < n; K++) {
        tmp(K, n) /= tmp(K, K);                     // нижнетреугольная матрица, на диагонали НЕ единица - делим столбец b
        tmp(K, K) = 1;                               // нижнетреугольная матрица, на диагонали НЕ единица - делим
        do_zero_under_diag(tmp, K);                 // обнуление ниже главной диагонали
        y.vec.push_back(tmp(K, n));
    }
    tmp = L_t + y;
    do_zero_upper_diag_without_ones(tmp);                               // обнуление выше главной диагонали, на диагонали НЕ нули
    for (size_t i = 0; i < n; ++i) x.vec.push_back(tmp(i, n));
    return x;
}
vctr SLAE_Tomas(Matrix& A, vctr& b) {
    if (!A.is_tridiagonal() || b.vec.size() != A.getRows()) {
        cout << "    Error: matrix is not tridiagonal for Tomas or wrong vector b\n\n";
        return NULL;
    }
    size_t n = A.getRows();
    vector<double> alpha = { -A(0,1) / A(0,0) }, beta = { b.vec[0] / A(0,0) }, gamma = { A(0,0) }, x(n);
    for (int i = 1; i < n - 1; ++i) {
        gamma.push_back(A(i, i) + A(i, i-1) * alpha[i-1]);
        alpha.push_back(-A(i, i + 1) / gamma[i]);
        beta.push_back((b.vec[i] - A(i, i-1) * beta[i-1]) / gamma[i]);
    }
    gamma.push_back(A(n-1, n-1) + A(n-1, n-2) * alpha[n-2]);
    beta.push_back((b.vec[n-1] - A(n-1, n-2) * beta[n-2]) / gamma[n-1]);

    x[n-1] = beta[n-1];
    for (int i = n - 2; i >= 0; --i)
        x[i] = alpha[i] * x[i+1] + beta[i];
    vctr answer(x, 0);
    return answer;
}
vctr SLAE_QR(Matrix& mtr, vctr& b) {
    if (mtr.getCols() != mtr.getRows() || b.vec.size() != mtr.getRows()) {
        cout << "    Error: wrong size for SLAE_QR without vector b or wrong vector b\n\n";
        return NULL;
    }
    Matrix A = mtr, Q, R;
    if (!A.find_QR_matrix(Q, R)) return NULL;
    else return SLAE_QR(Q, R, b);
}
vctr SLAE_QR(Matrix& Q, Matrix& R, vctr& b) {
    if (Q.getCols() != Q.getRows() || R.getCols() != R.getRows() || Q.getCols() != R.getCols() || b.vec.size() != Q.getRows()) {
        cout << "    Error: wrong size of Q or R or b for SLAE_QR\n\n";
        return NULL;
    }
    Matrix tmp_Q = transpos(Q);
    vctr tmp_b(0), x(0);
    size_t n = Q.getRows();
    tmp_Q = tmp_Q * b;
    for (size_t K = 0; K < n; K++) {
        tmp_b.vec.push_back(tmp_Q(K, 0));
    }
    Matrix tmp_R = R + tmp_b;
    do_zero_upper_diag_without_ones(tmp_R);                               // обнуление выше главной диагонали, на диагонали НЕ нули
    for (size_t i = 0; i < n; ++i) x.vec.push_back(tmp_R(i, n));
    return x;
}
vctr SLAE_Iteration(Matrix& A, vctr& b, double accuracy, vctr& initial_approximation) {
    if (A.getCols() != A.getRows() || b.vec.size() != A.getRows() || initial_approximation.vec.size() != A.getRows()) {
        cout << "    Error: wrong size of matrix A or vector b for SLAE_Iteration\n\n";
        return NULL;
    }
    int n = A.getRows(), k = 0;
    Matrix c(n,1);
    for (size_t i = 0; i < n; ++i) {
        if (!A(i, i)) {
            cout << "    Zero element on the main diagonal, SLAE_Iteration is not possible\n\n";
            return NULL;
        }
        else c(i, 0) = b.vec[i] / A(i, i);
    }
    Matrix B = B_matrix(A), x(initial_approximation), x_next = B * x;
    double e1 = (1 - B.norm_first()) * accuracy / B.norm_first();
    vctr ans;
    cout << "\n\n-----Нормы матрицы B (first, Evkl, Infinity):    " << B.norm_first() << "    " << B.norm_evkl() << "    " << B.norm_infinity();
    cout << "\n----------Начальное приближение:    " << initial_approximation << "\n";

    x_next = x_next + c;
    while ((x_next - x).norm_first() >= e1) {
        x = x_next;
        x_next = B * x;
        x_next = x_next + c;
        cout << "Итерация " << ++k << ":  " << (ans = x_next.to_vector());
    }
    return ans;
}
vctr SLAE_Seidel(Matrix& A, vctr& b, double accuracy, vctr& initial_approximation) {
    if (A.getCols() != A.getRows() || b.vec.size() != A.getRows() || initial_approximation.vec.size() != A.getRows()) {
        cout << "    Error: wrong size of matrix A or vector b for SLAE_Iteration\n\n";
        return NULL;
    }
    int n = A.getRows(), k = 0;
    Matrix c(n, 1);
    for (size_t i = 0; i < n; ++i) {
        if (!A(i, i)) {
            cout << "    Zero element on the main diagonal, SLAE_Iteration is not possible\n\n";
            return NULL;
        }
        else c(i, 0) = b.vec[i] / A(i, i);
    }
    Matrix B1, B2;
    B1_B2_matrix(A, B1, B2);
    double e2 = (1 - (B1 + B2).norm_first()) * accuracy / B2.norm_first();
    vctr x(initial_approximation), x_next(initial_approximation);
    cout << "\n\n-----Нормы матрицы B1 (first, Evkl, Infinity):    " << B1.norm_first() << "    " << B1.norm_evkl() << "    " << B1.norm_infinity();
    cout << "\n-----Нормы матрицы B2 (first, Evkl, Infinity):    " << B2.norm_first() << "    " << B2.norm_evkl() << "    " << B2.norm_infinity();
    cout << "\n----------Начальное приближение:    " << initial_approximation << "\n";

    for (int i = 0; i < n; ++i) {
        x_next.vec[i] = 0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (i < j) x_next.vec[i] += B2(i, j) * x.vec[j];
            if (i > j) x_next.vec[i] += B1(i, j) * x_next.vec[j];            
        }
        x_next.vec[i] += c(i, 0);
    }
    cout << "Итерация " << ++k << ":  " << x_next;
    while (norm_first(x = (x_next - x)) >= e2) {
        x = x_next;
        for (int i = 0; i < n; ++i) { //x(i,0)
            x_next.vec[i] = 0;
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                if (i < j) x_next.vec[i] += B2(i, j) * x.vec[j];
                if (i > j) x_next.vec[i] += B1(i, j) * x_next.vec[j];
            }
            x_next.vec[i] += c(i, 0);
        }
        cout << "Итерация " << ++k << ":  " << x_next;
    }
    return x_next;
}

double power_law_method(Matrix& A, double accuracy, vctr& x_next) {
    if (A.getCols() != A.getRows()) {
        cout << "    Error: wrong size of matrix A or vector b for SLAE_Iteration\n\n";
        return NULL;
    }
    vctr x(false);
    double lambda = 1, lambda_next;
    int k = 0;
    for (size_t i = 0; i < A.getRows(); ++i) x.vec.push_back(1);

    x_next = (A * x).to_vector();
    lambda_next = (x_next * x)/(x*x);
    cout << "Итерация " << ++k << "\nСобственное число : " << lambda_next << "\nВектор: " << x_next;
    while (fabs(lambda_next-lambda) > accuracy) {
        x = x_next;
        lambda = lambda_next;
        x_next = (A * x).to_vector();
        lambda_next = (x_next * x) / (x * x);
        cout << "Итерация " << ++k << "\nСобственное число : " << lambda_next << "\nВектор: " << x_next;
    }
    return lambda_next;
}
double power_law_method_with_normalization(Matrix& A, double accuracy, vctr& x_next) {
    if (A.getCols() != A.getRows()) {
        cout << "    Error: wrong size of matrix A or vector b for SLAE_Iteration\n\n";
        return NULL;
    }
    vctr x(false);
    double lambda = 1, lambda_next;
    int k = 0;
    for (size_t i = 0; i < A.getRows(); ++i) x.vec.push_back(1);

    x_next = (A * x).to_vector();
    lambda_next = x_next * x;
    x = x_next / norm_second(x_next);
    cout << "Итерация " << ++k << "\nСобственное число : " << lambda_next << "\nВектор: " << x_next;
    while (fabs(lambda_next - lambda) > accuracy) {
        lambda = lambda_next;
        x_next = (A * x).to_vector();
        lambda_next = x_next * x;
        x = x_next / norm_second(x_next);
        cout << "Итерация " << ++k << "\nСобственное число : " << lambda_next << "\nВектор: " << x_next;
    }
    x_next = x;
    return lambda_next;
}


//------------------------------Переопределение операторов------------------------------
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