#include "LinearAlgebra/LapackImpl.h"

#include <lapacke.h>

#include <iostream>

#include "FloatPointHelpers.h"

using cmplx = std::complex<double>;

void geev(const Matrix<double>& a, NumericArray<std::complex<double>>& w,
          Matrix<std::complex<double>>& v) {
  int n = a.rows();
  int lda = a.rows();
  int ldvr = v.rows();
  int ldvl = v.rows();
  int info = 0;
  NumericArray<double> vl(ldvr * n);
  NumericArray<double> vr(ldvl * n);
  NumericArray<double> wr(n);
  NumericArray<double> wi(n);

  info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, a.buffer().data(), lda, wr.buffer().data(),
                       wi.buffer().data(), vl.buffer().data(), ldvl, vr.buffer().data(), ldvr);

  if (info != 0) __builtin_trap();

  for (int i = 0; i < n; ++i) {
    w[i] = cmplx{wr[i], wi[i]};
  }

  for (int i = 0; i < n; ++i) {
    if (fp_eq(wi[i], 0.0)) {
      for (int j = 0; j < n; ++j) {
        v(j, i) = vr[j + i * ldvr];
      }
    } else {
      for (int j = 0; j < n; ++j) {
        v(j, i) = cmplx{vr[j + i * ldvr], vr[j + (i + 1) * ldvr]};
        v(j, i + 1) = cmplx{vr[j + i * ldvr], -vr[j + (i + 1) * ldvr]};
      }
      ++i;
    }
  }
}

void geev(const Matrix<std::complex<double>>& a, NumericArray<std::complex<double>>& w,
          Matrix<std::complex<double>>& v) {
  int n = a.rows();
  int lda = a.rows();
  int ldvr = v.rows();
  int ldvl = v.rows();
  NumericArray<std::complex<double>> vl(ldvr * n);
  int info = 0;

  info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, (_Complex double*)a.buffer().data(), lda,
                       (_Complex double*)w.buffer().data(), (_Complex double*)vl.buffer().data(),
                       ldvl, (_Complex double*)v.buffer().data(), ldvr);

  if (info != 0) __builtin_trap();
}

void syev(const Matrix<double>& a, NumericArray<double>& w, Matrix<double>& v) {
  int n = a.rows();
  int lda = a.rows();
  int ldz = v.rows();
  int info = 0;
  int m = a.rows();
  ;

  NumericArray<int> isuppz(2 * n);

  info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', n, a.buffer().data(), lda, 0.0, 0.0, 0, 0,
                        0.0, &m, w.buffer().data(), v.buffer().data(), ldz, isuppz.buffer().data());

  if (info != 0) __builtin_trap();
}

void diagonalize_hermitian(const Matrix<std::complex<double>>& a, NumericArray<double>& w,
                           Matrix<std::complex<double>>& v) {
  int n = a.rows();
  int lda = a.rows();
  int ldz = v.rows();
  int info = 0;
  int m;

  NumericArray<int> isuppz(2 * n);

  info = LAPACKE_zheevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', n, (_Complex double*)a.buffer().data(),
                        lda, 0.0, 0.0, 0, 0, 0.0, &m, w.buffer().data(),
                        (_Complex double*)v.buffer().data(), ldz, isuppz.buffer().data());

  if (info != 0) __builtin_trap();
}
