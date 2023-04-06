#include "LinearAlgebra/LapackImpl.h"

#include <lapacke.h>

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
