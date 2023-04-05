#include "LinearAlgebra/LapackImpl.h"

#if HAVE_LAPACK
#include <lapacke.h>
#endif

#define THROW() __builtin_unreachable()
#define NOT_IMPLEMENTED() __builtin_unreachable()

using cmplx = std::complex<double>;

void geev(const Matrix<double>& a, NumericArray<std::complex<double>>& w,
          Matrix<std::complex<double>>& v) {
#if HAVE_LAPACK
  int n = a.rows();
  int lda = a.rows();
  int ldvr = v.rows();
  int info = 0;
  NumericArray<double> vr(ldvr * n);
  NumericArray<double> wr(n);
  NumericArray<double> wi(n);

  info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, a.buffer().data(), lda, wr.buffer().data(),
                       wi.buffer().data(), nullptr, 1, vr.buffer().data(), ldvr);

  if (info != 0) THROW();

  for (int i = 0; i < n; ++i) {
    w[i] = cmplx{wr[i], wi[i]};
  }

  for (int i = 0; i < n; ++i) {
    if (std::abs(wi[i]) < 1e-6) {
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
#else
  NOT_IMPLEMENTED();
#endif
}
