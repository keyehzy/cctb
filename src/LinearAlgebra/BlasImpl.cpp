#include "LinearAlgebra/BlasImpl.h"

#if HAVE_BLAS
#include "cblas.h"
#endif

void scal(double alpha, NumericArray<double>& x) {
#if HAVE_BLAS
  cblas_dscal(x.size(), alpha, x.buffer().data(), 1);
#else
  for (int i = 0; i < x.size(); ++i) {
    x[i] *= alpha;
  }
#endif
}

void copy(const NumericArray<double>& x, NumericArray<double>& y) {
#if HAVE_BLAS
  cblas_dcopy(x.size(), x.buffer().data(), 1, y.buffer().data(), 1);
#else
  for (int i = 0; i < x.size(); ++i) {
    y[i] = x[i];
  }
#endif
}

void axpy(double alpha, const NumericArray<double>& x, NumericArray<double>& y) {
#if HAVE_BLAS
  cblas_daxpy(x.size(), alpha, x.buffer().data(), 1, y.buffer().data(), 1);
#else
  for (int i = 0; i < x.size(); ++i) {
    y[i] += alpha * x[i];
  }
#endif
}

void dot(const NumericArray<double>& x, const NumericArray<double>& y, double* result) {
#if HAVE_BLAS
  *result = cblas_ddot(x.size(), x.buffer().data(), 1, y.buffer().data(), 1);
#else
  *result = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    *result += x[i] * y[i];
  }
#endif
}

void nrm2(const NumericArray<double>& x, double* result) {
#if HAVE_BLAS
  *result = cblas_dnrm2(x.size(), x.buffer().data(), 1);
#else
  *result = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    *result += x[i] * x[i];
  }
  *result = sqrt(*result);
#endif
}

void asum(const NumericArray<double>& x, double* result) {
#if HAVE_BLAS
  *result = cblas_dasum(x.size(), x.buffer().data(), 1);
#else
  *result = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    *result += fabs(x[i]);
  }
#endif
}

void iamax(const NumericArray<double>& x, int* result) {
#if HAVE_BLAS
  *result = cblas_idamax(x.size(), x.buffer().data(), 1);
#else
  *result = 0;
  double max = std::abs(x[0]);
  for (int i = 1; i < x.size(); ++i) {
    if (std::abs(x[i]) > max) {
      max = std::abs(x[i]);
      *result = i;
    }
  }
#endif
}

void scal(std::complex<double> alpha, NumericArray<std::complex<double>>& x) {
#if HAVE_BLAS
  cblas_zscal(x.size(), &alpha, x.buffer().data(), 1);
#else
  for (int i = 0; i < x.size(); ++i) {
    x[i] *= alpha;
  }
#endif
}

void copy(const NumericArray<std::complex<double>>& x, NumericArray<std::complex<double>>& y) {
#if HAVE_BLAS
  cblas_zcopy(x.size(), x.buffer().data(), 1, y.buffer().data(), 1);
#else
  for (int i = 0; i < x.size(); ++i) {
    y[i] = x[i];
  }
#endif
}

void axpy(std::complex<double> alpha, const NumericArray<std::complex<double>>& x,
          NumericArray<std::complex<double>>& y) {
#if HAVE_BLAS
  cblas_zaxpy(x.size(), &alpha, x.buffer().data(), 1, y.buffer().data(), 1);
#else
  for (int i = 0; i < x.size(); ++i) {
    y[i] += alpha * x[i];
  }
#endif
}

void dot(const NumericArray<std::complex<double>>& x, const NumericArray<std::complex<double>>& y,
         std::complex<double>* result) {
#if HAVE_BLAS
  *result = cblas_zdotc(x.size(), x.buffer().data(), 1, y.buffer().data(), 1);
#else
  *result = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    *result += x[i] * std::conj(y[i]);
  }
#endif
}

void asum(const NumericArray<std::complex<double>>& x, std::complex<double>* result) {
#if HAVE_BLAS
  *result = cblas_dzasum(x.size(), x.buffer().data(), 1);
#else
  *result = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    *result += std::abs(x[i]);
  }
#endif
}

void iamax(const NumericArray<std::complex<double>>& x, int* result) {
#if HAVE_BLAS
  *result = cblas_izamax(x.size(), x.buffer().data(), 1);
#else
  *result = 0;
  std::complex<double> max = std::abs(x[0]);
  for (int i = 1; i < x.size(); ++i) {
    if (std::abs(x[i]) > max) {
      max = std::abs(x[i]);
      *result = i;
    }
  }
#endif
}

void gemv(double alpha, const Matrix<double>& A, const NumericArray<double>& x, double beta,
          NumericArray<double>& y) {
#if HAVE_BLAS
  cblas_dgemv(CblasRowMajor, CblasNoTrans, A.rows(), A.cols(), alpha, A.buffer().data(), A.cols(),
              x.buffer().data(), 1, beta, y.buffer().data(), 1);
#else
  for (int i = 0; i < A.rows(); ++i) {
    double result = 0;
    for (int j = 0; j < A.cols(); ++j) {
      result += A(i, j) * x[j];
    }
    y[i] = alpha * result + beta * y[i];
  }
#endif
}

void ger(double alpha, const NumericArray<double>& x, const NumericArray<double>& y,
         Matrix<double>& A) {
#if HAVE_BLAS
  cblas_dger(CblasRowMajor, A.rows(), A.cols(), alpha, x.buffer().data(), 1, y.buffer().data(), 1,
             A.buffer().data(), A.cols());
#else
  for (int i = 0; i < A.rows(); ++i) {
    for (int j = 0; j < A.cols(); ++j) {
      A(i, j) += alpha * x[i] * y[j];
    }
  }
#endif
}

void gemm(double alpha, const Matrix<double>& a, const Matrix<double>& b, double beta,
          Matrix<double>& c) {
#if HAVE_BLAS
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(), a.cols(), alpha,
              a.buffer().data(), a.cols(), b.buffer().data(), b.cols(), beta, c.buffer().data(),
              c.cols());
#else
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < b.cols(); ++j) {
      double result = 0;
      for (int k = 0; k < a.cols(); ++k) {
        result += a(i, k) * b(k, j);
      }
      c(i, j) = alpha * result + beta * c(i, j);
    }
  }
#endif
}

void gemv(std::complex<double> alpha, const Matrix<std::complex<double>>& A,
          const NumericArray<std::complex<double>>& x, std::complex<double> beta,
          NumericArray<std::complex<double>>& y) {
#if HAVE_BLAS
  cblas_zgemv(CblasRowMajor, CblasNoTrans, A.rows(), A.cols(), &alpha, A.buffer().data(), A.cols(),
              x.buffer().data(), 1, &beta, y.buffer().data(), 1);
#else
  for (int i = 0; i < A.rows(); ++i) {
    std::complex<double> result = 0;
    for (int j = 0; j < A.cols(); ++j) {
      result += A(i, j) * x[j];
    }
    y[i] = alpha * result + beta * y[i];
  }
#endif
}

void ger(std::complex<double> alpha, const NumericArray<std::complex<double>>& x,
         const NumericArray<std::complex<double>>& y, Matrix<std::complex<double>>& A) {
#if HAVE_BLAS
  cblas_zgerc(CblasRowMajor, A.rows(), A.cols(), &alpha, x.buffer().data(), 1, y.buffer().data(), 1,
              A.buffer().data(), A.cols());
#else
  for (int i = 0; i < A.rows(); ++i) {
    for (int j = 0; j < A.cols(); ++j) {
      A(i, j) += alpha * x[i] * std::conj(y[j]);
    }
  }
#endif
}

void gemm(std::complex<double> alpha, const Matrix<std::complex<double>>& a,
          const Matrix<std::complex<double>>& b, std::complex<double> beta,
          Matrix<std::complex<double>>& c) {
#if HAVE_BLAS
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(), a.cols(), &alpha,
              a.buffer().data(), a.cols(), b.buffer().data(), b.cols(), &beta, c.buffer().data(),
              c.cols());
#else
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < b.cols(); ++j) {
      std::complex<double> result = 0;
      for (int k = 0; k < a.cols(); ++k) {
        result += a(i, k) * b(k, j);
      }
      c(i, j) = alpha * result + beta * c(i, j);
    }
  }
#endif
}
