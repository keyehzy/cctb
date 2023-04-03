#include "LinearAlgebra/Matrix.h"

#if HAVE_BLAS
#include "cblas.h"
#endif
void sgemv(float alpha, const Matrix<float>& A, const NumericArray<float>& x, float beta,
           NumericArray<float>& y) {
#if HAVE_BLAS
  cblas_sgemv(CblasRowMajor, CblasNoTrans, A.rows(), A.cols(), alpha, A.buffer().data(), A.cols(),
              x.buffer().data(), 1, beta, y.buffer().data(), 1);
#else
  for (int i = 0; i < A.rows(); ++i) {
    float result = 0;
    for (int j = 0; j < A.cols(); ++j) {
      result += A(i, j) * x[j];
    }
    y[i] = alpha * result + beta * y[i];
  }
#endif
}

void sger(float alpha, const NumericArray<float>& x, const NumericArray<float>& y,
          Matrix<float>& A) {
#if HAVE_BLAS
  cblas_sger(CblasRowMajor, A.rows(), A.cols(), alpha, x.buffer().data(), 1, y.buffer().data(), 1,
             A.buffer().data(), A.cols());
#else
  for (int i = 0; i < A.rows(); ++i) {
    for (int j = 0; j < A.cols(); ++j) {
      A(i, j) += alpha * x[i] * y[j];
    }
  }
#endif
}

void sgemm(float alpha, const Matrix<float>& a, const Matrix<float>& b, float beta,
           Matrix<float>& c) {
#if HAVE_BLAS
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(), a.cols(), alpha,
              a.buffer().data(), a.cols(), b.buffer().data(), b.cols(), beta, c.buffer().data(),
              c.cols());
#else
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < b.cols(); ++j) {
      float result = 0;
      for (int k = 0; k < a.cols(); ++k) {
        result += a(i, k) * b(k, j);
      }
      c(i, j) = alpha * result + beta * c(i, j);
    }
  }
#endif
}

void dgemv(double alpha, const Matrix<double>& A, const NumericArray<double>& x, double beta,
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

void dger(double alpha, const NumericArray<double>& x, const NumericArray<double>& y,
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

void dgemm(double alpha, const Matrix<double>& a, const Matrix<double>& b, double beta,
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
