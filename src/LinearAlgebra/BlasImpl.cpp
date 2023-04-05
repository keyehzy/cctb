#include "LinearAlgebra/BlasImpl.h"

#include <cblas.h>

void scal(double alpha, NumericArray<double>& x) {
  cblas_dscal(x.size(), alpha, x.buffer().data(), 1);
}

void copy(const NumericArray<double>& x, NumericArray<double>& y) {
  cblas_dcopy(x.size(), x.buffer().data(), 1, y.buffer().data(), 1);
}

void axpy(double alpha, const NumericArray<double>& x, NumericArray<double>& y) {
  cblas_daxpy(x.size(), alpha, x.buffer().data(), 1, y.buffer().data(), 1);
}

void dot(const NumericArray<double>& x, const NumericArray<double>& y, double* result) {
  *result = cblas_ddot(x.size(), x.buffer().data(), 1, y.buffer().data(), 1);
}

void nrm2(const NumericArray<double>& x, double* result) {
  *result = cblas_dnrm2(x.size(), x.buffer().data(), 1);
}

void asum(const NumericArray<double>& x, double* result) {
  *result = cblas_dasum(x.size(), x.buffer().data(), 1);
}

void iamax(const NumericArray<double>& x, int* result) {
  *result = cblas_idamax(x.size(), x.buffer().data(), 1);
}

void scal(std::complex<double> alpha, NumericArray<std::complex<double>>& x) {
  cblas_zscal(x.size(), &alpha, x.buffer().data(), 1);
}

void copy(const NumericArray<std::complex<double>>& x, NumericArray<std::complex<double>>& y) {
  cblas_zcopy(x.size(), x.buffer().data(), 1, y.buffer().data(), 1);
}

void axpy(std::complex<double> alpha, const NumericArray<std::complex<double>>& x,
          NumericArray<std::complex<double>>& y) {
  cblas_zaxpy(x.size(), &alpha, x.buffer().data(), 1, y.buffer().data(), 1);
}

void dot(const NumericArray<std::complex<double>>& x, const NumericArray<std::complex<double>>& y,
         std::complex<double>* result) {
  *result = cblas_zdotc(x.size(), x.buffer().data(), 1, y.buffer().data(), 1);
}

void asum(const NumericArray<std::complex<double>>& x, std::complex<double>* result) {
  *result = cblas_dzasum(x.size(), x.buffer().data(), 1);
}

void iamax(const NumericArray<std::complex<double>>& x, int* result) {
  *result = cblas_izamax(x.size(), x.buffer().data(), 1);
}

void gemv(double alpha, const Matrix<double>& A, const NumericArray<double>& x, double beta,
          NumericArray<double>& y) {
  cblas_dgemv(CblasRowMajor, CblasNoTrans, A.rows(), A.cols(), alpha, A.buffer().data(), A.cols(),
              x.buffer().data(), 1, beta, y.buffer().data(), 1);
}

void ger(double alpha, const NumericArray<double>& x, const NumericArray<double>& y,
         Matrix<double>& A) {
  cblas_dger(CblasRowMajor, A.rows(), A.cols(), alpha, x.buffer().data(), 1, y.buffer().data(), 1,
             A.buffer().data(), A.cols());
}

void gemm(double alpha, const Matrix<double>& a, const Matrix<double>& b, double beta,
          Matrix<double>& c) {
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(), a.cols(), alpha,
              a.buffer().data(), a.cols(), b.buffer().data(), b.cols(), beta, c.buffer().data(),
              c.cols());
}

void gemv(std::complex<double> alpha, const Matrix<std::complex<double>>& A,
          const NumericArray<std::complex<double>>& x, std::complex<double> beta,
          NumericArray<std::complex<double>>& y) {
  cblas_zgemv(CblasRowMajor, CblasNoTrans, A.rows(), A.cols(), &alpha, A.buffer().data(), A.cols(),
              x.buffer().data(), 1, &beta, y.buffer().data(), 1);
}

void ger(std::complex<double> alpha, const NumericArray<std::complex<double>>& x,
         const NumericArray<std::complex<double>>& y, Matrix<std::complex<double>>& A) {
  cblas_zgerc(CblasRowMajor, A.rows(), A.cols(), &alpha, x.buffer().data(), 1, y.buffer().data(), 1,
              A.buffer().data(), A.cols());
}

void gemm(std::complex<double> alpha, const Matrix<std::complex<double>>& a,
          const Matrix<std::complex<double>>& b, std::complex<double> beta,
          Matrix<std::complex<double>>& c) {
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(), a.cols(), &alpha,
              a.buffer().data(), a.cols(), b.buffer().data(), b.cols(), &beta, c.buffer().data(),
              c.cols());
}
