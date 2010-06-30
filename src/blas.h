#ifndef VNL2_BLAS_H_
#define VNL2_BLAS_H_

/**
 * \file blas.h
 * \author Chuck Atkins
 *
 * \verbatim
 *   This provides an interface to the fortran routines in BLAS
 * \endverbatim
 */

#include <complex>

extern "C"
{

// Scalar multiplication
extern void sscal_( const int *n, const float *alpha, 
                    float *x, const int *incx);
extern void dscal_( const int *n, const double *alpha,
                    double *x, const int *incx);
extern void csscal_(const int *n, const float *alpha, 
                    std::complex<float> *x, const int *incx);
extern void cscal_( const int *n, const std::complex<float> *alpha, 
                    std::complex<float> *x, const int *incx);
extern void zdscal_(const int *n, const double *alpha, 
                    std::complex<double> *x, const int *incx);
extern void zscal_( const int *n, const std::complex<double> *alpha, 
                    std::complex<double> *x, const int *incx);
  
}

#endif
  
