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

#include "vnl2_types.h"

extern "C"
{

// alpha * x
extern void sscal_( const int *n, const float *alpha, 
                    float *x, const int *incx);
extern void dscal_( const int *n, const double *alpha,
                    double *x, const int *incx);
extern void csscal_(const int *n, const float *alpha, 
                    vnl2::complex_float *x, const int *incx);
extern void cscal_( const int *n, const vnl2::complex_float *alpha, 
                    vnl2::complex_float *x, const int *incx);
extern void zdscal_(const int *n, const double *alpha, 
                    vnl2::complex_double *x, const int *incx);
extern void zscal_( const int *n, const vnl2::complex_double *alpha, 
                    vnl2::complex_double *x, const int *incx);
 
// alpha * x + y 
extern void saxpy_(const int *n, const float *alpha,
                   const float *x, const int *incx,
                         float *y, const int *incy);
extern void daxpy_(const int *n, const double *alpha,
                   const double *x, const int *incx,
                         double *y, const int *incy);
extern void caxpy_(const int *n, const vnl2::complex_float *alpha,
                   const vnl2::complex_float *x, const int *incx,
                         vnl2::complex_float *y, const int *incy);
extern void zaxpy_(const int *n, const vnl2::complex_double *alpha,
                   const vnl2::complex_double *x, const int *incx,
                         vnl2::complex_double *y, const int *incy);
}

#endif
  
