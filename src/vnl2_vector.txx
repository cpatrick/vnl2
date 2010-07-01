#ifndef VNL2_VECTOR_TXX_
#define VNL2_VECTOR_TXX_

/**
 * \file vector.txx
 * \author Chuck Atkins
 *
 * \verbatim
 *   The implementation of the core data structure used for numerical vectors.
 * \endverbatim
 */
#include <cstring>
#include <stdexcept>

#include "vnl2_vector.h"
#include "vnl2_types.h"
#include "vnl2_blas.h"

namespace vnl2
{

/** Creates a zero filled vector 
 * \param len The number of elements
 */
template<typename T>
vector<T>
::vector(size_t len)
: m_shallow(false), m_len(len), m_data(NULL)
{
  if( this->m_len > 0 )
  {
    this->m_data = new T[this->m_len];
    std::memset(this->m_data, 0, sizeof(T)*this->m_len);
  }
}

  
/** Create a vector by filling every element with a specified value 
 * \param len   The number of elements
 * \param value Initialize all elements of the vector to value
 */
template<typename T>
vector<T>
::vector(size_t len, const T value)
: m_shallow(false), m_len(len), m_data(NULL)
{
  if( this->m_len > 0 )
  {
    this->m_data = new T[this->m_len];
    for( size_t i = 0; i < this->m_len; ++i )
    {
      this->m_data[i] = value;
    }
  }
}
 
 
/** Copy constructor from a fixed array 
 * \param v           The source array to copy from
 * \param shallowCopy If true, then only the underlying memory pointer
 *                    is copied, not the actual vector elements
 */
template<typename T>
vector<T>
::vector(const T v[], bool shallowCopy )
: m_shallow(shallowCopy), m_len(sizeof(v)/sizeof(T)), m_data(NULL)
{
  if( this->m_shallow )
  {
    //this->m_data = const_cast<T*>(v);
  }
  else if( this->m_len > 0 )
  {
    this->m_data = new T[this->m_len];
    std::memcpy(this->m_data, v, sizeof(v));
  }
}

  
/** Copy constructor from a dynamic array 
 * \param len         The number of elements
 * \param v           The source array to copy from
 * \param shallowCopy If true, then only the underlying memory pointer
 *                    is copied, not the actual vector elements
 */
template<typename T>
vector<T>
::vector(size_t len, const T* v, bool shallowCopy)
: m_shallow(shallowCopy), m_len(len), m_data(NULL)
{
  if( this->m_shallow )
  {
    this->m_data = const_cast<T*>(v);
  }
  else if( this->m_len > 0 )
  {
    m_data = new T[this->m_len];
    std::memcpy(this->m_data, v, sizeof(T)*this->m_len);
  }
}

  
/** Copy constructor from an vnl2::vector 
 * \param v           The source array to copy from
 * \param shallowCopy If true, then only the underlying memory pointer
 *                    is copied, not the actual vector elements
 */
template<typename T>
vector<T>
::vector(const vector<T>& v, bool shallowCopy)
: m_shallow(shallowCopy), m_len(v.m_len), m_data(NULL)
{
  if( this == &v )
  {
    return;
  }
  if( this->m_shallow )
  {
    this->m_data = v.m_data;
  }
  else if( this->m_len > 0 )
  {
    std::memcpy(this->m_data, v.m_data, sizeof(T)*this->m_len);
  }
}  


/** Copy constructor from an std::vector 
 * \param v           The source array to copy from
 * \param shallowCopy If true, then only the underlying memory pointer
 *                    is copied, not the actual vector elements
 */
template<typename T>
vector<T>
::vector(const std::vector<T>& v, bool shallowCopy) 
: m_shallow(shallowCopy), m_len(v.size()), m_data(NULL)
{
  if( this->m_shallow )
  {
    this->m_data = const_cast<T*>(&(*(v.begin())));
  }
  else if( this->m_len > 0 )
  {
    this->m_data = new T[this->m_len];
    std::memcpy(m_data, &(*(v.begin())), sizeof(T)*this->m_len);
  }
}  


/** Destructor 
 *  The shallow copy mechanism will only free memory if allocated by the
 *  vector class.  This allows various vector data sources to be thinly
 *  wrapped without copying the vector elements.
 */
template<typename T>
vector<T>
::~vector(void)
{
  if( !this->m_shallow && this->m_data )
  {
    delete[] this->m_data;
    this->m_data = NULL;
    this->m_len = 0;
  }
}


/** Assignment operator */
template<typename T>
const vector<T>& 
vector<T>
::operator=(const vector<T>& v)
{
  if( &v != this )
  {
    if( this->m_len != v.m_len )
    {
      if( !this->m_data )
      {
        delete[] this->m_data;
      }
      this->m_len = v.m_len;
      this->m_data = new T[this->m_len];
    }
    std::memcpy(this->m_data, v.m_data, sizeof(T)*this->m_len);
  }
  return *this;
}


/** Assignment operator */
template<typename T>
vector<T>& 
vector<T>
::operator=(vector<T>& v)
{
  if( &v != this )
  {
    if( this->m_len != v.m_len )
    {
      if( !this->m_data )
      {
        delete[] this->m_data;
      }
      this->m_len = v.m_len;
      this->m_data = new T[this->m_len];
    }
    std::memcpy(this->m_data, v.m_data, sizeof(T)*this->m_len);
  }
  return *this;
}


/** Indexing operator */
template<typename T>
const T& 
vector<T>
::operator[](size_t i) const
{
  if( i > this->m_len )
  {
    throw std::out_of_range("vnl2::vector: index out of range");
  }
  return this->m_data[i];
}


/** Indexing operator */
template<typename T>
T& 
vector<T>
::operator[](size_t i)
{
  if( i > this->m_len )
  {
    throw std::out_of_range("vnl2::vector: index out of range");
  }
  return this->m_data[i];
}


//------------------------------------------------------------
// Scalar Multiplication
//------------------------------------------------------------

template<>
template<>
vector<float> 
vector<float>
::operator  *(const float& alpha) const
{
  vector<float> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  sscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<double> 
vector<double>
::operator  *(const double& alpha) const
{
  vector<double> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  dscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<complex_float> 
vector<complex_float>
::operator  *(const float& alpha) const
{
  vector<complex_float> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  csscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<complex_float> 
vector<complex_float>
::operator  *(const complex_float& alpha) const
{
  vector<complex_float> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  cscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<complex_double> 
vector<complex_double>
::operator  *(const double& alpha) const
{
  vector<complex_double> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  zdscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<complex_double> 
vector<complex_double>
::operator  *(const complex_double& alpha) const
{
  vector<complex_double> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  zscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<float>& 
vector<float>
::operator *=(const float& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  sscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<double>& 
vector<double>
::operator *=(const double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  dscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<complex_float>& 
vector<complex_float>
::operator *=(const float& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  csscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<complex_float>& 
vector<complex_float>
::operator *=(const complex_float& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  cscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<complex_double>& 
vector<complex_double>
::operator *=(const complex_double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  zscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<complex_double>& 
vector<complex_double>
::operator *=(const double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  zdscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


//------------------------------------------------------------
// Scalar Division
//------------------------------------------------------------


template<>
template<>
vector<float> 
vector<float>
::operator  /(const float& alpha) const
{
  vector<float> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  float alpha_inv = 1.0f/alpha;
  sscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<double> 
vector<double>
::operator  /(const double& alpha) const
{
  vector<double> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  double alpha_inv = 1.0/alpha;
  dscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<complex_float> 
vector<complex_float>
::operator  /(const float& alpha) const
{
  vector<complex_float> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  const float alpha_inv = 1.0f/alpha;
  csscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<complex_float> 
vector<complex_float>
::operator  /(const complex_float& alpha) const
{
  vector<complex_float> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  const complex_float alpha_inv = 1.0f/alpha;
  cscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<complex_double> 
vector<complex_double>
::operator  /(const double& alpha) const
{
  vector<complex_double> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  const double alpha_inv = 1.0/alpha;
  zdscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<complex_double> 
vector<complex_double>
::operator  /(const complex_double& alpha) const
{
  vector<complex_double> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  const complex_double alpha_inv = 1.0/alpha;
  zscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<float>& 
vector<float>
::operator /=(const float& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const float alpha_inv = 1.0f/alpha;
  sscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<double>& 
vector<double>
::operator /=(const double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const double alpha_inv = 1.0/alpha;
  dscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<complex_float>& 
vector<complex_float>
::operator /=(const float& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const float alpha_inv = 1.0f/alpha;
  csscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<complex_float>& 
vector<complex_float>
::operator /=(const complex_float& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const complex_float alpha_inv = 1.0f/alpha;
  cscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<complex_double>& 
vector<complex_double>
::operator /=(const double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const double alpha_inv = 1.0/alpha;
  zdscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<complex_double>& 
vector<complex_double>
::operator /=(const complex_double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const complex_double alpha_inv = 1.0/alpha;
  zscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


//------------------------------------------------------------
//  Vector Addition
//------------------------------------------------------------


template<>
vector<float> 
vector<float>
::operator  +(const vector<float>& x) const
{
  vector<float> y(*this);
  const int inc = 1;
  const int n = static_cast<int>(y.m_len);
  const float alpha = 1.0f;
  saxpy_(&n, &alpha, x.m_data, &inc, y.m_data, &inc);
}


template<>
vector<double> 
vector<double>
::operator  +(const vector<double>& x) const
{
  vector<double> y(*this);
  const int inc = 1;
  const int n = static_cast<int>(y.m_len);
  const double alpha = 1.0;
  daxpy_(&n, &alpha, x.m_data, &inc, y.m_data, &inc);
}


template<>
vector<complex_float> 
vector<complex_float>
::operator  +(const vector<complex_float>& x) const
{
  vector<complex_float> y(*this);
  const int inc = 1;
  const int n = static_cast<int>(y.m_len);
  const complex_float alpha(1.0f, 0.0f);
  caxpy_(&n, &alpha, x.m_data, &inc, y.m_data, &inc);
}


template<>
vector<complex_double> 
vector<complex_double>
::operator  +(const vector<complex_double>& x) const
{
  vector<complex_double> y(*this);
  const int inc = 1;
  const int n = static_cast<int>(y.m_len);
  const complex_double alpha(1.0, 0.0);
  zaxpy_(&n, &alpha, x.m_data, &inc, y.m_data, &inc);
}


}

#endif // VNL2_VECTOR_TXX_
