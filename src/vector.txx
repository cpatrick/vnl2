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

#include "vector.h"
#include "blas.h"

namespace vnl2
{

/** Creates a zero filled vector 
 * \param len The number of elements
 */
template<typename T>
vector<T>
::vector(size_t len)
: m_len(len), m_data(NULL)
{
  if( this->m_len > 0 )
  {
    m_data = new T[this->m_len];
    std::memset(m_data, 0, sizeof(T)*this->m_len);
  }
}

  
/** Create a vector by filling every element with a specified value 
 * \param len   The number of elements
 * \param value Initialize all elements of the vector to value
 */
template<typename T>
vector<T>
::vector(size_t len, const T value)
: m_len(len), m_data(NULL)
{
  if( this->m_len > 0 )
  {
    m_data = new T[this->m_len];
    // FIX ME
    std::memset(m_data, 0, sizeof(T)*this->m_len);
  }
}
 
 
/** Copy constructor from a fixed array 
 * \param v   The source array to copy from
 */
template<typename T>
vector<T>
::vector(const T v[])
: m_len(sizeof(v)/sizeof(T)), m_data(NULL)
{
  if( this->m_len > 0 )
  {
    m_data = new T[this->m_len];
    std::memcpy(m_data, v, sizeof(v));
  }
}

  
/** Copy constructor from a dynamic array 
 * \param len The number of elements
 * \param v   The source array to copy from
 */
template<typename T>
vector<T>
::vector(size_t len, const T* v)
: m_len(len), m_data(NULL)
{
  if( this->m_len > 0 )
  {
    m_data = new T[this->m_len];
    std::memcpy(m_data, v, sizeof(T)*this->m_len);
  }
}

  
/** Copy constructor from an vnl2::vector 
 * \param v   The source array to copy from
 */
template<typename T>
vector<T>
::vector(const vector<T>& v)
: m_len(0), m_data(NULL)
{
  *this = v;
}  


/** Copy constructor from an std::vector 
 * \param v   The source array to copy from
 */
template<typename T>
vector<T>
::vector(const std::vector<T>& v) 
: m_len(v.size()), m_data(NULL)
{
  if( this->m_len > 0 )
  {
    m_data = new T[this->m_len];
    std::memcpy(m_data, &(*(v.begin())), sizeof(T)*this->m_len);
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
::operator*(const float& alpha)
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
::operator*(const double& alpha)
{
  vector<double> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  dscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<std::complex<float> > 
vector<std::complex<float> >
::operator*(const float& alpha)
{
  vector<std::complex<float> > r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  csscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<std::complex<float> > 
vector<std::complex<float> >
::operator*(const std::complex<float>& alpha)
{
  vector<std::complex<float> > r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  cscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<std::complex<double> > 
vector<std::complex<double> >
::operator*(const double& alpha)
{
  vector<std::complex<double> > r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  zdscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<std::complex<double> > 
vector<std::complex<double> >
::operator*(const std::complex<double>& alpha)
{
  vector<std::complex<double> > r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  zscal_(&n, &alpha, r.m_data, &incx);
}


template<>
template<>
vector<float>& 
vector<float>
::operator*=(const float& alpha)
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
::operator*=(const double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  dscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<std::complex<float> >& 
vector<std::complex<float> >
::operator*=(const float& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  csscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<std::complex<float> >& 
vector<std::complex<float> >
::operator*=(const std::complex<float>& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  cscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<std::complex<double> >& 
vector<std::complex<double> >
::operator*=(const std::complex<double>& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  zscal_(&n, &alpha, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<std::complex<double> >& 
vector<std::complex<double> >
::operator*=(const double& alpha)
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
::operator/(const float& alpha)
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
::operator/(const double& alpha)
{
  vector<double> r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  double alpha_inv = 1.0/alpha;
  dscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<std::complex<float> > 
vector<std::complex<float> >
::operator/(const float& alpha)
{
  vector<std::complex<float> > r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  const float alpha_inv = 1.0f/alpha;
  csscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<std::complex<float> > 
vector<std::complex<float> >
::operator/(const std::complex<float>& alpha)
{
  vector<std::complex<float> > r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  const std::complex<float> alpha_inv = 1.0f/alpha;
  cscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<std::complex<double> > 
vector<std::complex<double> >
::operator/(const double& alpha)
{
  vector<std::complex<double> > r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  const double alpha_inv = 1.0/alpha;
  zdscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<std::complex<double> > 
vector<std::complex<double> >
::operator/(const std::complex<double>& alpha)
{
  vector<std::complex<double> > r(*this);
  const int incx = 1;
  const int n = static_cast<int>(r.m_len);
  const std::complex<double> alpha_inv = 1.0/alpha;
  zscal_(&n, &alpha_inv, r.m_data, &incx);
}


template<>
template<>
vector<float>& 
vector<float>
::operator/=(const float& alpha)
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
::operator/=(const double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const double alpha_inv = 1.0/alpha;
  dscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<std::complex<float> >& 
vector<std::complex<float> >
::operator/=(const float& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const float alpha_inv = 1.0f/alpha;
  csscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<std::complex<float> >& 
vector<std::complex<float> >
::operator/=(const std::complex<float>& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const std::complex<float> alpha_inv = 1.0f/alpha;
  cscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<std::complex<double> >& 
vector<std::complex<double> >
::operator/=(const double& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const double alpha_inv = 1.0/alpha;
  zdscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


template<>
template<>
vector<std::complex<double> >& 
vector<std::complex<double> >
::operator/=(const std::complex<double>& alpha)
{
  const int incx = 1;
  const int n = static_cast<int>(this->m_len);
  const std::complex<double> alpha_inv = 1.0/alpha;
  zscal_(&n, &alpha_inv, this->m_data, &incx);
  return *this;
}


}

#endif // VNL2_VECTOR_TXX_
