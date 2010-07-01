#ifndef VNL2_VECTOR_H_
#define VNL2_VECTOR_H_

/**
 * \file vector.h
 * \author Chuck Atkins
 *
 * \verbatim
 *   The core data structure used for numerical vectors.
 * \endverbatim
 */

#include <cstddef>
#include <vector> 

namespace vnl2
{

// Forward declaration of the matrix class for "friendliness"
template<typename T>
class matrix;

/** 
 * \class vector
 * \brief A representation of numerical vectors
 *
 * The vector class is a fundabental data type for numerical methods.  In
 * this implementation, the underlying operators are to be performed by the
 * Level 1 BLAS routines.  There is also a shallow caopy mechanism which
 * allows data from other classes and sources to be thinly wrapped without
 * the need to copy unnecessarily.
 *
 * \tparam T The underlying vector type (double, float, complex, etc.)
 */
template<typename T>
class vector
{
public:

  // Allow matrix acces to the vector's data members for optimized operations
  friend class matrix<T>;

  //----------------------------------------
  // Constructors
  //----------------------------------------
  
  /** Creates a zero filled vector 
   * \param len The number of elements
   */
  vector(size_t len = 0);
  
  /** Create a vector by filling every element with a specified value 
   * \param len         The number of elements
   * \param value       Initialize all elements of the vector to value
   */
  vector(size_t len, const T value);
  
  /** Copy constructor from a fixed array 
   * \param v           The source array to copy from
   * \param shallowCopy If true, then only the underlying memory pointer
   *                    is copied, not the actual vector elements
   */
  vector(const T v[], bool shallowCopy = false);
  
  /** Copy constructor from a dynamic array 
   * \param len The number of elements
   * \param v           The source array to copy from.  An assumption is 
   *                    made that contains at least len elements.
   * \param shallowCopy If true, then only the underlying memory pointer
   *                    is copied, not the actual vector elements
   */
  vector(size_t len, const T* v, bool shallowCopy = false);
  
  /** Copy constructor from an vnl2::vector 
   * \param v           The source array to copy from
   * \param shallowCopy If true, then only the underlying memory pointer
   *                    is copied, not the actual vector elements
   */
  vector(const vector<T>& v, bool shallowCopy = false);

  /** Copy constructor from an std::vector 
   * \param v           The source array to copy from
   * \param shallowCopy If true, then only the underlying memory pointer
   *                    is copied, not the actual vector elements
   */
  vector(const std::vector<T>& v, bool shallowCopy = false);  

  /** Destructor 
   *  The shallow copy mechanism will only free memory if allocated by the
   *  vector class.  This allows various vector data sources to be thinly
   *  wrapped without copying the vector elements.
   */
  ~vector(void);

  //----------------------------------------
  // Accessors
  //----------------------------------------

  /** Get the number of elements in the vector */
  size_t len(void) const {  return this->m_len;  }

  //----------------------------------------
  // Operators
  //----------------------------------------
  
  /** Assignment operator */
  const vector<T>& operator=(const vector<T>& v);

  /** Assignment operator */
  vector<T>& operator=(vector<T>& v);

  /** Indexing operator */
  const T& operator[](size_t i) const;

  /** Indexing operator */
  T& operator[](size_t i);

  /** Scalar multiplication */
  template<typename T2> vector<T>  operator  *(const T2& alpha) const;  
  template<typename T2> vector<T>& operator *=(const T2& alpha);  

  /** Scalar division */
  template<typename T2> vector<T>  operator  /(const T2& alpha) const;  
  template<typename T2> vector<T>& operator /=(const T2& alpha);  

  /** Vector addition */
  vector<T>  operator  +(const vector<T>& x) const;  
  vector<T>& operator +=(const vector<T>& x);  

private:
  bool m_shallow;
  size_t m_len;
  T* m_data;
};

}

#endif // VNL2_VECTOR_H_

