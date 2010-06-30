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

/** 
 * \class vector
 * \brief A representation of numerical vectors
 *
 * For small sized vectors known at compile time (2, 3, and 4 elements) use
 * vector_fixed
 *
 * \tparam T The underlying vector type (double, float, compelx, etc.)
 */
template<typename T>
class vector
{
public:
  //----------------------------------------
  // Constructors
  //----------------------------------------
  
  /** Creates a zero filled vector 
   * \param len The number of elements
   */
  vector(size_t len = 0);
  
  /** Create a vector by filling every element with a specified value 
   * \param len   The number of elements
   * \param value Initialize all elements of the vector to value
   */
  vector(size_t len, const T value);
  
  /** Copy constructor from a fixed array 
   * \param v   The source array to copy from
   */
  vector(const T v[]);
  
  /** Copy constructor from a dynamic array 
   * \param len The number of elements
   * \param v   The source array to copy from.  An assumption is made that
   *            v contains at least len elements.
   */
  vector(size_t len, const T* v);
  
  /** Copy constructor from an vnl2::vector 
   * \param v   The source array to copy from
   */
  vector(const vector<T>& v);

  /** Copy constructor from an std::vector 
   * \param v   The source array to copy from
   */
  vector(const std::vector<T>& v);  

  /** Destructor */
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
  template<typename T2> vector<T>  operator  *(const T2& alpha);  
  template<typename T2> vector<T>& operator *=(const T2& alpha);  

  /** Scalar division */
  template<typename T2> vector<T>  operator  /(const T2& alpha);  
  template<typename T2> vector<T>& operator /=(const T2& alpha);  

private:
  size_t m_len;
  T* m_data;
};

}

#endif // VNL2_VECTOR_H_

