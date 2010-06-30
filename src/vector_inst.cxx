/**
 * \file vector_inst.cxx
 * \author Chuck Atkins
 *
 * \verbatim
 *   Instantiate the common instances of vnl2::vector
 * \endverbatim
 */

#include <complex>

#include "vector.h"
#include "vector.txx"

template class vnl2::vector<float>;
template class vnl2::vector<double>;
template class vnl2::vector<std::complex<float> >;
template class vnl2::vector<std::complex<double> >;

