/**
 * \file vector_inst.cxx
 * \author Chuck Atkins
 *
 * \verbatim
 *   Instantiate the common instances of vnl2::vector
 * \endverbatim
 */

#include "vnl2_types.h"
#include "vnl2_vector.h"
#include "vnl2_vector.txx"

template class vnl2::vector<float>;
template class vnl2::vector<double>;
template class vnl2::vector<vnl2::complex_float>;
template class vnl2::vector<vnl2::complex_double>;

