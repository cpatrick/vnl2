enable_language(Fortran)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library( vnl2
  vnl2_types.h
  vnl2_blas.h
  vnl2_vector.h  vnl2_vector.txx  vnl2_vector_inst.cxx
)
target_link_libraries( vnl2 
                       ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES}
                       ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
)

