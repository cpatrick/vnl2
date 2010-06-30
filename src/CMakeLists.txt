enable_language(Fortran)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library( vnl2
  vector.h  vector.txx  vector_inst.cxx
)
target_link_libraries( vnl2 
                       ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES}
                       ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
)

