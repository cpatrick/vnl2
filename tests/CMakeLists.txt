add_executable( test_vnl2_vector test_vnl2_vector.cxx )
target_link_libraries( test_vnl2_vector vnl2 )
add_test( test_vnl2_vector test_vnl2_vector )

