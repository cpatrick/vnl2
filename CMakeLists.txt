project(vnl2)

cmake_minimum_required( VERSION 2.8 )

if( NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY )
  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${vnl2_BINARY_DIR}/lib" )
endif()
if( NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY )
  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${vnl2_BINARY_DIR}/lib" )
endif()
if( NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY )
  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${vnl2_BINARY_DIR}/bin" )
endif()
set( CMAKE_MODULE_PATH "${vnl2_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )

include_directories(src)

add_subdirectory(src)

include(CTest)
if( BUILD_TESTING )
  add_subdirectory(tests)
endif( BUILD_TESTING )

