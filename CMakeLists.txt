cmake_minimum_required(VERSION 3.23)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project(cctb CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(compile_options)
include(extra_compile_options)

find_package(BLAS)
find_package(LAPACK)

if (BLAS_FOUND)
  add_definitions(-DHAVE_BLAS)
  find_path(BLAS_INCLUDE_DIRS cblas.h $ENV{CMAKE_PREFIX_PATH}/include)
  include_directories(${BLAS_INCLUDE_DIRS})
else()
  message(WARNING "BLAS not found")
endif()

if (LAPACK_FOUND)
  add_definitions(-DHAVE_LAPACK)
  find_path(LAPACK_INCLUDE_DIRS lapacke.h $ENV{CMAKE_PREFIX_PATH}/include)
  include_directories(${LAPACK_INCLUDE_DIRS})
else()
  message(WARNING "LAPACK not found")
endif()

add_subdirectory(src)
add_subdirectory(third-party)

enable_testing()

include(CTest)
add_subdirectory(test)
