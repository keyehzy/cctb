cmake_minimum_required(VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project(cctb CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(compile_options)
include(extra_compile_options)

find_package(BLAS)

if (BLAS_FOUND)
  add_definitions(-DHAVE_BLAS)
  find_path(BLAS_INCLUDE_DIRS cblas.h $ENV{CMAKE_PREFIX_PATH}/include)
  include_directories(${BLAS_INCLUDE_DIRS})
else()
  message(WARNING "BLAS not found")
endif()

add_subdirectory(src)
add_subdirectory(third-party)

enable_testing()

include(CTest)
add_subdirectory(test)
