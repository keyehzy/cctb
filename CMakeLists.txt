cmake_minimum_required(VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project(cctb CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(compile_options)
include(extra_compile_options)

add_subdirectory(src)
add_subdirectory(third-party)

enable_testing()

include(CTest)
add_subdirectory(test)
