cmake_minimum_required(VERSION 3.0)

add_executable(cctb-test
    test-lattice.cpp
    test-geometry.cpp
    test-vector.cpp
    test-matrix.cpp
)
target_include_directories(cctb-test PRIVATE .)
target_link_libraries(cctb-test
    PRIVATE
    Catch2::Catch2WithMain
    cctb-lib
)

add_test(
    NAME cctb-test
    COMMAND cctb-test
)

include(extra_compile_options)