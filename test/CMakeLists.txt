cmake_minimum_required(VERSION 3.23)

add_executable(cctb-test
    test-lattice.cpp
    test-geometry.cpp
)

target_compile_options(
    cctb-test
    PRIVATE
    "${CCTB_CXX_COMPILER_OPTIONS}"    
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
