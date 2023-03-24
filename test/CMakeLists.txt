cmake_minimum_required(VERSION 3.0)

add_executable(cctb-test basic-test.cpp)
target_include_directories(cctb-test PRIVATE .)
target_link_libraries(cctb-test
    PRIVATE
    gmock
    gmock_main
    gtest
    cctb-lib
)
add_test(
    NAME cctb-test
    COMMAND cctb-test
)