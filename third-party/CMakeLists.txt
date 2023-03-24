cmake_minimum_required(VERSION 3.0)

set(BUILD_GMOCK TRUE CACHE INTERNAL "")
set(INSTALL_GTEST FALSE CACHE INTERNAL "")
set(gmock_build_tests FALSE CACHE INTERNAL "")
set(gtest_build_samples FALSE CACHE INTERNAL "")
set(gtest_build_tests FALSE CACHE INTERNAL "")

add_subdirectory(googletest)