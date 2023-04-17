cmake_minimum_required(VERSION 3.23)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project(cctb CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CCTB_CXX_COMPILER_OPTIONS "" CACHE STRING "")
mark_as_advanced(CCTB_CXX_COMPILER_OPTIONS)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(third-party)
add_subdirectory(src)

if (BUILD_TESTING)
    include(CTest)
    add_subdirectory(test)

    if (CCTB_ENABLE_COVERAGE)
       include(CodeCoverage)
   
       target_compile_options(cctb-test PRIVATE --coverage)
       target_link_libraries(cctb-test PRIVATE gcov)
    
       target_compile_options(cctb-lib PRIVATE --coverage)
       target_link_libraries(cctb-lib PRIVATE gcov)
    
       setup_target_for_coverage_lcov(
          NAME coverage
          EXECUTABLE cctb-test
          DEPENDENCIES cctb-lib
          EXCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/third-party/*"
       )
    endif()
endif()
