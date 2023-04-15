cmake_minimum_required(VERSION 3.23)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project(cctb CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CCTB_CXX_COMPILER_OPTIONS "" CACHE STRING "")
mark_as_advanced(CCTB_CXX_COMPILER_OPTIONS)

include(FetchContent)
FetchContent_Declare(
    Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
)
FetchContent_MakeAvailable(Eigen3)
include_directories(${Eigen3_SOURCE_DIR})

add_subdirectory(third-party)
add_subdirectory(src)

if (BUILD_TESTING)
   include(CTest)
   add_subdirectory(test)
endif()

