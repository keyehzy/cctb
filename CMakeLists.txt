cmake_minimum_required(VERSION 3.23)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project(cctb CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CCTB_CXX_COMPILER_OPTIONS "" CACHE STRING "")
mark_as_advanced(CCTB_CXX_COMPILER_OPTIONS)

add_subdirectory(src)
add_subdirectory(third-party)

include(CTest)
if (BUILD_TESTING)
   add_subdirectory(test)
endif()
