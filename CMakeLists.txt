cmake_minimum_required(VERSION 3.0)
project(cctb CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -pedantic -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

add_subdirectory(src)
add_subdirectory(third-party)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(test)
endif()