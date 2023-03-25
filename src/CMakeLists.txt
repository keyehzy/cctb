cmake_minimum_required(VERSION 3.0)

add_executable(cctb main.cpp)
target_link_libraries(cctb PRIVATE cctb-lib)
add_library(
    cctb-lib
    cctb.cpp
)
target_include_directories(cctb-lib PUBLIC .)
target_link_libraries(cctb-lib PRIVATE -lm)