cmake_minimum_required(VERSION 3.0)

add_executable(cctb main.cpp)
target_link_libraries(cctb PRIVATE cctb-lib)
add_library(
    cctb-lib
    cctb.cpp
    lattice.cpp
    Painter/Tikz.cpp
    Painter/Asymptote.cpp
    Geometry/Line.cpp
    LinearAlgebra/NumericArray.cpp
    LinearAlgebra/Vector.cpp
    LinearAlgebra/Matrix.cpp
)
target_include_directories(cctb-lib PUBLIC .)
target_link_libraries(cctb-lib PRIVATE -lm)