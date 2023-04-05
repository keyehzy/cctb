cmake_minimum_required(VERSION 3.23)

add_executable(cctb main.cpp)
target_link_libraries(cctb PRIVATE cctb-lib)
add_library(
    cctb-lib
    cctb.cpp
    lattice.cpp
    Painter/Tikz.cpp
    Painter/Asymptote.cpp
    Geometry/Line.cpp
    LinearAlgebra/BlasImpl.cpp
    LinearAlgebra/LapackImpl.cpp
    LinearAlgebra/NumericArray.cpp
    LinearAlgebra/Matrix.cpp
)
target_include_directories(cctb-lib PUBLIC .)

if (BLAS_FOUND AND LAPACK_FOUND)
  target_link_libraries(cctb-lib PRIVATE -lm ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} -llapacke)
else()
  target_link_libraries(cctb-lib PRIVATE -lm)
endif()
