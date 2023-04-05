cmake_minimum_required(VERSION 3.23)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(
    cctb
    main.cpp
)

target_compile_options(
    cctb
    PRIVATE
    "${CCTB_CXX_COMPILER_OPTIONS}"    
)

set_target_properties(cctb PROPERTIES RUNTIME_OUTPUT_DIRECTORY ..)
target_link_libraries(cctb PRIVATE cctb-lib)

add_library(
    cctb-lib
    Geometry/Line.cpp 
    Geometry/Line.h
    Geometry/Point.h
    Geometry/Vector.h
    Lattice/lattice.cpp
    Lattice/lattice.h
    LinearAlgebra/BlasImpl.cpp
    LinearAlgebra/BlasImpl.h   
    LinearAlgebra/LapackImpl.cpp
    LinearAlgebra/LapackImpl.h
    LinearAlgebra/Matrix.cpp    
    LinearAlgebra/Matrix.h
    LinearAlgebra/NumericArray.cpp
    LinearAlgebra/NumericArray.h
    Painter/Asymptote.cpp
    Painter/Painter.h
    Painter/Tikz.cpp
)
target_include_directories(cctb-lib PUBLIC .)

target_link_libraries(
    cctb-lib
    PRIVATE
    BLAS::BLAS
)

