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
    Geometry/Region/Region.h
    Geometry/Region/Rectangle.h
    Geometry/Region/Rectangle.cpp
    Geometry/Region/Circle.h
    Geometry/Region/Circle.cpp    
    Geometry/Region/Parallelogram.h 
    Geometry/Region/Parallelogram.cpp
    Lattice/Lattice.h
    Lattice/Lattice.cpp
    Painter/Asymptote.cpp
    Painter/Painter.h
    Painter/Tikz.cpp
)
target_include_directories(cctb-lib PUBLIC .)

target_link_libraries(
    cctb-lib
    PRIVATE
    Eigen3::Eigen
)

