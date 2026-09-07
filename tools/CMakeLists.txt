cmake_minimum_required(VERSION 3.15)
project(brep2step CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/utf-8 /MD /EHsc)
endif()

set(OCC_ROOT "C:/Zsq/Develop/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0")
include_directories("${OCC_ROOT}/inc")
link_directories("${OCC_ROOT}/win64/vc14/lib")

add_executable(brep2step brep2step.cpp)

target_link_libraries(brep2step
    TKernel
    TKMath
    TKG2d
    TKG3d
    TKGeomBase
    TKBRep
    TKGeomAlgo
    TKTopAlgo
    TKXSBase
    TKSTEP
    TKSTEP209
    TKSTEPAttr
    TKSTEPBase
    TKXDESTEP
    TKIGES
    TKXCAF
    TKLCAF
)
