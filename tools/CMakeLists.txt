cmake_minimum_required(VERSION 3.15)
project(brep2step CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/utf-8 /MD /EHsc)
endif()

# OpenCASCADE 根目录（仅 brep2step 转换工具需要，SGK 迁移后为遗留）
# 用法：-DOCC_ROOT=<path> 或设置环境变量 OCC_ROOT
set(OCC_ROOT "" CACHE PATH "Path to OpenCASCADE root")
if(NOT OCC_ROOT AND DEFINED ENV{OCC_ROOT})
    set(OCC_ROOT $ENV{OCC_ROOT})
endif()
if(NOT OCC_ROOT)
    message(FATAL_ERROR
        "OCC_ROOT not set. Pass -DOCC_ROOT=<path> or set env var OCC_ROOT.")
endif()
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
