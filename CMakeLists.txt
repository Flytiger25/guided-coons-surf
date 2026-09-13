cmake_minimum_required(VERSION 3.10)
project(GuidedCoonsSurf)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================================
# 依赖查找 OpenCASCADE / Eigen
#   macOS : 默认 Homebrew 安装路径，可用 -DOpenCASCADE_DIR / -DEIGEN_INCLUDE_DIR 覆盖
#   Windows: 通过 -DOpenCASCADE_DIR=<OCC安装路径>/lib/cmake/opencascade、
#            -DEIGEN_INCLUDE_DIR=<eigen解压路径> 指定
# ============================================================================
if(APPLE)
    # Homebrew 默认路径
    set(OpenCASCADE_DIR "/opt/homebrew/Cellar/opencascade/7.9.1/lib/cmake/opencascade" CACHE PATH "OpenCASCADE CMake 配置目录")
    set(EIGEN_INCLUDE_DIR "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3" CACHE PATH "Eigen 头文件目录")
endif()

# 查找 OpenCASCADE
find_package(OpenCASCADE REQUIRED)
if(OpenCASCADE_INCLUDE_DIR)
    include_directories(${OpenCASCADE_INCLUDE_DIR})
endif()

# Eigen 头文件目录（用户未指定且未找到时给出提示）
if(NOT EIGEN_INCLUDE_DIR)
    message(FATAL_ERROR "未找到 Eigen 头文件目录，请通过 -DEIGEN_INCLUDE_DIR=<path> 指定")
endif()
include_directories(${EIGEN_INCLUDE_DIR})

# 项目头文件目录
include_directories("${CMAKE_SOURCE_DIR}/include")

# ============================================================================
# 源文件
#   Windows 暂不提供可视化，不编译 Viewer.mm（依赖 Cocoa）
# ============================================================================
set(SOURCES
    src/main.cpp
    src/GuidedCoonsSurfGenerator.cpp
    src/KnotUpdate.cpp
    src/CurveFair.cpp
)

if(APPLE)
    list(APPEND SOURCES src/Viewer.mm)
endif()

add_executable(guided_coons_surf ${SOURCES})

# OpenCASCADE 核心库（两平台通用）
set(OCC_CORE_LIBS
    TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo
    TKPrim TKBO TKBool TKHLR TKFillet TKOffset TKFeat TKMesh TKXSBase
    TKShHealing TKXCAF TKLCAF TKDESTEP
)

# 可视化库（仅 macOS，Windows 暂不需要）
if(APPLE)
    set(OCC_GUI_LIBS TKV3d TKService TKOpenGl "-framework Cocoa")
else()
    set(OCC_GUI_LIBS)
endif()

target_link_libraries(guided_coons_surf ${OCC_CORE_LIBS} ${OCC_GUI_LIBS})

# 输出信息
message(STATUS "OpenCASCADE found: ${OpenCASCADE_FOUND}")
message(STATUS "OpenCASCADE include dir: ${OpenCASCADE_INCLUDE_DIR}")
message(STATUS "OpenCASCADE library dir: ${OpenCASCADE_LIBRARY_DIR}")
message(STATUS "Platform: ${CMAKE_SYSTEM_NAME}")
