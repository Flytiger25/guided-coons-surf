cmake_minimum_required(VERSION 3.10)
project(GuidedCoonsSurf)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenCASCADE Windows 安装路径 (7.7.0, MSVC 2015+ x64)
set(OpenCASCADE_DIR "C:/Zsq/Develop/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/cmake")

# Eigen 为 header-only 库，只需包含头文件目录，无需 find_package
set(EIGEN_INCLUDE_DIR "C:/Zsq/Develop/eigen-3.4.0")

# 查找OpenCASCADE
find_package(OpenCASCADE REQUIRED)

# 包含OpenCASCADE头文件目录
include_directories("C:/Zsq/Develop/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/inc")

include_directories("${EIGEN_INCLUDE_DIR}")

# 项目头文件目录
include_directories("${CMAKE_SOURCE_DIR}/include")

# OpenCASCADE 7.7 自带 3rdparty 库的目录（链接 TKV3d 等模块时会用到）
# openvr_api.lib 属于 TKV3d，freeimage/tbb/freetype 等用于其他模块
set(OCC_3RDPARTY_ROOT "C:/Zsq/Develop/OpenCASCADE-7.7.0-vc14-64")
link_directories(
    "${OCC_3RDPARTY_ROOT}/openvr-1.14.15-64/lib/win64"
    "${OCC_3RDPARTY_ROOT}/freetype-2.5.5-vc14-64/lib/win64"
    "${OCC_3RDPARTY_ROOT}/tbb_2021.5-vc14-64/lib/win64"
    "${OCC_3RDPARTY_ROOT}/freeimage-3.17.0-vc14-64/lib/win64"
)

# 创建可执行文件
add_executable(guided_coons_surf
    src/main.cpp
    src/GuidedCoonsSurfGenerator.cpp
    src/KnotUpdate.cpp
    src/CurveFair.cpp
)

# 编译期注入数据目录，源码通过 GUIDED_COONS_DATA_DIR 定位输入/输出文件，
# 与运行时工作目录无关（Windows / macOS 通用）
target_compile_definitions(guided_coons_surf PRIVATE
    GUIDED_COONS_DATA_DIR="${CMAKE_SOURCE_DIR}/data")

# 链接OpenCASCADE库
target_link_libraries(guided_coons_surf 
    TKernel 
    TKMath 
    TKG2d 
    TKG3d 
    TKGeomBase 
    TKBRep 
    TKGeomAlgo 
    TKTopAlgo 
    TKPrim 
    TKBO 
    TKBool 
    TKHLR 
    TKFillet 
    TKOffset 
    TKFeat 
    TKMesh 
    TKXSBase 
    TKShHealing 
    TKXCAF 
    TKLCAF 
    TKV3d 
    TKService
    TKXDESTEP
)

# 输出信息
message(STATUS "OpenCASCADE found: ${OpenCASCADE_FOUND}")
message(STATUS "OpenCASCADE include dir: ${OpenCASCADE_INCLUDE_DIR}")
message(STATUS "OpenCASCADE library dir: ${OpenCASCADE_LIBRARY_DIR}")
