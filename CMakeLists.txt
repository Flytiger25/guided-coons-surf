cmake_minimum_required(VERSION 3.10)
project(GuidedCoonsSurf)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置OpenCASCADE路径 (Homebrew安装路径)
set(OpenCASCADE_DIR "C:/Zsq/Develop/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/cmake")

# 设置Eigen路径 (Homebrew安装路径)
set(Eigen_DIR "/opt/homebrew/Cellar/eigen/3.4.0_1/share/eigen3/cmake")

# 查找OpenCASCADE
find_package(OpenCASCADE REQUIRED)

# 包含OpenCASCADE头文件目录
include_directories("C:/Zsq/Develop/OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/inc")

include_directories("C:/Zsq/Develop/eigen-3.4.0")

# 项目头文件目录
include_directories("${CMAKE_SOURCE_DIR}/include")

# 创建可执行文件
add_executable(guided_coons_surf
    src/main.cpp
    src/GuidedCoonsSurfGenerator.cpp
    src/KnotUpdate.cpp
    src/CurveFair.cpp
)

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
