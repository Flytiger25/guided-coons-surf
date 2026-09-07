cmake_minimum_required(VERSION 3.10)
project(GuidedCoonsSurf)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 显式 Release + 动态运行时（必须与 SGK 的 Release DLL 匹配，否则跨 DLL 传 STL 容器会崩溃）
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# 源码含中文注释（/utf-8）；/MD 动态运行时；/EHsc 异常；/bigobj 大文件
if(MSVC)
    add_compile_options(/utf-8 /MD /EHsc /bigobj)
endif()

# Eigen 为 header-only 库
set(EIGEN_INCLUDE_DIR "C:/Zsq/Develop/eigen-3.4.0")
include_directories("${EIGEN_INCLUDE_DIR}")

# SGK 几何内核（替代 OpenCASCADE）
set(SGK_ROOT "C:/Zsq/Projects/conns-surf/occ_SGK/SGK1.0.24/SGK1.0.24/SGGK/x64-win")
set(SGK_INCLUDE "${SGK_ROOT}/include")
set(SGK_LIB "${SGK_ROOT}/lib")
set(SGK_BIN "${SGK_ROOT}/bin")

include_directories("${SGK_INCLUDE}")
link_directories("${SGK_LIB}")

# 项目头文件目录
include_directories("${CMAKE_SOURCE_DIR}/include")

# 创建可执行文件
add_executable(guided_coons_surf
    src/main.cpp
    src/GuidedCoonsSurfGenerator.cpp
    src/KnotUpdate.cpp
    src/CurveFair.cpp
)

# 编译期注入数据目录
target_compile_definitions(guided_coons_surf PRIVATE
    GUIDED_COONS_DATA_DIR="${CMAKE_SOURCE_DIR}/data")

# 链接 SGK 库（22 个模块，避免遗漏传递依赖）
target_link_libraries(guided_coons_surf
    SGGK_Foundation_ar
    SGGK_GeomBase_ar
    SGGK_Geometry_ar
    SGGK_GeomInt_ar
    SGGK_GeomProject_ar
    SGGK_GeomAlgo_ar
    SGGK_BSplineAlgo_ar
    SGGK_Topology_ar
    SGGK_StepExchange_ar
    SGGK_DataExchangeBase_ar
    SGGK_IgesExchange_ar
    SGGK_Boolean_ar
    SGGK_Offset_ar
    SGGK_ModelingBase_ar
    SGGK_ModelingPrim_ar
    SGGK_Heal_ar
    SGGK_HLR_ar
    SGGK_Lofting_ar
    SGGK_Math_ar
    SGGK_ModelAnalysis_ar
    SGGK_DataBase_ar
    SGGK_DxfReader_ar
)

# 编译后把 SGK 的 DLL 与许可证复制到 exe 目录（保证运行时能找到）
add_custom_command(TARGET guided_coons_surf POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${SGK_BIN}" "$<TARGET_FILE_DIR:guided_coons_surf>"
    COMMAND ${CMAKE_COMMAND} -E copy "${SGK_ROOT}/sggk.lic" "$<TARGET_FILE_DIR:guided_coons_surf>"
)

message(STATUS "SGK include dir: ${SGK_INCLUDE}")
message(STATUS "SGK library dir: ${SGK_LIB}")
