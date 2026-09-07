cmake_minimum_required(VERSION 3.15)
project(sgk_smoke CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 源码含中文注释，MSVC 需按 UTF-8 解析（否则 GBK 乱码导致语法错误）
# /MD：动态运行时，必须与 SGK DLL（依赖 msvcp140.dll）保持一致
# /EHsc：启用 C++ 异常
if(MSVC)
    add_compile_options(/utf-8 /MD /EHsc /bigobj)
endif()

# 显式 Release，确保 _ITERATOR_DEBUG_LEVEL=0 与 SGK 的 Release DLL 一致
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# SGK 根目录（x64-win 下含 include/lib/bin/sggk.lic）
set(SGK_ROOT "C:/Zsq/Projects/conns-surf/occ_SGK/SGK1.0.24/SGK1.0.24/SGGK/x64-win")
set(SGK_INCLUDE "${SGK_ROOT}/include")
set(SGK_LIB "${SGK_ROOT}/lib")
set(SGK_BIN "${SGK_ROOT}/bin")

include_directories("${SGK_INCLUDE}")
link_directories("${SGK_LIB}")

add_executable(sgk_smoke sgk_smoke.cpp)

add_executable(read_step read_step.cpp)
target_link_libraries(read_step
    SGGK_Foundation_ar
    SGGK_GeomBase_ar
    SGGK_Geometry_ar
    SGGK_Topology_ar
    SGGK_StepExchange_ar
    SGGK_DataExchangeBase_ar
)

# KnotUpdate 迁移验证：编译 src/KnotUpdate.cpp + 测试入口
add_executable(test_knotupdate test_knotupdate.cpp ../src/KnotUpdate.cpp)
target_include_directories(test_knotupdate PRIVATE ../include)
target_link_libraries(test_knotupdate
    SGGK_Foundation_ar
    SGGK_GeomBase_ar
    SGGK_Geometry_ar
)
add_custom_command(TARGET test_knotupdate POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${SGK_BIN}" "$<TARGET_FILE_DIR:test_knotupdate>"
    COMMAND ${CMAKE_COMMAND} -E copy "${SGK_ROOT}/sggk.lic" "$<TARGET_FILE_DIR:test_knotupdate>"
)

# CurveFair 迁移验证：编译 src/CurveFair.cpp + 测试入口（含 Eigen）
add_executable(test_curvefair test_curvefair.cpp ../src/CurveFair.cpp)
target_include_directories(test_curvefair PRIVATE ../include "C:/Zsq/Develop/eigen-3.4.0")
target_link_libraries(test_curvefair
    SGGK_Foundation_ar
    SGGK_GeomBase_ar
    SGGK_Geometry_ar
    SGGK_BSplineAlgo_ar
)
add_custom_command(TARGET test_curvefair POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${SGK_BIN}" "$<TARGET_FILE_DIR:test_curvefair>"
    COMMAND ${CMAKE_COMMAND} -E copy "${SGK_ROOT}/sggk.lic" "$<TARGET_FILE_DIR:test_curvefair>"
)

# GuidedCoonsSurfGenerator 迁移编译检查（仅编译，看错误清单）
add_library(gen_check OBJECT
    ../src/GuidedCoonsSurfGenerator.cpp
    ../src/KnotUpdate.cpp
    ../src/CurveFair.cpp
)
target_include_directories(gen_check PRIVATE ../include "C:/Zsq/Develop/eigen-3.4.0")
target_link_libraries(gen_check
    SGGK_Foundation_ar
    SGGK_GeomBase_ar
    SGGK_Geometry_ar
    SGGK_GeomProject_ar
    SGGK_BSplineAlgo_ar
    SGGK_Topology_ar
    SGGK_StepExchange_ar
    SGGK_DataExchangeBase_ar
)

target_link_libraries(sgk_smoke
    SGGK_Foundation_ar
    SGGK_GeomBase_ar
    SGGK_Geometry_ar
    SGGK_GeomProject_ar
    SGGK_BSplineAlgo_ar
    SGGK_Topology_ar
    SGGK_StepExchange_ar
    SGGK_DataExchangeBase_ar
)

# 编译后把 SGK 的 DLL 与许可证复制到 exe 目录
add_custom_command(TARGET sgk_smoke POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${SGK_BIN}" "$<TARGET_FILE_DIR:sgk_smoke>"
    COMMAND ${CMAKE_COMMAND} -E copy "${SGK_ROOT}/sggk.lic" "$<TARGET_FILE_DIR:sgk_smoke>"
)
