#pragma once
#include <TopoDS_Shape.hxx>

//! @brief 在 Cocoa 窗口中交互式显示 TopoDS_Shape
//! @param theShape  要显示的几何体 (Face, Shell, Solid 等)
//! @param theTitle  窗口标题
void DisplayShape(const TopoDS_Shape& theShape, const char* theTitle = "Guided Coons Surface");
