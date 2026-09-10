// 一次性工具：把 OCC 的 .brep 转为 .step（供 SGK 读取）。
// 仅用于数据准备，迁移完成后项目代码不再依赖 OCC。
#include <iostream>
#include <TopoDS_Shape.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <STEPControl_Writer.hxx>
#include <IFSelect_ReturnStatus.hxx>

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: brep2step <in.brep> <out.step>" << std::endl;
        return 1;
    }
    TopoDS_Shape shape;
    BRep_Builder builder;
    if (!BRepTools::Read(shape, argv[1], builder)) {
        std::cerr << "[ERROR] read failed: " << argv[1] << std::endl;
        return 1;
    }
    STEPControl_Writer writer;
    if (writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone) {
        std::cerr << "[ERROR] transfer failed" << std::endl;
        return 1;
    }
    if (writer.Write(argv[2]) != IFSelect_RetDone) {
        std::cerr << "[ERROR] write failed" << std::endl;
        return 1;
    }
    std::cout << "OK: " << argv[1] << " -> " << argv[2] << std::endl;
    return 0;
}
