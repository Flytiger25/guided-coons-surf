// 验证 SGK 读取 OCC 导出的 .step，并提取 B 样条曲线。
#include <iostream>
#include <Foundation/init.h>
#include <StepExchange/IStepReader.h>
#include <Topology/Brep/Body.h>
#include <Topology/Brep/Edge.h>
#include <Geometry/3D/Curve/Curve3D.h>
#include <Geometry/3D/Curve/BSplineCurve3D.h>

int main(int argc, char** argv)
{
    if (argc < 2) { std::cout << "usage: read_step <file.step>" << std::endl; return 1; }
    sggk::init();
    auto reader = sggk::IStepReader::Create();
    sggk::BodyPtr body = reader->ReadFromFile(argv[1]);
    if (!body) { std::cout << "read failed" << std::endl; sggk::fini(); return 1; }
    auto edges = body->QueryEdges();
    auto faces = body->QueryFaces();
    auto wires = body->QueryWires();
    std::cout << argv[1] << ": edges=" << edges.size()
              << " faces=" << faces.size() << " wires=" << wires.size() << std::endl;
    int idx = 0;
    for (auto& e : edges) {
        auto crv = e->GeomCurve();
        if (!crv) { std::cout << "  edge[" << idx++ << "] <no curve>" << std::endl; continue; }
        auto bs = crv->ToBSpline();
        std::cout << "  edge[" << idx++ << "] type=0x" << std::hex << (int)crv->CurveType() << std::dec;
        if (bs) std::cout << "  ->BSpline degree=" << bs->Degree() << " cp=" << bs->ControlPoints().size();
        std::cout << std::endl;
    }
    sggk::fini();
    return 0;
}
