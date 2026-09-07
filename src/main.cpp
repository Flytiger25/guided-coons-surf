#include <iostream>
#include <Foundation/init.h>
#include <Geometry/3D/Curve/BSplineCurve3D.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <StepExchange/IStepReader.h>
#include <StepExchange/IStepWriter.h>
#include <Topology/Tools/TopoBuilder.h>
#include <Topology/Brep/Body.h>
#include <Topology/Brep/Edge.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "GuidedCoonsSurfGenerator.h"

// 把曲面（带 4 条等参边界线 Loop）导出为 .step
void ExportSurfaceToStep(sggk::BSplineSurfacePtr surf, const std::string& filePath)
{
    if (!surf) return;

    auto msrf = sggk::TopoBuilder::MakeModelSurface(surf);
    auto face = sggk::TopoBuilder::MakeFace(msrf, true);

    // 4 条等参边界线
    sggk::Curve3DPtr c0 = surf->CalcUCurve(surf->MinParamV());  // v=min, u: min->max
    sggk::Curve3DPtr c1 = surf->CalcVCurve(surf->MaxParamU());  // u=max, v: min->max
    sggk::Curve3DPtr c2 = surf->CalcUCurve(surf->MaxParamV());  // v=max, u: min->max（需反向）
    sggk::Curve3DPtr c3 = surf->CalcVCurve(surf->MinParamU());  // u=min, v: min->max（需反向）
    c2->Reverse();
    c3->Reverse();

    sggk::EdgePtr e0 = sggk::TopoBuilder::MakeEdge(sggk::TopoBuilder::MakeModelCurve(c0), true);
    sggk::EdgePtr e1 = sggk::TopoBuilder::MakeEdge(sggk::TopoBuilder::MakeModelCurve(c1), true);
    sggk::EdgePtr e2 = sggk::TopoBuilder::MakeEdge(sggk::TopoBuilder::MakeModelCurve(c2), true);
    sggk::EdgePtr e3 = sggk::TopoBuilder::MakeEdge(sggk::TopoBuilder::MakeModelCurve(c3), true);

    sggk::CoedgeList coedges;
    coedges.push_back(sggk::TopoBuilder::MakeCoedge(e0, true));
    coedges.push_back(sggk::TopoBuilder::MakeCoedge(e1, true));
    coedges.push_back(sggk::TopoBuilder::MakeCoedge(e2, true));
    coedges.push_back(sggk::TopoBuilder::MakeCoedge(e3, true));
    auto loop = sggk::TopoBuilder::MakeLoop(coedges);

    sggk::TopoBuilder::FaceAddLoop(face, loop);
    auto body = sggk::TopoBuilder::MakeBody(face);

    auto writer = sggk::IStepWriter::Create();
    writer->WriteToFile(body, filePath.c_str());
    std::cout << "Exported: " << filePath << std::endl;
}

// 从 .step 读取所有边并提取为 B 样条曲线
void LoadBSplineCurves(const std::string& filePath, std::vector<sggk::BSplineCurve3DPtr>& curveArray)
{
    auto reader = sggk::IStepReader::Create();
    sggk::BodyPtr body = reader->ReadFromFile(filePath.c_str());
    if (!body)
    {
        std::cerr << "[ERROR] read step failed: " << filePath << std::endl;
        return;
    }

    auto edges = body->QueryEdges();
    for (auto& edge : edges)
    {
        auto crv = edge->GeomCurve();
        if (!crv) continue;
        sggk::BSplineCurve3DPtr bs = crv->ToBSpline();
        if (bs) curveArray.push_back(bs);
    }
}

int main()
{
try {
    sggk::init();

    // 获取 guideCurves（内部引导线）
    std::vector<sggk::BSplineCurve3DPtr> guideCurves;
    std::string fileName = std::string(GUIDED_COONS_DATA_DIR) + "/input/1_internal.step";
    LoadBSplineCurves(fileName, guideCurves);

    // 获取 boundary（边界线）
    std::vector<sggk::BSplineCurve3DPtr> boundary;
    fileName = std::string(GUIDED_COONS_DATA_DIR) + "/input/1_boundary.step";
    LoadBSplineCurves(fileName, boundary);

    // 生成带引导线的 Coons 曲面
    GuidedCoonsSurfGenerator msg = GuidedCoonsSurfGenerator(boundary, guideCurves);
    msg.Perform();
    sggk::BSplineSurfacePtr guidedSurf = msg.GuidedSurf();
    if (!msg.IsDone())
    {
        std::cout << "Case 1 Guided Failing!" << std::endl;
    }

    // 导出结果
    fileName = std::string(GUIDED_COONS_DATA_DIR) + "/output/1_guidedCoonsSurf.step";
    ExportSurfaceToStep(guidedSurf, fileName);

    sggk::fini();
    return 0;
} catch (const std::exception& e) {
    std::cerr << "[EXCEPTION] " << e.what() << std::endl;
    return 1;
} catch (...) {
    std::cerr << "[UNKNOWN EXCEPTION]" << std::endl;
    return 2;
}
}
