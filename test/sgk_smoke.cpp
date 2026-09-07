// SGK 冒烟测试：验证 SGK 内核在本机可编译、可运行，关键 API 可用。
#include <iostream>
#include <memory>
#include <cmath>
#include <cstdio>

#include <Foundation/init.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Point2D.h>
#include <Geometry/3D/Curve/BSplineCurve3D.h>
#include <Geometry/3D/Surface/BSplineSurface.h>
#include <GeomProject/GeomProject.h>
#include <StepExchange/IStepReader.h>
#include <StepExchange/IStepWriter.h>
#include <Topology/Tools/TopoBuilder.h>
#include <Topology/Brep/Body.h>

#define STEP(msg) std::cout << "[STEP] " << msg << std::endl

int main()
{
    std::cout << "== SGK smoke test ==" << std::endl;

    try {
        STEP("init");
        sggk::init();
        STEP("init done");

        STEP("BSplineCurve3D construct + CalcPoint");
        {
            sggk::Point3DArray ctrl = { sggk::Point3D(0, 0, 0), sggk::Point3D(1, 0, 0) };
            sggk::RealArray knots = { 0.0, 1.0 };
            sggk::UIntArray mults = { 2, 2 };
            sggk::BSplineCurve3D crv(1, ctrl, knots, mults);
            sggk::Point3D p = crv.CalcPoint(0.5);
            std::cout << "  CalcPoint(0.5) = (" << p.X() << "," << p.Y() << "," << p.Z() << ")" << std::endl;
        }

        STEP("BSplineSurface construct + CalcPoint");
        sggk::Point3DMatrix ctrl(2, sggk::Point3DArray(2));
        ctrl[0][0] = sggk::Point3D(0, 0, 0);
        ctrl[0][1] = sggk::Point3D(0, 1, 0);
        ctrl[1][0] = sggk::Point3D(1, 0, 0);
        ctrl[1][1] = sggk::Point3D(1, 1, 0);
        sggk::RealArray ku = { 0.0, 1.0 };
        sggk::RealArray kv = { 0.0, 1.0 };
        sggk::UIntArray mu = { 2, 2 };
        sggk::UIntArray mv = { 2, 2 };
        auto surf = std::make_shared<sggk::BSplineSurface>(1, 1, ctrl, ku, kv, mu, mv);
        {
            sggk::Point3D p = surf->CalcPoint(0.5, 0.5);
            std::cout << "  CalcPoint(0.5,0.5) = (" << p.X() << "," << p.Y() << "," << p.Z() << ")" << std::endl;
        }

        STEP("CalcNearestPnt 点投影（替代 PntSrfProject）");
        {
            sggk::Point3D q(0.4, 0.6, 1.0);
            sggk::Point2D uv;
            sggk::Point3D proj = surf->CalcNearestPnt(q, uv);
            std::cout << "  proj=(" << proj.X() << "," << proj.Y() << "," << proj.Z()
                      << ") uv=(" << uv.X() << "," << uv.Y() << ")" << std::endl;
        }

        STEP("构建完整 Face（4 边 Loop）+ STEP 写/读");
        {
            // 4 条线性 B 样条边界边
            auto makeLineEdge = [](const sggk::Point3D& a, const sggk::Point3D& b) {
                sggk::Point3DArray cp = { a, b };
                sggk::RealArray k = { 0.0, 1.0 };
                sggk::UIntArray m = { 2, 2 };
                sggk::BSplineCurve3D crv(1, cp, k, m);
                return sggk::TopoBuilder::MakeEdge(crv, false, true);
            };
            sggk::Point3D v0(0, 0, 0), v1(1, 0, 0), v2(1, 1, 0), v3(0, 1, 0);
            auto e0 = makeLineEdge(v0, v1);
            auto e1 = makeLineEdge(v1, v2);
            auto e2 = makeLineEdge(v2, v3);
            auto e3 = makeLineEdge(v3, v0);

            sggk::CoedgeList coedges;
            coedges.push_back(sggk::TopoBuilder::MakeCoedge(e0, true));
            coedges.push_back(sggk::TopoBuilder::MakeCoedge(e1, true));
            coedges.push_back(sggk::TopoBuilder::MakeCoedge(e2, true));
            coedges.push_back(sggk::TopoBuilder::MakeCoedge(e3, true));
            auto loop = sggk::TopoBuilder::MakeLoop(coedges);

            auto msrf = sggk::TopoBuilder::MakeModelSurface(surf);
            auto face = sggk::TopoBuilder::MakeFace(msrf, true);
            sggk::TopoBuilder::FaceAddLoop(face, loop);
            auto body = sggk::TopoBuilder::MakeBody(face);

            const char* path = "smoke_roundtrip.step";
            auto writer = sggk::IStepWriter::Create();
            writer->WriteToFile(body, path);
            STEP("STEP written: " + std::string(path));

            auto reader = sggk::IStepReader::Create();
            sggk::BodyPtr rbody = reader->ReadFromFile(path);
            if (rbody) {
                auto faces = rbody->QueryFaces();
                auto edges = rbody->QueryEdges();
                std::cout << "  read back faces = " << faces.size()
                          << ", edges = " << edges.size() << std::endl;
            } else {
                std::cout << "  read back body is null" << std::endl;
            }
            std::remove(path);
        }

        STEP("fini");
        sggk::fini();
        STEP("ALL DONE");
    } catch (const std::exception& e) {
        std::cout << "[EXCEPTION] " << e.what() << std::endl;
        return 2;
    } catch (...) {
        std::cout << "[EXCEPTION] unknown" << std::endl;
        return 3;
    }

    return 0;
}
