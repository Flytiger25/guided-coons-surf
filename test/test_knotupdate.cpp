// KnotUpdate 迁移验证：单独编译 + 运行，确认 SGK 版能正常工作。
#include "../include/KnotUpdate.h"
#include <Foundation/init.h>
#include <iostream>
#include <memory>

int main()
{
    sggk::init();

    // 构造一条三次 B 样条曲线（clamped，4 控制点）
    sggk::Point3DArray ctrl = { sggk::Point3D(0,0,0), sggk::Point3D(0.33,0,0),
                                sggk::Point3D(0.66,0,0), sggk::Point3D(1,0,0) };
    sggk::RealArray knots = { 0.0, 1.0 };
    sggk::UIntArray mults = { 4, 4 };
    auto bs = std::make_shared<sggk::BSplineCurve3D>(3, ctrl, knots, mults);

    std::vector<double> sequences = { 0,0,0,0,1,1,1,1 };
    std::vector<sggk::Point3D> pnts = { sggk::Point3D(0.1,0,0), sggk::Point3D(0.5,0,0), sggk::Point3D(0.9,0,0) };
    std::vector<double> params = { 0.1, 0.5, 0.9 };

    KnotUpdate ku(bs, sequences, pnts, params);

    double newKnot = ku.SelfSingleUpdate(MID_KNOT_BY_SINGLE_ERROR);
    std::cout << "KnotUpdate OK: newKnot=" << newKnot
              << " maxError=" << ku.getMaxError()
              << " sequences=" << ku.getSequences().size() << std::endl;

    sggk::fini();
    return 0;
}
