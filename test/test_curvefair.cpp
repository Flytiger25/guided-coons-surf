// CurveFair 迁移验证：单独编译 + 运行，确认 SGK 版 ComputeEnergyMatrix 正常。
#include "../include/CurveFair.h"
#include <Foundation/init.h>
#include <iostream>
#include <memory>

int main()
{
    sggk::init();

    // 三次 clamped B 样条曲线（4 控制点，节点 [0,0,0,0,1,1,1,1]）
    sggk::Point3DArray ctrl = { sggk::Point3D(0,0,0), sggk::Point3D(0.33,0,0),
                                sggk::Point3D(0.66,0,0), sggk::Point3D(1,0,0) };
    sggk::RealArray knots = { 0.0, 1.0 };
    sggk::UIntArray mults = { 4, 4 };
    auto bs = std::make_shared<sggk::BSplineCurve3D>(3, ctrl, knots, mults);

    CurveFair cf;
    Eigen::MatrixXd M = cf.ComputeEnergyMatrix(bs, 3);
    std::cout << "ComputeEnergyMatrix OK: size=" << M.rows() << "x" << M.cols()
              << " M(0,0)=" << M(0,0) << std::endl;

    // 再验证点投影（CalcNearestPoint 替代 GeomAPI_ProjectPointOnCurve）
    sggk::Point3D q(0.5, 0.3, 0.0);
    double param = -1;
    sggk::Point3D proj = bs->CalcNearestPoint(q, param);
    std::cout << "CalcNearestPoint: proj=(" << proj.X() << "," << proj.Y() << "," << proj.Z()
              << ") param=" << param << std::endl;

    sggk::fini();
    return 0;
}
