#include "GuidedCoonsSurfGenerator.h"

#include <GeomInt/GeomInt.h>
#include <limits>
#include <cmath>

namespace {

//! 查询 B 样条曲线在指定节点处的重数（节点不存在返回 0）
unsigned int CurveKnotMultiplicity(const sggk::BSplineCurve3DPtr& crv, double knot, double tol = 1e-7)
{
    const sggk::RealArray& K = crv->Knots();
    const sggk::UIntArray& M = crv->Mults();
    for (size_t i = 0; i < K.size(); ++i)
        if (std::fabs(K[i] - knot) <= tol)
            return M[i];
    return 0;
}

//! 查询 B 样条曲面 U 方向在指定节点处的重数（节点不存在返回 0）
unsigned int SurfaceUKnotMultiplicity(const sggk::BSplineSurfacePtr& surf, double knot, double tol = 1e-7)
{
    const sggk::RealArray& K = surf->KnotsU();
    const sggk::UIntArray& M = surf->MultsU();
    for (size_t i = 0; i < K.size(); ++i)
        if (std::fabs(K[i] - knot) <= tol)
            return M[i];
    return 0;
}

//! 查询 B 样条曲面 V 方向在指定节点处的重数（节点不存在返回 0）
unsigned int SurfaceVKnotMultiplicity(const sggk::BSplineSurfacePtr& surf, double knot, double tol = 1e-7)
{
    const sggk::RealArray& K = surf->KnotsV();
    const sggk::UIntArray& M = surf->MultsV();
    for (size_t i = 0; i < K.size(); ++i)
        if (std::fabs(K[i] - knot) <= tol)
            return M[i];
    return 0;
}

} // namespace

GuidedCoonsSurfGenerator::GuidedCoonsSurfGenerator(const std::vector<sggk::BSplineCurve3DPtr>& boundaryCurves, const std::vector<sggk::BSplineCurve3DPtr>& guideCurves, double theTol)
	: m_boundaryCurves(boundaryCurves), m_guideCurves(guideCurves), m_tol(theTol), m_isDone(false), m_iterateCount(0)
{
}

void GuidedCoonsSurfGenerator::Perform()
{
    // 三边情况，构造退化边
    if (m_boundaryCurves.size() == 3)
    {
        AddDegenerateCurve(m_boundaryCurves);
    }

    // 构造Coons
    ConstructCoonsSurf();

    // 无引导线，返回Coons曲面
    if (m_guideCurves.size() == 0)
    {
        m_guidedSurf = m_originalSurf;
        m_isDone = true;
        return;
    }

    // 调试用中间 STEP 导出（OCC 的 TopoDS/STEPControl 已移除，SGK 下暂不导出）
    //TopoDS_Face coonsFace = BRepBuilderAPI_MakeFace(m_originalSurf, 1e-7);
    //std::string filePath = std::string(GUIDED_COONS_DATA_DIR) + "/coons/";
    //filePath += "coons.step";
    //STEPControl_Writer stepWriter;
    //stepWriter.Transfer(coonsFace, STEPControl_AsIs);
    //stepWriter.Write(filePath.c_str());

    // 未达到容差要求，继续迭代
    while (!m_isDone && m_iterateCount < MAX_ITERATIONS)
    {
        std::cout << "m_iterateCount: " << m_iterateCount << std::endl;
        ConstructSurfWithGuideCrvs();
        m_originalSurf = m_guidedSurf;
        m_iterateCount++;

        //TopoDS_Face coonsFace = BRepBuilderAPI_MakeFace(m_originalSurf, 1e-7);
        //std::string filePath = std::string(GUIDED_COONS_DATA_DIR) + "/coons/";
        //filePath += "GuidedSurf_" + std::to_string(m_iterateCount) + ".step";
        //STEPControl_Writer stepWriter;
        //stepWriter.Transfer(coonsFace, STEPControl_AsIs);
        //stepWriter.Write(filePath.c_str());
    }
}

void GuidedCoonsSurfGenerator::ConstructCoonsSurf()
{
    // 重新拟合边界线
    ApproximateBoundaryCurves(m_boundaryCurves);
    // 构造Coons
    sggk::BSplineCurve3DPtr bsplineCurve1, bsplineCurve2, bsplineCurve3, bsplineCurve4;
    Arrange_Coons_G0(m_boundaryCurves, bsplineCurve1, bsplineCurve2, bsplineCurve3, bsplineCurve4);
    Coons_G0(bsplineCurve1, bsplineCurve2, bsplineCurve3, bsplineCurve4, m_originalSurf);
    m_coonsSurf = m_originalSurf;

    // 裁剪引导线
    TrimInternalCurves(m_guideCurves, m_boundaryCurves);
}

int GuidedCoonsSurfGenerator::Arrange_Coons_G0(std::vector<sggk::BSplineCurve3DPtr>& curveArray, sggk::BSplineCurve3DPtr& bslpineCurve1,
    sggk::BSplineCurve3DPtr& bslpineCurve2, sggk::BSplineCurve3DPtr& bslpineCurve3, sggk::BSplineCurve3DPtr& bslpineCurve4, double tol, int isModify)
{
    std::vector<sggk::BSplineCurve3DPtr> curveArraybak(curveArray);

    if (curveArraybak.size() != 4)
        return 0;

    bslpineCurve1 = curveArraybak[0];
    bslpineCurve2 = curveArraybak[1];
    bslpineCurve3 = curveArraybak[2];
    bslpineCurve4 = curveArraybak[3];

    double maxDis = -1;
    double dis;
    sggk::Point3D curve1startpoint = bslpineCurve1->CalcStart();
    sggk::Point3D curve1endpoint = bslpineCurve1->CalcEnd();

    curveArraybak.erase(curveArraybak.begin());
    sggk::Point3D curve2startpoint = bslpineCurve2->CalcStart();
    sggk::Point3D curve2endpoint = bslpineCurve2->CalcEnd();

    sggk::Point3D curve3startpoint = bslpineCurve3->CalcStart();
    sggk::Point3D curve3endpoint = bslpineCurve3->CalcEnd();

    sggk::Point3D curve4startpoint = bslpineCurve4->CalcStart();
    sggk::Point3D curve4endpoint = bslpineCurve4->CalcEnd();

    sggk::Point3D point;
    //fine curve2
    for (int i = 0; i < curveArraybak.size(); ++i)
    {
        bslpineCurve2 = curveArraybak[i];

        curve2startpoint = bslpineCurve2->CalcStart();
        curve2endpoint = bslpineCurve2->CalcEnd();

        if (curve1endpoint.DistanceTo(curve2startpoint) < tol)
        {
            dis = curve1endpoint.DistanceTo(curve2startpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve1endpoint = (curve1endpoint + curve2startpoint) / 2;
                bslpineCurve2->SetControlPoint(curve1endpoint, 1 - 1);
                bslpineCurve1->SetControlPoint(curve1endpoint, bslpineCurve1->ControlPoints().size() - 1);
            }

            curveArraybak.erase(curveArraybak.begin() + i);
            break;
        }
        if (curve1endpoint.DistanceTo(curve2endpoint) < tol)
        {
            dis = curve1endpoint.DistanceTo(curve2endpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve1endpoint = (curve1endpoint + curve2endpoint) / 2;
                bslpineCurve2->SetControlPoint(curve1endpoint, bslpineCurve2->ControlPoints().size() - 1);
                bslpineCurve1->SetControlPoint(curve1endpoint, bslpineCurve1->ControlPoints().size() - 1);
            }

            curveArraybak.erase(curveArraybak.begin() + i);
            bslpineCurve2->Reverse();

            curve2startpoint = bslpineCurve2->CalcStart();
            curve2endpoint = bslpineCurve2->CalcEnd();
            break;
        }
    }

    //fine curve3
    for (int i = 0; i < curveArraybak.size(); ++i)
    {
        bslpineCurve3 = curveArraybak[i];

        curve3startpoint = bslpineCurve3->CalcStart();
        curve3endpoint = bslpineCurve3->CalcEnd();

        if (curve2endpoint.DistanceTo(curve3endpoint) < tol)
        {
            dis = curve2endpoint.DistanceTo(curve3endpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve2endpoint = (curve2endpoint + curve3endpoint) / 2;
                bslpineCurve2->SetControlPoint(curve2endpoint, bslpineCurve2->ControlPoints().size() - 1);
                bslpineCurve3->SetControlPoint(curve2endpoint, bslpineCurve3->ControlPoints().size() - 1);
            }

            curveArraybak.erase(curveArraybak.begin() + i);
            break;
        }
        if (curve2endpoint.DistanceTo(curve3startpoint) < tol)
        {
            dis = curve2endpoint.DistanceTo(curve3startpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve2endpoint = (curve2endpoint + curve3startpoint) / 2;
                bslpineCurve2->SetControlPoint(curve2endpoint, bslpineCurve2->ControlPoints().size() - 1);
                bslpineCurve3->SetControlPoint(curve2endpoint, 1 - 1);
            }

            curveArraybak.erase(curveArraybak.begin() + i);
            bslpineCurve3->Reverse();

            curve3startpoint = bslpineCurve3->CalcStart();
            curve3endpoint = bslpineCurve3->CalcEnd();
            break;
        }
    }

    if (curveArraybak.size() != 1)
        return 0;

    bslpineCurve4 = curveArraybak[0];
    curve4startpoint = bslpineCurve4->CalcStart();
    curve4endpoint = bslpineCurve4->CalcEnd();
    if (curve3startpoint.DistanceTo(curve4endpoint) < tol)
    {
        dis = curve3startpoint.DistanceTo(curve4endpoint);
        if (dis > maxDis)
            maxDis = dis;

        if (isModify && IsLess(dis, tol))
        {
            curve3startpoint = (curve3startpoint + curve4endpoint) / 2;
            bslpineCurve3->SetControlPoint(curve3startpoint, 1 - 1);
            bslpineCurve4->SetControlPoint(curve3startpoint, bslpineCurve4->ControlPoints().size() - 1);
        }
    }
    else
        if (curve3startpoint.DistanceTo(curve4startpoint) < tol)
        {
            dis = curve3startpoint.DistanceTo(curve4startpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve3startpoint = (curve3startpoint + curve4startpoint) / 2;
                bslpineCurve3->SetControlPoint(curve3startpoint, 1 - 1);
                bslpineCurve4->SetControlPoint(curve3startpoint, 1 - 1);
            }

            bslpineCurve4->Reverse();
            curve4startpoint = bslpineCurve4->CalcStart();
            curve4endpoint = bslpineCurve4->CalcEnd();
        }
        else
            return 0;
    return 1;
}

void GuidedCoonsSurfGenerator::Coons_G0(sggk::BSplineCurve3DPtr& curve1, sggk::BSplineCurve3DPtr& curve2, sggk::BSplineCurve3DPtr& curve3, sggk::BSplineCurve3DPtr& curve4, sggk::BSplineSurfacePtr& mySurface_coons)
{
    // First make the opposite boundary curves compatible,
    // Second generate the ruled surface in the u/v direction and of the four corner points
    // Third generate the sum surface 

    //1: Make Curve Compatible
    //curve 1 and curve 3 are opposite curves, make them compatible
    SetSameDistribution(curve1, curve3);
    //curve 2 and curve 4 are opposite curves, make them compatible
    SetSameDistribution(curve2, curve4);

    //2: Generate the coons Bspline surface

    //Get the control point number in the u/v direction
    int NbUPoles = (int)curve1->ControlPoints().size();
    int NbVPoles = (int)curve2->ControlPoints().size();

    //control points of the coons surface
    sggk::Point3DMatrix Poles_result(NbUPoles, sggk::Point3DArray(NbVPoles));

    //control points of ruled surface in the v direction
    sggk::Point3DMatrix Poles_vruled(NbUPoles, sggk::Point3DArray(2));

    //control points of ruled surface in the u direction
    sggk::Point3DMatrix Poles_uruled(2, sggk::Point3DArray(NbVPoles));

    //control points of ruled surface of the four corner points
    sggk::Point3DMatrix Poles_uruled_vruled(2, sggk::Point3DArray(2));

    // Get the u knot vector
    int NbUKnot = (int)curve1->Knots().size();
    sggk::RealArray UKnots = curve1->Knots();
    sggk::UIntArray UMults = curve1->Mults();

    // Get the v knot vector
    int NbVKnot = (int)curve2->Knots().size();
    sggk::RealArray VKnots = curve2->Knots();
    sggk::UIntArray VMults = curve2->Mults();

    // Set v knots of the v ruled surface
    sggk::RealArray VKnots_RuledInVdirection(2);
    sggk::UIntArray VMults_RuledInVdirection(2);
    VKnots_RuledInVdirection[0] = curve2->MinParam();
    VKnots_RuledInVdirection[1] = curve2->MaxParam();

    VMults_RuledInVdirection[0] = 2;
    VMults_RuledInVdirection[1] = 2;

    // Set u knots of the u ruled surface
    sggk::RealArray UKnots_RuledInUdirection(2);
    sggk::UIntArray UMults_RuledInUdirection(2);
    UKnots_RuledInUdirection[0] = curve1->MinParam();
    UKnots_RuledInUdirection[1] = curve1->MaxParam();

    UMults_RuledInUdirection[0] = 2;
    UMults_RuledInUdirection[1] = 2;

    // 
    //2.1 generate the ruled surface in the v direction
    // 
    // 
    // Set control points of the v ruled surface
    for (int i = 0; i < NbUPoles; ++i)
    {
        Poles_vruled[i][0] = curve1->ControlPoints()[i];
        Poles_vruled[i][1] = curve3->ControlPoints()[i];
    }

    // Generate the bspline surface
    sggk::BSplineSurfacePtr mySurface_VRuled = std::make_shared<sggk::BSplineSurface>(
        curve1->Degree(), 1, Poles_vruled,
        UKnots, VKnots_RuledInVdirection,
        UMults, VMults_RuledInVdirection);

    // Make the surface compatible with the result coons surface
    if (curve2->Degree() > mySurface_VRuled->DegreeV())
        mySurface_VRuled->VDegreeElevation(curve2->Degree() - mySurface_VRuled->DegreeV());

    for (int i = 0; i < NbVKnot; ++i)
    {
        unsigned int curMult = SurfaceVKnotMultiplicity(mySurface_VRuled, VKnots[i]);
        if (VMults[i] > curMult)
            mySurface_VRuled->InsertVKnots(VKnots[i], VMults[i] - curMult);
    }

    // 
    //2.2 generate the ruled surface in the u direction
    // 

    // Set control points of the u ruled surface
    for (int i = 0; i < NbVPoles; ++i)
    {
        Poles_uruled[0][i] = curve4->ControlPoints()[i];
        Poles_uruled[1][i] = curve2->ControlPoints()[i];
    }

    sggk::BSplineSurfacePtr mySurface_URuled = std::make_shared<sggk::BSplineSurface>(
        1, curve2->Degree(), Poles_uruled,
        UKnots_RuledInUdirection, VKnots,
        UMults_RuledInUdirection, VMults);

    //increase degree
    if (curve1->Degree() > mySurface_URuled->DegreeU())
        mySurface_URuled->UDegreeElevation(curve1->Degree() - mySurface_URuled->DegreeU());

    for (int i = 0; i < NbUKnot; ++i)
    {
        unsigned int curMult = SurfaceUKnotMultiplicity(mySurface_URuled, UKnots[i]);
        if (UMults[i] > curMult)
            mySurface_URuled->InsertUKnots(UKnots[i], UMults[i] - curMult);
    }

    // 
    //2.3 generate the ruled surface of the four corner points
    // 
    // 
    Poles_uruled_vruled[0][0] = curve1->ControlPoints()[0];
    Poles_uruled_vruled[0][1] = curve3->ControlPoints()[0];
    Poles_uruled_vruled[1][0] = curve1->ControlPoints()[NbUPoles - 1];
    Poles_uruled_vruled[1][1] = curve3->ControlPoints()[NbUPoles - 1];

    sggk::BSplineSurfacePtr mySurface_RuledSurfaceof4CornerPoints = std::make_shared<sggk::BSplineSurface>(
        1, 1, Poles_uruled_vruled,
        UKnots_RuledInUdirection, VKnots_RuledInVdirection,
        UMults_RuledInUdirection, VMults_RuledInVdirection);

    if (curve1->Degree() > mySurface_RuledSurfaceof4CornerPoints->DegreeU())
        mySurface_RuledSurfaceof4CornerPoints->UDegreeElevation(curve1->Degree() - mySurface_RuledSurfaceof4CornerPoints->DegreeU());
    if (curve2->Degree() > mySurface_RuledSurfaceof4CornerPoints->DegreeV())
        mySurface_RuledSurfaceof4CornerPoints->VDegreeElevation(curve2->Degree() - mySurface_RuledSurfaceof4CornerPoints->DegreeV());

    for (int i = 0; i < NbVKnot; ++i)
    {
        unsigned int curMult = SurfaceVKnotMultiplicity(mySurface_RuledSurfaceof4CornerPoints, VKnots[i]);
        if (VMults[i] > curMult)
            mySurface_RuledSurfaceof4CornerPoints->InsertVKnots(VKnots[i], VMults[i] - curMult);
    }

    for (int i = 0; i < NbUKnot; ++i)
    {
        unsigned int curMult = SurfaceUKnotMultiplicity(mySurface_RuledSurfaceof4CornerPoints, UKnots[i]);
        if (UMults[i] > curMult)
            mySurface_RuledSurfaceof4CornerPoints->InsertUKnots(UKnots[i], UMults[i] - curMult);
    }

    //2.4 Generate the sum surface mySurface_VRuled + mySurface_URuled - mySurface_RuledSurfaceof4CornerPoints 
    //increase degree

    for (int i = 0; i < NbUPoles; ++i)
    {
        for (int j = 0; j < NbVPoles; ++j)
        {
            sggk::Point3D p1 = mySurface_VRuled->ControlPoints()[i][j];
            sggk::Point3D p2 = mySurface_URuled->ControlPoints()[i][j];
            sggk::Point3D p3 = mySurface_RuledSurfaceof4CornerPoints->ControlPoints()[i][j];
            Poles_result[i][j] = sggk::Point3D(
                p1.X() + p2.X() - p3.X(),
                p1.Y() + p2.Y() - p3.Y(),
                p1.Z() + p2.Z() - p3.Z());
        }
    }

    mySurface_coons = std::make_shared<sggk::BSplineSurface>(
        curve1->Degree(), curve2->Degree(), Poles_result,
        UKnots, VKnots,
        UMults, VMults);
}

void GuidedCoonsSurfGenerator::TrimInternalCurves(
    std::vector<sggk::BSplineCurve3DPtr>& theInternalBSplineCurves,
    const std::vector<sggk::BSplineCurve3DPtr>& theBoundaryCurveArray,
    double theToleranceDistance) // theToleranceDistance 用于判断交点是否足够接近零距离
{
    //std::vector<sggk::BSplineCurve3DPtr> aTrimmedCurvesResult; // 存储裁剪后的结果

    for (int i = 0; i < theInternalBSplineCurves.size(); i++)
    {
        sggk::BSplineCurve3DPtr aCurrentInternalCurve = theInternalBSplineCurves[i];
        std::vector<std::pair<double, double>> curTrimInterval; // 当前引导线线的裁剪区间

        // 存储所有交点在内部曲线上的参数
        std::vector<double> intersectionParams;
        intersectionParams.push_back(aCurrentInternalCurve->MinParam()); // 曲线起点
        intersectionParams.push_back(aCurrentInternalCurve->MaxParam());  // 曲线终点

        // 遍历所有边界曲线，查找交点
        for (auto& aBoundaryCurve : theBoundaryCurveArray)
        {
            // 注意：OCC 原版用 GeomAPI_ExtremaCurveCurve（极值）并以 theToleranceDistance 作为交点判据，
            // SGK 的线线极值未实现，改用 CrvCrvInt（求交）并显式传入相同容差，以容忍拟合带来的微小间隙。
            sggk::IntCrvCrvRet anExtrema = sggk::GeomInt::CrvCrvInt(*aCurrentInternalCurve, *aBoundaryCurve,
                sggk::CrvCrvIntOpts(sggk::Toler(theToleranceDistance)));

            // 遍历所有交点
            for (const auto& intPnt : anExtrema.IntInfos())
            {
                intersectionParams.push_back(intPnt.param1);
            }
        }

        // 对所有参数进行排序，并去重（避免重复的参数点）
        std::sort(intersectionParams.begin(), intersectionParams.end());
        intersectionParams.erase(std::unique(intersectionParams.begin(), intersectionParams.end(),
            [](double a, double b) { return std::fabs(a - b) < 1e-2; }), // 使用一个小的容差进行去重
            intersectionParams.end());

        // 如果只有起点、终点 + 两个内部交点
        if (intersectionParams.size() == 4)
        {
            double startParam = intersectionParams[1];
            double endParam = intersectionParams[2];

            // 排除两个点靠很近的情况
            sggk::Point3D p0 = aCurrentInternalCurve->CalcPoint(startParam);
            sggk::Point3D p1 = aCurrentInternalCurve->CalcPoint(endParam);
            if (!p0.DistanceTo(p1) < 10)
            {
                // 直接添加，不做包围盒判断
                curTrimInterval.emplace_back(startParam, endParam);
            }
        }
        else
        {
            // 根据交点参数划分曲线段并判断内部性
            for (int j = 0; j < intersectionParams.size() - 1; ++j)
            {
                double startParam = intersectionParams[j];
                double endParam = intersectionParams[j + 1];

                // 排除两个点靠很近的情况
                sggk::Point3D p0 = aCurrentInternalCurve->CalcPoint(startParam);
                sggk::Point3D p1 = aCurrentInternalCurve->CalcPoint(endParam);
                if (p0.DistanceTo(p1) < 10) continue;

                // 裁剪出当前区间段
                sggk::BSplineCurve3DPtr aBSplineSegment = aCurrentInternalCurve->TrimCurve(startParam, endParam)->ToBSpline();

                //std::vector<sggk::BSplineCurve3DPtr> aTempBoundaryCurveArray = theBoundaryCurveArray;
                //if (IsCurveInsideBoundaries(aBSplineSegment, aTempBoundaryCurveArray, theToleranceDistance))
                //{
                //    aTrimmedCurvesResult.push_back(aBSplineSegment);
                //}

                if (IsCurveInsideSurface(aBSplineSegment, m_originalSurf, 10))
                {
                    curTrimInterval.emplace_back(startParam, endParam);
                    //aTrimmedCurvesResult.push_back(aBSplineSegment);
                }
            }
        }

        m_guideCurvesTrimIntervals.emplace_back(curTrimInterval);
    }
    //theInternalBSplineCurves = aTrimmedCurvesResult; // 更新原始向量为裁剪后的曲线
}

// 将B样条曲线离散化为一系列三维点。
std::vector<sggk::Point3D> GuidedCoonsSurfGenerator::DiscretizeBSplineCurve(
    const sggk::BSplineCurve3DPtr& theCurve,
    int numSegments,
    bool theBoundaryFlag)
{
    std::vector<sggk::Point3D> discretizedPoints;
    if (numSegments < 1) return {};

    double firstParam = theCurve->MinParam();
    double lastParam = theCurve->MaxParam();
    double step = (lastParam - firstParam) / numSegments;

    for (int i = 0; i <= numSegments; ++i)
    {
        double currentParam = firstParam + i * step;
        if (i == numSegments) currentParam = lastParam;
        discretizedPoints.push_back(theCurve->CalcPoint(currentParam));
    }

    if (!theBoundaryFlag)
    {
        discretizedPoints.erase(discretizedPoints.end() - 1);
        discretizedPoints.erase(discretizedPoints.begin());
    }
    return discretizedPoints;
}

bool GuidedCoonsSurfGenerator::CurvesConnectedLoop(
    std::vector<sggk::BSplineCurve3DPtr>& theCurves,
    double theTolerance)
{
    if (theCurves.size() < 3)
        return false;

    std::vector<bool> used(theCurves.size(), false);
    std::vector<sggk::BSplineCurve3DPtr> orderedCurves;
    orderedCurves.reserve(theCurves.size());

    // 从第一条曲线开始
    sggk::BSplineCurve3DPtr current = theCurves[0];
    orderedCurves.push_back(current);
    used[0] = true;

    for (size_t i = 1; i < theCurves.size(); ++i)
    {
        sggk::Point3D currEnd = current->CalcPoint(current->MaxParam());
        bool foundNext = false;

        for (size_t j = 0; j < theCurves.size(); ++j)
        {
            if (used[j]) continue;

            sggk::BSplineCurve3DPtr candidate = theCurves[j];
            sggk::Point3D candStart = candidate->CalcPoint(candidate->MinParam());
            sggk::Point3D candEnd = candidate->CalcPoint(candidate->MaxParam());

            // 情况1：起点匹配，不用反转
            if (currEnd.DistanceTo(candStart) < theTolerance)
            {
                orderedCurves.push_back(candidate);
                used[j] = true;
                current = candidate;
                foundNext = true;
                break;
            }
            // 情况2：终点匹配，需要反转
            else if (currEnd.DistanceTo(candEnd) < theTolerance)
            {
                candidate->Reverse();
                orderedCurves.push_back(candidate);
                used[j] = true;
                current = candidate;
                foundNext = true;
                break;
            }
        }

        if (!foundNext)
        {
            // 找不到匹配项
            return false;
        }
    }

    // 最终闭合性检测（最后一条与第一条）
    sggk::Point3D lastEnd = orderedCurves.back()->CalcPoint(orderedCurves.back()->MaxParam());
    sggk::Point3D firstStart = orderedCurves.front()->CalcPoint(orderedCurves.front()->MinParam());

    if (lastEnd.DistanceTo(firstStart) > theTolerance)
    {
        return false; // 无法闭环
    }

    theCurves = orderedCurves; // 更新为首尾相接排序结果
    return true;
}


// 判断一个点是否在由一系列二维点构成的多边形内部（射线法）
bool GuidedCoonsSurfGenerator::IsPointInPolygon2D(
    const sggk::Point2D& theTestPoint,
    const std::vector<sggk::Point2D>& thePolygon2d,
    double theTolerance)
{
    if (thePolygon2d.empty()) return false;

    int intersectCount = 0;
    int n = (int)thePolygon2d.size();

    for (int i = 0; i < n; ++i)
    {
        const sggk::Point2D& p1 = thePolygon2d[i];
        const sggk::Point2D& p2 = thePolygon2d[(i + 1) % n]; // 形成闭合循环

        // 检查测试点是否在边界上（在容差范围内）
        // 可以根据需要增加更精确的点到线段距离判断
        if (std::fabs((p2.Y() - p1.Y()) * theTestPoint.X() - (p2.X() - p1.X()) * theTestPoint.Y() + p2.X() * p1.Y() - p2.Y() * p1.X()) < theTolerance &&
            ((theTestPoint.X() >= std::min(p1.X(), p2.X()) - theTolerance && theTestPoint.X() <= std::max(p1.X(), p2.X()) + theTolerance) ||
                (theTestPoint.Y() >= std::min(p1.Y(), p2.Y()) - theTolerance && theTestPoint.Y() <= std::max(p1.Y(), p2.Y()) + theTolerance)))
        {
            return true; // 点在边界线上，视为在内部
        }

        // 射线法：从 theTestPoint 向正X方向发射射线
        if (((p1.Y() > theTestPoint.Y() && p2.Y() <= theTestPoint.Y()) ||
            (p2.Y() > theTestPoint.Y() && p1.Y() <= theTestPoint.Y())))
        {
            double xIntersect = p1.X() + (theTestPoint.Y() - p1.Y()) / (p2.Y() - p1.Y()) * (p2.X() - p1.X());
            if (xIntersect > theTestPoint.X() - theTolerance) // 允许微小偏差
            {
                intersectCount++;
            }
        }
    }
    return (intersectCount % 2 == 1);
}

bool GuidedCoonsSurfGenerator::IsCurveInsideBoundaries(
    const sggk::BSplineCurve3DPtr& theCurve,
    std::vector<sggk::BSplineCurve3DPtr>& theBoundaryCurveArray,
    double theToleranceDistance)
{
    if (theBoundaryCurveArray.size() < 3)
    {
        return false;
    }

    CurvesConnectedLoop(theBoundaryCurveArray); // 假设这个函数确保了顺序和连接性

    // 获取边界的顶点
    std::vector<sggk::Point3D> boundaryVertices3D;
    for (const auto& boundaryCurve : theBoundaryCurveArray)
    {
        boundaryVertices3D.push_back(boundaryCurve->CalcPoint(boundaryCurve->MinParam()));
    }

    // 假设 boundaryVertices3D 至少有3个点，且它们定义了一个平面。
    // 在实际应用中，您可能需要更健壮的逻辑来处理共线或平面退化的情况。
    if (boundaryVertices3D.size() < 3)
    {
        return false; // 无法定义平面
    }

    sggk::Point3D p1 = boundaryVertices3D[0];
    sggk::Point3D p2 = boundaryVertices3D[1];
    sggk::Point3D p3 = boundaryVertices3D[2];

    // 计算平面的法向量
    sggk::Vector3D v12 = p2 - p1;
    sggk::Vector3D v13 = p3 - p1;
    sggk::Vector3D planeNormal = v12.Cross(v13);

    if (planeNormal.Length() < 1e-12)
    {
        return false;
    }
    planeNormal.Normalize(); // 归一化法向量

    // 建立平面局部二维坐标系基向量
    sggk::Vector3D uBasis = v12; // 沿 P1P2 方向
    if (uBasis.Length() < 1e-12)
    {
        // P1P2 距离太近，尝试用 P1P3
        uBasis = v13;
    }
    uBasis.Normalize();

    sggk::Vector3D vBasis = planeNormal.Cross(uBasis); // 确保与 uBasis 和法向量正交

    std::vector<sggk::Point2D> polygon2D;
    for (const auto& p3d : boundaryVertices3D)
    {
        sggk::Vector3D vecFromP1 = p3d - p1;
        double u = vecFromP1.Dot(uBasis);
        double v = vecFromP1.Dot(vBasis);
        polygon2D.emplace_back(u, v);
    }

    int numSamplePoints = 5;
    std::vector<sggk::Point3D> internalCurveSamplePoints = DiscretizeBSplineCurve(theCurve, numSamplePoints - 1, false);

    for (const sggk::Point3D& samplePoint3D : internalCurveSamplePoints)
    {
        // 将采样点投影到自定义平面
        sggk::Vector3D vecSP1 = samplePoint3D - p1;
        double distToPlane = vecSP1.Dot(planeNormal);
        sggk::Point3D projectedPoint3D = samplePoint3D + planeNormal * (-distToPlane); // 沿着法线方向移动到平面上

        // 将投影点转换到二维局部坐标系
        sggk::Vector3D vecProjP1 = projectedPoint3D - p1;
        double uProj = vecProjP1.Dot(uBasis);
        double vProj = vecProjP1.Dot(vBasis);
        sggk::Point2D samplePoint2D(uProj, vProj);

        // 判断投影点是否在二维多边形内部
        if (!IsPointInPolygon2D(samplePoint2D, polygon2D, theToleranceDistance))
        {
            // 只要有一个采样点不在内部，则整条曲线不在内部
            return false;
        }
    }

    return true; // 所有采样点都在自定义平面内的多边形区域内
}

bool GuidedCoonsSurfGenerator::IsCurveInsideSurface(
    const sggk::BSplineCurve3DPtr& theCurve,
    const sggk::SurfacePtr& theSurface,
    const double theTolerance)
{
    if (!theCurve || !theSurface)
        return false;

    // 构造 Coons 曲面的 AABB
    sggk::BndBox surfaceBox = theSurface->CalcBndBox(sggk::UVRange(theSurface->DomainU(), theSurface->DomainV()));
    surfaceBox.Expand(theTolerance);  // 加一点容差

    // 离散采样
    std::vector<sggk::Point3D> samplePoints = DiscretizeBSplineCurve(theCurve, 50, false);
    int totalCount = static_cast<int>(samplePoints.size());
    int insideCount = 0;

    for (const auto& p : samplePoints)
    {
        if (!surfaceBox.IsPointOut(p))
            ++insideCount;
    }

    double ratio = static_cast<double>(insideCount) / totalCount;
    return (ratio >= 0.7);
}

void GuidedCoonsSurfGenerator::TrimGuideCurves(std::vector<sggk::BSplineCurve3DPtr>& guideBSplineCurves, const std::vector<sggk::BSplineCurve3DPtr>& boundaryCurveArray,
    double toleranceDistance) 
{
    for (int i = 0; i < guideBSplineCurves.size(); ++i)
    {
        // 用于存储最近两条边界曲线的交点和裁剪参数
        sggk::Point3D replacePoints[2];
        double splitParams[2] = { 0 };
        int foundCount = 0;

        // 计算内部曲线与每条边界曲线的距离并排序
        std::vector<std::pair<double, sggk::BSplineCurve3DPtr>> aBoundaryCurves;
        for (auto& boundaryCurve : boundaryCurveArray)
        {
            double distance = ComputeCurveCurveDistance(guideBSplineCurves[i], boundaryCurve);
            aBoundaryCurves.emplace_back(distance, boundaryCurve);
        }

        // 按距离从小到大排序
        std::sort(aBoundaryCurves.begin(), aBoundaryCurves.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

        if (aBoundaryCurves[0].first >= toleranceDistance)
        {
            guideBSplineCurves.erase(guideBSplineCurves.begin() + i);
            i--;
            continue;
        }

        // 找到与内部曲线距离最近的两条边界曲线的交点
        for (size_t j = 0; j < 2 && j < aBoundaryCurves.size(); ++j)
        {
            sggk::IntCrvCrvRet extrema = sggk::GeomInt::CrvCrvInt(*guideBSplineCurves[i], *aBoundaryCurves[j].second);

            const auto& infos = extrema.IntInfos();
            if (!infos.empty())
            {
                splitParams[foundCount] = infos.front().param1;
                replacePoints[foundCount] = infos.front().point;
                foundCount++;
            }
        }

        // 如果没有找到两个有效的交点，则跳过
        if (foundCount != 2)
            continue;

        // 确保裁剪参数按升序排列
        if (splitParams[0] > splitParams[1])
        {
            std::swap(splitParams[0], splitParams[1]);
            std::swap(replacePoints[0], replacePoints[1]);
        }

        if (splitParams[0] == splitParams[1])
        {
            guideBSplineCurves.erase(guideBSplineCurves.begin() + i);
            i--;
            continue;
        }

        // 裁剪内部曲线并更新为新的B样条曲线
        sggk::BSplineCurve3DPtr modifiedCurve = guideBSplineCurves[i]->TrimCurve(splitParams[0], splitParams[1])->ToBSpline();

        // 设置裁剪后的端点
        modifiedCurve->SetControlPoint(replacePoints[0], 1 - 1);
        modifiedCurve->SetControlPoint(replacePoints[1], modifiedCurve->ControlPoints().size() - 1);

        // 更新内部曲线
        guideBSplineCurves[i] = modifiedCurve;
    }
}

double GuidedCoonsSurfGenerator::ComputeCurveCurveDistance(const sggk::BSplineCurve3DPtr& guideCurve, const sggk::BSplineCurve3DPtr& boundaryCurve)
{
    sggk::IntCrvCrvRet aExtrema = sggk::GeomInt::CrvCrvInt(*guideCurve, *boundaryCurve);
    // 有交点则认为距离为 0
    if (!aExtrema.IntInfos().empty())
    {
        return 0.0;
    }

    // 无交点时，通过离散采样计算最小距离
    double aMinDistance = std::numeric_limits<double>::max();
    int numSamples = 100;
    for (int i = 0; i <= numSamples; ++i)
    {
        double t = guideCurve->MinParam() + (guideCurve->MaxParam() - guideCurve->MinParam()) * i / numSamples;
        sggk::Point3D p = guideCurve->CalcPoint(t);
        double aDist = boundaryCurve->CalcDistanceToPt(p);
        if (aDist < aMinDistance)
        {
            aMinDistance = aDist;
        }
    }
    return aMinDistance;
}

void GuidedCoonsSurfGenerator::ApproximateBoundaryCurves(std::vector<sggk::BSplineCurve3DPtr>& curves, int samplingNum)
{
    for (auto& curve : curves)
    {
        const sggk::RealArray& curveKnots = curve->Knots();

        // 重新参数化曲线的节点
        if (!(curveKnots.front() == 0 && curveKnots.back() == 1))
        {
            curve->AdjustKnots(sggk::Interval(0.0, 1.0));
        }

        
        std::vector<sggk::Point3D> samplingPnts;
        std::vector<double> samplingParams;
        // 改为弧长采样
        double totalLength = CurveFair::ComputeCurveLength(curve);
        if (totalLength > 0)
        {
            for (int i = 1; i <= samplingNum; ++i) {
                sggk::Point3D pnt;
                // 首末采样点直接用曲线端点，避免 CalcParaByLength 在弧长==总长时反求失败返回 (0,0,0)
                if (i == 1) {
                    pnt = curve->CalcStart();
                }
                else if (i == samplingNum) {
                    pnt = curve->CalcEnd();
                }
                else {
                    double targetLen = totalLength * (i - 1.0) / (samplingNum - 1.0);
                    curve->CalcParaByLength(targetLen, curve->MinParam(), pnt);
                }
                samplingParams.push_back((i - 1.0) / (samplingNum - 1.0));
                samplingPnts.push_back(pnt);
            }
        }
        else
        {
            double vMin = curve->MinParam();
            double vMax = curve->MaxParam();
            for (int j = 1; j <= samplingNum; ++j)
            {
                double param = vMin + (vMax - vMin) * (j - 1) / (samplingNum - 1);
                sggk::Point3D pnt = curve->CalcPoint(param);
                samplingParams.push_back((j - 1.0) / (samplingNum - 1.0));
                samplingPnts.push_back(pnt);
            }
        }

        //GeomAPI_PointsToBSpline approx(points);
        //curve = approx.Curve();
        
        /*
        // 采样点与参数生成
        sggk::Point3DArray points(1, samplingNum);
        sggk::RealArray params(1, samplingNum);
        std::vector<sggk::Point3D> samplingPnts;
        std::vector<double> samplingParams;
        double vMin = curve->MinParam();
        double vMax = curve->MaxParam();

        for (int j = 1; j <= samplingNum; ++j)
        {
            double param = vMin + (vMax - vMin) * (j - 1) / (samplingNum - 1);
            sggk::Point3D pnt = curve->CalcPoint(param);
            samplingParams.push_back((j - 1.0) / (samplingNum - 1.0));
            samplingPnts.push_back(pnt);
        }
        */
        // 初始化节点并进行拟合
        std::vector<double> init_knots = KnotGernerationByParams(samplingParams, APPROXIMATE_KNOTS_NUM, 3);
        std::vector<double> insertKnots;
        curve = IterateApproximate(insertKnots, samplingPnts, samplingParams, init_knots, 3);
    }
}

sggk::BSplineCurve3DPtr GuidedCoonsSurfGenerator::IterateApproximate(std::vector<double>& insertKnots, const std::vector<sggk::Point3D>& pnts, std::vector<double>& pntsParams,
    std::vector<double>& initKnots, int degree, int maxIterNum, double toler)
{
    int itNum = 1;
    double currentMaxError = 100;
    sggk::BSplineCurve3DPtr IterBspineCurve;
    std::vector<double> CurrentKnots = initKnots;

    while (currentMaxError > toler && itNum <= maxIterNum)
    {
        IterBspineCurve = ApproximateCurve(pnts, pntsParams, CurrentKnots, degree);

        KnotUpdate knotUpdate(IterBspineCurve, CurrentKnots, pnts, pntsParams);
        auto newKnot = knotUpdate.SelfSingleUpdate(PARAM_BASED_BY_INTERVAL_ERROR);
        currentMaxError = knotUpdate.getMaxError();

        if (currentMaxError > toler)
        {
            //update knot vector
            insertKnots.push_back(newKnot);
            CurrentKnots = knotUpdate.getSequences();
        }
        else
        {  
            return IterBspineCurve;
        }
        itNum++;
    }

    return IterBspineCurve;
}

sggk::BSplineCurve3DPtr GuidedCoonsSurfGenerator::ApproximateCurve(const std::vector<sggk::Point3D>& pnts, std::vector<double>& params, std::vector<double>& fKnots, int degree)
{
    std::vector<double> Knots;
    std::vector<int> Mutis;
    SequenceToKnots(fKnots, Knots, Mutis);
    sggk::RealArray KnotsArr(Knots.begin(), Knots.end());
    sggk::UIntArray MutisArr(Mutis.begin(), Mutis.end());
    return ApproximateCurve(pnts, params, KnotsArr, MutisArr, fKnots, degree);
}

sggk::BSplineCurve3DPtr GuidedCoonsSurfGenerator::ApproximateCurve(const std::vector<sggk::Point3D>& pnts, std::vector<double>& pntsParams, sggk::RealArray& knots,
    sggk::UIntArray& mutis, std::vector<double>& fKnots, int degree)
{
    int n = (int)fKnots.size() - degree - 2;
    int m = (int)pnts.size() - 1;

    //1.Construct matrix N
    Eigen::MatrixXd matN(m - 1, n - 1);
    for (int i = 0; i < m - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            double value = CalBasicFunction(pntsParams[i + 1], j + 1, degree, fKnots);
            matN(i, j) = value;
        }
    }
    Eigen::MatrixXd matTransposeN = matN.transpose();
    Eigen::MatrixXd NTN = matTransposeN * matN;
    //2.Construct matrix R
    //2.1 Compute Ri(2-(m-1))
    Eigen::VectorXd VRx(m - 1);
    Eigen::VectorXd VRy(m - 1);
    Eigen::VectorXd VRz(m - 1);
    for (int i = 1; i <= m - 1; ++i) {
        sggk::Vector3D VecTemp = CalResPnt(i, pnts, pntsParams, degree, fKnots, n);
        double x = VecTemp.X(), y = VecTemp.Y(), z = VecTemp.Z();
        VRx(i - 1) = x;
        VRy(i - 1) = y;
        VRz(i - 1) = z;
    }
    //2.2 Compute the component of R
    Eigen::VectorXd Rx = matTransposeN * VRx;
    Eigen::VectorXd Ry = matTransposeN * VRy;
    Eigen::VectorXd Rz = matTransposeN * VRz;

    //3.solve NtN P = R
    Eigen::VectorXd Sx = NTN.lu().solve(Rx);
    Eigen::VectorXd Sy = NTN.lu().solve(Ry);
    Eigen::VectorXd Sz = NTN.lu().solve(Rz);

    //4.construct bspline curve
    sggk::Point3DArray ctrlPnts(n + 1);
    ctrlPnts[0] = pnts[0];
    ctrlPnts[n] = pnts[m];
    for (int i = 1; i < n; i++) {
        sggk::Point3D pntTemp(Sx(i - 1), Sy(i - 1), Sz(i - 1));
        ctrlPnts[i] = pntTemp;
    }
    sggk::BSplineCurve3DPtr bspline = std::make_shared<sggk::BSplineCurve3D>(degree, ctrlPnts, knots, mutis);
    return bspline;
}

sggk::Vector3D GuidedCoonsSurfGenerator::CalResPnt(int k, const std::vector<sggk::Point3D>& dataPoints, const std::vector<double>& parameters, int p, std::vector<double>& knots, int ctrlPntNum)
{
    double aCoeff1 = CalBasicFunction(parameters[k], 0, p, knots);
    double aCoeff2 = CalBasicFunction(parameters[k], ctrlPntNum, p, knots);
    sggk::Vector3D vecTemp0 = dataPoints[0];
    sggk::Vector3D vecTempm = dataPoints[dataPoints.size() - 1];
    sggk::Vector3D vecTempk = dataPoints[k];
    sggk::Vector3D vectemp = vecTempk - aCoeff1 * vecTemp0 - aCoeff2 * vecTempm;
    return vectemp;
}

void GuidedCoonsSurfGenerator::SequenceToKnots(const std::vector<double>& sequence, std::vector<double>& knots, std::vector<int>& multiplicities)
{
    if (sequence.empty()) return;

    std::map<double, int> knotMap;

    // 使用map来统计每个节点的重复次数
    for (double value : sequence) {
        bool found = false;
        for (auto& knot : knotMap) {
            if (IsEqual(value, knot.first)) {
                knot.second++;
                found = true;
                break;
            }
        }
        if (!found) {
            knotMap[value] = 1;
        }
    }

    // 将map的内容转移到knots和multiplicities向量
    for (const auto& knot : knotMap) {
        knots.push_back(knot.first);
        multiplicities.push_back(knot.second);
    }
}

std::vector<double> GuidedCoonsSurfGenerator::KnotGernerationByParams(const std::vector<double>& params, int n, int p)
{
    int m = (int)params.size() - 1;
    double d = (m + 1) / (n - p + 1);
    std::vector<double> Knots(n + p + 2);
    int temp;
    double alpha;
    for (int i = 0; i <= p; ++i)
    {
        Knots[i] = 0.0;
    }
    for (int j = 1; j <= n - p; ++j)
    {
        temp = int(j * d);
        alpha = j * d - temp;
        Knots[p + j] = (1 - alpha) * params[temp - 1] + alpha * params[temp];
    }
    for (int i = n + 1; i <= n + p + 1; ++i)
    {
        Knots[i] = 1;
    }
    return Knots;
}

void GuidedCoonsSurfGenerator::AddDegenerateCurve(std::vector<sggk::BSplineCurve3DPtr>& boundaryCurves)
{
    if (boundaryCurves.size() == 3)
    {
        // 初始化一个向量用于存储每条曲线的交点计数以及对应的样条曲线
        std::vector<std::pair<int, sggk::BSplineCurve3DPtr>> anInterCount =
        {
            {0, boundaryCurves[0]},
            {0, boundaryCurves[1]},
            {0, boundaryCurves[2]}
        };

        // 按交点计数从大到小对曲线进行排序
        std::sort(anInterCount.begin(), anInterCount.end(),
            [](const std::pair<int, sggk::BSplineCurve3DPtr>& curve1,
                const std::pair<int, sggk::BSplineCurve3DPtr>& curve2)
        {
            return curve1.first < curve2.first;
        });

        if (anInterCount[2].first >= 100)
        {
            // 清空并重新调整tempArray，将交点最多的两条边作为u方向的边
            boundaryCurves.clear();
            boundaryCurves.resize(4);
            boundaryCurves[0] = anInterCount[0].second;
            boundaryCurves[2] = anInterCount[1].second;
            boundaryCurves[1] = anInterCount[2].second;

            // 将交点最多的两条边的交点作为退化边
            sggk::Point3D aDegeneratePoint(0, 0, 0);
            if (boundaryCurves[0]->CalcStart().DistanceTo(boundaryCurves[2]->CalcStart()) > 10
                && boundaryCurves[0]->CalcStart().DistanceTo(boundaryCurves[2]->CalcEnd()) > 10)
            {
                aDegeneratePoint = boundaryCurves[0]->CalcEnd();
            }
            else if (boundaryCurves[0]->CalcEnd().DistanceTo(boundaryCurves[2]->CalcStart()) > 10
                && boundaryCurves[0]->CalcEnd().DistanceTo(boundaryCurves[2]->CalcEnd()) > 10)
            {
                aDegeneratePoint = boundaryCurves[0]->CalcStart();
            }

            // 构建退化边
            sggk::Point3DArray poles(2);
            poles[0] = aDegeneratePoint;
            poles[1] = aDegeneratePoint;

            sggk::RealArray knots(2);
            knots[0] = 0.0;
            knots[1] = 1.0;

            sggk::UIntArray multiplicities(2);
            multiplicities[0] = 2;
            multiplicities[1] = 2;

            boundaryCurves[3] = std::make_shared<sggk::BSplineCurve3D>(1, poles, knots, multiplicities);

            // 调整次序，要求首尾相接
            double tol = boundaryCurves[0]->CalcEnd().DistanceTo(boundaryCurves[0]->CalcStart()) / 1000;
            if (boundaryCurves[0]->CalcStart().DistanceTo(boundaryCurves[1]->CalcStart()) < tol)
            {
                boundaryCurves[0]->Reverse();
            }
            if (boundaryCurves[2]->CalcEnd().DistanceTo(boundaryCurves[1]->CalcEnd()) < tol)
            {
                boundaryCurves[2]->Reverse();
            }
        }
        else
        {
            sggk::Point3D pnt1 = boundaryCurves[0]->CalcStart(), pnt2 = boundaryCurves[0]->CalcEnd(), pnt3 = boundaryCurves[1]->CalcStart();
            sggk::Point3D pnt4 = boundaryCurves[1]->CalcEnd(), pnt5 = boundaryCurves[2]->CalcStart(), pnt6 = boundaryCurves[2]->CalcEnd();

            double tol = boundaryCurves[0]->CalcEnd().DistanceTo(boundaryCurves[0]->CalcStart()) / 1000;
            if (boundaryCurves[1]->CalcEnd().DistanceTo(boundaryCurves[0]->CalcEnd()) < tol)
            {
                boundaryCurves[1]->Reverse();
            }
            else if (boundaryCurves[2]->CalcStart().DistanceTo(boundaryCurves[0]->CalcEnd()) < tol)
            {
                std::swap(boundaryCurves[1], boundaryCurves[2]);
            }
            else if (boundaryCurves[2]->CalcEnd().DistanceTo(boundaryCurves[0]->CalcEnd()) < tol)
            {
                std::swap(boundaryCurves[1], boundaryCurves[2]);
                boundaryCurves[1]->Reverse();
            }

            if (boundaryCurves[2]->CalcEnd().DistanceTo(boundaryCurves[1]->CalcEnd()) < tol)
            {
                boundaryCurves[2]->Reverse();
            }

            pnt1 = boundaryCurves[0]->CalcStart(); pnt2 = boundaryCurves[0]->CalcEnd(); pnt3 = boundaryCurves[1]->CalcStart();
            pnt4 = boundaryCurves[1]->CalcEnd(); pnt5 = boundaryCurves[2]->CalcStart(); pnt6 = boundaryCurves[2]->CalcEnd();

            // 三边情况，创建退化边构成四边
            std::vector<sggk::Point3D> boundaryPoints = { boundaryCurves[0]->CalcStart(), boundaryCurves[1]->CalcStart(), boundaryCurves[2]->CalcStart() };

            // 定义边
            sggk::Vector3D line_01(boundaryPoints[1] - boundaryPoints[0]);
            sggk::Vector3D line_12(boundaryPoints[2] - boundaryPoints[1]);
            sggk::Vector3D line_20(boundaryPoints[0] - boundaryPoints[2]);

            auto calculateAngle = [](const sggk::Vector3D& v1, const sggk::Vector3D& v2)
            {
                double dotProduct = v1.Dot(v2);
                double magnitudes = v1.Length() * v2.Length();
                return std::acos(dotProduct / magnitudes);  // 返回角度
            };

            // 计算三个角的夹角
            double angleAtPoint0 = calculateAngle(-line_20, line_01);  // 点0的夹角
            double angleAtPoint1 = calculateAngle(-line_01, line_12);  // 点1的夹角
            double angleAtPoint2 = calculateAngle(-line_12, line_20);  // 点2的夹角

            // 找出最大角度
            double maxAngle = std::max({ angleAtPoint0, angleAtPoint1, angleAtPoint2 });
            int maxAngleIndex = 0;
            if (maxAngle == angleAtPoint1) maxAngleIndex = 1;
            else if (maxAngle == angleAtPoint2) maxAngleIndex = 2;

            // 构造退化边
            auto CreateDegenerateEdge = [](const sggk::Point3D& p1, const sggk::Point3D& p2)
            {
                sggk::Point3DArray poles(2);
                poles[0] = p1;
                poles[1] = p2;

                sggk::RealArray knots(2);
                knots[0] = 0.0;
                knots[1] = 1.0;

                sggk::UIntArray multiplicities(2);
                multiplicities[0] = 2;
                multiplicities[1] = 2;

                return std::make_shared<sggk::BSplineCurve3D>(1, poles, knots, multiplicities);
            };
            if (maxAngleIndex == 0)
            {
                //boundaryCurves.push_back(CreateDegenerateEdge(boundaryPoints[maxAngleIndex]));
                boundaryCurves.push_back(CreateDegenerateEdge(boundaryCurves[2]->CalcEnd(), boundaryCurves[0]->CalcStart()));
            }
            else
            {
                boundaryCurves.insert(boundaryCurves.begin() + maxAngleIndex, CreateDegenerateEdge(boundaryCurves[maxAngleIndex - 1]->CalcEnd(), boundaryCurves[maxAngleIndex]->CalcStart()));
            }
            while (boundaryCurves[3]->CalcStart().DistanceTo(boundaryCurves[3]->CalcEnd()) > 10.0)
            {
                boundaryCurves.insert(boundaryCurves.begin(), boundaryCurves[3]);
                boundaryCurves.pop_back();  // 使用 pop_back 替代 erase
            }
        }
    }
}

void GuidedCoonsSurfGenerator::ConstructSurfWithGuideCrvs()
{
    // 获取采样点
    if (m_iterateCount == 0)
    {
        GetGuideSamples(); // **等分参数的采样（点的个数取决于弧长和步长）
    }

    // 获取偏移量和投影点参数
    std::vector<sggk::Point2D> aPntParams;
    std::vector<sggk::Point3D> anOffsets;
    GetSamplesOffset(aPntParams, anOffsets, true);

    // 若初始已满足容差，则直接返回Coons
    double maxDis = 0;
    for (const sggk::Point3D& aPnt : anOffsets)
    {
        double x = aPnt.X();
        double y = aPnt.Y();
        double z = aPnt.Z();
        double dis = std::sqrt(x * x + y * y + z * z);
        maxDis = std::max(maxDis, dis);
    }
    if (IsLess(maxDis, m_tol))
    {
        m_isDone = true;
        m_guidedSurf = m_originalSurf;
        return;
    }

    std::vector<Eigen::Vector3d> anOffsetsEigen;
    std::vector<double> aPntParamsU, aPntParamsV;
    anOffsetsEigen.reserve(anOffsets.size());
    aPntParamsU.reserve(aPntParams.size());
    aPntParamsV.reserve(aPntParams.size());
    for (const sggk::Point3D& offset : anOffsets)
    {
        anOffsetsEigen.push_back(Eigen::Vector3d(offset.X(), offset.Y(), offset.Z()));
    }
    for (const sggk::Point2D& param : aPntParams)
    {
        aPntParamsU.push_back(param.X());
        aPntParamsV.push_back(param.Y());
    }

    // 获取初始曲面的节点和次数
    const sggk::RealArray& uKnotsOCC = m_originalSurf->KnotsU();
    const sggk::RealArray& vKnotsOCC = m_originalSurf->KnotsV();
    const sggk::UIntArray& uMultsOCC = m_originalSurf->MultsU();
    const sggk::UIntArray& vMultsOCC = m_originalSurf->MultsV();
    const int uDeg = (int)m_originalSurf->DegreeU();
    const int vDeg = (int)m_originalSurf->DegreeV();

    std::vector<double> aUKnots, aVKnots;
    aUKnots.reserve(uKnotsOCC.size());
    aVKnots.reserve(vKnotsOCC.size());
    for (size_t i = 0; i < uKnotsOCC.size(); ++i)
    {
        for (unsigned j = 0; j < uMultsOCC[i]; ++j)
        {
            aUKnots.push_back(uKnotsOCC[i]);
        }
    }
    for (size_t i = 0; i < vKnotsOCC.size(); ++i)
    {
        for (unsigned j = 0; j < vMultsOCC[i]; ++j)
        {
            aVKnots.push_back(vKnotsOCC[i]);
        }
    }

    const int aCtrlPntsUNum = (int)aUKnots.size() - uDeg - 1;
    const int aCtrlPntsVNum = (int)aVKnots.size() - vDeg - 1;
    std::vector<Eigen::Vector3d> aCtrlPnts(aCtrlPntsUNum * aCtrlPntsVNum);

    // 求解偏移曲面控制点
    FitOffsetSurface(anOffsetsEigen, aPntParamsU, aPntParamsV, aUKnots, aVKnots, uDeg, vDeg, aCtrlPnts);

    sggk::Point3DMatrix aCtrlPntsOCC(aCtrlPntsUNum, sggk::Point3DArray(aCtrlPntsVNum));

    for (int i = 0; i < aCtrlPntsUNum; ++i)
    {
        for (int j = 0; j < aCtrlPntsVNum; ++j)
        {
            aCtrlPntsOCC[i][j].SetX(aCtrlPnts[i * aCtrlPntsVNum + j](0));
            aCtrlPntsOCC[i][j].SetY(aCtrlPnts[i * aCtrlPntsVNum + j](1));
            aCtrlPntsOCC[i][j].SetZ(aCtrlPnts[i * aCtrlPntsVNum + j](2));
        }
    }

    // 创建新的引导后的曲面
    sggk::BSplineSurfacePtr offsetSurf = std::make_shared<sggk::BSplineSurface>(
        uDeg, vDeg, aCtrlPntsOCC, uKnotsOCC, vKnotsOCC, uMultsOCC, vMultsOCC);

    // 叠加曲面
    // 共同节点和次数
    const sggk::RealArray& knotsUCommon = m_originalSurf->KnotsU();
    const sggk::RealArray& knotsVCommon = m_originalSurf->KnotsV();
    const sggk::UIntArray& multsUCommon = m_originalSurf->MultsU();
    const sggk::UIntArray& multsVCommon = m_originalSurf->MultsV();
    const int degUCommon = (int)m_originalSurf->DegreeU();
    const int degVCommon = (int)m_originalSurf->DegreeV();

    // 计算控制点 
    const sggk::Point3DMatrix& originalCtrlPnts = m_originalSurf->ControlPoints();
    const int ctrlPntsUNumCom = (int)m_originalSurf->ControlPoints().size();
    const int ctrlPntsVNumCom = (int)m_originalSurf->ControlPoints()[0].size();

    sggk::Point3DMatrix ctrlPntsAdd(ctrlPntsUNumCom, sggk::Point3DArray(ctrlPntsVNumCom));

    for (int i = 0; i < ctrlPntsUNumCom; ++i)
    {
        for (int j = 0; j < ctrlPntsVNumCom; ++j)
        {
            sggk::Point3D aCoord = originalCtrlPnts[i][j] + STEP_LENGTH * (sggk::Vector3D)aCtrlPntsOCC[i][j];
            ctrlPntsAdd[i][j].SetXYZ(aCoord.X(), aCoord.Y(), aCoord.Z());
        }
    }

    m_guidedSurf = std::make_shared<sggk::BSplineSurface>(
        degUCommon, degVCommon, ctrlPntsAdd, knotsUCommon,
        knotsVCommon, multsUCommon, multsVCommon);
    
    /*
    // 0.9 0.1
    const sggk::Point3DMatrix& guidedCtrlPnts = m_guidedSurf->ControlPoints();
    sggk::Point3DMatrix ctrlPntsAdd2(1, ctrlPntsUNumCom, 1, ctrlPntsVNumCom);

    for (int i = 1; i <= ctrlPntsUNumCom; ++i)
    {
        for (int j = 1; j <= ctrlPntsVNumCom; ++j)
        {
            sggk::Point3D aCoord = 0.9 * originalCtrlPnts(i, j).Coord() + 0.1 * guidedCtrlPnts(i, j).Coord();
            ctrlPntsAdd2(i, j).SetCoord(aCoord.X(), aCoord.Y(), aCoord.Z());
        }
    }

    m_guidedSurf = new Geom_BSplineSurface(ctrlPntsAdd2, knotsUCommon,
        knotsVCommon, multsUCommon, multsVCommon, degUCommon, degVCommon);
    */

    // 获取新曲面的采样点对应投影参数和误差
    std::vector<sggk::Point2D> aPntParams2;
    std::vector<sggk::Point3D> anOffsets2;
    GetSamplesOffset(aPntParams2, anOffsets2, false);

    std::vector<double> errorU, errorV;

    double maxErrDis = 0;
    int pos = 0;
    for (const sggk::Point3D& aPnt : anOffsets2)
    {
        double x = aPnt.X();
        double y = aPnt.Y();
        double z = aPnt.Z();
        double dis = std::sqrt(x * x + y * y + z * z);
        maxErrDis = std::max(maxErrDis, dis);
        if (IsGreater(dis, m_tol))
        {
            errorU.push_back(aPntParams2[pos].X());
            errorV.push_back(aPntParams2[pos].Y());
        }
        pos++;
    }
    std::cout << "error params number: " << errorU.size() << std::endl;
    std::cout << "max distance: " << maxErrDis << std::endl;

    // 分区间筛选要插入的节点
    std::sort(errorU.begin(), errorU.end());
    std::sort(errorV.begin(), errorV.end());
    std::vector<double> knotsUToInsert, knotsVToInsert;
    int left = 0, right = 0;
    while (left < errorU.size() && right < errorU.size())
    {
        for (right = left; right < errorU.size(); ++right)
        {
            if (errorU[right] - errorU[left] >= KNOT_INTERVAL_THRESHOLD)
            {
                knotsUToInsert.push_back((errorU[right] + errorU[left]) / 2.0);
                left = right;
                break;
            }
            if (right == errorU.size() - 1)
            {
                knotsUToInsert.push_back((errorU[right] + errorU[left]) / 2.0);
            }
        }
    }
    left = 0;
    right = 0;
    while (left < errorV.size() && right < errorV.size())
    {
        for (right = left; right < errorV.size(); ++right)
        {
            if (errorV[right] - errorV[left] >= KNOT_INTERVAL_THRESHOLD)
            {
                knotsVToInsert.push_back((errorV[right] + errorV[left]) / 2.0);
                left = right;
                break;
            }
            if (right == errorV.size() - 1)
            {
                knotsVToInsert.push_back((errorV[right] + errorV[left]) / 2.0);
            }
        }
    }

    // 插入节点
    if (errorU.size() > 0 || errorV.size() > 0)
    {
        for (int i = 0; i < knotsUToInsert.size(); ++i)
        {
            m_guidedSurf->InsertUKnots(knotsUToInsert[i], 1);
        }
        for (int i = 0; i < knotsVToInsert.size(); ++i)
        {
            m_guidedSurf->InsertVKnots(knotsVToInsert[i], 1);
        }

        m_isDone = false;
    }
    else 
    {
        m_isDone = true; // 所有点都容差内
    }
}

void GuidedCoonsSurfGenerator::GetGuideSamples()
{
    int crvIdx = 0;
    for (const sggk::BSplineCurve3DPtr& guideCurve : m_guideCurves)
    {
        // 得到采样点
        std::vector<sggk::Point3D> tempSamples;
        //samples = SampleGuideCurve(guideCurve, 0, 1, theSamplesNum);
        for (const auto& trimInterval : m_guideCurvesTrimIntervals[crvIdx])
        {
            // 根据长度自适应采样
            double length = CurveFair::ComputeCurveLengthBetweenParameters(guideCurve, trimInterval.first, trimInterval.second);
            int sampleNum = length / GUIDE_SAMPLING_INTERVAL + 1; // 至少采一个点

            tempSamples = SampleGuideCurve(guideCurve, trimInterval.first, trimInterval.second, sampleNum);
            for (const auto& sample : tempSamples)
            {
                m_samples.push_back(sample);
            }
        }
        crvIdx++;
    }
}

void GuidedCoonsSurfGenerator::GetSamplesOffset(std::vector<sggk::Point2D>& thePntParams, std::vector<sggk::Point3D>& theOffsets, bool isOriginal)
{
    // 得到投影点和参数
    std::vector<sggk::Point2D> aPntParams;
    std::vector<sggk::Point3D> projectionPoints;
    if (isOriginal)
    {
        // 初始曲面
        theOffsets = ProjectPntsToSurf(m_samples, projectionPoints, m_originalSurf, thePntParams);
    }
    else
    {
        // 迭代后曲面
        theOffsets = ProjectPntsToSurf(m_samples, projectionPoints, m_guidedSurf, thePntParams);
    }
}

/*
void GuidedCoonsSurfGenerator::GetSamplesOffset(std::vector<sggk::Point2D>& thePntParams, std::vector<sggk::Point3D>& theOffsets, bool isOriginal, int theSamplesNum)
{
    int crvIdx = 0;
    std::vector<sggk::Point3D> allSamples, allProjectionPoints;

    for (const sggk::BSplineCurve3DPtr& guideCurve : m_guideCurves) 
    {
        // 得到采样点
        std::vector<sggk::Point3D> tempSamples, samples, anOffsets;
        // 根据长度自适应采样
        //samples = SampleGuideCurve(guideCurve, 0, 1, theSamplesNum);
        for (const auto& trimInterval : m_guideCurvesTrimIntervals[crvIdx])
        {
            tempSamples = SampleGuideCurve(guideCurve, trimInterval.first, trimInterval.second, theSamplesNum);
            for (const auto& sample : tempSamples)
            {
                samples.push_back(sample);
                allSamples.push_back(sample);
            }
        }
        crvIdx++;

        // 得到投影点和参数
        std::vector<sggk::Point2D> aPntParams;
        std::vector<sggk::Point3D> projectionPoints;
        if (isOriginal)
        {
            // 初始曲面
            anOffsets = ProjectPntsToSurf(samples, projectionPoints, m_originalSurf, aPntParams);
        }
        else
        {
            // 迭代后曲面
            anOffsets = ProjectPntsToSurf(samples, projectionPoints, m_guidedSurf, aPntParams);
        }

        for (const auto& projectionPoint : projectionPoints)
        {
            allProjectionPoints.push_back(projectionPoint);
        }

        for (const sggk::Point2D& aParam : aPntParams) 
        {
            // reserve
            thePntParams.push_back(aParam);
        }

        for (const sggk::Point3D& aPnt : anOffsets)
        {
            theOffsets.push_back(aPnt);
        }
    }

    
    // 导出采样点
    if (isOriginal)
    {
        TopoDS_Compound cp;
        BRep_Builder builder;
        builder.MakeCompound(cp);
        for (const auto& sample : allSamples)
        {
            builder.Add(cp, BRepBuilderAPI_MakeVertex(sample));
        }
        std::string cpName = "data\\Coons\\samples_points.step";
        STEPControl_Writer stepWriter;
        stepWriter.Transfer(cp, STEPControl_AsIs);
        IFSelect_ReturnStatus status = stepWriter.Write(cpName.c_str());
    }

    // 导出投影点
    if (isOriginal)
    {
        TopoDS_Compound cp;
        BRep_Builder builder;
        builder.MakeCompound(cp);
        for (const auto& sample : allProjectionPoints)
        {
            builder.Add(cp, BRepBuilderAPI_MakeVertex(sample));
        }
        std::string cpName = "data\\Coons\\projection_points.step";
        STEPControl_Writer stepWriter;
        stepWriter.Transfer(cp, STEPControl_AsIs);
        IFSelect_ReturnStatus status = stepWriter.Write(cpName.c_str());
    }
    
}
*/

std::vector<sggk::Point3D> GuidedCoonsSurfGenerator::SampleGuideCurve(const sggk::BSplineCurve3DPtr& theCurve, double startParam, double endParam, int theSamplesNum)
{
    std::vector<sggk::Point3D> samples;
    for (int i = 1; i <= theSamplesNum; ++i)
    {
        //double t = theCurve->MinParam() + i * (theCurve->MaxParam() - theCurve->MinParam()) / (theSamplesNum + 1);
        double t = startParam + i * (endParam - startParam) / (theSamplesNum + 1);
        sggk::Point3D sample = theCurve->CalcPoint(t);
        samples.push_back(sample);
    }
    return samples;
}

std::vector<sggk::Point3D> GuidedCoonsSurfGenerator::ProjectPntsToSurf(const std::vector<sggk::Point3D>& thePoints, std::vector<sggk::Point3D>& theProjectionPoints, const sggk::BSplineSurfacePtr& theSurface, std::vector<sggk::Point2D>& thePntParams)
{
    std::vector<sggk::Point3D> offsets;
    for (const auto& point : thePoints)
    {
        // 求点到曲面的最近点及参数
        sggk::Point2D uv;
        sggk::Point3D projectedPoint = theSurface->CalcNearestPnt(point, uv);
        double u = uv.X();
        double v = uv.Y();

        if (IsGreater(u, 0) && IsLess(u, 1) && IsGreater(v, 0) && IsLess(v, 1))
        {
            sggk::Vector3D offset = point - projectedPoint;
            if (std::abs(offset.X()) > MAX_OFFSETDISTANCE ||
                std::abs(offset.Y()) > MAX_OFFSETDISTANCE ||
                std::abs(offset.Z()) > MAX_OFFSETDISTANCE)
            {
                continue;
            }
            offsets.push_back(offset);
            theProjectionPoints.push_back(projectedPoint);
            thePntParams.push_back(uv);
        }
    }

    return offsets;
}

std::vector<sggk::Point3D> GuidedCoonsSurfGenerator::CalOffsets(const std::vector<sggk::Point3D>& theSamples, const std::vector<sggk::Point3D>& theProjections)
{
    std::vector<sggk::Point3D> offsets;
    for (int i = 0; i < theSamples.size(); ++i)
    {
        sggk::Vector3D offset = theSamples[i] - theProjections[i];
        offsets.push_back(offset);
    }
    return offsets;
}

Eigen::MatrixXd ConstructConvMat(int theRow, int theCol, const Eigen::Matrix3d& kernel)
{
    // 输出矩阵的尺寸：展平后的像素向量长度为 theRow * theCol
    // 卷积后输出的尺寸为 (theRow-2) * (theCol-2)
    int outputSize = (theRow - 2) * (theCol - 2);
    int inputSize = theRow * theCol;
    Eigen::MatrixXd convMat = Eigen::MatrixXd::Zero(outputSize, inputSize);

    for (int i = 1; i < theRow - 1; ++i)
    {
        for (int j = 1; j < theCol - 1; ++j)
        {
            int outputIdx = (i - 1) * (theCol - 2) + (j - 1);
            for (int ki = -1; ki <= 1; ++ki)
            {
                for (int kj = -1; kj <= 1; ++kj)
                {
                    int inputRow = i + ki;
                    int inputCol = j + kj;
                    int inputIdx = inputRow * theCol + inputCol;
                    convMat(outputIdx, inputIdx) = kernel(ki + 1, kj + 1);
                }
            }
        }
    }

    return convMat;
}

Eigen::MatrixXd ConstructVariableConvMat(int ctrlPtsUNum, int ctrlPtsVNum, const std::vector<double>& lDiffs, const std::vector<double>& sDiffs)
{
    //
    if (ctrlPtsUNum < 3 || ctrlPtsVNum < 3)
    {
        return Eigen::MatrixXd(0, ctrlPtsUNum * ctrlPtsVNum);
    }

    // U和V方向内部点的数量
    int internalUpoints = ctrlPtsUNum - 2;
    int internalVpoints = ctrlPtsVNum - 2;

    int outputRows = internalUpoints * internalVpoints;
    int inputCols = ctrlPtsUNum * ctrlPtsVNum;

    Eigen::MatrixXd convMat = Eigen::MatrixXd::Zero(outputRows, inputCols);

    for (int i = 1; i <= internalUpoints; ++i)
    {
        for (int j = 1; j <= internalVpoints; ++j)
        {

            // 当前行在输出矩阵中的索引
            int outputIdx = (i - 1) * internalVpoints + (j - 1);

            double lPrev = lDiffs[i - 1];
            double lCurr = lDiffs[i];


            double sPrev = sDiffs[j - 1];
            double sCurr = sDiffs[j];

            double denL = lPrev + lCurr;
            double denS = sPrev + sCurr;

            if (std::abs(denL) < 1e-7) denL = 1.0; // 避免除零，如果lPrev=lCurr=0，则权重为0
            if (std::abs(denS) < 1e-7) denS = 1.0;

            // P_{i-1,j} 的权重
            double wIm1J = lCurr / (2.0 * denL);
            // P_{i+1,j} 的权重
            double wIp1J = lPrev / (2.0 * denL);
            // P_{i,j-1} 的权重
            double wIJm1 = sCurr / (2.0 * denS);
            // P_{i,j+1} 的权重
            double wIJp1 = sPrev / (2.0 * denS);

            // P_center (P_i,j) 在扁平化输入向量中的列索引
            // 在原始的 ctrlPtsUNum x ctrlPtsVNum 网格中，P_i,j 的索引就是 (i,j)。
            int idxCenter = i * ctrlPtsVNum + j;
            convMat(outputIdx, idxCenter) = -1.0;

            // P_{i-1,j} 的列索引
            int idxIm1J = (i - 1) * ctrlPtsVNum + j;
            convMat(outputIdx, idxIm1J) = wIm1J;

            // P_{i+1,j} 的列索引
            int idxIp1J = (i + 1) * ctrlPtsVNum + j;
            convMat(outputIdx, idxIp1J) = wIp1J;

            // P_{i,j-1} 的列索引
            int idxIJm1 = i * ctrlPtsVNum + (j - 1);
            convMat(outputIdx, idxIJm1) = wIJm1;

            // P_{i,j+1} 的列索引
            int idxIJp1 = i * ctrlPtsVNum + (j + 1);
            convMat(outputIdx, idxIJp1) = wIJp1;
        }
    }
    return convMat;
}

std::vector<double> CalGrevilleAbscissae1D(const sggk::RealArray& knots, int ctrlPtsNum, int degree)
{
    std::vector<double> greville_coords(ctrlPtsNum);

    for (int i = 0; i < ctrlPtsNum; ++i)
    {
        double sum = 0.0;
        for (int j = 1; j <= degree; ++j)
        {
            int knotIdx = i + j;
            sum += knots[knotIdx];
        }
        greville_coords[i] = sum / static_cast<double>(degree);
    }
    return greville_coords;
}

std::pair<std::vector<double>, std::vector<double>> CalGrevilleCoordDiffs(const sggk::BSplineSurfacePtr& surface)
{
    // U 方向
    int ctrlPtsUNum = (int)surface->ControlPoints().size();
    const sggk::RealArray& uKnotsComp = surface->KnotsU();
    const sggk::UIntArray& uMults = surface->MultsU();
    std::vector<double> uKnots;
    for (size_t i = 0; i < uKnotsComp.size(); ++i)
        for (unsigned j = 0; j < uMults[i]; ++j)
            uKnots.push_back(uKnotsComp[i]);

    // 计算U方向的Greville坐标
    std::vector<double> grevilleU = CalGrevilleAbscissae1D(uKnots, ctrlPtsUNum, (int)surface->DegreeU());

    std::vector<double> lDiffs; // 存储U方向的差值
    if (ctrlPtsUNum > 1) {
        lDiffs.resize(ctrlPtsUNum - 1);
        for (int i = 0; i < ctrlPtsUNum - 1; ++i) {
            lDiffs[i] = grevilleU[i + 1] - grevilleU[i];
        }
    }

    // V 方向
    int ctrlPtsVNum = (int)surface->ControlPoints()[0].size();
    const sggk::RealArray& vKnotsComp = surface->KnotsV();
    const sggk::UIntArray& vMults = surface->MultsV();
    std::vector<double> vKnots;
    for (size_t i = 0; i < vKnotsComp.size(); ++i)
        for (unsigned j = 0; j < vMults[i]; ++j)
            vKnots.push_back(vKnotsComp[i]);

    // 计算V方向的Greville坐标
    std::vector<double> grevilleV = CalGrevilleAbscissae1D(vKnots, ctrlPtsVNum, (int)surface->DegreeV());

    std::vector<double> sDiffs; // 存储V方向的差值
    if (ctrlPtsVNum > 1) {
        sDiffs.resize(ctrlPtsVNum - 1);
        for (int j = 0; j < ctrlPtsVNum - 1; ++j) {
            sDiffs[j] = grevilleV[j + 1] - grevilleV[j];
        }
    }

    return { lDiffs, sDiffs }; // 返回包含 l 和 s 差值向量的pair
}





Eigen::MatrixXd GuidedCoonsSurfGenerator::ConstructCurveSmoothingMatrix(
    const sggk::BSplineCurve3DPtr& theBSplineCurve,
    int derivative_order,
    double tolerance
) {
    // 创建CurveFair实例并调用ComputeEnergyMatrix
    CurveFair curveFair;

    // 获取曲线的次数
    int degree = theBSplineCurve->Degree();

    // 调用CurveFair的ComputeEnergyMatrix方法
    Eigen::MatrixXd energyMatrix = curveFair.ComputeEnergyMatrix(
        theBSplineCurve,
        degree,
        tolerance
    );

    return energyMatrix;
}

void GuidedCoonsSurfGenerator::ConstructUIsoparam(
    int n_u, int n_v,
    int p_u,
    const std::vector<double>& knots_u,
    double u_value,
    std::vector<Eigen::Triplet<double>>& triplets,
    Eigen::MatrixXd* dense_matrix
) {
    // u-等参线：Q_j = sum_{i=0}^{n_u} N_{i,p_u}(u_value) * P_{ij}
    // 矩阵C的第j行对应等参线控制点Q_j的计算

    // 预先计算所有u方向基函数值
    std::vector<double> basis_values_u(n_u + 1, 0.0);
    for (int i = 0; i <= n_u; ++i) {
        basis_values_u[i] = CalBasicFunction(u_value, i, p_u, knots_u);
    }

    // 构造矩阵C
    for (int j = 0; j <= n_v; ++j) {  // 等参线控制点索引
        for (int i = 0; i <= n_u; ++i) {  // 曲面控制点u方向索引
            double basis_value = basis_values_u[i];

            if (std::abs(basis_value) > 1e-15) {  // 只处理非零基函数值
                // 曲面控制点P_{ij}在P_vec中的索引
                int col_index = GetControlPointIndex(i, j, n_v);

                // 添加到三元组
                triplets.emplace_back(j, col_index, basis_value);

                // 如果提供了稠密矩阵，同时填充
                if (dense_matrix) {
                    (*dense_matrix)(j, col_index) = basis_value;
                }
            }
        }
    }
}

void GuidedCoonsSurfGenerator::ConstructVIsoparam(
    int n_u, int n_v,
    int p_v,
    const std::vector<double>& knots_v,
    double v_value,
    std::vector<Eigen::Triplet<double>>& triplets,
    Eigen::MatrixXd* dense_matrix
) {
    // v-等参线：R_i = sum_{j=0}^{n_v} N_{j,p_v}(v_value) * P_{ij}
    // 矩阵C的第i行对应等参线控制点R_i的计算

    // 预先计算所有v方向基函数值
    std::vector<double> basis_values_v(n_v + 1, 0.0);
    for (int j = 0; j <= n_v; ++j) {
        basis_values_v[j] = CalBasicFunction(v_value, j, p_v, knots_v);
    }

    // 构造矩阵C
    for (int i = 0; i <= n_u; ++i) {  // 等参线控制点索引
        for (int j = 0; j <= n_v; ++j) {  // 曲面控制点v方向索引
            double basis_value = basis_values_v[j];

            if (std::abs(basis_value) > 1e-15) {  // 只处理非零基函数值
                // 曲面控制点P_{ij}在P_vec中的索引
                int col_index = GetControlPointIndex(i, j, n_v);

                // 添加到三元组
                triplets.emplace_back(i, col_index, basis_value);

                // 如果提供了稠密矩阵，同时填充
                if (dense_matrix) {
                    (*dense_matrix)(i, col_index) = basis_value;
                }
            }
        }
    }
}

Eigen::MatrixXd GuidedCoonsSurfGenerator::ConstructDenseMatrix(
    int n_u, int n_v,
    int p_u, int p_v,
    const std::vector<double>& knots_u,
    const std::vector<double>& knots_v,
    int direction,
    double param_value
) {
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::MatrixXd dense_matrix;

    if (direction == 0) {
        // u-等参线：Q = C * P_vec，其中Q是(n_v+1)维，P_vec是(n_u+1)*(n_v+1)维
        int rows = n_v + 1;
        int cols = (n_u + 1) * (n_v + 1);
        dense_matrix = Eigen::MatrixXd::Zero(rows, cols);

        ConstructUIsoparam(n_u, n_v, p_u, knots_u, param_value, triplets, &dense_matrix);
    }
    else {
        // v-等参线：R = C * P_vec，其中R是(n_u+1)维，P_vec是(n_u+1)*(n_v+1)维
        int rows = n_u + 1;
        int cols = (n_u + 1) * (n_v + 1);
        dense_matrix = Eigen::MatrixXd::Zero(rows, cols);

        ConstructVIsoparam(n_u, n_v, p_v, knots_v, param_value, triplets, &dense_matrix);
    }

    return dense_matrix;
}

Eigen::MatrixXd GuidedCoonsSurfGenerator::ConstructBidirectionalSmoothingMatrixImpl(
    const sggk::BSplineSurfacePtr& theBSplineSurface,
    int n_u, int n_v,
    int p_u, int p_v,
    const std::vector<double>& knots_u,
    const std::vector<double>& knots_v,
    const std::vector<double>& u_params,
    const std::vector<double>& v_params,
    int derivative_order,
    double tolerance
) {
    // 计算总的控制点数量和矩阵维度
    int total_control_points = (n_u + 1) * (n_v + 1);

    // 获取Lobatto权重数据（如果需要）
    std::vector<double> u_weights, v_weights;
    if (0) {
        // 获取u方向的权重
        int u_sample_num = u_params.size();
        const QuadratureData* u_lobatto_data = getQuadratureData(u_sample_num);
        if (u_lobatto_data != nullptr) {
            u_weights = u_lobatto_data->weights_01;
        }
        else {
            // 如果获取失败，使用均匀权重
            u_weights.resize(u_sample_num, 1.0);
        }

        // 获取v方向的权重
        int v_sample_num = v_params.size();
        const QuadratureData* v_lobatto_data = getQuadratureData(v_sample_num);
        if (v_lobatto_data != nullptr) {
            v_weights = v_lobatto_data->weights_01;
        }
        else {
            // 如果获取失败，使用均匀权重
            v_weights.resize(v_sample_num, 1.0);
        }
    }
    else {
        // 非Lobatto模式，使用均匀权重
        u_weights.resize(u_params.size(), 1.0);
        v_weights.resize(v_params.size(), 1.0);
    }

    // 初始化总能量矩阵
    Eigen::MatrixXd M_total;
    M_total = Eigen::MatrixXd::Zero(total_control_points, total_control_points);
    std::vector<double> utraces;
    std::vector<double> vtraces;

    // ========================================================================
    // 1. 处理u-等参线（固定u值，v方向变化）
    // ========================================================================

    // 对每条u-等参线进行处理
    for (size_t u_idx = 0; u_idx < u_params.size(); ++u_idx) {
        double u_val = u_params[u_idx];
        double u_weight = u_weights[u_idx];
        // 从曲面中提取u-等参线（固定u，v变化）
        sggk::BSplineCurve3DPtr u_isoline = theBSplineSurface->CalcVCurve(u_val)->ToBSpline();

        // 构造该等参线的光顺能量矩阵 M_v
        Eigen::MatrixXd M_v;
        M_v = ConstructCurveSmoothingMatrix(u_isoline, derivative_order, tolerance);

        // 构造u-等参线的提取矩阵 C_u
        Eigen::MatrixXd C_u;
        C_u = ConstructDenseMatrix(n_u, n_v, p_u, p_v, knots_u, knots_v, 0, u_val);

        // 计算该u-等参线对总能量矩阵的贡献: C_u^T * M_v * C_u
        Eigen::MatrixXd contribution = C_u.transpose() * M_v * C_u;
        double normal_factor = 1.0 / contribution.trace();
        utraces.push_back(contribution.trace());
        M_total += normal_factor * u_weight * contribution;
    }

    // ========================================================================
    // 2. 处理v-等参线（固定v值，u方向变化）
    // ========================================================================

    // 对每条v-等参线进行处理
    for (size_t v_idx = 0; v_idx < v_params.size(); ++v_idx) {
        double v_val = v_params[v_idx];
        double v_weight = v_weights[v_idx];
        // 从曲面中提取v-等参线（固定v，u变化）
        sggk::BSplineCurve3DPtr v_isoline = theBSplineSurface->CalcUCurve(v_val)->ToBSpline();

        // 构造该等参线的光顺能量矩阵 M_u
        Eigen::MatrixXd M_u;
        M_u = ConstructCurveSmoothingMatrix(v_isoline, derivative_order, tolerance);

        // 构造v-等参线的提取矩阵 D_v
        Eigen::MatrixXd D_v;
        D_v = ConstructDenseMatrix(n_u, n_v, p_u, p_v, knots_u, knots_v, 1, v_val);

        // 计算该v-等参线对总能量矩阵的贡献: D_v^T * M_u * D_v
        Eigen::MatrixXd contribution = D_v.transpose() * M_u * D_v;
        double normal_factor = 1.0 / contribution.trace();
        vtraces.push_back(contribution.trace());
        M_total += normal_factor * v_weight * contribution;
    }

     return M_total;
}

Eigen::MatrixXd GuidedCoonsSurfGenerator::ConstructBidirectionalSmoothingMatrix(
    const sggk::BSplineSurfacePtr& theBSplineSurface,
    const std::vector<double>& u_params,
    const std::vector<double>& v_params,
    int derivative_order,
    double tolerance
) {
    // 从曲面对象中提取参数
    int n_u = (int)theBSplineSurface->ControlPoints().size() - 1;
    int n_v = (int)theBSplineSurface->ControlPoints()[0].size() - 1;
    int p_u = (int)theBSplineSurface->DegreeU();
    int p_v = (int)theBSplineSurface->DegreeV();

    // 提取节点向量
    std::vector<double> knots_u, knots_v;

    // 提取U方向节点向量
    const sggk::RealArray& uKnots = theBSplineSurface->KnotsU();
    const sggk::UIntArray& uMults = theBSplineSurface->MultsU();
    for (size_t i = 0; i < uKnots.size(); i++) {
        for (unsigned j = 0; j < uMults[i]; j++) {
            knots_u.push_back(uKnots[i]);
        }
    }

    // 提取V方向节点向量
    const sggk::RealArray& vKnots = theBSplineSurface->KnotsV();
    const sggk::UIntArray& vMults = theBSplineSurface->MultsV();
    for (size_t i = 0; i < vKnots.size(); i++) {
        for (unsigned j = 0; j < vMults[i]; j++) {
            knots_v.push_back(vKnots[i]);
        }
    }

    return ConstructBidirectionalSmoothingMatrixImpl(
        theBSplineSurface, n_u, n_v, p_u, p_v, knots_u, knots_v, u_params, v_params, derivative_order, tolerance
    );
}

//! @brief 在容差意义下比较 x 是否小于 y
//! @param [In] x 第一个数
//! @param [In] y 第二个数
//! @return x 小于 y 则返回true， 否则返回false
bool IsLess(double x, double y, double tol = 1e-12)
{
    return (y - x) > tol;
}

//! @brief 在容差意义下比较 x 是否大于等于 y
//! @param [In] x 第一个数
//! @param [In] y 第二个数
//! @return x 大于等于 y 则返回true， 否则返回false
bool IsGreaterOrEqual(double x, double y, double tol = 1e-12)
{
    return (x - y) > -tol;
}

double CalBasicFunction(double param, int index, int deg, const std::vector<double>& knots)
{
    double nip, uleft, uright, saved, temp;
    int m = (int)knots.size() - 1;
    std::vector<double> N(deg + 1);

    if ((index == 0 && std::fabs(param - knots[0]) < 1e-12) || (index == m - deg - 1 && std::fabs(param - knots[m]) < 1e-12))
    {
        return 1.0;
    }
    if (IsLess(param, knots[index]) || IsGreaterOrEqual(param, knots[index + deg + 1]))
    {
        return 0.0;
    }
    for (int j = 0; j <= deg; ++j)
    {
        if (IsGreaterOrEqual(param, knots[index + j]) && IsLess(param, knots[index + j + 1]))
        {
            N[j] = 1.0;
        }
        else
        {
            N[j] = 0.0;
        }
    }
    for (int k = 1; k <= deg; ++k)
    {
        if (N[0] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = ((param - knots[index]) * N[0]) / (knots[index + k] - knots[index]);
        }
        for (int j = 0; j < deg - k + 1; ++j)
        {
            uleft = knots[index + j + 1];
            uright = knots[index + j + k + 1];
            if (N[j + 1] == 0.0)
            {
                N[j] = saved;
                saved = 0.0;
            }
            else
            {
                temp = N[j + 1] / (uright - uleft);
                N[j] = saved + (uright - param) * temp;
                saved = (param - uleft) * temp;
            }
        }
    }
    nip = N[0];
    return nip;
}

std::vector<double> GenerateSamplePoints(
    int theNumSamples, const std::vector<double>& theKnots, int theDegree)
{
    std::vector<double> samples;
    samples.reserve(theNumSamples);

    // 在节点区间内均匀采样，避开重复节点
    double start = theKnots[theDegree];
    double end = theKnots[theKnots.size() - theDegree - 1];

    for (int i = 0; i < theNumSamples; ++i) {
        double t = start + (end - start) * i / (theNumSamples - 1);
        samples.push_back(t);
    }

    return samples;
}

double GuidedCoonsSurfGeneratorComputeIntegrationWeight(
    double theU, double theV,
    const std::vector<double>& theUSamples,
    const std::vector<double>& theVSamples,
    int theUIndex, int theVIndex)
{
    // 简单的梯形积分权重
    double du = (theUIndex == 0 || theUIndex == (int)theUSamples.size() - 1) ? 0.5 : 1.0;
    double dv = (theVIndex == 0 || theVIndex == (int)theVSamples.size() - 1) ? 0.5 : 1.0;

    return du * dv;
}

double CalBasicFunctionDerivative(
    double theParam, int theIndex,
    int theDegree, const std::vector<double>& theKnots,
    int theDerivOrder)
{
    // B样条基函数导数计算 - 使用递归公式
    if (theDerivOrder == 0) {
        return CalBasicFunction(theParam, theIndex, theDegree, theKnots);
    }

    // 一阶导数公式
    if (theDegree == 0) {
        return 0.0; // 0次B样条导数恒为0
    }

    double left = 0.0, right = 0.0;
    double denomLeft = theKnots[theIndex + theDegree] - theKnots[theIndex];
    double denomRight = theKnots[theIndex + theDegree + 1] - theKnots[theIndex + 1];

    if (denomLeft > 1e-10) {
        left = CalBasicFunctionDerivative(theParam, theIndex, theDegree - 1, theKnots, theDerivOrder - 1);
        left *= theDegree / denomLeft;
    }

    if (denomRight > 1e-10) {
        right = CalBasicFunctionDerivative(theParam, theIndex + 1, theDegree - 1, theKnots, theDerivOrder - 1);
        right *= theDegree / denomRight;
    }

    return left - right;
}

Eigen::VectorXd ComputeBasisFunctions(
    double theParam, int theCtrlPtsNum,
    int theDegree, const std::vector<double>& theKnots,
    int theDerivOrder)
{
    Eigen::VectorXd result = Eigen::VectorXd::Zero(theCtrlPtsNum);

    for (int i = 0; i < theCtrlPtsNum; ++i) {
        if (theDerivOrder == 0) {
            // 零阶导数 - 基函数值
            result(i) = CalBasicFunction(theParam, i, theDegree, theKnots);
        }
        else if (theDerivOrder == 1) {
            // 一阶导数
            result(i) = CalBasicFunctionDerivative(theParam, i, theDegree, theKnots, 1);
        }
        else if (theDerivOrder == 2) {
            // 二阶导数  
            result(i) = CalBasicFunctionDerivative(theParam, i, theDegree, theKnots, 2);
        }
    }

    return result;
}

double ComputeIntegrationWeight(
    double theU, double theV,
    const std::vector<double>& theUSamples,
    const std::vector<double>& theVSamples,
    int theUIndex, int theVIndex)
{
    int numU = (int)theUSamples.size();
    int numV = (int)theVSamples.size();

    // 计算U方向的有效步长
    double uStep = 0.0;
    if (numU == 1) {
        uStep = 1.0; // 单点情况
    }
    else if (theUIndex == 0) {
        // 第一个点：使用右半区间
        uStep = (theUSamples[1] - theUSamples[0]) / 2.0;
    }
    else if (theUIndex == numU - 1) {
        // 最后一个点：使用左半区间  
        uStep = (theUSamples[numU - 1] - theUSamples[numU - 2]) / 2.0;
    }
    else {
        // 内部点：使用左右邻接区间的平均值
        uStep = (theUSamples[theUIndex + 1] - theUSamples[theUIndex - 1]) / 2.0;
    }

    // 计算V方向的有效步长
    double vStep = 0.0;
    if (numV == 1) {
        vStep = 1.0;
    }
    else if (theVIndex == 0) {
        vStep = (theVSamples[1] - theVSamples[0]) / 2.0;
    }
    else if (theVIndex == numV - 1) {
        vStep = (theVSamples[numV - 1] - theVSamples[numV - 2]) / 2.0;
    }
    else {
        vStep = (theVSamples[theVIndex + 1] - theVSamples[theVIndex - 1]) / 2.0;
    }

    // 返回面积微元
    return uStep * vStep;
}

Eigen::MatrixXd ConstructThinPlateMatrix(
    int theCtrlPtsUNum, int theCtrlPtsVNum,
    const std::vector<double>& theUKnots, const std::vector<double>& theVKnots,
    int theDegU, int theDegV)
{
    // 在参数域采样用于数值积分
    int numSamplesU = 200; // 可根据需要调整
    int numSamplesV = 200;

    std::vector<double> uSamples = GenerateSamplePoints(numSamplesU, theUKnots, theDegU);
    std::vector<double> vSamples = GenerateSamplePoints(numSamplesV, theVKnots, theDegV);

    int totalSamples = numSamplesU * numSamplesV;
    int totalCtrlPts = theCtrlPtsUNum * theCtrlPtsVNum;

    Eigen::MatrixXd D_tp = Eigen::MatrixXd::Zero(totalSamples, totalCtrlPts);

    int sampleIndex = 0;

    // 数值积分计算薄板能量
    for (int i = 0; i < numSamplesU; ++i) {
        double u = uSamples[i];

        for (int j = 0; j < numSamplesV; ++j) {
            double v = vSamples[j];

            // 计算基函数在各方向的一阶和二阶导数
            Eigen::VectorXd basisU_0 = ComputeBasisFunctions(u, theCtrlPtsUNum, theDegU, theUKnots, 0);
            Eigen::VectorXd basisU_1 = ComputeBasisFunctions(u, theCtrlPtsUNum, theDegU, theUKnots, 1);
            Eigen::VectorXd basisU_2 = ComputeBasisFunctions(u, theCtrlPtsUNum, theDegU, theUKnots, 2);

            Eigen::VectorXd basisV_0 = ComputeBasisFunctions(v, theCtrlPtsVNum, theDegV, theVKnots, 0);
            Eigen::VectorXd basisV_1 = ComputeBasisFunctions(v, theCtrlPtsVNum, theDegV, theVKnots, 1);
            Eigen::VectorXd basisV_2 = ComputeBasisFunctions(v, theCtrlPtsVNum, theDegV, theVKnots, 2);

            // 构造当前采样点的薄板能量行向量
            Eigen::RowVectorXd row_uu = Eigen::KroneckerProduct(basisU_2.transpose(), basisV_0.transpose());
            Eigen::RowVectorXd row_vv = Eigen::KroneckerProduct(basisU_0.transpose(), basisV_2.transpose());
            Eigen::RowVectorXd row_uv = Eigen::KroneckerProduct(basisU_1.transpose(), basisV_1.transpose());

            // 薄板能量: uu + 2*uv + vv
            Eigen::RowVectorXd thinPlateRow = row_uu + 2.0 * row_uv + row_vv;

            // 乘以积分权重
            double weight = ComputeIntegrationWeight(u, v, uSamples, vSamples, i, j);
            D_tp.row(sampleIndex) = thinPlateRow * weight;

            sampleIndex++;
        }
    }

    return D_tp;
}



// 偏移曲面拟合
void GuidedCoonsSurfGenerator::FitOffsetSurface(const std::vector<Eigen::Vector3d>& theSamplePntOffsets, const std::vector<double>& thePntParamsU, const std::vector<double>& thePntParamsV,
    const std::vector<double>& theUKnots, const std::vector<double>& theVKnots, int theDegU, int theDegV, std::vector<Eigen::Vector3d>& theCtrlPoints)
{
    int aCtrlPtsUNum = (int)theUKnots.size() - theDegU - 1;
    int aCtrlPtsVNum = (int)theVKnots.size() - theDegV - 1;

    // 构建M矩阵: 边界控制点约束
    Eigen::MatrixXd M;
    BuildMatrixConstraint(aCtrlPtsUNum, aCtrlPtsVNum, M); // 约束项为边界控制点
    Eigen::MatrixXd MT = M.transpose();

    // 构造N矩阵：偏移向量能量矩阵
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(thePntParamsU.size(), aCtrlPtsUNum * aCtrlPtsVNum);
    Eigen::MatrixXd Ni(thePntParamsU.size(), aCtrlPtsUNum);
    Eigen::MatrixXd Nj(thePntParamsV.size(), aCtrlPtsVNum);
    double valueTemp = 0.0;
    for (int i = 0; i < thePntParamsU.size(); ++i)
    {
        for (int j = 0; j < aCtrlPtsUNum; ++j)
        {
            valueTemp = CalBasicFunction(thePntParamsU[i], j, theDegU, theUKnots);
            Ni(i, j) = valueTemp;
        }
        for (int j = 0; j < aCtrlPtsVNum; ++j)
        {
            valueTemp = CalBasicFunction(thePntParamsV[i], j, theDegV, theVKnots);
            Nj(i, j) = valueTemp;
        }
    }
    //N = CalKroneckerProduct(Ni, Nj);

    for (int i = 0; i < thePntParamsU.size(); i++) {
        N.row(i) = Eigen::KroneckerProduct(Ni.row(i), Nj.row(i)).eval();
    }


    Eigen::MatrixXd NT = N.transpose();

    Eigen::MatrixXd O = Eigen::MatrixXd::Zero(N.rows(), 3);
    for (int i = 0; i < (int)O.rows(); ++i)
    {
        O(i, 0) = theSamplePntOffsets[i](0);
        O(i, 1) = theSamplePntOffsets[i](1);
        O(i, 2) = theSamplePntOffsets[i](2);
    }

    sggk::Point3DMatrix originalCtrlPnts = m_originalSurf->ControlPoints();
    Eigen::MatrixXd P0 = Eigen::MatrixXd::Zero(aCtrlPtsUNum * aCtrlPtsVNum, 3);
    int t = 0;
    for (int i = 0; i < aCtrlPtsUNum; ++i)
    {
        for (int j = 0; j < aCtrlPtsVNum; ++j)
        {
            P0(t, 0) = originalCtrlPnts[i][j].X();
            P0(t, 1) = originalCtrlPnts[i][j].Y();
            P0(t, 2) = originalCtrlPnts[i][j].Z();
            t++;
        }
    }
    /*
	// 构造Z矩阵：选择矩阵，只对内部控制点进行调整
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(originalCtrlPnts.Size(), originalCtrlPnts.Size());

    for (int i = 0; i < aCtrlPtsUNum; ++i)
    {
        for (int j = 0; j < aCtrlPtsVNum; ++j)
        {
            int k = i * aCtrlPtsVNum + j;  // 展平索引 (i,j) → k

            bool isBoundary = (i == 0) || (i == aCtrlPtsUNum - 1) || (j == 0) || (j == aCtrlPtsVNum - 1);
            if (!isBoundary)
                Z(k, k) = 1.0;  // 内部点
            else
                Z(k, k) = 0.0;  // 边界点
        }
    }

    // 构造D矩阵：光顺能量矩阵（用格雷维尔坐标的laplace平滑）
    auto diffs = CalGrevilleCoordDiffs(m_originalSurf);
    Eigen::MatrixXd D = ConstructVariableConvMat(aCtrlPtsUNum, aCtrlPtsVNum, diffs.first, diffs.second);

    Eigen::MatrixXd A = (1 - FAIRNESS_WEIGHT) * Z * (NT * N) * Z
        + FAIRNESS_WEIGHT * Z * (D.transpose() * D) * Z;
    Eigen::MatrixXd B = (1 - FAIRNESS_WEIGHT) * Z * (NT * O)
        - FAIRNESS_WEIGHT * Z * (D.transpose() * D) * P0;
    Eigen::MatrixXd P = A.ldlt().solve(B);

    for (int i = 0; i < aCtrlPtsUNum * aCtrlPtsVNum; ++i)
    {
        theCtrlPoints[i] = Eigen::Vector3d(P(i, 0), P(i, 1), P(i, 2));
        //std::cout << P(i, 0) << " " << P(i, 1) << " " << P(i, 2) << std::endl;
    }
    */
    
    // 构建光顺能量矩阵
    Eigen::MatrixXd H, b;
    Eigen::SparseMatrix<double> H_sparse, b_sparse;
    bool isCurveFair = true;

    if (!isCurveFair)
    {

        // 构造D矩阵：光顺能量矩阵（用格雷维尔坐标的laplace平滑）
        auto diffs = CalGrevilleCoordDiffs(m_originalSurf);
        Eigen::MatrixXd D = ConstructVariableConvMat(aCtrlPtsUNum, aCtrlPtsVNum, diffs.first, diffs.second);
        H = (1 - FAIRNESS_WEIGHT) * (NT * N) + FAIRNESS_WEIGHT * D.transpose() * D;
        b = (1 - FAIRNESS_WEIGHT) * (N.transpose() * O) - FAIRNESS_WEIGHT * (D.transpose() * D * P0);

        H_sparse = H.sparseView();
        b_sparse = b.sparseView();

        /*
        // 替换原来的 ConstructVariableConvMat 调用
        Eigen::MatrixXd D_tp = ConstructThinPlateMatrix(
            aCtrlPtsUNum, aCtrlPtsVNum, theUKnots, theVKnots, theDegU, theDegV);


        // 添加正则化项 - 关键修复
        double regularization = 1e-8;

        H = (1 - FAIRNESS_WEIGHT) * (NT * N) + FAIRNESS_WEIGHT * (D_tp.transpose() * D_tp);
        Eigen::MatrixXd regMatrix = regularization * Eigen::MatrixXd::Identity(H.rows(), H.cols());
        H += regMatrix;
        b = (1 - FAIRNESS_WEIGHT) * (N.transpose() * O);
        */
    }
    else 
    {
        // 在0-1之间选取10个参数值（不包含0和1）
        std::vector<double> u_params, v_params;
        //for (int i = 1; i <= 100; i++) {
        //    double param = static_cast<double>(i) / 101.0; // 生成1/11, 2/11, ..., 10/11
        //    u_params.push_back(param);
        //    v_params.push_back(param);
        //}
        if (0) {
            // 使用Lobatto节点进行采样
            int sample_num = 50; // 最多50个节点

            // 确保sample_num在有效范围内（Lobatto数据支持3-50个节点）
            sample_num = std::max(3, sample_num);

            // 获取Lobatto积分数据
            const QuadratureData* lobatto_data = getQuadratureData(sample_num);

            if (lobatto_data != nullptr) {
                // 使用[0,1]区间的Lobatto节点
                u_params = lobatto_data->nodes_01;
                v_params = lobatto_data->nodes_01;

                std::cout << "使用Lobatto采样，节点数量: " << sample_num << std::endl;
            }
            else {
                // 如果获取Lobatto数据失败，回退到均匀采样
                std::cout << "警告：无法获取Lobatto数据（节点数=" << sample_num << "），使用均匀采样" << std::endl;
                for (int i = 1; i <= 100; i++) {
                    double param = static_cast<double>(i) / 101.0; // 生成1/101, 2/101, ..., 100/101
                    u_params.push_back(param);
                    v_params.push_back(param);
                }
            }
        }
        else {
            // 使用均匀采样（原始方法）
            for (int i = 1; i <= 100; i++) {
                double param = static_cast<double>(i) / 101.0; // 生成1/101, 2/101, ..., 100/101
                u_params.push_back(param);
                v_params.push_back(param);
            }
        }

        auto start = std::chrono::high_resolution_clock::now();
        // 构造D矩阵：光顺能量矩阵（等参线光顺）
        auto D = ConstructBidirectionalSmoothingMatrix(
            m_originalSurf,
            u_params,
            v_params,
            2,      // derivative_order = 2 (默认值)
            1e-6    // tolerance = 1e-6 (默认值)
        );
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "smoothing matrix calculation time: " << elapsed.count() / 1000 << " s\n";
        //double trace_D = D.trace();
        //Eigen::MatrixXd NTN = N.transpose() * N;
        //double trace_NTN = NTN.trace();

        // 归一化
        //double lambda_trace = 1.0; // 默认值
        //if (IsGreater(std::abs(trace_D), 1e-12))
        //{
        //    lambda_trace = trace_NTN / trace_D;
        //}
        //D = D * lambda_trace;

        //Eigen::MatrixXd H = (1 - FAIRNESS_WEIGHT) * (NT * N) + FAIRNESS_WEIGHT * D;
        //Eigen::MatrixXd b = (1 - FAIRNESS_WEIGHT) * (N.transpose() * O) - FAIRNESS_WEIGHT * (D * P0);

        //Eigen::MatrixXd I = Eigen::MatrixXd::Identity(H.rows(), H.cols());
        //H = H + 0.00001 * I;


        Eigen::SparseMatrix<double> N_sparse = N.sparseView();
        Eigen::SparseMatrix<double> NT_sparse = N.transpose().sparseView();
        Eigen::SparseMatrix<double> D_sparse = D.sparseView();
        Eigen::SparseMatrix<double> DT_sparse = D.transpose().sparseView();
        Eigen::SparseMatrix<double> O_sparse = O.sparseView();
        Eigen::SparseMatrix<double> OT_sparse = O.transpose().sparseView();
        Eigen::SparseMatrix<double> P0_sparse = P0.sparseView();
        Eigen::SparseMatrix<double> P0T_sparse = P0.transpose().sparseView();

        H_sparse = (1 - FAIRNESS_WEIGHT) * (NT_sparse * N_sparse) + FAIRNESS_WEIGHT * D_sparse;
        b_sparse = (1 - FAIRNESS_WEIGHT) * (NT_sparse * O_sparse) - FAIRNESS_WEIGHT * (D_sparse * P0_sparse);

        H = H_sparse.toDense();
        b = b_sparse.toDense();
    }
    

    // 对曲面光顺
    //Eigen::MatrixXd I = Eigen::MatrixXd::Identity(D.rows(), D.cols());
    //Eigen::MatrixXd H = (1 - FAIRNESS_WEIGHT) * I + FAIRNESS_WEIGHT * D;
    //Eigen::MatrixXd b = - FAIRNESS_WEIGHT * (D * P0);

    // -------------------------------------------------------------------------
    // 优化：通过消除边界控制点 (Static Condensation) 求解内部控制点的位移
    // -------------------------------------------------------------------------
    // M 矩阵的每一行对应一个边界控制点的约束。
    // 在这里，我们的约束是所有边界控制点的位移为 0，即 delta_P_boundary = 0
    // 因此，我们可以将系统分解为：
    // [ H_ii H_ib ] [ P_i ] = [ b_i ]
    // [ H_bi H_bb ] [ P_b ]   [ b_b ]
    // 其中 i 代表内部点(free)，b 代表边界点(fixed=0)
    // 既然 P_b = 0，那么第一行等式直接变成 H_ii * P_i = b_i
    // 这就是一个对称正定（SPD）的系统，使用 SimplicialLDLT 或 SparseLU 求解非常快

    int totalPts = aCtrlPtsUNum * aCtrlPtsVNum;

    // 1. 找出所有内部控制点的索引
    std::vector<int> innerIndices;
    std::vector<int> boundaryIndices;
    for (int u = 0; u < aCtrlPtsUNum; ++u) {
        for (int v = 0; v < aCtrlPtsVNum; ++v) {
            int idx = u * aCtrlPtsVNum + v;
            bool isBoundary = (u == 0 || u == aCtrlPtsUNum - 1 || v == 0 || v == aCtrlPtsVNum - 1);
            if (!isBoundary) {
                innerIndices.push_back(idx);
            } else {
                boundaryIndices.push_back(idx);
            }
        }
    }

    int numInner = innerIndices.size();

    if (numInner > 0) {
        // 2. 提取子矩阵 H_ii 和 子向量 b_i
        Eigen::SparseMatrix<double> H_ii(numInner, numInner);
        Eigen::MatrixXd b_i(numInner, 3);

        std::vector<Eigen::Triplet<double>> H_ii_triplets;

        // 遍历 H_sparse，提取内部点到内部点的块
        for (int k=0; k < H_sparse.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(H_sparse,k); it; ++it) {
                int row = it.row();
                int col = it.col();

                // 查找原索引在 innerIndices 中的位置
                auto itRow = std::lower_bound(innerIndices.begin(), innerIndices.end(), row);
                auto itCol = std::lower_bound(innerIndices.begin(), innerIndices.end(), col);

                if (itRow != innerIndices.end() && *itRow == row &&
                    itCol != innerIndices.end() && *itCol == col)
                {
                    int newRow = std::distance(innerIndices.begin(), itRow);
                    int newCol = std::distance(innerIndices.begin(), itCol);
                    H_ii_triplets.push_back(Eigen::Triplet<double>(newRow, newCol, it.value()));
                }
            }
        }
        H_ii.setFromTriplets(H_ii_triplets.begin(), H_ii_triplets.end());

        // 提取 b_i
        Eigen::MatrixXd b_dense = b_sparse.toDense();
        for (int i = 0; i < numInner; ++i) {
            b_i.row(i) = b_dense.row(innerIndices[i]);
        }

        // 3. 求解 H_ii * P_i = b_i (因为 H_ii 是 SPD 的，使用 SimplicialLDLT)
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(H_ii);

        if (solver.info() != Eigen::Success) {
            std::cout << "Warning: SimplicialLDLT failed, falling back to SparseLU" << std::endl;
            Eigen::SparseLU<Eigen::SparseMatrix<double>> luSolver;
            luSolver.compute(H_ii);

            Eigen::MatrixXd P_i_x = luSolver.solve(b_i.col(0));
            Eigen::MatrixXd P_i_y = luSolver.solve(b_i.col(1));
            Eigen::MatrixXd P_i_z = luSolver.solve(b_i.col(2));

            // 4. 将解写回总的控制点数组 (边界点保持为0)
            for (int i = 0; i < numInner; ++i) {
                int originalIdx = innerIndices[i];
                theCtrlPoints[originalIdx] = Eigen::Vector3d(P_i_x(i,0), P_i_y(i,0), P_i_z(i,0));
            }
        } else {
            Eigen::MatrixXd P_i_x = solver.solve(b_i.col(0));
            Eigen::MatrixXd P_i_y = solver.solve(b_i.col(1));
            Eigen::MatrixXd P_i_z = solver.solve(b_i.col(2));

            // 4. 将解写回总的控制点数组 (边界点保持为0)
            for (int i = 0; i < numInner; ++i) {
                int originalIdx = innerIndices[i];
                theCtrlPoints[originalIdx] = Eigen::Vector3d(P_i_x(i,0), P_i_y(i,0), P_i_z(i,0));
            }
        }

        // 边界点位移强制设为0 (实际上它们初始化时就是0了，为了保险再写一遍)
        for(int idx : boundaryIndices) {
            theCtrlPoints[idx] = Eigen::Vector3d(0,0,0);
        }

    } else {
        // 如果没有内部控制点（全为边界），则位移全为0
        for (int i = 0; i < totalPts; ++i) {
            theCtrlPoints[i] = Eigen::Vector3d(0, 0, 0);
        }
    }


/*
// 构建W矩阵：设置光顺能量的参数alpha
Eigen::MatrixXd aMatrixW;
BuildMatrixWeight((int)thePntParamsU.size(), aCtrlPtsUNum, aCtrlPtsVNum, FAIRNESS_WEIGHT, aMatrixW);

    // 计算N矩阵：投影点参数的基函数系数矩阵
    Eigen::MatrixXd aMatrixN;
    BuildMatrixUnconstraint(thePntParamsU, thePntParamsV, theUKnots, theVKnots, theDegU, theDegV, aCtrlPtsUNum, aCtrlPtsVNum, aMatrixN);

    // 构建O矩阵：偏移向量矩阵
    Eigen::MatrixXd theMatrixO = Eigen::MatrixXd::Zero(aMatrixN.rows(), 3);
    for (int i = (int)theMatrixO.rows() - (int)theSamplePntOffsets.size(); i < (int)theMatrixO.rows(); ++i)
    {
        theMatrixO(i, 0) = theSamplePntOffsets[i - theMatrixO.rows() + theSamplePntOffsets.size()](0);
        theMatrixO(i, 1) = theSamplePntOffsets[i - theMatrixO.rows() + theSamplePntOffsets.size()](1);
        theMatrixO(i, 2) = theSamplePntOffsets[i - theMatrixO.rows() + theSamplePntOffsets.size()](2);
    }

    // 改为叠加曲面光顺

    auto diffs = CalGrevilleCoordDiffs(m_originalSurf);
    Eigen::MatrixXd I = ConstructVariableConvMat(aCtrlPtsUNum, aCtrlPtsVNum, diffs.first, diffs.second);
    Eigen::MatrixXd O0 = Eigen::MatrixXd::Zero(theMatrixO.rows(), I.cols());
    O0.block(0, 0, I.rows(), I.cols()) = I;
    
    sggk::Point3DMatrix originalCtrlPnts = m_originalSurf->ControlPoints();
    Eigen::MatrixXd theCtrlPoints0 = Eigen::MatrixXd::Zero(originalCtrlPnts.Size(), 3);
    int t = 0;
    for (int i = 1; i <= aCtrlPtsUNum; ++i)
    {
        for (int j = 1; j <= aCtrlPtsVNum; ++j)
        {
            theCtrlPoints0(t, 0) = originalCtrlPnts(i, j).X();
            theCtrlPoints0(t, 1) = originalCtrlPnts(i, j).Y();
            theCtrlPoints0(t, 2) = originalCtrlPnts(i, j).Z();
            t++;
        }
    }

    theMatrixO = theMatrixO - O0 * theCtrlPoints0;
    
    Eigen::MatrixXd MT = aMatrixM.transpose();
    //Eigen::MatrixXd NTWN = aMatrixN.transpose() * aMatrixW * aMatrixN;
    //Eigen::MatrixXd NTWO = aMatrixN.transpose() * aMatrixW * theMatrixO;
    // 改为稀疏矩阵
    Eigen::SparseMatrix<double> N_sparse = aMatrixN.sparseView();
    Eigen::SparseMatrix<double> NT_sparse = aMatrixN.transpose().sparseView();
    Eigen::SparseMatrix<double> W_sparse = aMatrixW.sparseView();
    Eigen::SparseMatrix<double> O_sparse = theMatrixO.sparseView();

    Eigen::SparseMatrix<double> NTWN_sparse = NT_sparse * W_sparse * N_sparse;
    Eigen::SparseMatrix<double> NTWO_sparse = NT_sparse * W_sparse * O_sparse;

    Eigen::MatrixXd NTWN = NTWN_sparse.toDense();
    Eigen::MatrixXd NTWO = NTWO_sparse.toDense();

    // 写成分块矩阵形式
    //[ NTWN MT ] [ P ] = [ NTWO ]
    //[  M   0  ] [ L ] = [   0  ]
    // 
    // 即 Ax = B

    int rowsA = (int)NTWN.rows() + (int)aMatrixM.rows();
    int colsA = (int)NTWN.cols() + (int)MT.cols();
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rowsA, colsA);
    A.block(0, 0, NTWN.rows(), NTWN.cols()) = NTWN;
    A.block(0, NTWN.cols(), MT.rows(), MT.cols()) = MT;
    A.block(NTWN.rows(), 0, aMatrixM.rows(), aMatrixM.cols()) = aMatrixM;

    int rowsB = (int)A.rows();
    int colsB = (int)NTWO.cols();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(rowsB, colsB);
    B.block(0, 0, NTWO.rows(), NTWO.cols()) = NTWO;

    // QR分解
    Eigen::VectorXd VBx(B.rows());
    Eigen::VectorXd VBy(B.rows());
    Eigen::VectorXd VBz(B.rows());
    for (int i = 0; i < (int)B.rows(); ++i) {
        VBx(i) = B(i, 0);
        VBy(i) = B(i, 1);
        VBz(i) = B(i, 2);
    }

    //Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(A);
    //Eigen::VectorXd Sx = qr.solve(VBx);
    //Eigen::VectorXd Sy = qr.solve(VBy);
    //Eigen::VectorXd Sz = qr.solve(VBz);

    // 改为稀疏矩阵求解：最小二乘共轭梯度法
    //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(A.sparseView());
    Eigen::VectorXd Sx = solver.solve(VBx);
    Eigen::VectorXd Sy = solver.solve(VBy);
    Eigen::VectorXd Sz = solver.solve(VBz);

    for (int i = 0; i < aCtrlPtsUNum * aCtrlPtsVNum; ++i)
    {
        theCtrlPoints[i] = Eigen::Vector3d(Sx(i), Sy(i), Sz(i));
    }
    */
}

// 构建约束项矩阵
void GuidedCoonsSurfGenerator::BuildMatrixConstraint(int theCtrlPtsUNum, int theCtrlPtsVNum, Eigen::MatrixXd& theMatrixM)
{
    // 计算矩阵 M 的尺寸
    int numRows = 2 * (theCtrlPtsUNum + theCtrlPtsVNum - 2);
    int numCols = theCtrlPtsUNum * theCtrlPtsVNum;

    // 初始化矩阵 M
    theMatrixM = Eigen::MatrixXd::Zero(numRows, numCols);

    // 定义n+1阶单位矩阵 I
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(theCtrlPtsVNum, theCtrlPtsVNum);

    // 定义J
    // 1 0 · · · 0 0
    // 0 0 · · · 0 1
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2, theCtrlPtsVNum);
    J(0, 0) = 1;
    J(1, theCtrlPtsVNum - 1) = 1;

    // 填充矩阵 M 中 I 的左上部分
    theMatrixM.block(0, 0, theCtrlPtsVNum, theCtrlPtsVNum) = I;
    // 填充矩阵 M 中 I 的右下部分
    int rowstart = theCtrlPtsVNum + 2 * theCtrlPtsUNum - 4;
    int colstart = (theCtrlPtsUNum - 1) * theCtrlPtsVNum;
    int blockrows = theCtrlPtsVNum;
    int blockcols = theCtrlPtsVNum;
    theMatrixM.block(rowstart, colstart, blockrows, blockcols) = I;

    // 填充矩阵 M 中 J 的部分
    for (int i = 0; i < theCtrlPtsUNum - 2; i++)
    {
        rowstart = theCtrlPtsVNum + i * 2;
        colstart = theCtrlPtsVNum + i * theCtrlPtsVNum;
        blockrows = 2;
        blockcols = theCtrlPtsVNum;
        theMatrixM.block(rowstart, colstart, blockrows, blockcols) = J;
    }

}

// 构建非约束项矩阵
void GuidedCoonsSurfGenerator::BuildMatrixUnconstraint(const std::vector<double>& thePntParamsU, const std::vector<double>& thePntParamsV, const std::vector<double>& theUKnots, const std::vector<double>& theVKnots,
    int theDegU, int theDegV, int theCtrlPtsUNum, int theCtrlPtsVNum, Eigen::MatrixXd& theMatrixN)
{
    // 构造矩阵 N^
    Eigen::MatrixXd Nhat;
    Eigen::MatrixXd Ni(thePntParamsU.size(), theCtrlPtsUNum);
    Eigen::MatrixXd Nj(thePntParamsV.size(), theCtrlPtsVNum);
    double valueTemp = 0.0;
    for (int i = 0; i < thePntParamsU.size(); ++i)
    {
        for (int j = 0; j < theCtrlPtsUNum; ++j)
        {
            valueTemp = CalBasicFunction(thePntParamsU[i], j, theDegU, theUKnots);
            Ni(i, j) = valueTemp;
        }
        for (int j = 0; j < theCtrlPtsVNum; ++j)
        {
            valueTemp = CalBasicFunction(thePntParamsV[i], j, theDegV, theVKnots);
            Nj(i, j) = valueTemp;
        }
    }
    Nhat = CalKroneckerProduct(Ni, Nj);

    // 以下光顺能量可选其一

    // 论文中光顺能量
    //Eigen::Matrix3d kernel;
    //kernel << 0, -1, 0,
    //         -1, 4, -1,
    //          0, -1, 0;
    //Eigen::MatrixXd I = ConstructConvMat(theCtrlPtsUNum, theCtrlPtsVNum, kernel);

    // 结合格雷维尔横坐标的光顺能量
    auto diffs = CalGrevilleCoordDiffs(m_originalSurf);
    Eigen::MatrixXd I = ConstructVariableConvMat(theCtrlPtsUNum, theCtrlPtsVNum, diffs.first, diffs.second);

    /*
    // 构造矩阵 N
    theMatrixN = Eigen::MatrixXd::Zero(I.rows() + Nhat.rows(), I.cols());
    theMatrixN.block(0, 0, I.rows(), I.cols()) = I;
    theMatrixN.block(I.rows(), 0, Nhat.rows(), Nhat.cols()) = Nhat;
    */

}

// 构建能量权重矩阵 
void GuidedCoonsSurfGenerator::BuildMatrixWeight(int thePntParamsSize, int theCtrlPtsUNum, int theCtrlPtsVNum, double alpha, Eigen::MatrixXd& theMatrixW)
{
    int size = (theCtrlPtsUNum - 2) * (theCtrlPtsVNum - 2);
    Eigen::MatrixXd Alpha = Eigen::MatrixXd::Zero(size, size);
    for (int i = 0; i < size; ++i)
    {
        Alpha(i, i) = alpha;
    }

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(thePntParamsSize, thePntParamsSize);
    for (int i = 0; i < thePntParamsSize; ++i)
    {
        I(i, i) = 1 - alpha;
    }

    theMatrixW = Eigen::MatrixXd::Zero(Alpha.rows() + I.rows(), Alpha.rows() + I.rows());
    theMatrixW.block(0, 0, Alpha.rows(), Alpha.cols()) = Alpha;
    theMatrixW.block(Alpha.rows(), Alpha.cols(), I.rows(), I.cols()) = I;

}

double GuidedCoonsSurfGenerator::CalBasicFunction(double param, int index, int deg, const std::vector<double>& knots)
{
    double nip, uleft, uright, saved, temp; 
    int m = (int)knots.size() - 1;
    std::vector<double> N(deg + 1);

    if ((index == 0 && IsEqual(param, knots[0])) || (index == m - deg - 1 && IsEqual(param, knots[m])))
    {
        return 1.0;
    }
    if (IsLess(param, knots[index]) || IsGreaterOrEqual(param, knots[index + deg + 1]))
    {
        return 0.0;
    }
    for (int j = 0; j <= deg; ++j)
    {
        if (IsGreaterOrEqual(param, knots[index + j]) && IsLess(param, knots[index + j + 1]))
        {
            N[j] = 1.0;
        }
        else
        {
            N[j] = 0.0;
        }
    }
    for (int k = 1; k <= deg; ++k)
    {
        if (N[0] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = ((param - knots[index]) * N[0]) / (knots[index + k] - knots[index]);
        }
        for (int j = 0; j < deg - k + 1; ++j)
        {
            uleft = knots[index + j + 1];
            uright = knots[index + j + k + 1];
            if (N[j + 1] == 0.0)
            {
                N[j] = saved;
                saved = 0.0;
            }
            else
            {
                temp = N[j + 1] / (uright - uleft);
                N[j] = saved + (uright - param) * temp;
                saved = (param - uleft) * temp;
            }
        }
    }
    nip = N[0];
    return nip;
}

Eigen::MatrixXd GuidedCoonsSurfGenerator::CalKroneckerProduct(const Eigen::MatrixXd& theMatA, const Eigen::MatrixXd& theMatB)
{
    int col_num1, col_num2, row_num1, row_num2;
    col_num1 = (int)theMatA.cols();
    row_num1 = (int)theMatA.rows();
    col_num2 = (int)theMatB.cols();
    row_num2 = (int)theMatB.rows();
    Eigen::MatrixXd resMat(row_num1, col_num1 * col_num2);
    int index = 0;
    for (int rowLoopIndex = 0; rowLoopIndex < row_num1; ++rowLoopIndex)
    {
        index = 0;
        for (int i = 0; i < col_num1; ++i)
        {
            for (int j = 0; j < col_num2; ++j)
            {
                resMat(rowLoopIndex, index++) = theMatA(rowLoopIndex, i) * theMatB(rowLoopIndex, j);
            }
        }
    }
    return resMat;
}

int GuidedCoonsSurfGenerator::SetSameDistribution(sggk::BSplineCurve3DPtr& C1, sggk::BSplineCurve3DPtr& C2)
{
    const sggk::RealArray& K1 = C1->Knots();
    const sggk::RealArray& K2 = C2->Knots();
    double K11 = K1.front();
    double K12 = K1.back();
    double K21 = K2.front();
    double K22 = K2.back();

    // 统一两条曲线的参数域
    if ((K12 - K11) > (K22 - K21)) {
        C2->AdjustKnots(sggk::Interval(K11, K12));
    }
    else if ((K12 - K11) < (K22 - K21)) {
        C1->AdjustKnots(sggk::Interval(K21, K22));
    }
    else if (std::fabs(K12 - K11) > 1e-7) {
        C2->AdjustKnots(sggk::Interval(K11, K12));
    }

    // 重新获取节点与重数（AdjustKnots 后可能变化）
    sggk::RealArray K1c = C1->Knots();
    sggk::UIntArray M1c = C1->Mults();
    sggk::RealArray K2c = C2->Knots();
    sggk::UIntArray M2c = C2->Mults();

    // 互相插入节点，使两条曲线具有一致的节点向量
    // 注意：OCC 的 InsertKnots(k, M) 中 M 是目标重数，SGK 的 InsertKnots(t, times) 中 times 是插入次数，
    // 因此需换算为「目标重数 - 当前重数」后再插入
    for (size_t i = 0; i < K2c.size(); ++i) {
        unsigned int curMult = CurveKnotMultiplicity(C1, K2c[i]);
        if (M2c[i] > curMult)
            C1->InsertKnots(K2c[i], M2c[i] - curMult);
    }
    for (size_t i = 0; i < K1c.size(); ++i) {
        unsigned int curMult = CurveKnotMultiplicity(C2, K1c[i]);
        if (M1c[i] > curMult)
            C2->InsertKnots(K1c[i], M1c[i] - curMult);
    }

    return (int)C1->ControlPoints().size();
}