#include "GuidedCoonsSurfGenerator.h"

GuidedCoonsSurfGenerator::GuidedCoonsSurfGenerator(const std::vector<Handle(Geom_BSplineCurve)>& boundaryCurves, const std::vector<Handle(Geom_BSplineCurve)>& guideCurves, Standard_Real theTol)
	: m_boundaryCurves(boundaryCurves), m_guideCurves(guideCurves), m_tol(theTol), m_isDone(Standard_False), m_iterateCount(0)
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
        m_isDone = Standard_True;
        return;
    }

    TopoDS_Face coonsFace = BRepBuilderAPI_MakeFace(m_originalSurf, Precision::Confusion());
    std::string filePath = "/Users/flytiger25/work/occ/data/coons/";
    filePath += "coons.step";
    STEPControl_Writer stepWriter;
    stepWriter.Transfer(coonsFace, STEPControl_AsIs);
    stepWriter.Write(filePath.c_str());

    // 未达到容差要求，继续迭代
    while (!m_isDone && m_iterateCount < MAX_ITERATIONS)
    {
        std::cout << "m_iterateCount: " << m_iterateCount << std::endl;
        ConstructSurfWithGuideCrvs();
        m_originalSurf = m_guidedSurf;
        m_iterateCount++;

        TopoDS_Face coonsFace = BRepBuilderAPI_MakeFace(m_originalSurf, Precision::Confusion());
        std::string filePath = "/Users/flytiger25/work/occ/data/coons/";
        filePath += "GuidedSurf_" + std::to_string(m_iterateCount) + ".step";
        STEPControl_Writer stepWriter;
        stepWriter.Transfer(coonsFace, STEPControl_AsIs);
        stepWriter.Write(filePath.c_str());
    }
}

void GuidedCoonsSurfGenerator::ConstructCoonsSurf()
{
    // 重新拟合边界线
    ApproximateBoundaryCurves(m_boundaryCurves);
    // 构造Coons
    Handle(Geom_BSplineCurve) bsplineCurve1, bsplineCurve2, bsplineCurve3, bsplineCurve4;
    Arrange_Coons_G0(m_boundaryCurves, bsplineCurve1, bsplineCurve2, bsplineCurve3, bsplineCurve4);
    Coons_G0(bsplineCurve1, bsplineCurve2, bsplineCurve3, bsplineCurve4, m_originalSurf);
    m_coonsSurf = m_originalSurf;

    // 裁剪引导线
    TrimInternalCurves(m_guideCurves, m_boundaryCurves);
}

Standard_Integer GuidedCoonsSurfGenerator::Arrange_Coons_G0(std::vector<Handle(Geom_BSplineCurve)>& curveArray, Handle(Geom_BSplineCurve)& bslpineCurve1,
    Handle(Geom_BSplineCurve)& bslpineCurve2, Handle(Geom_BSplineCurve)& bslpineCurve3, Handle(Geom_BSplineCurve)& bslpineCurve4, Standard_Real tol, Standard_Integer isModify)
{
    std::vector<Handle(Geom_BSplineCurve)> curveArraybak(curveArray);

    if (curveArraybak.size() != 4)
        return 0;

    bslpineCurve1 = curveArraybak[0];
    bslpineCurve2 = curveArraybak[1];
    bslpineCurve3 = curveArraybak[2];
    bslpineCurve4 = curveArraybak[3];

    Standard_Real maxDis = -1;
    Standard_Real dis;
    gp_Pnt curve1startpoint = bslpineCurve1->StartPoint();
    gp_Pnt curve1endpoint = bslpineCurve1->EndPoint();

    curveArraybak.erase(curveArraybak.begin());
    gp_Pnt curve2startpoint = bslpineCurve2->StartPoint();
    gp_Pnt curve2endpoint = bslpineCurve2->EndPoint();

    gp_Pnt curve3startpoint = bslpineCurve3->StartPoint();
    gp_Pnt curve3endpoint = bslpineCurve3->EndPoint();

    gp_Pnt curve4startpoint = bslpineCurve4->StartPoint();
    gp_Pnt curve4endpoint = bslpineCurve4->EndPoint();

    gp_Pnt point;
    //fine curve2
    for (Standard_Integer i = 0; i < curveArraybak.size(); ++i)
    {
        bslpineCurve2 = curveArraybak[i];

        curve2startpoint = bslpineCurve2->StartPoint();
        curve2endpoint = bslpineCurve2->EndPoint();

        if (curve1endpoint.IsEqual(curve2startpoint, tol))
        {
            dis = curve1endpoint.Distance(curve2startpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve1endpoint.BaryCenter(0.5, curve2startpoint, 0.5);
                bslpineCurve2->SetPole(1, curve1endpoint);
                bslpineCurve1->SetPole(bslpineCurve1->NbPoles(), curve1endpoint);
            }

            curveArraybak.erase(curveArraybak.begin() + i);
            break;
        }
        if (curve1endpoint.IsEqual(curve2endpoint, tol))
        {
            dis = curve1endpoint.Distance(curve2endpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve1endpoint.BaryCenter(0.5, curve2endpoint, 0.5);
                bslpineCurve2->SetPole(bslpineCurve2->NbPoles(), curve1endpoint);
                bslpineCurve1->SetPole(bslpineCurve1->NbPoles(), curve1endpoint);
            }

            curveArraybak.erase(curveArraybak.begin() + i);
            bslpineCurve2->Reverse();

            curve2startpoint = bslpineCurve2->StartPoint();
            curve2endpoint = bslpineCurve2->EndPoint();
            break;
        }
    }

    //fine curve3
    for (Standard_Integer i = 0; i < curveArraybak.size(); ++i)
    {
        bslpineCurve3 = curveArraybak[i];

        curve3startpoint = bslpineCurve3->StartPoint();
        curve3endpoint = bslpineCurve3->EndPoint();

        if (curve2endpoint.IsEqual(curve3endpoint, tol))
        {
            dis = curve2endpoint.Distance(curve3endpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve2endpoint.BaryCenter(0.5, curve3endpoint, 0.5);
                bslpineCurve2->SetPole(bslpineCurve2->NbPoles(), curve2endpoint);
                bslpineCurve3->SetPole(bslpineCurve3->NbPoles(), curve2endpoint);
            }

            curveArraybak.erase(curveArraybak.begin() + i);
            break;
        }
        if (curve2endpoint.IsEqual(curve3startpoint, tol))
        {
            dis = curve2endpoint.Distance(curve3startpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve2endpoint.BaryCenter(0.5, curve3startpoint, 0.5);
                bslpineCurve2->SetPole(bslpineCurve2->NbPoles(), curve2endpoint);
                bslpineCurve3->SetPole(1, curve2endpoint);
            }

            curveArraybak.erase(curveArraybak.begin() + i);
            bslpineCurve3->Reverse();

            curve3startpoint = bslpineCurve3->StartPoint();
            curve3endpoint = bslpineCurve3->EndPoint();
            break;
        }
    }

    if (curveArraybak.size() != 1)
        return 0;

    bslpineCurve4 = curveArraybak[0];
    curve4startpoint = bslpineCurve4->StartPoint();
    curve4endpoint = bslpineCurve4->EndPoint();
    if (curve3startpoint.IsEqual(curve4endpoint, tol))
    {
        dis = curve3startpoint.Distance(curve4endpoint);
        if (dis > maxDis)
            maxDis = dis;

        if (isModify && IsLess(dis, tol))
        {
            curve3startpoint.BaryCenter(0.5, curve4endpoint, 0.5);
            bslpineCurve3->SetPole(1, curve3startpoint);
            bslpineCurve4->SetPole(bslpineCurve4->NbPoles(), curve3startpoint);
        }
    }
    else
        if (curve3startpoint.IsEqual(curve4startpoint, tol))
        {
            dis = curve3startpoint.Distance(curve4startpoint);
            if (dis > maxDis)
                maxDis = dis;

            if (isModify && IsLess(dis, tol))
            {
                curve3startpoint.BaryCenter(0.5, curve4startpoint, 0.5);
                bslpineCurve3->SetPole(1, curve3startpoint);
                bslpineCurve4->SetPole(1, curve3startpoint);
            }

            bslpineCurve4->Reverse();
            curve4startpoint = bslpineCurve4->StartPoint();
            curve4endpoint = bslpineCurve4->EndPoint();
        }
        else
            return 0;
    return 1;
}

void GuidedCoonsSurfGenerator::Coons_G0(Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4, Handle(Geom_BSplineSurface)& mySurface_coons)
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
    Standard_Integer NbUPoles = curve1->NbPoles();
    Standard_Integer NbVPoles = curve2->NbPoles();

    //control points of the coons surface
    TColgp_Array2OfPnt Poles_result(1, NbUPoles, 1, NbVPoles);

    //control points of ruled surface in the v direction
    TColgp_Array2OfPnt Poles_vruled(1, NbUPoles, 1, 2);

    //control points of ruled surface in the u direction
    TColgp_Array2OfPnt Poles_uruled(1, 2, 1, NbVPoles);

    //control points of ruled surface of the four corner points
    TColgp_Array2OfPnt Poles_uruled_vruled(1, 2, 1, 2);

    // Get the u knot vector
    Standard_Integer NbUKnot = curve1->NbKnots();
    TColStd_Array1OfReal    UKnots(1, NbUKnot);
    TColStd_Array1OfInteger UMults(1, NbUKnot);
    curve1->Knots(UKnots);
    curve1->Multiplicities(UMults);

    // Get the v knot vector
    Standard_Integer NbVKnot = curve2->NbKnots();
    TColStd_Array1OfReal    VKnots(1, NbVKnot);
    TColStd_Array1OfInteger VMults(1, NbVKnot);
    curve2->Knots(VKnots);
    curve2->Multiplicities(VMults);

    // Set v knots of the v ruled surface
    TColStd_Array1OfReal    VKnots_RuledInVdirection(1, 2);
    TColStd_Array1OfInteger VMults_RuledInVdirection(1, 2);
    VKnots_RuledInVdirection(1) = curve2->FirstParameter();
    VKnots_RuledInVdirection(2) = curve2->LastParameter();

    VMults_RuledInVdirection(1) = 2;
    VMults_RuledInVdirection(2) = 2;

    // Set u knots of the u ruled surface
    TColStd_Array1OfReal    UKnots_RuledInUdirection(1, 2);
    TColStd_Array1OfInteger UMults_RuledInUdirection(1, 2);
    UKnots_RuledInUdirection(1) = curve1->FirstParameter();
    UKnots_RuledInUdirection(2) = curve1->LastParameter();

    UMults_RuledInUdirection(1) = 2;
    UMults_RuledInUdirection(2) = 2;

    // 
    //2.1 generate the ruled surface in the v direction
    // 
    // 
    // Set control points of the v ruled surface
    for (Standard_Integer i = 1; i <= NbUPoles; ++i)
    {
        Poles_vruled(i, 1) = curve1->Pole(i);
        Poles_vruled(i, 2) = curve3->Pole(i);
    }

    // Generate the bspline surface
    Handle(Geom_BSplineSurface) mySurface_VRuled = new Geom_BSplineSurface(Poles_vruled,
        UKnots, VKnots_RuledInVdirection,
        UMults, VMults_RuledInVdirection,
        curve1->Degree(), 1);

    // Make the surface compatible with the result coons surface
    mySurface_VRuled->IncreaseDegree(curve1->Degree(), curve2->Degree());

    for (Standard_Integer i = 1; i <= NbVKnot; ++i)
        mySurface_VRuled->InsertVKnot(VKnots(i), VMults(i), Precision::Confusion());

    // 
    //2.2 generate the ruled surface in the u direction
    // 

    // Set control points of the u ruled surface
    for (Standard_Integer i = 1; i <= NbVPoles; ++i)
    {
        Poles_uruled(1, i) = curve4->Pole(i);
        Poles_uruled(2, i) = curve2->Pole(i);
    }

    Handle(Geom_BSplineSurface) mySurface_URuled = new Geom_BSplineSurface(Poles_uruled,
        UKnots_RuledInUdirection, VKnots,
        UMults_RuledInUdirection, VMults,
        1, curve2->Degree());

    //increase degree
    mySurface_URuled->IncreaseDegree(curve1->Degree(), curve2->Degree());

    for (Standard_Integer i = 1; i <= NbUKnot; ++i)
        mySurface_URuled->InsertUKnot(UKnots(i), UMults(i), Precision::Confusion());

    // 
    //2.3 generate the ruled surface of the four corner points
    // 
    // 
    Poles_uruled_vruled(1, 1) = curve1->Pole(1);
    Poles_uruled_vruled(1, 2) = curve3->Pole(1);
    Poles_uruled_vruled(2, 1) = curve1->Pole(NbUPoles);
    Poles_uruled_vruled(2, 2) = curve3->Pole(NbUPoles);

    Handle(Geom_BSplineSurface) mySurface_RuledSurfaceof4CornerPoints = new Geom_BSplineSurface(Poles_uruled_vruled,
        UKnots_RuledInUdirection, VKnots_RuledInVdirection,
        UMults_RuledInUdirection, VMults_RuledInVdirection,
        1, 1);

    mySurface_RuledSurfaceof4CornerPoints->IncreaseDegree(curve1->Degree(), curve2->Degree());

    for (Standard_Integer i = 1; i <= NbVKnot; ++i)
        mySurface_RuledSurfaceof4CornerPoints->InsertVKnot(VKnots(i), VMults(i), Precision::Confusion());

    for (Standard_Integer i = 1; i <= NbUKnot; ++i)
        mySurface_RuledSurfaceof4CornerPoints->InsertUKnot(UKnots(i), UMults(i), Precision::Confusion());

    //2.4 Generate the sum surface mySurface_VRuled + mySurface_URuled - mySurface_RuledSurfaceof4CornerPoints 
    //increase degree

    for (Standard_Integer i = 1; i <= NbUPoles; ++i)
    {
        for (Standard_Integer j = 1; j <= NbVPoles; ++j)
        {
            Poles_result(i, j).SetXYZ(mySurface_VRuled->Pole(i, j).XYZ() + mySurface_URuled->Pole(i, j).XYZ() - mySurface_RuledSurfaceof4CornerPoints->Pole(i, j).XYZ());
        }
    }

    mySurface_coons = new Geom_BSplineSurface(Poles_result,
        UKnots, VKnots,
        UMults, VMults,
        curve1->Degree(), curve2->Degree());
}

void GuidedCoonsSurfGenerator::TrimInternalCurves(
    std::vector<Handle(Geom_BSplineCurve)>& theInternalBSplineCurves,
    const std::vector<Handle(Geom_BSplineCurve)>& theBoundaryCurveArray,
    Standard_Real theToleranceDistance) // theToleranceDistance 用于判断交点是否足够接近零距离
{
    //std::vector<Handle(Geom_BSplineCurve)> aTrimmedCurvesResult; // 存储裁剪后的结果

    for (Standard_Integer i = 0; i < theInternalBSplineCurves.size(); i++)
    {
        Handle(Geom_BSplineCurve) aCurrentInternalCurve = theInternalBSplineCurves[i];
        std::vector<std::pair<Standard_Real, Standard_Real>> curTrimInterval; // 当前引导线线的裁剪区间

        // 存储所有交点在内部曲线上的参数
        std::vector<Standard_Real> intersectionParams;
        intersectionParams.push_back(aCurrentInternalCurve->FirstParameter()); // 曲线起点
        intersectionParams.push_back(aCurrentInternalCurve->LastParameter());  // 曲线终点

        // 遍历所有边界曲线，查找零距离的极值点（即交点）
        for (auto& aBoundaryCurve : theBoundaryCurveArray)
        {
            GeomAPI_ExtremaCurveCurve anExtrema(aCurrentInternalCurve, aBoundaryCurve);

            // 遍历所有的极值点
            for (Standard_Integer k = 1; k <= anExtrema.NbExtrema(); ++k)
            {
                // 获取第 k 个极值点对之间的距离**
                Standard_Real dist = anExtrema.Distance(k);

                // 如果距离足够小（视为交点），则获取其参数
                if (dist < theToleranceDistance)
                {
                    Standard_Real paramOnInternalCurve, paramOnBoundaryCurve;
                    anExtrema.Parameters(k, paramOnInternalCurve, paramOnBoundaryCurve);
                    intersectionParams.push_back(paramOnInternalCurve);
                }
            }
        }

        // 对所有参数进行排序，并去重（避免重复的参数点）
        std::sort(intersectionParams.begin(), intersectionParams.end());
        intersectionParams.erase(std::unique(intersectionParams.begin(), intersectionParams.end(),
            [](Standard_Real a, Standard_Real b) { return Abs(a - b) < 1e-2; }), // 使用一个小的容差进行去重
            intersectionParams.end());

        // 如果只有起点、终点 + 两个内部交点
        if (intersectionParams.size() == 4)
        {
            Standard_Real startParam = intersectionParams[1];
            Standard_Real endParam = intersectionParams[2];

            // 排除两个点靠很近的情况
            gp_Pnt p0, p1;
            aCurrentInternalCurve->D0(startParam, p0);
            aCurrentInternalCurve->D0(endParam, p1);
            if (!p0.IsEqual(p1, 10))
            {
                // 直接添加，不做包围盒判断
                curTrimInterval.emplace_back(startParam, endParam);
            }
        }
        else
        {
            // 根据交点参数划分曲线段并判断内部性
            for (Standard_Integer j = 0; j < intersectionParams.size() - 1; ++j)
            {
                Standard_Real startParam = intersectionParams[j];
                Standard_Real endParam = intersectionParams[j + 1];

                // 排除两个点靠很近的情况
                gp_Pnt p0, p1;
                aCurrentInternalCurve->D0(startParam, p0);
                aCurrentInternalCurve->D0(endParam, p1);
                if (p0.IsEqual(p1, 10)) continue;

                // 裁剪出当前区间段
                Handle(Geom_TrimmedCurve) aSegment = new Geom_TrimmedCurve(aCurrentInternalCurve, startParam, endParam);
                Handle(Geom_BSplineCurve) aBSplineSegment = GeomConvert::CurveToBSplineCurve(aSegment, Convert_TgtThetaOver2);

                //std::vector<Handle(Geom_BSplineCurve)> aTempBoundaryCurveArray = theBoundaryCurveArray;
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
std::vector<gp_Pnt> GuidedCoonsSurfGenerator::DiscretizeBSplineCurve(
    const Handle(Geom_BSplineCurve)& theCurve,
    Standard_Integer numSegments,
    Standard_Boolean theBoundaryFlag)
{
    std::vector<gp_Pnt> discretizedPoints;
    if (numSegments < 1) return {};

    Standard_Real firstParam = theCurve->FirstParameter();
    Standard_Real lastParam = theCurve->LastParameter();
    Standard_Real step = (lastParam - firstParam) / numSegments;

    for (Standard_Integer i = 0; i <= numSegments; ++i)
    {
        Standard_Real currentParam = firstParam + i * step;
        if (i == numSegments) currentParam = lastParam;
        discretizedPoints.push_back(theCurve->Value(currentParam));
    }

    if (!theBoundaryFlag)
    {
        discretizedPoints.erase(discretizedPoints.end() - 1);
        discretizedPoints.erase(discretizedPoints.begin());
    }
    return discretizedPoints;
}

Standard_Boolean GuidedCoonsSurfGenerator::CurvesConnectedLoop(
    std::vector<Handle(Geom_BSplineCurve)>& theCurves,
    Standard_Real theTolerance)
{
    if (theCurves.size() < 3)
        return Standard_False;

    std::vector<bool> used(theCurves.size(), false);
    std::vector<Handle(Geom_BSplineCurve)> orderedCurves;
    orderedCurves.reserve(theCurves.size());

    // 从第一条曲线开始
    Handle(Geom_BSplineCurve) current = theCurves[0];
    orderedCurves.push_back(current);
    used[0] = true;

    for (size_t i = 1; i < theCurves.size(); ++i)
    {
        gp_Pnt currEnd = current->Value(current->LastParameter());
        Standard_Boolean foundNext = Standard_False;

        for (size_t j = 0; j < theCurves.size(); ++j)
        {
            if (used[j]) continue;

            Handle(Geom_BSplineCurve) candidate = theCurves[j];
            gp_Pnt candStart = candidate->Value(candidate->FirstParameter());
            gp_Pnt candEnd = candidate->Value(candidate->LastParameter());

            // 情况1：起点匹配，不用反转
            if (currEnd.Distance(candStart) < theTolerance)
            {
                orderedCurves.push_back(candidate);
                used[j] = true;
                current = candidate;
                foundNext = Standard_True;
                break;
            }
            // 情况2：终点匹配，需要反转
            else if (currEnd.Distance(candEnd) < theTolerance)
            {
                candidate->Reverse();
                orderedCurves.push_back(candidate);
                used[j] = true;
                current = candidate;
                foundNext = Standard_True;
                break;
            }
        }

        if (!foundNext)
        {
            // 找不到匹配项
            return Standard_False;
        }
    }

    // 最终闭合性检测（最后一条与第一条）
    gp_Pnt lastEnd = orderedCurves.back()->Value(orderedCurves.back()->LastParameter());
    gp_Pnt firstStart = orderedCurves.front()->Value(orderedCurves.front()->FirstParameter());

    if (lastEnd.Distance(firstStart) > theTolerance)
    {
        return Standard_False; // 无法闭环
    }

    theCurves = orderedCurves; // 更新为首尾相接排序结果
    return Standard_True;
}


// 判断一个点是否在由一系列二维点构成的多边形内部（射线法）
Standard_Boolean GuidedCoonsSurfGenerator::IsPointInPolygon2D(
    const gp_Pnt2d& theTestPoint,
    const std::vector<gp_Pnt2d>& thePolygon2d,
    Standard_Real theTolerance)
{
    if (thePolygon2d.empty()) return Standard_False;

    Standard_Integer intersectCount = 0;
    Standard_Integer n = (Standard_Integer)thePolygon2d.size();

    for (Standard_Integer i = 0; i < n; ++i)
    {
        const gp_Pnt2d& p1 = thePolygon2d[i];
        const gp_Pnt2d& p2 = thePolygon2d[(i + 1) % n]; // 形成闭合循环

        // 检查测试点是否在边界上（在容差范围内）
        // 可以根据需要增加更精确的点到线段距离判断
        if (Abs((p2.Y() - p1.Y()) * theTestPoint.X() - (p2.X() - p1.X()) * theTestPoint.Y() + p2.X() * p1.Y() - p2.Y() * p1.X()) < theTolerance &&
            ((theTestPoint.X() >= std::min(p1.X(), p2.X()) - theTolerance && theTestPoint.X() <= std::max(p1.X(), p2.X()) + theTolerance) ||
                (theTestPoint.Y() >= std::min(p1.Y(), p2.Y()) - theTolerance && theTestPoint.Y() <= std::max(p1.Y(), p2.Y()) + theTolerance)))
        {
            return Standard_True; // 点在边界线上，视为在内部
        }

        // 射线法：从 theTestPoint 向正X方向发射射线
        if (((p1.Y() > theTestPoint.Y() && p2.Y() <= theTestPoint.Y()) ||
            (p2.Y() > theTestPoint.Y() && p1.Y() <= theTestPoint.Y())))
        {
            Standard_Real xIntersect = p1.X() + (theTestPoint.Y() - p1.Y()) / (p2.Y() - p1.Y()) * (p2.X() - p1.X());
            if (xIntersect > theTestPoint.X() - theTolerance) // 允许微小偏差
            {
                intersectCount++;
            }
        }
    }
    return (intersectCount % 2 == 1);
}

Standard_Boolean GuidedCoonsSurfGenerator::IsCurveInsideBoundaries(
    const Handle(Geom_BSplineCurve)& theCurve,
    std::vector<Handle(Geom_BSplineCurve)>& theBoundaryCurveArray,
    Standard_Real theToleranceDistance)
{
    if (theBoundaryCurveArray.size() < 3)
    {
        return Standard_False;
    }

    CurvesConnectedLoop(theBoundaryCurveArray); // 假设这个函数确保了顺序和连接性

    // 获取边界的顶点
    std::vector<gp_Pnt> boundaryVertices3D;
    for (const auto& boundaryCurve : theBoundaryCurveArray)
    {
        boundaryVertices3D.push_back(boundaryCurve->Value(boundaryCurve->FirstParameter()));
    }

    // 假设 boundaryVertices3D 至少有3个点，且它们定义了一个平面。
    // 在实际应用中，您可能需要更健壮的逻辑来处理共线或平面退化的情况。
    if (boundaryVertices3D.size() < 3)
    {
        return Standard_False; // 无法定义平面
    }

    gp_Pnt p1 = boundaryVertices3D[0];
    gp_Pnt p2 = boundaryVertices3D[1];
    gp_Pnt p3 = boundaryVertices3D[2];

    // 计算平面的法向量
    gp_Vec v12(p1, p2);
    gp_Vec v13(p1, p3);
    gp_Vec planeNormal = v12.Crossed(v13);

    if (planeNormal.Magnitude() < gp::Resolution())
    {
        return Standard_False;
    }
    planeNormal.Normalize(); // 归一化法向量

    // 建立平面局部二维坐标系基向量
    gp_Vec uBasis = v12; // 沿 P1P2 方向
    if (uBasis.Magnitude() < gp::Resolution())
    {
        // P1P2 距离太近，尝试用 P1P3
        uBasis = v13;
    }
    uBasis.Normalize();

    gp_Vec vBasis = planeNormal.Crossed(uBasis); // 确保与 uBasis 和法向量正交

    std::vector<gp_Pnt2d> polygon2D;
    for (const auto& p3d : boundaryVertices3D)
    {
        gp_Vec vecFromP1(p1, p3d);
        Standard_Real u = vecFromP1.Dot(uBasis);
        Standard_Real v = vecFromP1.Dot(vBasis);
        polygon2D.emplace_back(u, v);
    }

    Standard_Integer numSamplePoints = 5;
    std::vector<gp_Pnt> internalCurveSamplePoints = DiscretizeBSplineCurve(theCurve, numSamplePoints - 1, Standard_False);

    for (const gp_Pnt& samplePoint3D : internalCurveSamplePoints)
    {
        // 将采样点投影到自定义平面
        gp_Vec vecSP1(p1, samplePoint3D);
        Standard_Real distToPlane = vecSP1.Dot(planeNormal);
        gp_Pnt projectedPoint3D = samplePoint3D.Translated(gp_Vec(planeNormal).Multiplied(-distToPlane)); // 沿着法线方向移动到平面上

        // 将投影点转换到二维局部坐标系
        gp_Vec vecProjP1(p1, projectedPoint3D);
        Standard_Real uProj = vecProjP1.Dot(uBasis);
        Standard_Real vProj = vecProjP1.Dot(vBasis);
        gp_Pnt2d samplePoint2D(uProj, vProj);

        // 判断投影点是否在二维多边形内部
        if (!IsPointInPolygon2D(samplePoint2D, polygon2D, theToleranceDistance))
        {
            // 只要有一个采样点不在内部，则整条曲线不在内部
            return Standard_False;
        }
    }

    return Standard_True; // 所有采样点都在自定义平面内的多边形区域内
}

Standard_Boolean GuidedCoonsSurfGenerator::IsCurveInsideSurface(
    const Handle(Geom_BSplineCurve)& theCurve,
    const Handle(Geom_Surface)& theSurface,
    const Standard_Real theTolerance)
{
    if (theCurve.IsNull() || theSurface.IsNull())
        return Standard_False;

    // 构造 Coons 曲面的 AABB
    Bnd_Box surfaceBox;
    BRepBndLib::Add(BRepLib_MakeFace(theSurface, 1e-6), surfaceBox);
    surfaceBox.Enlarge(theTolerance);  // 加一点容差

    // 离散采样
    std::vector<gp_Pnt> samplePoints = DiscretizeBSplineCurve(theCurve, 50, Standard_False);
    Standard_Integer totalCount = static_cast<Standard_Integer>(samplePoints.size());
    Standard_Integer insideCount = 0;

    for (const auto& p : samplePoints)
    {
        if (!surfaceBox.IsOut(p))
            ++insideCount;
    }

    Standard_Real ratio = static_cast<Standard_Real>(insideCount) / totalCount;
    return (ratio >= 0.7);
}

Standard_Real GuidedCoonsSurfGenerator::ComputeCurveCurveDistance(const Handle(Geom_BSplineCurve)& guideCurve, const Handle(Geom_BSplineCurve)& boundaryCurve)
{
    GeomAPI_ExtremaCurveCurve aExtrema(guideCurve, boundaryCurve);
    // 检查是否找到了极值点
    if (aExtrema.NbExtrema() > 0)
    {
        // 遍历所有极值点，找到最小距离
        Standard_Real aMinDistance = RealLast();
        for (Standard_Integer i = 1; i <= aExtrema.NbExtrema(); ++i)
        {
            Standard_Real aDist = aExtrema.Distance(i);
            if (aDist < aMinDistance)
            {
                aMinDistance = aDist;
            }
        }
        return aMinDistance;
    }

    return INT_MAX;
}

void GuidedCoonsSurfGenerator::ApproximateBoundaryCurves(std::vector<Handle(Geom_BSplineCurve)>& curves, Standard_Integer samplingNum)
{
    for (auto& curve : curves)
    {
        TColStd_Array1OfReal curveKnots(1, curve->NbKnots());
        curve->Knots(curveKnots);

        // 重新参数化曲线的节点
        if (!(curveKnots(curveKnots.Lower()) == 0 && curveKnots(curveKnots.Upper()) == 1))
        {
            BSplCLib::Reparametrize(0, 1, curveKnots);
            curve->SetKnots(curveKnots);
        }

        
        std::vector<gp_Pnt> samplingPnts;
        std::vector<Standard_Real> samplingParams;
        //TColgp_Array1OfPnt points(1, samplingNum);
        //TColStd_Array1OfReal params(1, samplingNum);
        GeomAdaptor_Curve adaptor(curve);
        // 改为弧长采样
        GCPnts_UniformAbscissa  uniformAbscissa(adaptor, samplingNum);
        if (uniformAbscissa.IsDone())
        {
            for (int i = 1; i <= samplingNum; ++i) {
                Standard_Real param = uniformAbscissa.Parameter(i);
                gp_Pnt pnt = curve->Value(param);
                samplingParams.push_back((i - 1.0) / (samplingNum - 1.0));
                samplingPnts.push_back(pnt);
            }
        }
        else
        {
            Standard_Real vMin = curve->FirstParameter();
            Standard_Real vMax = curve->LastParameter();
            for (Standard_Integer j = 1; j <= samplingNum; ++j)
            {
                Standard_Real param = vMin + (vMax - vMin) * (j - 1) / (samplingNum - 1);
                gp_Pnt pnt = curve->Value(param);
                samplingParams.push_back((j - 1.0) / (samplingNum - 1.0));
                samplingPnts.push_back(pnt);
            }
        }

        //GeomAPI_PointsToBSpline approx(points);
        //curve = approx.Curve();
        
        /*
        // 采样点与参数生成
        TColgp_Array1OfPnt points(1, samplingNum);
        TColStd_Array1OfReal params(1, samplingNum);
        std::vector<gp_Pnt> samplingPnts;
        std::vector<Standard_Real> samplingParams;
        Standard_Real vMin = curve->FirstParameter();
        Standard_Real vMax = curve->LastParameter();

        for (Standard_Integer j = 1; j <= samplingNum; ++j)
        {
            Standard_Real param = vMin + (vMax - vMin) * (j - 1) / (samplingNum - 1);
            gp_Pnt pnt = curve->Value(param);
            samplingParams.push_back((j - 1.0) / (samplingNum - 1.0));
            samplingPnts.push_back(pnt);
        }
        */
        // 初始化节点并进行拟合
        std::vector<Standard_Real> init_knots = KnotGernerationByParams(samplingParams, APPROXIMATE_KNOTS_NUM, 3);
        std::vector<Standard_Real> insertKnots;
        curve = IterateApproximate(insertKnots, samplingPnts, samplingParams, init_knots, 3);
    }
}

Handle(Geom_BSplineCurve) GuidedCoonsSurfGenerator::IterateApproximate(std::vector<Standard_Real>& insertKnots, const std::vector<gp_Pnt>& pnts, std::vector<Standard_Real>& pntsParams,
    std::vector<Standard_Real>& initKnots, Standard_Integer degree, Standard_Integer maxIterNum, Standard_Real toler)
{
    Standard_Integer itNum = 1;
    Standard_Real currentMaxError = 100;
    Handle(Geom_BSplineCurve) IterBspineCurve;
    std::vector<Standard_Real> CurrentKnots = initKnots;

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

Handle(Geom_BSplineCurve) GuidedCoonsSurfGenerator::ApproximateCurve(const std::vector<gp_Pnt>& pnts, std::vector<Standard_Real>& params, std::vector<Standard_Real>& fKnots, Standard_Integer degree)
{
    std::vector<Standard_Real> Knots;
    std::vector<Standard_Integer> Mutis;
    SequenceToKnots(fKnots, Knots, Mutis);
    TColStd_Array1OfReal Knots_OCC(1, (Standard_Integer)Knots.size());
    TColStd_Array1OfInteger Mutis_OCC(1, (Standard_Integer)Mutis.size());
    for (Standard_Integer i = 0; i < (Standard_Integer)Knots.size(); ++i) {
        Knots_OCC[i + 1] = Knots[i];
        Mutis_OCC[i + 1] = Mutis[i];
    }
    return ApproximateCurve(pnts, params, Knots_OCC, Mutis_OCC, fKnots, degree);
}

Handle(Geom_BSplineCurve) GuidedCoonsSurfGenerator::ApproximateCurve(const std::vector<gp_Pnt>& pnts, std::vector<Standard_Real>& pntsParams, TColStd_Array1OfReal& knots,
    TColStd_Array1OfInteger& mutis, std::vector<Standard_Real>& fKnots, Standard_Integer degree)
{
    Standard_Integer n = (Standard_Integer)fKnots.size() - degree - 2;
    Standard_Integer m = (Standard_Integer)pnts.size() - 1;

    //1.Construct matrix N
    Eigen::MatrixXd matN(m - 1, n - 1);
    for (Standard_Integer i = 0; i < m - 1; ++i) {
        for (Standard_Integer j = 0; j < n - 1; ++j) {
            Standard_Real value = CalBasicFunction(pntsParams[i + 1], j + 1, degree, fKnots);
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
    for (Standard_Integer i = 1; i <= m - 1; ++i) {
        gp_Vec VecTemp = CalResPnt(i, pnts, pntsParams, degree, fKnots, n);
        Standard_Real x, y, z;
        VecTemp.Coord(x, y, z);
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
    TColgp_Array1OfPnt ctrlPnts(1, n + 1);
    ctrlPnts[1] = pnts[0];
    ctrlPnts[n + 1] = pnts[m];
    for (Standard_Integer i = 2; i <= n; i++) {
        gp_Pnt pntTemp(Sx(i - 2), Sy(i - 2), Sz(i - 2));
        ctrlPnts[i] = pntTemp;
    }
    Handle(Geom_BSplineCurve) bspline = new Geom_BSplineCurve(ctrlPnts, knots, mutis, degree);
    return bspline;
}

gp_Vec GuidedCoonsSurfGenerator::CalResPnt(Standard_Integer k, const std::vector<gp_Pnt>& dataPoints, const std::vector<Standard_Real>& parameters, Standard_Integer p, std::vector<Standard_Real>& knots, Standard_Integer ctrlPntNum)
{
    Standard_Real aCoeff1 = CalBasicFunction(parameters[k], 0, p, knots);
    Standard_Real aCoeff2 = CalBasicFunction(parameters[k], ctrlPntNum, p, knots);
    gp_Vec vecTemp0(dataPoints[0].Coord());
    gp_Vec vecTempm(dataPoints[dataPoints.size() - 1].Coord());
    gp_Vec vecTempk(dataPoints[k].Coord());
    gp_Vec vectemp = vecTempk - aCoeff1 * vecTemp0 - aCoeff2 * vecTempm;
    return vectemp;
}

void GuidedCoonsSurfGenerator::SequenceToKnots(const std::vector<Standard_Real>& sequence, std::vector<Standard_Real>& knots, std::vector<Standard_Integer>& multiplicities)
{
    if (sequence.empty()) return;

    std::map<Standard_Real, Standard_Integer> knotMap;

    // 使用map来统计每个节点的重复次数
    for (Standard_Real value : sequence) {
        Standard_Boolean found = false;
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

std::vector<Standard_Real> GuidedCoonsSurfGenerator::KnotGernerationByParams(const std::vector<Standard_Real>& params, Standard_Integer n, Standard_Integer p)
{
    Standard_Integer m = (Standard_Integer)params.size() - 1;
    Standard_Real d = (m + 1) / (n - p + 1);
    std::vector<Standard_Real> Knots(n + p + 2);
    Standard_Integer temp;
    Standard_Real alpha;
    for (Standard_Integer i = 0; i <= p; ++i)
    {
        Knots[i] = 0.0;
    }
    for (Standard_Integer j = 1; j <= n - p; ++j)
    {
        temp = Standard_Integer(j * d);
        alpha = j * d - temp;
        Knots[p + j] = (1 - alpha) * params[temp - 1] + alpha * params[temp];
    }
    for (Standard_Integer i = n + 1; i <= n + p + 1; ++i)
    {
        Knots[i] = 1;
    }
    return Knots;
}

void GuidedCoonsSurfGenerator::AddDegenerateCurve(std::vector<Handle(Geom_BSplineCurve)>& boundaryCurves)
{
    if (boundaryCurves.size() == 3)
    {
        // 初始化一个向量用于存储每条曲线的交点计数以及对应的样条曲线
        std::vector<std::pair<Standard_Integer, Handle(Geom_BSplineCurve)>> anInterCount =
        {
            {0, boundaryCurves[0]},
            {0, boundaryCurves[1]},
            {0, boundaryCurves[2]}
        };

        // 按交点计数从大到小对曲线进行排序
        std::sort(anInterCount.begin(), anInterCount.end(),
            [](const std::pair<Standard_Integer, Handle(Geom_BSplineCurve)>& curve1,
                const std::pair<Standard_Integer, Handle(Geom_BSplineCurve)>& curve2)
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
            gp_Pnt aDegeneratePoint(0, 0, 0);
            if (boundaryCurves[0]->StartPoint().Distance(boundaryCurves[2]->StartPoint()) > 10
                && boundaryCurves[0]->StartPoint().Distance(boundaryCurves[2]->EndPoint()) > 10)
            {
                aDegeneratePoint = boundaryCurves[0]->EndPoint();
            }
            else if (boundaryCurves[0]->EndPoint().Distance(boundaryCurves[2]->StartPoint()) > 10
                && boundaryCurves[0]->EndPoint().Distance(boundaryCurves[2]->EndPoint()) > 10)
            {
                aDegeneratePoint = boundaryCurves[0]->StartPoint();
            }

            // 构建退化边
            TColgp_Array1OfPnt poles(1, 2);
            poles.SetValue(1, aDegeneratePoint);
            poles.SetValue(2, aDegeneratePoint);

            TColStd_Array1OfReal knots(1, 2);
            knots.SetValue(1, 0.0);
            knots.SetValue(2, 1.0);

            TColStd_Array1OfInteger multiplicities(1, 2);
            multiplicities.SetValue(1, 2);
            multiplicities.SetValue(2, 2);

            boundaryCurves[3] = new Geom_BSplineCurve(poles, knots, multiplicities, 1);

            // 调整次序，要求首尾相接
            Standard_Real tol = boundaryCurves[0]->EndPoint().Distance(boundaryCurves[0]->StartPoint()) / 1000;
            if (boundaryCurves[0]->StartPoint().Distance(boundaryCurves[1]->StartPoint()) < tol)
            {
                boundaryCurves[0]->Reverse();
            }
            if (boundaryCurves[2]->EndPoint().Distance(boundaryCurves[1]->EndPoint()) < tol)
            {
                boundaryCurves[2]->Reverse();
            }
        }
        else
        {
            gp_Pnt pnt1 = boundaryCurves[0]->StartPoint(), pnt2 = boundaryCurves[0]->EndPoint(), pnt3 = boundaryCurves[1]->StartPoint();
            gp_Pnt pnt4 = boundaryCurves[1]->EndPoint(), pnt5 = boundaryCurves[2]->StartPoint(), pnt6 = boundaryCurves[2]->EndPoint();

            Standard_Real tol = boundaryCurves[0]->EndPoint().Distance(boundaryCurves[0]->StartPoint()) / 1000;
            if (boundaryCurves[1]->EndPoint().Distance(boundaryCurves[0]->EndPoint()) < tol)
            {
                boundaryCurves[1]->Reverse();
            }
            else if (boundaryCurves[2]->StartPoint().Distance(boundaryCurves[0]->EndPoint()) < tol)
            {
                std::swap(boundaryCurves[1], boundaryCurves[2]);
            }
            else if (boundaryCurves[2]->EndPoint().Distance(boundaryCurves[0]->EndPoint()) < tol)
            {
                std::swap(boundaryCurves[1], boundaryCurves[2]);
                boundaryCurves[1]->Reverse();
            }

            if (boundaryCurves[2]->EndPoint().Distance(boundaryCurves[1]->EndPoint()) < tol)
            {
                boundaryCurves[2]->Reverse();
            }

            pnt1 = boundaryCurves[0]->StartPoint(); pnt2 = boundaryCurves[0]->EndPoint(); pnt3 = boundaryCurves[1]->StartPoint();
            pnt4 = boundaryCurves[1]->EndPoint(); pnt5 = boundaryCurves[2]->StartPoint(); pnt6 = boundaryCurves[2]->EndPoint();

            // 三边情况，创建退化边构成四边
            std::vector<gp_Pnt> boundaryPoints = { boundaryCurves[0]->StartPoint(), boundaryCurves[1]->StartPoint(), boundaryCurves[2]->StartPoint() };

            // 定义边
            gp_Vec line_01(boundaryPoints[1].XYZ() - boundaryPoints[0].XYZ());
            gp_Vec line_12(boundaryPoints[2].XYZ() - boundaryPoints[1].XYZ());
            gp_Vec line_20(boundaryPoints[0].XYZ() - boundaryPoints[2].XYZ());

            auto calculateAngle = [](const gp_Vec& v1, const gp_Vec& v2)
            {
                Standard_Real dotProduct = v1.Dot(v2);
                Standard_Real magnitudes = v1.Magnitude() * v2.Magnitude();
                return std::acos(dotProduct / magnitudes);  // 返回角度
            };

            // 计算三个角的夹角
            Standard_Real angleAtPoint0 = calculateAngle(-line_20, line_01);  // 点0的夹角
            Standard_Real angleAtPoint1 = calculateAngle(-line_01, line_12);  // 点1的夹角
            Standard_Real angleAtPoint2 = calculateAngle(-line_12, line_20);  // 点2的夹角

            // 找出最大角度
            Standard_Real maxAngle = std::max({ angleAtPoint0, angleAtPoint1, angleAtPoint2 });
            Standard_Integer maxAngleIndex = 0;
            if (maxAngle == angleAtPoint1) maxAngleIndex = 1;
            else if (maxAngle == angleAtPoint2) maxAngleIndex = 2;

            // 构造退化边
            auto CreateDegenerateEdge = [](const gp_Pnt& p1, const gp_Pnt& p2)
            {
                TColgp_Array1OfPnt poles(1, 2);
                poles.SetValue(1, p1);
                poles.SetValue(2, p2);

                TColStd_Array1OfReal knots(1, 2);
                knots.SetValue(1, 0.0);
                knots.SetValue(2, 1.0);

                TColStd_Array1OfInteger multiplicities(1, 2);
                multiplicities.SetValue(1, 2);
                multiplicities.SetValue(2, 2);

                return new Geom_BSplineCurve(poles, knots, multiplicities, 1);
            };
            if (maxAngleIndex == 0)
            {
                //boundaryCurves.push_back(CreateDegenerateEdge(boundaryPoints[maxAngleIndex]));
                boundaryCurves.push_back(CreateDegenerateEdge(boundaryCurves[2]->EndPoint(), boundaryCurves[0]->StartPoint()));
            }
            else
            {
                boundaryCurves.insert(boundaryCurves.begin() + maxAngleIndex, CreateDegenerateEdge(boundaryCurves[maxAngleIndex - 1]->EndPoint(), boundaryCurves[maxAngleIndex]->StartPoint()));
            }
            while (boundaryCurves[3]->StartPoint().Distance(boundaryCurves[3]->EndPoint()) > 10.0)
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
        GetGuideSamples();
    }

    // 获取偏移量和投影点参数
    std::vector<gp_Pnt2d> aPntParams;
    std::vector<gp_Pnt> anOffsets;
    GetSamplesOffset(aPntParams, anOffsets, Standard_True);

    // 若初始已满足容差，则直接返回Coons
    Standard_Real maxDis = 0;
    for (const gp_Pnt& aPnt : anOffsets)
    {
        Standard_Real x = aPnt.X();
        Standard_Real y = aPnt.Y();
        Standard_Real z = aPnt.Z();
        Standard_Real dis = std::sqrt(x * x + y * y + z * z);
        maxDis = std::max(maxDis, dis);
    }
    if (IsLess(maxDis, m_tol))
    {
        m_isDone = Standard_True;
        m_guidedSurf = m_originalSurf;
        return;
    }

    std::vector<Eigen::Vector3d> anOffsetsEigen;
    std::vector<Standard_Real> aPntParamsU, aPntParamsV;
    anOffsetsEigen.reserve(anOffsets.size());
    aPntParamsU.reserve(aPntParams.size());
    aPntParamsV.reserve(aPntParams.size());
    for (const gp_Pnt& offset : anOffsets)
    {
        anOffsetsEigen.push_back(Eigen::Vector3d(offset.X(), offset.Y(), offset.Z()));
    }
    for (const gp_Pnt2d& param : aPntParams)
    {
        aPntParamsU.push_back(param.X());
        aPntParamsV.push_back(param.Y());
    }

    // 获取初始曲面的节点和次数
    const TColStd_Array1OfReal& uKnotsOCC = m_originalSurf->UKnots();
    const TColStd_Array1OfReal& vKnotsOCC = m_originalSurf->VKnots();
    const TColStd_Array1OfInteger& uMultsOCC = m_originalSurf->UMultiplicities();
    const TColStd_Array1OfInteger& vMultsOCC = m_originalSurf->VMultiplicities();
    const Standard_Integer uDeg = m_originalSurf->UDegree();
    const Standard_Integer vDeg = m_originalSurf->VDegree();

    std::vector<Standard_Real> aUKnots, aVKnots;
    aUKnots.reserve(uKnotsOCC.Length());
    aVKnots.reserve(vKnotsOCC.Length());
    for (Standard_Integer i = 1; i <= uKnotsOCC.Length(); ++i)
    {
        for (Standard_Integer j = 1; j <= uMultsOCC[i]; ++j)
        {
            aUKnots.push_back(uKnotsOCC[i]);
        }
    }
    for (Standard_Integer i = 1; i <= vKnotsOCC.Length(); ++i)
    {
        for (Standard_Integer j = 1; j <= vMultsOCC[i]; ++j)
        {
            aVKnots.push_back(vKnotsOCC[i]);
        }
    }

    const Standard_Integer aCtrlPntsUNum = (Standard_Integer)aUKnots.size() - uDeg - 1;
    const Standard_Integer aCtrlPntsVNum = (Standard_Integer)aVKnots.size() - vDeg - 1;
    std::vector<Eigen::Vector3d> aCtrlPnts(aCtrlPntsUNum * aCtrlPntsVNum);

    // 求解偏移曲面控制点
    FitOffsetSurface(anOffsetsEigen, aPntParamsU, aPntParamsV, aUKnots, aVKnots, uDeg, vDeg, aCtrlPnts);

    TColgp_Array2OfPnt aCtrlPntsOCC(1, aCtrlPntsUNum, 1, aCtrlPntsVNum);

    for (Standard_Integer i = 1; i <= aCtrlPntsUNum; ++i)
    {
        for (Standard_Integer j = 1; j <= aCtrlPntsVNum; ++j)
        {
            aCtrlPntsOCC(i, j).SetX(aCtrlPnts[(i - 1) * aCtrlPntsVNum + (j - 1)](0));
            aCtrlPntsOCC(i, j).SetY(aCtrlPnts[(i - 1) * aCtrlPntsVNum + (j - 1)](1));
            aCtrlPntsOCC(i, j).SetZ(aCtrlPnts[(i - 1) * aCtrlPntsVNum + (j - 1)](2));
        }
    }

    // 创建新的引导后的曲面
    Handle(Geom_BSplineSurface) offsetSurf = new Geom_BSplineSurface(aCtrlPntsOCC, uKnotsOCC, vKnotsOCC, uMultsOCC, vMultsOCC, uDeg, vDeg);

    // 叠加曲面
    // 共同节点和次数
    const TColStd_Array1OfReal& knotsUCommon = m_originalSurf->UKnots();
    const TColStd_Array1OfReal& knotsVCommon = m_originalSurf->VKnots();
    const TColStd_Array1OfInteger& multsUCommon = m_originalSurf->UMultiplicities();
    const TColStd_Array1OfInteger& multsVCommon = m_originalSurf->VMultiplicities();
    const Standard_Integer degUCommon = m_originalSurf->UDegree();
    const Standard_Integer degVCommon = m_originalSurf->VDegree();

    // 计算控制点 
    const TColgp_Array2OfPnt& originalCtrlPnts = m_originalSurf->Poles();
    const TColgp_Array2OfPnt& offsetCtrlPnts = offsetSurf->Poles();
    const Standard_Integer ctrlPntsUNumCom = m_originalSurf->NbUPoles();
    const Standard_Integer ctrlPntsVNumCom = m_originalSurf->NbVPoles();

    TColgp_Array2OfPnt ctrlPntsAdd(1, ctrlPntsUNumCom, 1, ctrlPntsVNumCom);

    for (Standard_Integer i = 1; i <= ctrlPntsUNumCom; ++i)
    {
        for (Standard_Integer j = 1; j <= ctrlPntsVNumCom; ++j)
        {
            gp_XYZ aCoord = originalCtrlPnts(i, j).Coord() + STEP_LENGTH * aCtrlPntsOCC(i, j).Coord();
            ctrlPntsAdd(i, j).SetCoord(aCoord.X(), aCoord.Y(), aCoord.Z());
        }
    }

    m_guidedSurf = new Geom_BSplineSurface(ctrlPntsAdd, knotsUCommon,
        knotsVCommon, multsUCommon, multsVCommon, degUCommon, degVCommon);
    
    /*
    // 0.9 0.1
    const TColgp_Array2OfPnt& guidedCtrlPnts = m_guidedSurf->Poles();
    TColgp_Array2OfPnt ctrlPntsAdd2(1, ctrlPntsUNumCom, 1, ctrlPntsVNumCom);

    for (Standard_Integer i = 1; i <= ctrlPntsUNumCom; ++i)
    {
        for (Standard_Integer j = 1; j <= ctrlPntsVNumCom; ++j)
        {
            gp_XYZ aCoord = 0.9 * originalCtrlPnts(i, j).Coord() + 0.1 * guidedCtrlPnts(i, j).Coord();
            ctrlPntsAdd2(i, j).SetCoord(aCoord.X(), aCoord.Y(), aCoord.Z());
        }
    }

    m_guidedSurf = new Geom_BSplineSurface(ctrlPntsAdd2, knotsUCommon,
        knotsVCommon, multsUCommon, multsVCommon, degUCommon, degVCommon);
    */

    // 获取新曲面的采样点对应投影参数和误差
    std::vector<gp_Pnt2d> aPntParams2;
    std::vector<gp_Pnt> anOffsets2;
    GetSamplesOffset(aPntParams2, anOffsets2, Standard_False);

    std::vector<Standard_Real> errorU, errorV;

    Standard_Real maxErrDis = 0;
    Standard_Integer pos = 0;
    for (const gp_Pnt& aPnt : anOffsets2)
    {
        Standard_Real x = aPnt.X();
        Standard_Real y = aPnt.Y();
        Standard_Real z = aPnt.Z();
        Standard_Real dis = std::sqrt(x * x + y * y + z * z);
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
    std::vector<Standard_Real> knotsUToInsert, knotsVToInsert;
    Standard_Integer left = 0, right = 0;
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
        for (Standard_Integer i = 0; i < knotsUToInsert.size(); ++i)
        {
            m_guidedSurf->InsertUKnot(knotsUToInsert[i], 1, Precision::Confusion());
        }
        for (Standard_Integer i = 0; i < knotsVToInsert.size(); ++i)
        {
            m_guidedSurf->InsertVKnot(knotsVToInsert[i], 1, Precision::Confusion());
        }

        m_isDone = Standard_False;
    }
    else 
    {
        m_isDone = Standard_True; // 所有点都容差内
    }
}

void GuidedCoonsSurfGenerator::GetGuideSamples()
{
    Standard_Integer crvIdx = 0;
    for (const Handle(Geom_BSplineCurve)& guideCurve : m_guideCurves)
    {
        // 得到采样点
        std::vector<gp_Pnt> tempSamples;
        //samples = SampleGuideCurve(guideCurve, 0, 1, theSamplesNum);
        for (const auto& trimInterval : m_guideCurvesTrimIntervals[crvIdx])
        {
            // 根据长度自适应采样
            GeomAdaptor_Curve curveAdaptor(guideCurve, trimInterval.first, trimInterval.second);
            Standard_Real length = GCPnts_AbscissaPoint::Length(curveAdaptor, curveAdaptor.FirstParameter(), 
                curveAdaptor.LastParameter(), Precision::Confusion());
            Standard_Integer sampleNum = length / GUIDE_SAMPLING_INTERVAL + 1; // 至少采一个点

            tempSamples = SampleGuideCurve(guideCurve, trimInterval.first, trimInterval.second, sampleNum);
            for (const auto& sample : tempSamples)
            {
                m_samples.push_back(sample);
            }
        }
        crvIdx++;
    }
}

void GuidedCoonsSurfGenerator::GetSamplesOffset(std::vector<gp_Pnt2d>& thePntParams, std::vector<gp_Pnt>& theOffsets, Standard_Boolean isOriginal)
{
    // 得到投影点和参数
    std::vector<gp_Pnt2d> aPntParams;
    std::vector<gp_Pnt> projectionPoints;
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

std::vector<gp_Pnt> GuidedCoonsSurfGenerator::SampleGuideCurve(const Handle(Geom_BSplineCurve)& theCurve, Standard_Real startParam, Standard_Real endParam, Standard_Integer theSamplesNum)
{
    std::vector<gp_Pnt> samples;
    for (Standard_Integer i = 1; i <= theSamplesNum; ++i)
    {
        //Standard_Real t = theCurve->FirstParameter() + i * (theCurve->LastParameter() - theCurve->FirstParameter()) / (theSamplesNum + 1);
        Standard_Real t = startParam + i * (endParam - startParam) / (theSamplesNum + 1);
        gp_Pnt sample = theCurve->Value(t);
        samples.push_back(sample);
    }
    return samples;
}

std::vector<gp_Pnt> GuidedCoonsSurfGenerator::ProjectPntsToSurf(const std::vector<gp_Pnt>& thePoints, std::vector<gp_Pnt>& theProjectionPoints, const Handle(Geom_BSplineSurface)& theSurface, std::vector<gp_Pnt2d>& thePntParams)
{
    std::vector<gp_Pnt> offsets;
    for (const auto& point : thePoints)
    {
        // 创建投影对象，不指定参数域
        GeomAPI_ProjectPointOnSurf projector(point, theSurface);
        Standard_Real u, v;

        if (projector.NbPoints() > 0)
        {
            projector.LowerDistanceParameters(u, v);
            if (IsGreater(u, 0) && IsLess(u, 1) && IsGreater(v, 0) && IsLess(v, 1))
            {
                gp_Pnt projectedPoint = projector.NearestPoint();
                gp_Pnt offset = point.XYZ() - projectedPoint.XYZ();
                if (std::abs(offset.X()) > MAX_OFFSETDISTANCE ||
                    std::abs(offset.Y()) > MAX_OFFSETDISTANCE ||
                    std::abs(offset.Z()) > MAX_OFFSETDISTANCE)
                {
                    continue;
                }
                offsets.push_back(offset);
                theProjectionPoints.push_back(projectedPoint);
                gp_Pnt2d uv(u, v);
                thePntParams.push_back(uv);
            }
        }
    }

    return offsets;
}

std::vector<gp_Pnt> GuidedCoonsSurfGenerator::CalOffsets(const std::vector<gp_Pnt>& theSamples, const std::vector<gp_Pnt>& theProjections)
{
    std::vector<gp_Pnt> offsets;
    for (Standard_Integer i = 0; i < theSamples.size(); ++i)
    {
        gp_Pnt offset = theSamples[i].XYZ() - theProjections[i].XYZ();
        offsets.push_back(offset);
    }
    return offsets;
}

Eigen::MatrixXd ConstructConvMat(Standard_Integer theRow, Standard_Integer theCol, const Eigen::Matrix3d& kernel)
{
    // 输出矩阵的尺寸：展平后的像素向量长度为 theRow * theCol
    // 卷积后输出的尺寸为 (theRow-2) * (theCol-2)
    Standard_Integer outputSize = (theRow - 2) * (theCol - 2);
    Standard_Integer inputSize = theRow * theCol;
    Eigen::MatrixXd convMat = Eigen::MatrixXd::Zero(outputSize, inputSize);

    for (Standard_Integer i = 1; i < theRow - 1; ++i)
    {
        for (Standard_Integer j = 1; j < theCol - 1; ++j)
        {
            Standard_Integer outputIdx = (i - 1) * (theCol - 2) + (j - 1);
            for (Standard_Integer ki = -1; ki <= 1; ++ki)
            {
                for (Standard_Integer kj = -1; kj <= 1; ++kj)
                {
                    Standard_Integer inputRow = i + ki;
                    Standard_Integer inputCol = j + kj;
                    Standard_Integer inputIdx = inputRow * theCol + inputCol;
                    convMat(outputIdx, inputIdx) = kernel(ki + 1, kj + 1);
                }
            }
        }
    }

    return convMat;
}

Eigen::MatrixXd ConstructVariableConvMat(Standard_Integer ctrlPtsUNum, Standard_Integer ctrlPtsVNum, const std::vector<Standard_Real>& lDiffs, const std::vector<Standard_Real>& sDiffs)
{
    //
    if (ctrlPtsUNum < 3 || ctrlPtsVNum < 3)
    {
        return Eigen::MatrixXd(0, ctrlPtsUNum * ctrlPtsVNum);
    }

    // U和V方向内部点的数量
    Standard_Integer internalUpoints = ctrlPtsUNum - 2;
    Standard_Integer internalVpoints = ctrlPtsVNum - 2;

    Standard_Integer outputRows = internalUpoints * internalVpoints;
    Standard_Integer inputCols = ctrlPtsUNum * ctrlPtsVNum;

    Eigen::MatrixXd convMat = Eigen::MatrixXd::Zero(outputRows, inputCols);

    for (Standard_Integer i = 1; i <= internalUpoints; ++i)
    {
        for (Standard_Integer j = 1; j <= internalVpoints; ++j)
        {

            // 当前行在输出矩阵中的索引
            Standard_Integer outputIdx = (i - 1) * internalVpoints + (j - 1);

            Standard_Real lPrev = lDiffs[i - 1];
            Standard_Real lCurr = lDiffs[i];


            Standard_Real sPrev = sDiffs[j - 1];
            Standard_Real sCurr = sDiffs[j];

            Standard_Real denL = lPrev + lCurr;
            Standard_Real denS = sPrev + sCurr;

            if (std::abs(denL) < Precision::Confusion()) denL = 1.0; // 避免除零，如果lPrev=lCurr=0，则权重为0
            if (std::abs(denS) < Precision::Confusion()) denS = 1.0;

            // P_{i-1,j} 的权重
            Standard_Real wIm1J = lCurr / (2.0 * denL);
            // P_{i+1,j} 的权重
            Standard_Real wIp1J = lPrev / (2.0 * denL);
            // P_{i,j-1} 的权重
            Standard_Real wIJm1 = sCurr / (2.0 * denS);
            // P_{i,j+1} 的权重
            Standard_Real wIJp1 = sPrev / (2.0 * denS);

            // P_center (P_i,j) 在扁平化输入向量中的列索引
            // 在原始的 ctrlPtsUNum x ctrlPtsVNum 网格中，P_i,j 的索引就是 (i,j)。
            Standard_Integer idxCenter = i * ctrlPtsVNum + j;
            convMat(outputIdx, idxCenter) = -1.0;

            // P_{i-1,j} 的列索引
            Standard_Integer idxIm1J = (i - 1) * ctrlPtsVNum + j;
            convMat(outputIdx, idxIm1J) = wIm1J;

            // P_{i+1,j} 的列索引
            Standard_Integer idxIp1J = (i + 1) * ctrlPtsVNum + j;
            convMat(outputIdx, idxIp1J) = wIp1J;

            // P_{i,j-1} 的列索引
            Standard_Integer idxIJm1 = i * ctrlPtsVNum + (j - 1);
            convMat(outputIdx, idxIJm1) = wIJm1;

            // P_{i,j+1} 的列索引
            Standard_Integer idxIJp1 = i * ctrlPtsVNum + (j + 1);
            convMat(outputIdx, idxIJp1) = wIJp1;
        }
    }
    return convMat;
}

std::vector<Standard_Real> CalGrevilleAbscissae1D(const TColStd_Array1OfReal& knots, Standard_Integer ctrlPtsNum, Standard_Integer degree)
{
    std::vector<Standard_Real> greville_coords(ctrlPtsNum);

    for (Standard_Integer i = 0; i < ctrlPtsNum; ++i)
    {
        Standard_Real sum = 0.0;
        for (Standard_Integer j = 1; j <= degree; ++j)
        {
            Standard_Integer knotIdx = i + j + 1;
            sum += knots.Value(knotIdx);
        }
        greville_coords[i] = sum / static_cast<Standard_Real>(degree);
    }
    return greville_coords;
}

std::pair<std::vector<Standard_Real>, std::vector<Standard_Real>> CalGrevilleCoordDiffs(const Handle(Geom_BSplineSurface)& surface)
{
    // U 方向
    Standard_Integer ctrlPtsUNum = surface->NbUPoles();
    const TColStd_Array1OfReal& uKnots = surface->UKnotSequence();

    // 计算U方向的Greville坐标
    std::vector<Standard_Real> grevilleU = CalGrevilleAbscissae1D(uKnots, ctrlPtsUNum, surface->UDegree());

    std::vector<Standard_Real> lDiffs; // 存储U方向的差值
    if (ctrlPtsUNum > 1) {
        lDiffs.resize(ctrlPtsUNum - 1);
        for (Standard_Integer i = 0; i < ctrlPtsUNum - 1; ++i) {
            lDiffs[i] = grevilleU[i + 1] - grevilleU[i];
        }
    }

    // V 方向
    Standard_Integer ctrlPtsVNum = surface->NbVPoles();
    const TColStd_Array1OfReal& vKnots = surface->VKnotSequence();

    // 计算V方向的Greville坐标
    std::vector<Standard_Real> grevilleV = CalGrevilleAbscissae1D(vKnots, ctrlPtsVNum, surface->VDegree());

    std::vector<Standard_Real> sDiffs; // 存储V方向的差值
    if (ctrlPtsVNum > 1) {
        sDiffs.resize(ctrlPtsVNum - 1);
        for (Standard_Integer j = 0; j < ctrlPtsVNum - 1; ++j) {
            sDiffs[j] = grevilleV[j + 1] - grevilleV[j];
        }
    }

    return { lDiffs, sDiffs }; // 返回包含 l 和 s 差值向量的pair
}





Eigen::MatrixXd GuidedCoonsSurfGenerator::ConstructCurveSmoothingMatrix(
    const Handle(Geom_BSplineCurve)& theBSplineCurve,
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
    const Handle(Geom_BSplineSurface)& theBSplineSurface,
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
        const Lobatto::QuadratureData* u_lobatto_data = getQuadratureData(u_sample_num);
        if (u_lobatto_data != nullptr) {
            u_weights = u_lobatto_data->weights_01;
        }
        else {
            // 如果获取失败，使用均匀权重
            u_weights.resize(u_sample_num, 1.0);
        }

        // 获取v方向的权重
        int v_sample_num = v_params.size();
        const Lobatto::QuadratureData* v_lobatto_data = getQuadratureData(v_sample_num);
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
    std::vector<Standard_Real> utraces;
    std::vector<Standard_Real> vtraces;

    // ========================================================================
    // 1. 处理u-等参线（固定u值，v方向变化）
    // ========================================================================

    // 对每条u-等参线进行处理
    for (size_t u_idx = 0; u_idx < u_params.size(); ++u_idx) {
        double u_val = u_params[u_idx];
        double u_weight = u_weights[u_idx];
        // 从曲面中提取u-等参线（固定u，v变化）
        Handle(Geom_BSplineCurve) u_isoline = Handle(Geom_BSplineCurve)::DownCast(
            theBSplineSurface->UIso(u_val)
        );

        // 构造该等参线的光顺能量矩阵 M_v
        Eigen::MatrixXd M_v;
        M_v = ConstructCurveSmoothingMatrix(u_isoline, derivative_order, tolerance);

        // 构造u-等参线的提取矩阵 C_u
        Eigen::MatrixXd C_u;
        C_u = ConstructDenseMatrix(n_u, n_v, p_u, p_v, knots_u, knots_v, 0, u_val);

        // 计算该u-等参线对总能量矩阵的贡献: C_u^T * M_v * C_u
        Eigen::MatrixXd contribution = C_u.transpose() * M_v * C_u;
        Standard_Real normal_factor = 1.0 / contribution.trace();
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
        Handle(Geom_BSplineCurve) v_isoline = Handle(Geom_BSplineCurve)::DownCast(
            theBSplineSurface->VIso(v_val)
        );

        // 构造该等参线的光顺能量矩阵 M_u
        Eigen::MatrixXd M_u;
        M_u = ConstructCurveSmoothingMatrix(v_isoline, derivative_order, tolerance);

        // 构造v-等参线的提取矩阵 D_v
        Eigen::MatrixXd D_v;
        D_v = ConstructDenseMatrix(n_u, n_v, p_u, p_v, knots_u, knots_v, 1, v_val);

        // 计算该v-等参线对总能量矩阵的贡献: D_v^T * M_u * D_v
        Eigen::MatrixXd contribution = D_v.transpose() * M_u * D_v;
        Standard_Real normal_factor = 1.0 / contribution.trace();
        vtraces.push_back(contribution.trace());
        M_total += normal_factor * v_weight * contribution;
    }

     return M_total;
}

Eigen::MatrixXd GuidedCoonsSurfGenerator::ConstructBidirectionalSmoothingMatrix(
    const Handle(Geom_BSplineSurface)& theBSplineSurface,
    const std::vector<double>& u_params,
    const std::vector<double>& v_params,
    int derivative_order,
    double tolerance
) {
    // 从曲面对象中提取参数
    int n_u = theBSplineSurface->NbUPoles() - 1;
    int n_v = theBSplineSurface->NbVPoles() - 1;
    int p_u = theBSplineSurface->UDegree();
    int p_v = theBSplineSurface->VDegree();

    // 提取节点向量
    std::vector<double> knots_u, knots_v;

    // 提取U方向节点向量
    for (int i = 1; i <= theBSplineSurface->NbUKnots(); i++) {
        double knot = theBSplineSurface->UKnot(i);
        int mult = theBSplineSurface->UMultiplicity(i);
        for (int j = 0; j < mult; j++) {
            knots_u.push_back(knot);
        }
    }

    // 提取V方向节点向量
    for (int i = 1; i <= theBSplineSurface->NbVKnots(); i++) {
        double knot = theBSplineSurface->VKnot(i);
        int mult = theBSplineSurface->VMultiplicity(i);
        for (int j = 0; j < mult; j++) {
            knots_v.push_back(knot);
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
Standard_Boolean IsLess(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
{
    return (y - x) > tol;
}

//! @brief 在容差意义下比较 x 是否大于等于 y
//! @param [In] x 第一个数
//! @param [In] y 第二个数
//! @return x 大于等于 y 则返回true， 否则返回false
Standard_Boolean IsGreaterOrEqual(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
{
    return (x - y) > -tol;
}

Standard_Real CalBasicFunction(Standard_Real param, Standard_Integer index, Standard_Integer deg, const std::vector<Standard_Real>& knots)
{
    Standard_Real nip, uleft, uright, saved, temp;
    Standard_Integer m = (Standard_Integer)knots.size() - 1;
    std::vector<Standard_Real> N(deg + 1);

    if ((index == 0 && IsEqual(param, knots[0])) || (index == m - deg - 1 && IsEqual(param, knots[m])))
    {
        return 1.0;
    }
    if (IsLess(param, knots[index]) || IsGreaterOrEqual(param, knots[index + deg + 1]))
    {
        return 0.0;
    }
    for (Standard_Integer j = 0; j <= deg; ++j)
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
    for (Standard_Integer k = 1; k <= deg; ++k)
    {
        if (N[0] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = ((param - knots[index]) * N[0]) / (knots[index + k] - knots[index]);
        }
        for (Standard_Integer j = 0; j < deg - k + 1; ++j)
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

std::vector<Standard_Real> GenerateSamplePoints(
    Standard_Integer theNumSamples, const std::vector<Standard_Real>& theKnots, Standard_Integer theDegree)
{
    std::vector<Standard_Real> samples;
    samples.reserve(theNumSamples);

    // 在节点区间内均匀采样，避开重复节点
    Standard_Real start = theKnots[theDegree];
    Standard_Real end = theKnots[theKnots.size() - theDegree - 1];

    for (Standard_Integer i = 0; i < theNumSamples; ++i) {
        Standard_Real t = start + (end - start) * i / (theNumSamples - 1);
        samples.push_back(t);
    }

    return samples;
}


Standard_Real CalBasicFunctionDerivative(
    Standard_Real theParam, Standard_Integer theIndex,
    Standard_Integer theDegree, const std::vector<Standard_Real>& theKnots,
    Standard_Integer theDerivOrder)
{
    // B样条基函数导数计算 - 使用递归公式
    if (theDerivOrder == 0) {
        return CalBasicFunction(theParam, theIndex, theDegree, theKnots);
    }

    // 一阶导数公式
    if (theDegree == 0) {
        return 0.0; // 0次B样条导数恒为0
    }

    Standard_Real left = 0.0, right = 0.0;
    Standard_Real denomLeft = theKnots[theIndex + theDegree] - theKnots[theIndex];
    Standard_Real denomRight = theKnots[theIndex + theDegree + 1] - theKnots[theIndex + 1];

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
    Standard_Real theParam, Standard_Integer theCtrlPtsNum,
    Standard_Integer theDegree, const std::vector<Standard_Real>& theKnots,
    Standard_Integer theDerivOrder)
{
    Eigen::VectorXd result = Eigen::VectorXd::Zero(theCtrlPtsNum);

    for (Standard_Integer i = 0; i < theCtrlPtsNum; ++i) {
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


// 偏移曲面拟合
void GuidedCoonsSurfGenerator::FitOffsetSurface(const std::vector<Eigen::Vector3d>& theSamplePntOffsets, const std::vector<Standard_Real>& thePntParamsU, const std::vector<Standard_Real>& thePntParamsV,
    const std::vector<Standard_Real>& theUKnots, const std::vector<Standard_Real>& theVKnots, Standard_Integer theDegU, Standard_Integer theDegV, std::vector<Eigen::Vector3d>& theCtrlPoints)
{
    Standard_Integer aCtrlPtsUNum = (Standard_Integer)theUKnots.size() - theDegU - 1;
    Standard_Integer aCtrlPtsVNum = (Standard_Integer)theVKnots.size() - theDegV - 1;

    // 构建M矩阵: 边界控制点约束
    Eigen::MatrixXd M;
    BuildMatrixConstraint(aCtrlPtsUNum, aCtrlPtsVNum, M); // 约束项为边界控制点
    Eigen::MatrixXd MT = M.transpose();

    // 构造N矩阵：偏移向量能量矩阵
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(thePntParamsU.size(), aCtrlPtsUNum * aCtrlPtsVNum);
    Eigen::MatrixXd Ni(thePntParamsU.size(), aCtrlPtsUNum);
    Eigen::MatrixXd Nj(thePntParamsV.size(), aCtrlPtsVNum);
    Standard_Real valueTemp = 0.0;
    for (Standard_Integer i = 0; i < thePntParamsU.size(); ++i)
    {
        for (Standard_Integer j = 0; j < aCtrlPtsUNum; ++j)
        {
            valueTemp = CalBasicFunction(thePntParamsU[i], j, theDegU, theUKnots);
            Ni(i, j) = valueTemp;
        }
        for (Standard_Integer j = 0; j < aCtrlPtsVNum; ++j)
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
    for (Standard_Integer i = 0; i < (Standard_Integer)O.rows(); ++i)
    {
        O(i, 0) = theSamplePntOffsets[i](0);
        O(i, 1) = theSamplePntOffsets[i](1);
        O(i, 2) = theSamplePntOffsets[i](2);
    }

    TColgp_Array2OfPnt originalCtrlPnts = m_originalSurf->Poles();
    Eigen::MatrixXd P0 = Eigen::MatrixXd::Zero(originalCtrlPnts.Size(), 3);
    int t = 0;
    for (Standard_Integer i = 1; i <= aCtrlPtsUNum; ++i)
    {
        for (Standard_Integer j = 1; j <= aCtrlPtsVNum; ++j)
        {
            P0(t, 0) = originalCtrlPnts(i, j).X();
            P0(t, 1) = originalCtrlPnts(i, j).Y();
            P0(t, 2) = originalCtrlPnts(i, j).Z();
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

    for (Standard_Integer i = 0; i < aCtrlPtsUNum * aCtrlPtsVNum; ++i)
    {
        theCtrlPoints[i] = Eigen::Vector3d(P(i, 0), P(i, 1), P(i, 2));
        //std::cout << P(i, 0) << " " << P(i, 1) << " " << P(i, 2) << std::endl;
    }
    */
    
    // 构建光顺能量矩阵
    Eigen::MatrixXd H, b;
    Eigen::SparseMatrix<Standard_Real> H_sparse, b_sparse;
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
            const Lobatto::QuadratureData* lobatto_data = getQuadratureData(sample_num);

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
        //Standard_Real trace_D = D.trace();
        //Eigen::MatrixXd NTN = N.transpose() * N;
        //Standard_Real trace_NTN = NTN.trace();

        // 归一化
        //Standard_Real lambda_trace = 1.0; // 默认值
        //if (IsGreater(std::abs(trace_D), 1e-12))
        //{
        //    lambda_trace = trace_NTN / trace_D;
        //}
        //D = D * lambda_trace;

        //Eigen::MatrixXd H = (1 - FAIRNESS_WEIGHT) * (NT * N) + FAIRNESS_WEIGHT * D;
        //Eigen::MatrixXd b = (1 - FAIRNESS_WEIGHT) * (N.transpose() * O) - FAIRNESS_WEIGHT * (D * P0);

        //Eigen::MatrixXd I = Eigen::MatrixXd::Identity(H.rows(), H.cols());
        //H = H + 0.00001 * I;


        Eigen::SparseMatrix<Standard_Real> N_sparse = N.sparseView();
        Eigen::SparseMatrix<Standard_Real> NT_sparse = N.transpose().sparseView();
        Eigen::SparseMatrix<Standard_Real> D_sparse = D.sparseView();
        Eigen::SparseMatrix<Standard_Real> DT_sparse = D.transpose().sparseView();
        Eigen::SparseMatrix<Standard_Real> O_sparse = O.sparseView();
        Eigen::SparseMatrix<Standard_Real> OT_sparse = O.transpose().sparseView();
        Eigen::SparseMatrix<Standard_Real> P0_sparse = P0.sparseView();
        Eigen::SparseMatrix<Standard_Real> P0T_sparse = P0.transpose().sparseView();

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

    Standard_Integer totalPts = aCtrlPtsUNum * aCtrlPtsVNum;

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

    Standard_Integer numInner = innerIndices.size();

    if (numInner > 0) {
        // 2. 提取子矩阵 H_ii 和 子向量 b_i
        Eigen::SparseMatrix<Standard_Real> H_ii(numInner, numInner);
        Eigen::MatrixXd b_i(numInner, 3);

        std::vector<Eigen::Triplet<Standard_Real>> H_ii_triplets;

        // 遍历 H_sparse，提取内部点到内部点的块
        for (int k=0; k < H_sparse.outerSize(); ++k) {
            for (Eigen::SparseMatrix<Standard_Real>::InnerIterator it(H_sparse,k); it; ++it) {
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
                    H_ii_triplets.push_back(Eigen::Triplet<Standard_Real>(newRow, newCol, it.value()));
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
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<Standard_Real>> solver;
        solver.compute(H_ii);

        if (solver.info() != Eigen::Success) {
            std::cout << "Warning: SimplicialLDLT failed, falling back to SparseLU" << std::endl;
            Eigen::SparseLU<Eigen::SparseMatrix<Standard_Real>> luSolver;
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

}

// 构建约束项矩阵
void GuidedCoonsSurfGenerator::BuildMatrixConstraint(Standard_Integer theCtrlPtsUNum, Standard_Integer theCtrlPtsVNum, Eigen::MatrixXd& theMatrixM)
{
    // 计算矩阵 M 的尺寸
    Standard_Integer numRows = 2 * (theCtrlPtsUNum + theCtrlPtsVNum - 2);
    Standard_Integer numCols = theCtrlPtsUNum * theCtrlPtsVNum;

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
    Standard_Integer rowstart = theCtrlPtsVNum + 2 * theCtrlPtsUNum - 4;
    Standard_Integer colstart = (theCtrlPtsUNum - 1) * theCtrlPtsVNum;
    Standard_Integer blockrows = theCtrlPtsVNum;
    Standard_Integer blockcols = theCtrlPtsVNum;
    theMatrixM.block(rowstart, colstart, blockrows, blockcols) = I;

    // 填充矩阵 M 中 J 的部分
    for (Standard_Integer i = 0; i < theCtrlPtsUNum - 2; i++)
    {
        rowstart = theCtrlPtsVNum + i * 2;
        colstart = theCtrlPtsVNum + i * theCtrlPtsVNum;
        blockrows = 2;
        blockcols = theCtrlPtsVNum;
        theMatrixM.block(rowstart, colstart, blockrows, blockcols) = J;
    }

}

// 构建非约束项矩阵
void GuidedCoonsSurfGenerator::BuildMatrixUnconstraint(const std::vector<Standard_Real>& thePntParamsU, const std::vector<Standard_Real>& thePntParamsV, const std::vector<Standard_Real>& theUKnots, const std::vector<Standard_Real>& theVKnots,
    Standard_Integer theDegU, Standard_Integer theDegV, Standard_Integer theCtrlPtsUNum, Standard_Integer theCtrlPtsVNum, Eigen::MatrixXd& theMatrixN)
{
    // 构造矩阵 N^
    Eigen::MatrixXd Nhat;
    Eigen::MatrixXd Ni(thePntParamsU.size(), theCtrlPtsUNum);
    Eigen::MatrixXd Nj(thePntParamsV.size(), theCtrlPtsVNum);
    Standard_Real valueTemp = 0.0;
    for (Standard_Integer i = 0; i < thePntParamsU.size(); ++i)
    {
        for (Standard_Integer j = 0; j < theCtrlPtsUNum; ++j)
        {
            valueTemp = CalBasicFunction(thePntParamsU[i], j, theDegU, theUKnots);
            Ni(i, j) = valueTemp;
        }
        for (Standard_Integer j = 0; j < theCtrlPtsVNum; ++j)
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

}

// 构建能量权重矩阵 
void GuidedCoonsSurfGenerator::BuildMatrixWeight(Standard_Integer thePntParamsSize, Standard_Integer theCtrlPtsUNum, Standard_Integer theCtrlPtsVNum, Standard_Real alpha, Eigen::MatrixXd& theMatrixW)
{
    Standard_Integer size = (theCtrlPtsUNum - 2) * (theCtrlPtsVNum - 2);
    Eigen::MatrixXd Alpha = Eigen::MatrixXd::Zero(size, size);
    for (Standard_Integer i = 0; i < size; ++i)
    {
        Alpha(i, i) = alpha;
    }

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(thePntParamsSize, thePntParamsSize);
    for (Standard_Integer i = 0; i < thePntParamsSize; ++i)
    {
        I(i, i) = 1 - alpha;
    }

    theMatrixW = Eigen::MatrixXd::Zero(Alpha.rows() + I.rows(), Alpha.rows() + I.rows());
    theMatrixW.block(0, 0, Alpha.rows(), Alpha.cols()) = Alpha;
    theMatrixW.block(Alpha.rows(), Alpha.cols(), I.rows(), I.cols()) = I;

}

Standard_Real GuidedCoonsSurfGenerator::CalBasicFunction(Standard_Real param, Standard_Integer index, Standard_Integer deg, const std::vector<Standard_Real>& knots)
{
    Standard_Real nip, uleft, uright, saved, temp; 
    Standard_Integer m = (Standard_Integer)knots.size() - 1;
    std::vector<Standard_Real> N(deg + 1);

    if ((index == 0 && IsEqual(param, knots[0])) || (index == m - deg - 1 && IsEqual(param, knots[m])))
    {
        return 1.0;
    }
    if (IsLess(param, knots[index]) || IsGreaterOrEqual(param, knots[index + deg + 1]))
    {
        return 0.0;
    }
    for (Standard_Integer j = 0; j <= deg; ++j)
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
    for (Standard_Integer k = 1; k <= deg; ++k)
    {
        if (N[0] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = ((param - knots[index]) * N[0]) / (knots[index + k] - knots[index]);
        }
        for (Standard_Integer j = 0; j < deg - k + 1; ++j)
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
    Standard_Integer col_num1, col_num2, row_num1, row_num2;
    col_num1 = (Standard_Integer)theMatA.cols();
    row_num1 = (Standard_Integer)theMatA.rows();
    col_num2 = (Standard_Integer)theMatB.cols();
    row_num2 = (Standard_Integer)theMatB.rows();
    Eigen::MatrixXd resMat(row_num1, col_num1 * col_num2);
    Standard_Integer index = 0;
    for (Standard_Integer rowLoopIndex = 0; rowLoopIndex < row_num1; ++rowLoopIndex)
    {
        index = 0;
        for (Standard_Integer i = 0; i < col_num1; ++i)
        {
            for (Standard_Integer j = 0; j < col_num2; ++j)
            {
                resMat(rowLoopIndex, index++) = theMatA(rowLoopIndex, i) * theMatB(rowLoopIndex, j);
            }
        }
    }
    return resMat;
}

Standard_Integer GuidedCoonsSurfGenerator::SetSameDistribution(Handle(Geom_BSplineCurve)& C1, Handle(Geom_BSplineCurve)& C2)
{
    Standard_Integer nbp1 = C1->NbPoles();
    Standard_Integer nbk1 = C1->NbKnots();
    TColgp_Array1OfPnt      P1(1, nbp1);
    TColStd_Array1OfReal    W1(1, nbp1);
    W1.Init(1.);
    TColStd_Array1OfReal    K1(1, nbk1);
    TColStd_Array1OfInteger M1(1, nbk1);

    C1->Poles(P1);
    if (C1->IsRational())
        C1->Weights(W1);
    C1->Knots(K1);
    C1->Multiplicities(M1);

    Standard_Integer nbp2 = C2->NbPoles();
    Standard_Integer nbk2 = C2->NbKnots();
    TColgp_Array1OfPnt      P2(1, nbp2);
    TColStd_Array1OfReal    W2(1, nbp2);
    W2.Init(1.);
    TColStd_Array1OfReal    K2(1, nbk2);
    TColStd_Array1OfInteger M2(1, nbk2);

    C2->Poles(P2);
    if (C2->IsRational())
        C2->Weights(W2);
    C2->Knots(K2);
    C2->Multiplicities(M2);

    Standard_Real K11 = K1(1);
    Standard_Real K12 = K1(nbk1);
    Standard_Real K21 = K2(1);
    Standard_Real K22 = K2(nbk2);

    if ((K12 - K11) > (K22 - K21)) {
        BSplCLib::Reparametrize(K11, K12, K2);
        C2->SetKnots(K2);
    }
    else if ((K12 - K11) < (K22 - K21)) {
        BSplCLib::Reparametrize(K21, K22, K1);
        C1->SetKnots(K1);
    }
    else if (Abs(K12 - K11) > Precision::PConfusion()) {
        BSplCLib::Reparametrize(K11, K12, K2);
        C2->SetKnots(K2);
    }

    Standard_Integer NP, NK;
    if (BSplCLib::PrepareInsertKnots(C1->Degree(), Standard_False,
        K1, M1, K2, &M2, NP, NK, Precision::PConfusion(),
        Standard_False)) {
        TColgp_Array1OfPnt      NewP(1, NP);
        TColStd_Array1OfReal    NewW(1, NP);
        TColStd_Array1OfReal    NewK(1, NK);
        TColStd_Array1OfInteger NewM(1, NK);
        BSplCLib::InsertKnots(C1->Degree(), Standard_False,
            P1, &W1, K1, M1, K2, &M2,
            NewP, &NewW, NewK, NewM, Precision::PConfusion(),
            Standard_False);
        if (C1->IsRational()) {
            C1 = new Geom_BSplineCurve(NewP, NewW, NewK, NewM, C1->Degree());
        }
        else {
            C1 = new Geom_BSplineCurve(NewP, NewK, NewM, C1->Degree());
        }
        BSplCLib::InsertKnots(C2->Degree(), Standard_False,
            P2, &W2, K2, M2, K1, &M1,
            NewP, &NewW, NewK, NewM, Precision::PConfusion(),
            Standard_False);
        if (C2->IsRational()) {
            C2 = new Geom_BSplineCurve(NewP, NewW, NewK, NewM, C2->Degree());
        }
        else {
            C2 = new Geom_BSplineCurve(NewP, NewK, NewM, C2->Degree());
        }
    }
    else {
        throw Standard_ConstructionError(" ");
    }

    return C1->NbPoles();
}