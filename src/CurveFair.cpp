#include "CurveFair.h"
//#include "CPnts/CPnts_AbscissaPoint.cxx"
// 定义静态成员变量
Handle(Geom_BSplineCurve) CurveFair::m_ArcLengthMappingFunction;

bool isEqual(Standard_Real x, Standard_Real y, Standard_Real epsilon = 1e-10)
{
    return std::fabs(x - y) < epsilon;
}
bool isGreaterThan(Standard_Real x, Standard_Real y, Standard_Real epsilon = 1e-10)
{
    return (x - y) > epsilon;
}
bool isLessThan(Standard_Real x, Standard_Real y, Standard_Real epsilon = 1e-10)
{
    return (y - x) > epsilon;
}
bool isGreaterThanOrEqual(Standard_Real x, Standard_Real y, Standard_Real epsilon = 1e-10)
{
    return (x - y) > -epsilon;
}
bool isLessThanOrEqual(Standard_Real x, Standard_Real y, Standard_Real epsilon = 1e-10)
{
    return (y - x) > -epsilon;
}


void CurveFair::UniformCurve(Handle(Geom_BSplineCurve)& curve)
{
    TColStd_Array1OfReal curveKnots(1, curve->NbKnots());
    curve->Knots(curveKnots);
    if (!(curveKnots(curveKnots.Lower()) == 0 && curveKnots(curveKnots.Upper()) == 1))
    {
        BSplCLib::Reparametrize(0, 1, curveKnots);
        curve->SetKnots(curveKnots);
    }
}

// To compute the value of a b-spline basic function value
// write by: xin.c
static Standard_Real OneBasicFun(
    const Standard_Real u,
    const Standard_Integer i,
    const Standard_Integer p,
    const std::vector<Standard_Real>& Knots)
{
    Standard_Real Nip, uleft, uright, saved, temp;
    Standard_Integer m = Knots.size() - 1;
    std::vector<Standard_Real>N(p + 1);
    if ((i == 0 && isEqual(u, Knots[0])) || (i == m - p - 1 && isEqual(u, Knots[m])))
    {
        return 1.0;
    }

    if (isLessThan(u, Knots[i]) || isGreaterThanOrEqual(u, Knots[i + p + 1]))
    {
        return 0.0;
    }


    for (size_t j = 0; j <= p; j++)
    {
        if (isGreaterThanOrEqual(u, Knots[i + j]) && isLessThan(u, Knots[i + j + 1]))
        {
            N[j] = 1.0;
        }
        else
        {
            N[j] = 0.0;
        }
    }
    for (size_t k = 1; k <= p; k++)
    {
        if (N[0] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = ((u - Knots[i]) * N[0]) / (Knots[i + k] - Knots[i]);
        }
        for (size_t j = 0; j < p - k + 1; j++)
        {
            uleft = Knots[i + j + 1];
            uright = Knots[i + j + k + 1];
            if (N[j + 1] == 0.0)
            {
                N[j] = saved;
                saved = 0.0;
            }
            else
            {
                temp = N[j + 1] / (uright - uleft);
                N[j] = saved + (uright - u) * temp;
                saved = (u - uleft) * temp;
            }
        }
    }
    Nip = N[0];
    return Nip;
}
//std::string format_as(Eigen::MatrixXd M)
//{
//    std::stringstream ss;
//    for (Standard_Integer i = 0; i < M.rows(); i++)
//    {
//        for (Standard_Integer j = 0; j < M.cols(); j++)
//        {
//            ss << M(i, j) << " ";
//        }
//        ss << std::endl;
//    }
//    return ss.str();
//}
//std::string format_as(TColgp_Array1OfPnt Poles)
//{
//    std::stringstream ss;
//    ss << "[";
//    for (Standard_Integer i = Poles.Lower(); i <= Poles.Upper(); i++)
//    {
//        ss << "(" << Poles.Value(i).X() << ", " << Poles.Value(i).Y() << ", " << Poles.Value(i).Z() << ")";
//    }
//    ss << "]" << std::endl;
//    return ss.str();
//}
//std::string format_as(std::vector<Standard_Real> Knots)
//{
//    std::stringstream ss;
//    ss << "[";
//    for (Standard_Integer i = 0; i < Knots.size(); i++)
//    {
//        ss << Knots[i];
//        if (i != Knots.size() - 1)
//        {
//            ss << ", ";
//        }
//    }
//    ss << "]" << std::endl;
//    return ss.str();
//}
//std::string format_as(const Handle(Geom_BSplineCurve) theCurve)
//{
//    std::stringstream ss;
//    Standard_Real aFirstParameter = theCurve->FirstParameter();
//    Standard_Real aLastParameter = theCurve->LastParameter();
//    TColStd_Array1OfReal aKnotArray = theCurve->Knots();
//    ss << "FirstParameter:" << aFirstParameter << ", " << std::endl;
//    ss << "LastParameter:" << aLastParameter << std::endl;
//    ss << "[";
//    for (Standard_Integer i = 1; i <= aKnotArray.Size(); i++)
//    {
//        ss << aKnotArray.Value(i);
//        if (i == aKnotArray.Size())
//        {
//            ss << "]";
//        }
//        else
//        {
//            ss << ",";
//        }
//    }
//    return ss.str();
//}
//std::string format_as(const std::vector<Handle(Geom_BSplineCurve)> theCurveArray)
//{
//    std::stringstream ss;
//
//    for (auto curve : theCurveArray)
//    {
//        ss << format_as(curve) << std::endl;
//    }
//    return ss.str();
//}
//std::string format_as(gp_Vec Vec)
//{
//    std::stringstream ss;
//    ss << "(" << Vec.X() << ", " << Vec.Y() << ", " << Vec.Z() << ")" << " 模长 : " << Vec.Magnitude();
//    return ss.str();
//}
//std::string to_string(Standard_Real value, Standard_Integer precision)
//{
//    std::ostringstream oss;
//    oss << std::fixed << std::setprecision(precision) << value;
//    return oss.str();
//}

Handle(Geom_BSplineCurve) CurveFair::GetArclengthParameterMapping(const Handle(Geom_BSplineCurve)& theCurve, const Standard_Real theTolerance)
{
    // bspline属性
    Standard_Real aFirstParam = theCurve->FirstParameter();
    Standard_Real aLastParam = theCurve->LastParameter();

    auto aKnotArray = theCurve->Knots();
    const Standard_Integer aLowerKnot = aKnotArray.Lower();
    const Standard_Integer aUpperKnot = aKnotArray.Upper();

    // bspline: 计算弧长
    GeomAdaptor_Curve aCurveAdapter(theCurve);
    Standard_Real aBsplineLen = 0;
    Standard_Real aAvgLen = 0;
    Standard_Real aParamStep = 0;
    try
    {
        aBsplineLen = GCPnts_AbscissaPoint::Length(aCurveAdapter, aFirstParam, aLastParam, theTolerance);
        aAvgLen = aBsplineLen / (m_ArcLengthMappingSampleNum - 1.0);
        aParamStep = (aLastParam - aFirstParam) / (m_ArcLengthMappingSampleNum - 1.0);
    }
    catch (const Standard_Failure& e)
    {
        m_errorCode.push_back(e.GetMessageString());
        return nullptr;
    }
    catch (const std::exception& e)
    {
        m_errorCode.push_back(e.what());
        return nullptr;
    }
    catch (...)
    {
        m_errorCode.push_back("Unknown Exception occurred!");
        return nullptr;
    }

    // 弧长参数化: 计算新的节点向量
    TColStd_Array1OfReal reparamKnots(aLowerKnot, aLowerKnot + m_ArcLengthMappingSampleNum - 1);
    TColgp_Array1OfPnt reparamPoints(aLowerKnot, aLowerKnot + m_ArcLengthMappingSampleNum - 1);
    reparamKnots.SetValue(aLowerKnot, aFirstParam);
    reparamPoints.SetValue(aLowerKnot, gp_Pnt(aFirstParam, aFirstParam, aFirstParam));

    Standard_Real curLen = 0.0;
    Standard_Real curParam = aFirstParam;
    for (Standard_Integer i = aLowerKnot + 1; i < m_ArcLengthMappingSampleNum; i++)
    {
        curLen += aAvgLen;
        curParam += aParamStep;

        try
        {
            GCPnts_AbscissaPoint gap(aCurveAdapter, curLen, aFirstParam);
            Standard_Real tparam = gap.Parameter();
            reparamKnots.SetValue(i, curParam);
            reparamPoints.SetValue(i, gp_Pnt(tparam, tparam, tparam));
        }
        catch (const Standard_Failure& e)
        {
            m_errorCode.push_back(e.GetMessageString());
            return nullptr;
        }
        catch (const std::exception& e)
        {
            m_errorCode.push_back(e.what());
            return nullptr;
        }
        catch (...)
        {
            m_errorCode.push_back("Unknown Exception occurred!");
            return nullptr;
        }
    }

    reparamKnots.SetValue(m_ArcLengthMappingSampleNum, aLastParam);
    reparamPoints.SetValue(m_ArcLengthMappingSampleNum, gp_Pnt(aLastParam, aLastParam, aLastParam));

    return GeomAPI_PointsToBSpline
    (
        reparamPoints,
        reparamKnots,
        theCurve->Degree(),
        theCurve->Degree(),
        GeomAbs_C0,
        1e-11
    );
}

Standard_Real CurveFair::f(
    const Standard_Real theParameter,
    const Standard_Integer k)
{
    if (k == 0) return m_ArcLengthMappingFunction->Value(theParameter).X();
    gp_Vec DerivateVector = m_ArcLengthMappingFunction->DN(theParameter, k);
    return DerivateVector.X();
}

gp_Vec CurveFair::Ds(Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real sParameter, Standard_Integer k)
{
    if (theBSplineCurve.IsNull())
        return gp_Vec(0, 0, 0);

    // 弧长参数 s 转换为曲线参数 u
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve);
    Standard_Real u = f(sParameter, 0);

    // f'(s), f''(s), f'''(s)
    Standard_Real f1 = f(sParameter, 1);
    Standard_Real f2 = f(sParameter, 2);
    Standard_Real f3 = f(sParameter, 3);

    // 调用内置接口获取高阶导数
    gp_Vec D1 = theBSplineCurve->DN(u, 1);
    if (k == 1)
        return D1.Multiplied(f1);

    gp_Vec D2 = theBSplineCurve->DN(u, 2);
    if (k == 2)
        return D2.Multiplied(f1 * f1).Added(D1.Multiplied(f2));

    gp_Vec D3 = theBSplineCurve->DN(u, 3);
    if (k == 3)
        return D3.Multiplied(f1 * f1 * f1)
        .Added(D2.Multiplied(3.0 * f1 * f2))
        .Added(D1.Multiplied(f3));

    if (k == 0)
        return gp_Vec(theBSplineCurve->Value(u).X(), theBSplineCurve->Value(u).Y(), theBSplineCurve->Value(u).Z());

    return gp_Vec(0, 0, 0);
}

// CurveFair::D 函数的实现
gp_Vec CurveFair::Dt(Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real tParameter, Standard_Integer k)
{
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve);
    return Ds(theBSplineCurve, f_inverse(theBSplineCurve, tParameter), k);
}


Standard_Real CurveFair::BasisFunctionDerivative(
    const Standard_Real u,
    const Standard_Integer i,
    const Standard_Integer p,
    const Standard_Integer k,
    const std::vector<Standard_Real>& Knots)
{
    // 0阶导，直接返回基函数值
    if (k == 0)
    {
        return OneBasicFun(u, i, p, Knots);
    }

    // 当k > p时，导数为0, 当p = 0时，导数为0
    if (k > p || p == 0)
    {
        return 0.0;
    }

    // 计算第一项: (p / (u_{i+p} - u_i)) * N^{(k-1)}_{i,p-1}
    Standard_Real Term1 = 0.0;
    Standard_Real Denominator1 = Knots[i + p] - Knots[i];
    if (std::abs(Denominator1) > std::numeric_limits<Standard_Real>::epsilon())
    {
        Term1 = (p / Denominator1) * BasisFunctionDerivative(u, i, p - 1, k - 1, Knots);
    }

    // 计算第二项: (p / (u_{i+p+1} - u_{i+1})) * N^{(k-1)}_{i+1,p-1}
    Standard_Real Term2 = 0.0;
    Standard_Real Denominator2 = Knots[i + p + 1] - Knots[i + 1];
    if (std::abs(Denominator2) > std::numeric_limits<Standard_Real>::epsilon())
    {
        Term2 = (p / Denominator2) * BasisFunctionDerivative(u, i + 1, p - 1, k - 1, Knots);
    }

    // 根据公式计算k阶导数
    return Term1 - Term2;
}

Handle(Geom_BSplineCurve) CurveFair::CreateNewBSplineCurve(
    const Handle(Geom_BSplineCurve)& theOriginalCurve,
    const Eigen::MatrixXd& newD)
{
    // 获取原始曲线的 Degree
    Standard_Integer degree = theOriginalCurve->Degree();
    // 获取原始曲线的 节点向量和次数
    TColStd_Array1OfReal knots(1, theOriginalCurve->NbKnots());
    theOriginalCurve->Knots(knots);

    TColStd_Array1OfInteger multiplicities(1, theOriginalCurve->NbKnots());
    theOriginalCurve->Multiplicities(multiplicities);

    // 将 Eigen::Matrix3d 转换为 OCC 控制点数组
    TColgp_Array1OfPnt Poles(1, newD.rows()); // OCC 控制点数组
    for (Standard_Integer i = 0; i < newD.rows(); i++)
    {
        Poles.SetValue(i + 1, gp_Pnt(newD(i, 0), newD(i, 1), newD(i, 2)));
    }
    // 创建新的 B 样条曲线
    Handle(Geom_BSplineCurve) newCurve;
    try
    {
        newCurve = new Geom_BSplineCurve(
            Poles, knots, multiplicities, degree);
    }
    catch (const Standard_Failure& e)
    {
        m_errorCode.push_back(e.GetMessageString());
        return nullptr;
    }
    catch (const std::exception& e)
    {
        m_errorCode.push_back(e.what());
        return nullptr;
    }
    catch (...)
    {
        m_errorCode.push_back("Unknown Exception occurred!");
        return nullptr;
    }

    UniformCurve(newCurve);
    return newCurve;
}

Standard_Real CurveFair::GetCurveCurveHausdorffDistance(
    const Handle(Geom_BSplineCurve) theOriginalCurve,
    const Handle(Geom_BSplineCurve) theOperateCurve)
{
    Standard_Real firstParameter = theOriginalCurve->FirstParameter();
    Standard_Real lastParameter = theOriginalCurve->LastParameter();
    const Standard_Integer sampleNum = 100;
    const Standard_Real step = (lastParameter - firstParameter) / sampleNum;
    Standard_Real hausdorffDistanceResult = INT_MIN;
    for (Standard_Real t = firstParameter; t <= lastParameter; t += step)
    {
        gp_Pnt aPntOnOriginalCurve = theOriginalCurve->Value(t);
        GeomAPI_ProjectPointOnCurve aProjector(aPntOnOriginalCurve, theOperateCurve);
        // 检查是否找到投影点
        if (!aProjector.NbPoints())
        {
            continue;
        }

        Standard_Real minDistance = INT_MAX;
        for (Standard_Integer i = 1; i <= aProjector.NbPoints(); i++)
        {
            minDistance = std::min(minDistance, aProjector.Distance(i));
        }

        hausdorffDistanceResult = std::max(hausdorffDistanceResult, minDistance);
    }
    return hausdorffDistanceResult;
}

Standard_Real CurveFair::GetFitPointsCurveHausdorffDistance(
    const std::vector<std::pair<gp_Pnt, Standard_Real>> theFitPointParams,
    const Handle(Geom_BSplineCurve)& theOperateCurve)
{
    Handle(Geom_BSplineCurve) mappingCurve = GetArclengthParameterMapping(theOperateCurve);
    Standard_Real hausdorffDistance = -1.0;
    std::vector<gp_Pnt> pointsOnCurve;
    for (auto p : theFitPointParams)
    {
        gp_Pnt fitPoint = p.first;
        Standard_Real s = p.second;
        Standard_Real t = mappingCurve->Value(s).X();
        gp_Pnt pointOnCurve = theOperateCurve->Value(t);
        pointsOnCurve.push_back(pointOnCurve);
        Standard_Real distance = fitPoint.Distance(pointOnCurve);
        hausdorffDistance = std::max(hausdorffDistance, distance);
    }

    //SurfaceModelingTool::ExportPoints(pointsOnCurve, CurveFair::ExportFilePath + "/Test/pointsOnCurve.step", false);
    return hausdorffDistance;
}

Standard_Real CurveFair::GetPointCurveHausdorffDistance(
    const std::vector<gp_Pnt>& thePoints,
    const Handle(Geom_BSplineCurve)& theCurve)
{
    if (thePoints.empty() || theCurve.IsNull())
        return -1.0;

    Standard_Real hausdorffDistance = -1.0;

    for (const auto& pt : thePoints)
    {
        // 使用 GeomAPI_ProjectPointOnCurve 投影点到曲线，获取最近点
        GeomAPI_ProjectPointOnCurve projector(pt, theCurve);
        if (!projector.NbPoints()) continue;

        gp_Pnt projected = projector.NearestPoint();
        Standard_Real dist = pt.Distance(projected);

        hausdorffDistance = std::max(hausdorffDistance, dist);
    }

    return hausdorffDistance;
}


Standard_Real CurveFair::f_inverse(const Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real t)
{
    Standard_Real leftParam = theBSplineCurve->FirstParameter();
    Standard_Real rightParam = theBSplineCurve->LastParameter();
    if (t == leftParam) return leftParam;
    if (t == rightParam) return rightParam;
    while (true)
    {
        Standard_Real mid = leftParam + (rightParam - leftParam) / 2;
        Standard_Real t_mid = f(mid);
        if (std::abs(t_mid - t) < 1e-10) return mid;
        if (t_mid > t)
        {
            rightParam = mid;
        }
        else if (t_mid < t)
        {
            leftParam = mid;
        }
    }
    return 0;
}

std::vector<std::pair<gp_Pnt, Standard_Real>> CurveFair::ReCalculateFitPointParameters(
    const std::vector<gp_Pnt>& theFitPoints,
    const Handle(Geom_BSplineCurve)& theCurve
)
{
    std::vector<std::pair<gp_Pnt, Standard_Real>> result;
    if (theFitPoints.empty())
    {
        return result;
    }

    if (!theCurve.IsNull())
    {
        // 有曲线，进行投影和弧长映射
        GCPnts_AbscissaPoint abscissa;
        GeomAdaptor_Curve adaptor(theCurve);
        Standard_Real first = adaptor.FirstParameter();
        Standard_Real last = adaptor.LastParameter();
        Standard_Real totalLength = GCPnts_AbscissaPoint::Length(adaptor, first, last);

        for (const auto& pt : theFitPoints)
        {
            GeomAPI_ProjectPointOnCurve projector(pt, theCurve);
            if (!projector.NbPoints())
            {
                result.emplace_back(pt, first); // 若无法投影，使用起始参数
                continue;
            }

            Standard_Real u = projector.LowerDistanceParameter();
            Standard_Real partialLength = GCPnts_AbscissaPoint::Length(adaptor, first, u);
            Standard_Real normalized = (partialLength / totalLength) * (last - first) + first;
            result.emplace_back(pt, normalized);
        }
    }
    else
    {
        // 无曲线，使用折线方式估算
        Standard_Real totalLength = 0.0;
        std::vector<Standard_Real> cumulativeLengths(theFitPoints.size(), 0.0);
        for (size_t i = 1; i < theFitPoints.size(); ++i)
        {
            Standard_Real d = theFitPoints[i].Distance(theFitPoints[i - 1]);
            totalLength += d;
            cumulativeLengths[i] = cumulativeLengths[i - 1] + d;
        }

        for (size_t i = 0; i < theFitPoints.size(); ++i)
        {
            Standard_Real ratio = (totalLength > 0) ? cumulativeLengths[i] / totalLength : 0.0;
            result.emplace_back(theFitPoints[i], ratio);
        }
    }

    return result;
}

Standard_Real CurveFair::ComputeCurveLengthBetweenParameters(const Handle(Geom_BSplineCurve)& theCurve, Standard_Real theParameter1, Standard_Real theParameter2)
{
    GeomAdaptor_Curve adaptor(theCurve);
    try
    {
        Standard_Real length = CPnts_AbscissaPoint::Length(adaptor, theParameter1, theParameter2);
        return length;
    }
    catch (Standard_ConstructionError& e)
    {
        std::cerr << "Error while calculating length: " << e.GetMessageString() << std::endl;
        return -1.0; // 错误标识
    }
}

Standard_Real CurveFair::ComputeCurveLengthBetweenParameters(const Handle(Geom_BSplineCurve)& theCurve, gp_Pnt thePnt1, gp_Pnt thePnt2)
{
    Standard_Real theParameter1 = GetPntParameterOnCurve(theCurve, thePnt1);
    Standard_Real theParameter2 = GetPntParameterOnCurve(theCurve, thePnt2);

    if (theParameter1 > theParameter2)
    {
        std::swap(theParameter1, theParameter2);
    }
    return ComputeCurveLengthBetweenParameters(theCurve, theParameter1, theParameter2);
}

// 函数：根据 gp_Pnt 找到 B-Spline 曲线上的最近点参数
Standard_Real CurveFair::GetPntParameterOnCurve(const Handle(Geom_BSplineCurve)& theCurve, const gp_Pnt& thePoint)
{
    // 创建投影工具
    GeomAPI_ProjectPointOnCurve projector(thePoint, theCurve);

    // 检查是否找到投影点
    if (!projector.NbPoints())
    {
        std::cerr << "No projection found for the point on the curve." << std::endl;
        return -1.0; // 错误标识
    }

    // 返回第一个最近点的参数
    return projector.LowerDistanceParameter();
}
std::pair<std::vector<gp_Pnt>, std::vector <gp_Pnt>> CurveFair::SampleCurveWithArclengthMapping(const Handle(Geom_BSplineCurve)& theCurve,
    const Standard_Integer nSamples)
{
    std::vector<gp_Pnt> originPnts;
    std::vector<gp_Pnt> mappingPnts;

    // 获取变换后的曲线
    Handle(Geom_BSplineCurve) mappedCurve = CurveFair::GetArclengthParameterMapping(theCurve);
    if (mappedCurve.IsNull())
    {
        std::cerr << "Mapping failed." << std::endl;
    }

    // 原曲线参数区间
    Standard_Real uStart = theCurve->FirstParameter();
    Standard_Real uEnd = theCurve->LastParameter();
    std::vector<Standard_Real> originLengths;
    std::vector<Standard_Real> mappingLengths;

    for (int i = 0; i < nSamples; ++i)
    {
        Standard_Real u = uStart + (uEnd - uStart) * i / (nSamples - 1);

        originLengths.push_back(ComputeCurveLengthBetweenParameters(m_OriginalCurve, uStart, u));
        gp_Pnt pOriginal, pMapped;
        theCurve->D0(u, pOriginal);  // 原始点

        // 获取映射后的参数值
        Standard_Real uMapped = mappedCurve->Value(u).X();  // 假设映射为 X 坐标存储
        mappingLengths.push_back(ComputeCurveLengthBetweenParameters(m_OriginalCurve, uStart, uMapped));
        theCurve->D0(uMapped, pMapped);                     // 映射点仍在原曲线上取点

        originPnts.push_back(pOriginal);
        mappingPnts.push_back(pMapped);

    }
    return std::make_pair(originPnts, mappingPnts);
}




Handle(Geom_BSplineCurve) CurveFair::GetTempFairCurve(
    const Handle(Geom_BSplineCurve)& theCurve,
    Eigen::MatrixXd M,
    Eigen::MatrixXd V,
    Eigen::MatrixXd& D0,
    Eigen::MatrixXd& D,
    Standard_Real theAlpha)
{
    // 控制点的总数
    Standard_Integer n = D0.rows();

    // 固定控制点
    Eigen::Vector3d D0_start = D0.row(0);
    Eigen::Vector3d D0_end = D0.row(n - 1);

    // 获取需要优化的控制点原始值
    Eigen::MatrixXd D0_internal = D0.block(1, 0, n - 2, 3);
    Eigen::VectorXd D0_internal_x = D0.col(0).segment(1, n - 2);  // x分量
    Eigen::VectorXd D0_internal_y = D0.col(1).segment(1, n - 2);  // y分量
    Eigen::VectorXd D0_internal_z = D0.col(2).segment(1, n - 2);  // z分量

    // A = α * M_internal + V_internal
    // 获取需要优化的控制点对应的 M 矩阵和 V 矩阵
    Eigen::MatrixXd M_internal = M.block(1, 1, n - 2, n - 2);
    Eigen::MatrixXd V_internal = V.block(1, 1, n - 2, n - 2);

    Standard_Integer aCompouteTimes = 0;
    Handle(Geom_BSplineCurve) aResultCurve = nullptr;

    while (Standard_True)
    {
        aCompouteTimes++;
        Eigen::MatrixXd A = theAlpha * M_internal + V_internal;

        // 分别求解x、y、z分量
        Eigen::VectorXd b_x = Eigen::VectorXd::Zero(n - 2);
        Eigen::VectorXd b_y = Eigen::VectorXd::Zero(n - 2);
        Eigen::VectorXd b_z = Eigen::VectorXd::Zero(n - 2);

        // 减去首尾控制点的影响，分别计算每个分量的b
        for (Standard_Integer k = 0; k < n - 2; k++)
        {
            Standard_Real D0_term_x = M(0, k + 1) * D0_start.x();
            Standard_Real D0_term_y = M(0, k + 1) * D0_start.y();
            Standard_Real D0_term_z = M(0, k + 1) * D0_start.z();
            Standard_Real Dn_term_x = M(n - 1, k + 1) * D0_end.x();
            Standard_Real Dn_term_y = M(n - 1, k + 1) * D0_end.y();
            Standard_Real Dn_term_z = M(n - 1, k + 1) * D0_end.z();

            Standard_Real reg_term_x = V_internal(k, k) * D0_internal_x(k);
            Standard_Real reg_term_y = V_internal(k, k) * D0_internal_y(k);
            Standard_Real reg_term_z = V_internal(k, k) * D0_internal_z(k);

            b_x(k) = -theAlpha * (D0_term_x + Dn_term_x) + reg_term_x;
            b_y(k) = -theAlpha * (D0_term_y + Dn_term_y) + reg_term_y;
            b_z(k) = -theAlpha * (D0_term_z + Dn_term_z) + reg_term_z;
        }

        Eigen::VectorXd D_internal_x;
        Eigen::VectorXd D_internal_y;
        Eigen::VectorXd D_internal_z;
        try
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
            Eigen::MatrixXd A_pinv = svd.solve(Eigen::MatrixXd::Identity(A.rows(), A.cols()));
            D_internal_x = A_pinv * b_x;
            D_internal_y = A_pinv * b_y;
            D_internal_z = A_pinv * b_z;
        }
        catch (const Standard_Failure& e)
        {
            m_errorCode.push_back(e.GetMessageString());
            return nullptr;
        }
        catch (const std::exception& e)
        {
            m_errorCode.push_back(e.what());
            return nullptr;
        }
        catch (...)
        {
            m_errorCode.push_back("Unknown Exception occurred!");
            return nullptr;
        }

        // 合并结果
        Eigen::MatrixXd D_internal(n - 2, 3);
        D_internal.col(0) = D_internal_x;
        D_internal.col(1) = D_internal_y;
        D_internal.col(2) = D_internal_z;

        // 构建最终控制点矩阵
        D = Eigen::MatrixXd::Zero(n, 3);
        D.row(0) = D0_start;                    // 起始点
        D.block(1, 0, n - 2, 3) = D_internal;   // 内部点
        D.row(n - 1) = D0_end;                  // 终止点

        // 创建新的B样条曲线
        aResultCurve = CreateNewBSplineCurve(m_OriginalCurve, D);

        // 计算新旧曲线的控制点偏差
        Standard_Real aHausdorffDistance = 0;
        if (m_OriginalFitPoints.size() > 0)
        {
            aHausdorffDistance = CurveFair::GetFitPointsCurveHausdorffDistance(m_FitPointParameters, aResultCurve);
            m_FitPointParameters = ReCalculateFitPointParameters(m_OriginalFitPoints, m_OriginalCurve);
        }
        else
        {
            aHausdorffDistance = CurveFair::GetCurveCurveHausdorffDistance(m_OriginalCurve, aResultCurve);
        }

        m_HausdorffDistanceResult = aHausdorffDistance;
        if (aHausdorffDistance <= m_HausdorffDistanceTol || aCompouteTimes >= 100)
        {
            break;
        }
        else
        {
            theAlpha /= 2;
        }

    }
    return aResultCurve;
}

Eigen::MatrixXd CurveFair::ComputeContinuityMatrix(
    const Handle(Geom_BSplineCurve)& theBSplineCurve,
    const Standard_Integer p)
{
    const Standard_Integer maxDerivate = 3; // C³连续
    const Standard_Integer aNum = theBSplineCurve->NbPoles();
    TColStd_Array1OfReal aKnots = theBSplineCurve->Knots();
    Standard_Integer aDeg = theBSplineCurve->Degree();
    TColStd_Array1OfReal aKnotSeq = theBSplineCurve->KnotSequence();
    std::vector<Standard_Real> KnotSeqVector = OccArrayConvertoVector(aKnotSeq);
    std::vector<Standard_Real> KnotVector = OccArrayConvertoVector(aKnots);
    // 获取内部节点（排除首尾以及重复节点）
    std::vector<Standard_Real> internalKnots(KnotVector.begin() + 1, KnotVector.end() - 1);

    Standard_Integer m = internalKnots.size(); // 约束条件数量
    Eigen::MatrixXd C(m, aNum); // 约束条件数量 , 控制点数量
    C.setZero();

    for (Standard_Integer j = 0; j < m; ++j)
    {
        Standard_Real u_j = internalKnots[j];
        Standard_Real s_j = f_inverse(theBSplineCurve, u_j);
        Standard_Real Left_uj = u_j - Precision::PConfusion();
        Standard_Real Right_uj = u_j + Precision::PConfusion();
        // 计算 f(s) = u_j 及其导数
        Standard_Real f0 = f(s_j, 0);   // f(s)
        Standard_Real f1 = f(s_j, 1);   // f'(s)
        Standard_Real f2 = f(s_j, 2);   // f''(s)
        Standard_Real f3 = f(s_j, 3);   // f'''(s)
        // 计算三阶导数的左右极限
        math_Matrix leftBasis(1, maxDerivate + 1, 1, aDeg + 1);
        math_Matrix rightBasis(1, maxDerivate + 1, 1, aDeg + 1);
        Standard_Integer leftFirstIndex, rightFirstIndex;


        // 左极限
        BSplCLib::EvalBsplineBasis(maxDerivate, aDeg + 1, aKnotSeq, Left_uj, leftFirstIndex, leftBasis);

        // 右极限
        BSplCLib::EvalBsplineBasis(maxDerivate, aDeg + 1, aKnotSeq, Right_uj, rightFirstIndex, rightBasis);

        // 计算跳跃值
        for (Standard_Integer localIdx = 0; localIdx < aDeg + 1; ++localIdx)
        {
            Standard_Integer leftGlobalIdx = leftFirstIndex + localIdx - 1;
            Standard_Integer rightGlobalIdx = rightFirstIndex + localIdx - 1;

            // 确保索引有效
            if (leftGlobalIdx >= 0 && leftGlobalIdx < aNum)
            {
                //Standard_Real d1Ni_df1 = BasisFunctionDerivative(Left_uj, leftGlobalIdx, p, 1, KnotSeqVector);
                //Standard_Real d2Ni_df2 = BasisFunctionDerivative(Left_uj, leftGlobalIdx, p, 2, KnotSeqVector);
                //Standard_Real d3Ni_df3 = BasisFunctionDerivative(Left_uj, leftGlobalIdx, p, 3, KnotSeqVector);
                Standard_Real d1Ni_df1 = leftBasis(2, 1 + localIdx);
                Standard_Real d2Ni_df2 = leftBasis(3, 1 + localIdx);
                Standard_Real d3Ni_df3 = leftBasis(4, 1 + localIdx);

                // 计算弧长参数化链式求导的三阶导数
                Standard_Real D1Ni_Dt1 = d1Ni_df1 * f1;
                Standard_Real D2Ni_Dt2 = d2Ni_df2 * pow(f1, 2) + d1Ni_df1 * f2;
                Standard_Real D3Ni_Left = d3Ni_df3 * pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;
                C(j, leftGlobalIdx) += D3Ni_Left;
            }
            if (rightGlobalIdx >= 0 && rightGlobalIdx < aNum)
            {
                //Standard_Real d1Ni_df1 = BasisFunctionDerivative(Right_uj, rightGlobalIdx, p, 1, KnotSeqVector);
                //Standard_Real d2Ni_df2 = BasisFunctionDerivative(Right_uj, rightGlobalIdx, p, 2, KnotSeqVector);
                //Standard_Real d3Ni_df3 = BasisFunctionDerivative(Right_uj, rightGlobalIdx, p, 3, KnotSeqVector);

                Standard_Real d1Ni_df1 = rightBasis(2, 1 + localIdx);
                Standard_Real d2Ni_df2 = rightBasis(3, 1 + localIdx);
                Standard_Real d3Ni_df3 = rightBasis(4, 1 + localIdx);


                // 计算弧长参数化链式求导的三阶导数
                Standard_Real D1Ni_Dt1 = d1Ni_df1 * f1;
                Standard_Real D2Ni_Dt2 = d2Ni_df2 * pow(f1, 2) + d1Ni_df1 * f2;
                Standard_Real D3Ni_Right = d3Ni_df3 * pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;
                C(j, rightGlobalIdx) -= D3Ni_Right;
            }
        }
    }
    return C;
}

// 非并行处理
Eigen::MatrixXd CurveFair::ComputeEnergyMatrix(
    const Handle(Geom_BSplineCurve)& theBSplineCurve,
    const Standard_Integer p,
    const Standard_Real tol)
{
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve); // 计算f(t)
    const Standard_Integer maxDerivate = 3;
    const Standard_Integer aGaussNum = 30;
    Standard_Integer aNum = theBSplineCurve->NbPoles();
    TColStd_Array1OfReal aKnots = theBSplineCurve->Knots();
    Standard_Integer aDeg = theBSplineCurve->Degree();
    TColStd_Array1OfReal aKnotSeq = theBSplineCurve->KnotSequence();
    math_Vector aGaussPnts(1, aGaussNum);
    math_Vector aGausssWeights(1, aGaussNum);
    math::GaussPoints(aGaussNum, aGaussPnts);
    math::GaussWeights(aGaussNum, aGausssWeights);
    std::vector<Standard_Real> aKnotVector = OccArrayConvertoVector(aKnots);
    std::vector<Standard_Real> aKnotSeqVector = OccArrayConvertoVector(aKnotSeq);
    Eigen::MatrixXd M(aNum, aNum);
    M.setZero();

    for (Standard_Integer i = aKnots.Lower(); i < aKnots.Upper(); ++i)
    {
        Standard_Real sStart = f_inverse(theBSplineCurve, aKnots[i]);
        Standard_Real sEnd = f_inverse(theBSplineCurve, aKnots[i + 1]);
        //将各段的两个高斯积分点相加
        for (Standard_Integer GaussIndex = 0; GaussIndex < aGaussNum; ++GaussIndex)
        {
            Standard_Real s = (sEnd - sStart) * aGaussPnts(aGaussNum - GaussIndex) / 2.0 + (sStart + sEnd) / 2.0;
            Standard_Integer aFirstIndex;
            math_Matrix aBsplineBasis(1, maxDerivate + 1, 1, aDeg + 1);
            BSplCLib::EvalBsplineBasis(maxDerivate, aDeg + 1, aKnotSeq, f(s), aFirstIndex, aBsplineBasis);

            // 计算f(t)及其导数
            Standard_Real f0 = f(s, 0);   // f(t)
            Standard_Real f1 = f(s, 1);   // f'(t)
            Standard_Real f2 = f(s, 2);   // f''(t)
            Standard_Real f3 = f(s, 3);   // f'''(t)
            for (Standard_Integer m = 0; m < aDeg + 1; ++m)
            {
                Standard_Integer globalI = m + aFirstIndex - 1;
                // 弧长参数化 k 阶导数 参数 f(t)
                //Standard_Real d1Ni_df1 = BasisFunctionDerivative(f0, globalI, p, 1, aKnotSeqVector);
                //Standard_Real d2Ni_df2 = BasisFunctionDerivative(f0, globalI, p, 2, aKnotSeqVector);
                //Standard_Real d3Ni_df3 = BasisFunctionDerivative(f0, globalI, p, 3, aKnotSeqVector);
                Standard_Real d1Ni_df1 = aBsplineBasis(2, 1 + m);
                Standard_Real d2Ni_df2 = aBsplineBasis(3, 1 + m);
                Standard_Real d3Ni_df3 = aBsplineBasis(4, 1 + m);

                // 计算弧长参数化链式求导的三阶导数
                Standard_Real D1Ni_Dt1 = d1Ni_df1 * f1;
                Standard_Real D2Ni_Dt2 = d2Ni_df2 * pow(f1, 2) + d1Ni_df1 * f2;
                Standard_Real D3Ni_Dt3 = d3Ni_df3 * pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;

                for (Standard_Integer n = 0; n < aDeg + 1; ++n)
                {
                    Standard_Integer globalJ = n + aFirstIndex - 1;
                    // 弧长参数化 k 阶导数 参数 f(t)
                    //Standard_Real d1Nj_df1 = BasisFunctionDerivative(f0, globalJ, p, 1, aKnotSeqVector);
                    //Standard_Real d2Nj_df2 = BasisFunctionDerivative(f0, globalJ, p, 2, aKnotSeqVector);
                    //Standard_Real d3Nj_df3 = BasisFunctionDerivative(f0, globalJ, p, 3, aKnotSeqVector);
                    Standard_Real d1Nj_df1 = aBsplineBasis(2, 1 + n);
                    Standard_Real d2Nj_df2 = aBsplineBasis(3, 1 + n);
                    Standard_Real d3Nj_df3 = aBsplineBasis(4, 1 + n);
                    // 计算弧长参数化链式求导的三阶导数
                    Standard_Real D1Nj_Dt1 = d1Nj_df1 * f1;
                    Standard_Real D2Nj_Dt2 = d2Nj_df2 * pow(f1, 2) + d1Nj_df1 * f2;
                    Standard_Real D3Nj_Dt3 = d3Nj_df3 * pow(f1, 3) + 3 * d2Nj_df2 * f1 * f2 + d1Nj_df1 * f3;

                    // 弧长参数化三阶导数
                    Standard_Real aElement =
                        D3Ni_Dt3 *
                        D3Nj_Dt3 *
                        aGausssWeights(aGaussNum - GaussIndex) * (sEnd - sStart) / 2.0;
                    M(globalI, globalJ) += aElement;

                }
            }
        }
    }
    return M;
}

/*
Eigen::MatrixXd CurveFair::ComputeEnergyMatrix(
    const Handle(Geom_BSplineCurve)& theBSplineCurve,
    const Standard_Integer p,
    const Standard_Real theTolerance)
{
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve);
    const Standard_Integer maxDerivate = 3;
    const Standard_Integer aGaussNum = 30;
    Standard_Integer aNum = theBSplineCurve->NbPoles();
    TColStd_Array1OfReal aKnots = theBSplineCurve->Knots();
    Standard_Integer aDeg = theBSplineCurve->Degree();
    TColStd_Array1OfReal aKnotSeq = theBSplineCurve->KnotSequence();

    math_Vector aGaussPnts(1, aGaussNum);
    math_Vector aGausssWeights(1, aGaussNum);
    math::GaussPoints(aGaussNum, aGaussPnts);
    math::GaussWeights(aGaussNum, aGausssWeights);

    std::vector<Standard_Real> aKnotSeqVector = OccArrayConvertoVector(aKnotSeq);

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(aNum, aNum);

    // 并行处理每一段
#pragma omp parallel for schedule(dynamic)
    for (Standard_Integer i = aKnots.Lower(); i < aKnots.Upper(); ++i)
    {
        Eigen::MatrixXd M_local = Eigen::MatrixXd::Zero(aNum, aNum); // 每线程局部累加

        Standard_Real sStart = f_inverse(theBSplineCurve, aKnots[i]);
        Standard_Real sEnd = f_inverse(theBSplineCurve, aKnots[i + 1]);

        for (Standard_Integer GaussIndex = 0; GaussIndex < aGaussNum; ++GaussIndex)
        {
            Standard_Real s = (sEnd - sStart) * aGaussPnts(aGaussNum - GaussIndex) / 2.0
                + (sStart + sEnd) / 2.0;

            Standard_Integer aFirstIndex;
            math_Matrix aBsplineBasis(1, maxDerivate + 1, 1, aDeg + 1);
            BSplCLib::EvalBsplineBasis(maxDerivate, aDeg + 1, aKnotSeq, f(s), aFirstIndex, aBsplineBasis);

            Standard_Real f0 = f(s, 0);
            Standard_Real f1 = f(s, 1);
            Standard_Real f2 = f(s, 2);
            Standard_Real f3 = f(s, 3);

            for (Standard_Integer m = 0; m < aDeg + 1; ++m)
            {
                Standard_Integer globalI = m + aFirstIndex - 1;

                Standard_Real d1Ni_df1 = aBsplineBasis(2, 1 + m);
                Standard_Real d2Ni_df2 = aBsplineBasis(3, 1 + m);
                Standard_Real d3Ni_df3 = aBsplineBasis(4, 1 + m);

                Standard_Real D1Ni_Dt1 = d1Ni_df1 * f1;
                Standard_Real D2Ni_Dt2 = d2Ni_df2 * pow(f1, 2) + d1Ni_df1 * f2;
                Standard_Real D3Ni_Dt3 = d3Ni_df3 * pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;

                for (Standard_Integer n = 0; n < aDeg + 1; ++n)
                {
                    Standard_Integer globalJ = n + aFirstIndex - 1;

                    Standard_Real d1Nj_df1 = aBsplineBasis(2, 1 + n);
                    Standard_Real d2Nj_df2 = aBsplineBasis(3, 1 + n);
                    Standard_Real d3Nj_df3 = aBsplineBasis(4, 1 + n);

                    Standard_Real D1Nj_Dt1 = d1Nj_df1 * f1;
                    Standard_Real D2Nj_Dt2 = d2Nj_df2 * pow(f1, 2) + d1Nj_df1 * f2;
                    Standard_Real D3Nj_Dt3 = d3Nj_df3 * pow(f1, 3) + 3 * d2Nj_df2 * f1 * f2 + d1Nj_df1 * f3;

                    Standard_Real aElement = D3Ni_Dt3 * D3Nj_Dt3 * aGausssWeights(aGaussNum - GaussIndex) * (sEnd - sStart) / 2.0;
                    M_local(globalI, globalJ) += aElement;
                }
            }
        }

        // 合并线程局部矩阵
#pragma omp critical
        M += M_local;
    }

    return M;
}
*/

Eigen::SparseMatrix<Standard_Real> CurveFair::ComputeSparseEnergyMatrix(
    const Handle(Geom_BSplineCurve)& theBSplineCurve,
    const Standard_Integer p,
    const Standard_Real theTolerance)
{
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve);
    const Standard_Integer maxDerivate = 3;
    const Standard_Integer aGaussNum = 30;
    Standard_Integer aNum = theBSplineCurve->NbPoles();
    TColStd_Array1OfReal aKnots = theBSplineCurve->Knots();
    Standard_Integer aDeg = theBSplineCurve->Degree();
    TColStd_Array1OfReal aKnotSeq = theBSplineCurve->KnotSequence();

    math_Vector aGaussPnts(1, aGaussNum);
    math_Vector aGausssWeights(1, aGaussNum);
    math::GaussPoints(aGaussNum, aGaussPnts);
    math::GaussWeights(aGaussNum, aGausssWeights);

    std::vector<Standard_Real> aKnotSeqVector = OccArrayConvertoVector(aKnotSeq);

    // 最终收集所有非零项
    std::vector<Eigen::Triplet<Standard_Real>> triplets;

    {
        std::vector<Eigen::Triplet<Standard_Real>> local_triplets;


        for (Standard_Integer i = aKnots.Lower(); i < aKnots.Upper(); ++i)
        {
            Standard_Real sStart = f_inverse(theBSplineCurve, aKnots[i]);
            Standard_Real sEnd = f_inverse(theBSplineCurve, aKnots[i + 1]);

            for (Standard_Integer GaussIndex = 0; GaussIndex < aGaussNum; ++GaussIndex)
            {
                Standard_Real s = (sEnd - sStart) * aGaussPnts(aGaussNum - GaussIndex) / 2.0
                    + (sStart + sEnd) / 2.0;

                Standard_Integer aFirstIndex;
                math_Matrix aBsplineBasis(1, maxDerivate + 1, 1, aDeg + 1);
                BSplCLib::EvalBsplineBasis(maxDerivate, aDeg + 1, aKnotSeq, f(s), aFirstIndex, aBsplineBasis);

                Standard_Real f1 = f(s, 1);
                Standard_Real f2 = f(s, 2);
                Standard_Real f3 = f(s, 3);

                for (Standard_Integer m = 0; m < aDeg + 1; ++m)
                {
                    Standard_Integer globalI = m + aFirstIndex - 1;
                    Standard_Real d1Ni_df1 = aBsplineBasis(2, 1 + m);
                    Standard_Real d2Ni_df2 = aBsplineBasis(3, 1 + m);
                    Standard_Real d3Ni_df3 = aBsplineBasis(4, 1 + m);

                    Standard_Real D3Ni_Dt3 = d3Ni_df3 * pow(f1, 3)
                        + 3 * d2Ni_df2 * f1 * f2
                        + d1Ni_df1 * f3;

                    for (Standard_Integer n = 0; n < aDeg + 1; ++n)
                    {
                        Standard_Integer globalJ = n + aFirstIndex - 1;
                        Standard_Real d1Nj_df1 = aBsplineBasis(2, 1 + n);
                        Standard_Real d2Nj_df2 = aBsplineBasis(3, 1 + n);
                        Standard_Real d3Nj_df3 = aBsplineBasis(4, 1 + n);

                        Standard_Real D3Nj_Dt3 = d3Nj_df3 * pow(f1, 3)
                            + 3 * d2Nj_df2 * f1 * f2
                            + d1Nj_df1 * f3;

                        Standard_Real aElement = D3Ni_Dt3 * D3Nj_Dt3
                            * aGausssWeights(aGaussNum - GaussIndex)
                            * (sEnd - sStart) / 2.0;

                        if (std::abs(aElement) > theTolerance)
                            local_triplets.emplace_back(globalI, globalJ, aElement);
                    }
                }
            }
        }


        triplets.insert(triplets.end(), local_triplets.begin(), local_triplets.end());
    }

    Eigen::SparseMatrix<Standard_Real> M(aNum, aNum);
    M.setFromTriplets(triplets.begin(), triplets.end());

    return M;
}

Standard_Real CurveFair::ComputeCurveFairEnergy(const Handle(Geom_BSplineCurve)& theCurve)
{
    const TColgp_Array1OfPnt& Poles = theCurve->Poles();
    const Standard_Integer aNumPoles = theCurve->NbPoles();

    Eigen::MatrixXd M = ComputeEnergyMatrix(theCurve, 3); // 返回 M = CᵀC
    std::cout << M << std::endl;
    Eigen::MatrixXd P(aNumPoles, 3);
    for (Standard_Integer i = Poles.Lower(); i <= Poles.Upper(); ++i)
    {
        const gp_Pnt& pt = Poles(i);
        P(i - Poles.Lower(), 0) = pt.X();
        P(i - Poles.Lower(), 1) = pt.Y();
        P(i - Poles.Lower(), 2) = pt.Z();
    }

    // 能量 = trace(P^T * M * P)
    return (P.transpose() * M * P).trace();
}


Eigen::MatrixXd CurveFair::ComputeConstraintMatrix(
    const Eigen::MatrixXd& theD0,
    Eigen::VectorXd& theH)
{
    Standard_Integer n = theD0.rows();
    Standard_Integer num_unknowns = n - 2;

    // C 矩阵的尺寸为 4 x (3 * (n-2))，因为有4个约束方程，
    // 未知量为 n-2 个点的 x,y,z 共 3*(n-2) 个变量。
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(4, 3 * num_unknowns);
    theH = Eigen::VectorXd::Zero(4);

    if (num_unknowns <= 0)
    {
        return C;
    }

    // 提取首尾控制点
    Eigen::Vector3d P0 = theD0.row(0);
    Eigen::Vector3d P1 = theD0.row(1);
    Eigen::Vector3d Pn_1 = theD0.row(n - 2);
    Eigen::Vector3d Pn = theD0.row(n - 1);

    // --- 起始点切向约束 (P1 - P0) · u = 0 ---
    Eigen::Vector3d T_start_0 = P1 - P0;
    if (T_start_0.norm() > 1e-9) // 避免切向为零向量
    {
        // 寻找与 T_start_0 正交的两个单位向量 u1, u2
        Eigen::Vector3d v_arb = (abs(T_start_0.x()) < 0.9) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d u1 = T_start_0.cross(v_arb).normalized();
        Eigen::Vector3d u2 = T_start_0.cross(u1).normalized();

        // 填充 C 矩阵的第一行和第二行 (对应 P1)
        // P1 是未知点的第 0 个
        C(0, 0) = u1.x(); C(0, num_unknowns) = u1.y(); C(0, 2 * num_unknowns) = u1.z();
        C(1, 0) = u2.x(); C(1, num_unknowns) = u2.y(); C(1, 2 * num_unknowns) = u2.z();

        // 填充 h 向量的前两个元素
        theH(0) = P0.dot(u1);
        theH(1) = P0.dot(u2);
    }

    // --- 终止点切向约束 (Pn - Pn-1) · v = 0 ---
    Eigen::Vector3d T_end_0 = Pn - Pn_1;
    if (T_end_0.norm() > 1e-9) // 避免切向为零向量
    {
        // 寻找与 T_end_0 正交的两个单位向量 v1, v2
        Eigen::Vector3d v_arb = (abs(T_end_0.x()) < 0.9) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d v1 = T_end_0.cross(v_arb).normalized();
        Eigen::Vector3d v2 = T_end_0.cross(v1).normalized();

        // 填充 C 矩阵的第三行和第四行 (对应 Pn-1)
        // Pn-1 是未知点的第 (num_unknowns - 1) 个
        Standard_Integer last_idx = num_unknowns - 1;
        C(2, last_idx) = -v1.x(); C(2, last_idx + num_unknowns) = -v1.y(); C(2, last_idx + 2 * num_unknowns) = -v1.z();
        C(3, last_idx) = -v2.x(); C(3, last_idx + num_unknowns) = -v2.y(); C(3, last_idx + 2 * num_unknowns) = -v2.z();

        // 填充 h 向量的后两个元素
        theH(2) = -Pn.dot(v1);
        theH(3) = -Pn.dot(v2);
    }

    return C;
}

Eigen::SparseMatrix<Standard_Real> CurveFair::ComputeSparseConstraintMatrix(
    const Eigen::MatrixXd& theD0,
    Eigen::VectorXd& theH)
{
    Standard_Integer n = theD0.rows();
    Standard_Integer num_unknowns = n - 2;

    Eigen::SparseMatrix<Standard_Real> C(4, 3 * num_unknowns);
    theH = Eigen::VectorXd::Zero(4);

    if (num_unknowns <= 0)
    {
        return C;
    }

    std::vector<Eigen::Triplet<Standard_Real>> triplets;
    triplets.reserve(12); // 每行最多 3 个非零，总共 4 行

    // 提取首尾控制点
    Eigen::Vector3d P0 = theD0.row(0);
    Eigen::Vector3d P1 = theD0.row(1);
    Eigen::Vector3d Pn_1 = theD0.row(n - 2);
    Eigen::Vector3d Pn = theD0.row(n - 1);

    // 起始点切向约束
    Eigen::Vector3d T_start_0 = P1 - P0;
    if (T_start_0.norm() > 1e-9)
    {
        Eigen::Vector3d v_arb = (std::abs(T_start_0.x()) < 0.9)
            ? Eigen::Vector3d(1, 0, 0)
            : Eigen::Vector3d(0, 1, 0);

        Eigen::Vector3d u1 = T_start_0.cross(v_arb).normalized();
        Eigen::Vector3d u2 = T_start_0.cross(u1).normalized();

        // 行 0
        triplets.emplace_back(0, 0, u1.x());
        triplets.emplace_back(0, num_unknowns, u1.y());
        triplets.emplace_back(0, 2 * num_unknowns, u1.z());

        // 行 1
        triplets.emplace_back(1, 0, u2.x());
        triplets.emplace_back(1, num_unknowns, u2.y());
        triplets.emplace_back(1, 2 * num_unknowns, u2.z());

        theH(0) = P0.dot(u1);
        theH(1) = P0.dot(u2);
    }

    // 终止点切向约束
    Eigen::Vector3d T_end_0 = Pn - Pn_1;
    if (T_end_0.norm() > 1e-9)
    {
        Eigen::Vector3d v_arb = (std::abs(T_end_0.x()) < 0.9)
            ? Eigen::Vector3d(1, 0, 0)
            : Eigen::Vector3d(0, 1, 0);

        Eigen::Vector3d v1 = T_end_0.cross(v_arb).normalized();
        Eigen::Vector3d v2 = T_end_0.cross(v1).normalized();

        Standard_Integer last_idx = num_unknowns - 1;

        // 行 2
        triplets.emplace_back(2, last_idx, -v1.x());
        triplets.emplace_back(2, last_idx + num_unknowns, -v1.y());
        triplets.emplace_back(2, last_idx + 2 * num_unknowns, -v1.z());

        // 行 3
        triplets.emplace_back(3, last_idx, -v2.x());
        triplets.emplace_back(3, last_idx + num_unknowns, -v2.y());
        triplets.emplace_back(3, last_idx + 2 * num_unknowns, -v2.z());

        theH(2) = -Pn.dot(v1);
        theH(3) = -Pn.dot(v2);
    }

    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}

//// 普通矩阵求解
//Handle(Geom_BSplineCurve) CurveFair::GetTempFairCurveWithTangentConstraint(
//    const Handle(Geom_BSplineCurve)& theCurve,
//    Eigen::MatrixXd M,
//    Eigen::MatrixXd V,
//    Eigen::MatrixXd& D0,
//    Eigen::MatrixXd& D,
//    Standard_Real theAlpha)
//{
//    // 控制点的总数
//    Standard_Integer n = D0.rows();
//    if (n <= 2) return Handle(Geom_BSplineCurve)(); // 无法优化
//
//    Standard_Integer num_unknowns = n - 2;
//
//    // 固定控制点
//    Eigen::Vector3d D0_start = D0.row(0);
//    Eigen::Vector3d D0_end = D0.row(n - 1);
//
//    // 获取需要优化的控制点原始值 (x, y, z 分量)
//    Eigen::VectorXd D0_internal_x = D0.col(0).segment(1, num_unknowns);
//    Eigen::VectorXd D0_internal_y = D0.col(1).segment(1, num_unknowns);
//    Eigen::VectorXd D0_internal_z = D0.col(2).segment(1, num_unknowns);
//
//    // 获取 M 和 V 的内部块
//    Eigen::MatrixXd M_internal = M.block(1, 1, num_unknowns, num_unknowns);
//    Eigen::MatrixXd V_internal = V.block(1, 1, num_unknowns, num_unknowns);
//
//    Standard_Real average = M_internal.sum() / (M_internal.rows() * M_internal.cols());
//    //std::cout << M_internal << std::endl;
//    //std::cout << V_internal << std::endl;
//    // --- 新增：调用函数计算约束矩阵 C 和右端项 h ---
//    Eigen::VectorXd h;
//    Eigen::MatrixXd C_mat = ComputeConstraintMatrix(D0, h);
//
//
//    Standard_Integer aCompouteTimes = 0;
//    Handle(Geom_BSplineCurve) aResultCurve = nullptr;
//
//    while (Standard_True)
//    {
//        aCompouteTimes++;
//
//        // --- 构建大型线性系统 [BigA  C^T] [ X ] = [ BigB ] ---
//        //                       [ C     0  ] [ λ ]   [   h  ]
//
//        // 1. 构建系数矩阵 A (对x,y,z分量相同)
//        Eigen::MatrixXd A = theAlpha * M_internal + V_internal;
//
//        // 2. 构建右端项 b (分x,y,z)
//        Eigen::VectorXd b_x = Eigen::VectorXd::Zero(num_unknowns);
//        Eigen::VectorXd b_y = Eigen::VectorXd::Zero(num_unknowns);
//        Eigen::VectorXd b_z = Eigen::VectorXd::Zero(num_unknowns);
//        for (Standard_Integer k = 0; k < num_unknowns; k++)
//        {
//            Standard_Real D0_term_x = M(0, k + 1) * D0_start.x();
//            Standard_Real Dn_term_x = M(n - 1, k + 1) * D0_end.x();
//            Standard_Real reg_term_x = V_internal(k, k) * D0_internal_x(k);
//            b_x(k) = -theAlpha * (D0_term_x + Dn_term_x) + reg_term_x;
//
//            Standard_Real D0_term_y = M(0, k + 1) * D0_start.y();
//            Standard_Real Dn_term_y = M(n - 1, k + 1) * D0_end.y();
//            Standard_Real reg_term_y = V_internal(k, k) * D0_internal_y(k);
//            b_y(k) = -theAlpha * (D0_term_y + Dn_term_y) + reg_term_y;
//
//            Standard_Real D0_term_z = M(0, k + 1) * D0_start.z();
//            Standard_Real Dn_term_z = M(n - 1, k + 1) * D0_end.z();
//            Standard_Real reg_term_z = V_internal(k, k) * D0_internal_z(k);
//            b_z(k) = -theAlpha * (D0_term_z + Dn_term_z) + reg_term_z;
//        }
//
//        // 3. 组装大型系统的左端矩阵 (LHS) 和右端向量 (RHS)
//        Standard_Integer total_vars = 3 * num_unknowns;
//        Standard_Integer system_size = total_vars + 4; // 4个拉格朗日乘子
//        Eigen::MatrixXd BigLHS = Eigen::MatrixXd::Zero(system_size, system_size);
//        Eigen::VectorXd BigRHS = Eigen::VectorXd::Zero(system_size);
//
//        // 填充对角块 A
//        BigLHS.block(0, 0, num_unknowns, num_unknowns) = A;
//        BigLHS.block(num_unknowns, num_unknowns, num_unknowns, num_unknowns) = A;
//        BigLHS.block(2 * num_unknowns, 2 * num_unknowns, num_unknowns, num_unknowns) = A;
//
//        // 填充约束块 C 和 C^T
//        BigLHS.block(total_vars, 0, 4, total_vars) = C_mat;
//        BigLHS.block(0, total_vars, total_vars, 4) = C_mat.transpose();
//
//        // 填充右端项
//        BigRHS.segment(0, num_unknowns) = b_x;
//        BigRHS.segment(num_unknowns, num_unknowns) = b_y;
//        BigRHS.segment(2 * num_unknowns, num_unknowns) = b_z;
//        BigRHS.segment(total_vars, 4) = h;
//
//        Eigen::VectorXd D_internal_x, D_internal_y, D_internal_z;
//
//        try
//        {
//            // 4. 求解大型线性系统
//            Eigen::JacobiSVD<Eigen::MatrixXd> svd(BigLHS, Eigen::ComputeThinU | Eigen::ComputeThinV);
//            Eigen::VectorXd solution = svd.solve(BigRHS);
//
//            // 5. 从解中提取各分量
//            D_internal_x = solution.segment(0, num_unknowns);
//            D_internal_y = solution.segment(num_unknowns, num_unknowns);
//            D_internal_z = solution.segment(2 * num_unknowns, num_unknowns);
//        }
//        catch (...)
//        {
//            m_errorCode.push_back("Exception occurred during solving constrained system!");
//            return nullptr;
//        }
//
//        // 合并结果
//        Eigen::MatrixXd D_internal(num_unknowns, 3);
//        D_internal.col(0) = D_internal_x;
//        D_internal.col(1) = D_internal_y;
//        D_internal.col(2) = D_internal_z;
//
//        // 构建最终控制点矩阵
//        D = Eigen::MatrixXd::Zero(n, 3);
//        D.row(0) = D0_start;
//        D.block(1, 0, num_unknowns, 3) = D_internal;
//        D.row(n - 1) = D0_end;
//
//        // 创建新的B样条曲线
//        aResultCurve = CreateNewBSplineCurve(m_OriginalCurve, D);
//
//        // 计算 Hausdorff 距离
//        Standard_Real aHausdorffDistance = 0;
//        if (m_OriginalFitPoints.size() > 0)
//        {
//            aHausdorffDistance = CurveFair::GetFitPointsCurveHausdorffDistance(m_FitPointParameters, aResultCurve);
//            m_FitPointParameters = ReCalculateFitPointParameters(m_OriginalFitPoints, m_OriginalCurve);
//        }
//        else
//        {
//            aHausdorffDistance = CurveFair::GetCurveCurveHausdorffDistance(m_OriginalCurve, aResultCurve);
//        }
//
//        m_HausdorffDistanceResult = aHausdorffDistance;
//        if (aHausdorffDistance <= m_HausdorffDistanceTol || aCompouteTimes >= 100)
//        {
//            break;
//        }
//        else
//        {
//            theAlpha /= 2;
//        }
//    }
//    return aResultCurve;
//}

// 稀疏矩阵求解
Handle(Geom_BSplineCurve) CurveFair::GetTempFairCurveWithTangentConstraint(
    const Handle(Geom_BSplineCurve)& theCurve,
    Eigen::MatrixXd M,
    Eigen::MatrixXd V,
    Eigen::MatrixXd& D0,
    Eigen::MatrixXd& D,
    Standard_Real theAlpha)
{
    Standard_Integer n = D0.rows();
    if (n <= 2) return Handle(Geom_BSplineCurve)();

    Standard_Integer num_unknowns = n - 2;
    Eigen::Vector3d D0_start = D0.row(0);
    Eigen::Vector3d D0_end = D0.row(n - 1);

    Eigen::VectorXd D0_internal_x = D0.col(0).segment(1, num_unknowns);
    Eigen::VectorXd D0_internal_y = D0.col(1).segment(1, num_unknowns);
    Eigen::VectorXd D0_internal_z = D0.col(2).segment(1, num_unknowns);

    Eigen::MatrixXd M_internal = M.block(1, 1, num_unknowns, num_unknowns);
    Eigen::MatrixXd V_internal = V.block(1, 1, num_unknowns, num_unknowns);

    Eigen::VectorXd h;
    Eigen::MatrixXd C_mat = ComputeConstraintMatrix(D0, h);

    Standard_Integer aCompouteTimes = 0;
    Handle(Geom_BSplineCurve) aResultCurve = nullptr;

    while (Standard_True)
    {
        aCompouteTimes++;

        Eigen::MatrixXd A = theAlpha * M_internal + V_internal;

        Eigen::VectorXd b_x(num_unknowns), b_y(num_unknowns), b_z(num_unknowns);
        for (Standard_Integer k = 0; k < num_unknowns; k++)
        {
            Standard_Real D0_term_x = M(0, k + 1) * D0_start.x();
            Standard_Real Dn_term_x = M(n - 1, k + 1) * D0_end.x();
            Standard_Real reg_term_x = V_internal(k, k) * D0_internal_x(k);
            b_x(k) = -theAlpha * (D0_term_x + Dn_term_x) + reg_term_x;

            Standard_Real D0_term_y = M(0, k + 1) * D0_start.y();
            Standard_Real Dn_term_y = M(n - 1, k + 1) * D0_end.y();
            Standard_Real reg_term_y = V_internal(k, k) * D0_internal_y(k);
            b_y(k) = -theAlpha * (D0_term_y + Dn_term_y) + reg_term_y;

            Standard_Real D0_term_z = M(0, k + 1) * D0_start.z();
            Standard_Real Dn_term_z = M(n - 1, k + 1) * D0_end.z();
            Standard_Real reg_term_z = V_internal(k, k) * D0_internal_z(k);
            b_z(k) = -theAlpha * (D0_term_z + Dn_term_z) + reg_term_z;
        }

        Standard_Integer total_vars = 3 * num_unknowns;
        Standard_Integer system_size = total_vars + 4;

        Eigen::SparseMatrix<Standard_Real> BigLHS_sparse(system_size, system_size);
        std::vector<Eigen::Triplet<Standard_Real>> tripletList;
        tripletList.reserve(system_size * 10);

        // 对角块 A
        for (int i = 0; i < num_unknowns; i++)
        {
            for (int j = 0; j < num_unknowns; j++)
            {
                if (std::abs(A(i, j)) > 1e-12)
                {
                    tripletList.push_back(Eigen::Triplet<Standard_Real>(i, j, A(i, j)));
                    tripletList.push_back(Eigen::Triplet<Standard_Real>(i + num_unknowns, j + num_unknowns, A(i, j)));
                    tripletList.push_back(Eigen::Triplet<Standard_Real>(i + 2 * num_unknowns, j + 2 * num_unknowns, A(i, j)));
                }
            }
        }

        // 约束块 C_mat 和 C_mat^T
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < total_vars; j++)
            {
                if (std::abs(C_mat(i, j)) > 1e-12)
                {
                    tripletList.push_back(Eigen::Triplet<Standard_Real>(total_vars + i, j, C_mat(i, j)));
                    tripletList.push_back(Eigen::Triplet<Standard_Real>(j, total_vars + i, C_mat(i, j)));
                }
            }
        }

        BigLHS_sparse.setFromTriplets(tripletList.begin(), tripletList.end());

        Eigen::VectorXd BigRHS(system_size);
        BigRHS.segment(0, num_unknowns) = b_x;
        BigRHS.segment(num_unknowns, num_unknowns) = b_y;
        BigRHS.segment(2 * num_unknowns, num_unknowns) = b_z;
        BigRHS.segment(total_vars, 4) = h;

        Eigen::VectorXd solution(system_size);
        try
        {
            Eigen::SparseLU<Eigen::SparseMatrix<Standard_Real>> solver;
            solver.analyzePattern(BigLHS_sparse);
            solver.factorize(BigLHS_sparse);
            solution = solver.solve(BigRHS);
        }
        catch (...)
        {
            m_errorCode.push_back("Exception occurred during solving constrained system!");
            return nullptr;
        }

        Eigen::VectorXd D_internal_x = solution.segment(0, num_unknowns);
        Eigen::VectorXd D_internal_y = solution.segment(num_unknowns, num_unknowns);
        Eigen::VectorXd D_internal_z = solution.segment(2 * num_unknowns, num_unknowns);

        Eigen::MatrixXd D_internal(num_unknowns, 3);
        D_internal.col(0) = D_internal_x;
        D_internal.col(1) = D_internal_y;
        D_internal.col(2) = D_internal_z;

        D = Eigen::MatrixXd::Zero(n, 3);
        D.row(0) = D0_start;
        D.block(1, 0, num_unknowns, 3) = D_internal;
        D.row(n - 1) = D0_end;

        aResultCurve = CreateNewBSplineCurve(m_OriginalCurve, D);

        Standard_Real aHausdorffDistance = 0;
        if (m_OriginalFitPoints.size() > 0)
        {
            aHausdorffDistance = CurveFair::GetFitPointsCurveHausdorffDistance(m_FitPointParameters, aResultCurve);
            m_FitPointParameters = ReCalculateFitPointParameters(m_OriginalFitPoints, m_OriginalCurve);
        }
        else
        {
            aHausdorffDistance = CurveFair::GetCurveCurveHausdorffDistance(m_OriginalCurve, aResultCurve);
        }

        m_HausdorffDistanceResult = aHausdorffDistance;
        if (aHausdorffDistance <= m_HausdorffDistanceTol || aCompouteTimes >= 100)
            break;
        else
            theAlpha /= 2;
    }

    return aResultCurve;
}


Handle(Geom_BSplineCurve) CurveFair::GetTempContinuniousFairCurve(
    const Handle(Geom_BSplineCurve)& theCurve,
    Eigen::MatrixXd M,          // n×n 矩阵
    Eigen::MatrixXd C,          // m×n 矩阵（包含首尾控制点）
    Eigen::MatrixXd V,          // n×n 矩阵
    Eigen::MatrixXd& D0,        // n×3 初始控制点
    Eigen::MatrixXd& D,         // n×3 输出控制点
    Standard_Real theAlpha)
{
    // 控制点的总数
    Standard_Integer n = D0.rows();
    Standard_Integer m = C.rows();       // 约束条件数量
    Standard_Integer internal_n = n - 2;
    // 固定控制点
    Eigen::Vector3d D0_start = D0.row(0);
    Eigen::Vector3d D0_end = D0.row(n - 1);

    // 获取需要优化的控制点原始值
    Eigen::MatrixXd D0_internal = D0.block(1, 0, internal_n, 3);
    Eigen::VectorXd D0_internal_x = D0.col(0).segment(1, internal_n);  // x分量
    Eigen::VectorXd D0_internal_y = D0.col(1).segment(1, internal_n);  // y分量
    Eigen::VectorXd D0_internal_z = D0.col(2).segment(1, internal_n);  // z分量

    // A = α * M_internal + V_internal
    // 获取需要优化的控制点对应的 M 矩阵和 V 矩阵
    Eigen::MatrixXd M_internal = M.block(1, 1, internal_n, internal_n); // (n - 2) × (n - 2)
    Eigen::MatrixXd V_internal = V.block(1, 1, internal_n, internal_n); // (n - 2) × (n - 2)
    Eigen::MatrixXd C_internal = C.block(0, 1, m, internal_n); // 仅保留内部控制点对应的列 // (m) × (n - 2)

    // 计算右端项 h（首尾控制点对约束的贡献）
    Eigen::VectorXd h_x(m);
    Eigen::VectorXd h_y(m);
    Eigen::VectorXd h_z(m);
    for (Standard_Integer j = 0; j < m; ++j)
    {
        h_x(j) = -C(j, 0) * D0_start.x() - C(j, n - 1) * D0_end.x();
        h_y(j) = -C(j, 0) * D0_start.y() - C(j, n - 1) * D0_end.y();
        h_z(j) = -C(j, 0) * D0_start.z() - C(j, n - 1) * D0_end.z();
    }

    // 构造扩展系统 [ A  C_internal^T ]
    //              [ C_internal  0    ]
    Eigen::MatrixXd extendedA(internal_n + m, internal_n + m);
    extendedA.topRightCorner(internal_n, m) = C_internal.transpose();
    extendedA.bottomLeftCorner(m, internal_n) = C_internal;
    extendedA.bottomRightCorner(m, m).setZero();
    Standard_Integer aCompouteTimes = 0;
    Handle(Geom_BSplineCurve) aResultCurve;
    while (++aCompouteTimes <= 100)
    {
        Eigen::MatrixXd A = theAlpha * M_internal + V_internal;
        extendedA.topLeftCorner(internal_n, internal_n) = A;

        //     分别求解x、y、z分量
        Eigen::VectorXd b_x = Eigen::VectorXd::Zero(internal_n);
        Eigen::VectorXd b_y = Eigen::VectorXd::Zero(internal_n);
        Eigen::VectorXd b_z = Eigen::VectorXd::Zero(internal_n);

        //    减去首尾控制点的影响，分别计算每个分量的b
        for (Standard_Integer k = 0; k < n - 2; k++)
        {
            Standard_Real D0_term_x = M(0, k + 1) * D0_start.x();
            Standard_Real D0_term_y = M(0, k + 1) * D0_start.y();
            Standard_Real D0_term_z = M(0, k + 1) * D0_start.z();
            Standard_Real Dn_term_x = M(n - 1, k + 1) * D0_end.x();
            Standard_Real Dn_term_y = M(n - 1, k + 1) * D0_end.y();
            Standard_Real Dn_term_z = M(n - 1, k + 1) * D0_end.z();

            Standard_Real reg_term_x = V_internal(k, k) * D0_internal_x(k);
            Standard_Real reg_term_y = V_internal(k, k) * D0_internal_y(k);
            Standard_Real reg_term_z = V_internal(k, k) * D0_internal_z(k);

            b_x(k) = -theAlpha * (D0_term_x + Dn_term_x) + reg_term_x;
            b_y(k) = -theAlpha * (D0_term_y + Dn_term_y) + reg_term_y;
            b_z(k) = -theAlpha * (D0_term_z + Dn_term_z) + reg_term_z;
        }

        Eigen::VectorXd rhs_x(internal_n + m);
        rhs_x << b_x, h_x;
        Eigen::VectorXd rhs_y(internal_n + m);
        rhs_y << b_y, h_y;
        Eigen::VectorXd rhs_z(internal_n + m);
        rhs_z << b_z, h_z;

        Eigen::VectorXd solution_x;
        Eigen::VectorXd solution_y;
        Eigen::VectorXd solution_z;
        try
        {
            Eigen::BDCSVD<Eigen::MatrixXd> svd(extendedA, Eigen::ComputeThinU | Eigen::ComputeThinV);
            solution_x = svd.solve(rhs_x);
            solution_y = svd.solve(rhs_y);
            solution_z = svd.solve(rhs_z);
        }
        catch (const Standard_Failure& e)
        {
            m_errorCode.push_back(e.GetMessageString());
            return nullptr;
        }
        catch (const std::exception& e)
        {
            m_errorCode.push_back(e.what());
            return nullptr;
        }
        catch (...)
        {
            m_errorCode.push_back("Unknown Exception occurred!");
            return nullptr;
        }


        Eigen::VectorXd D_internal_x = solution_x.head(internal_n);
        Eigen::VectorXd D_internal_y = solution_y.head(internal_n);
        Eigen::VectorXd D_internal_z = solution_z.head(internal_n);

        //     合并结果
        Eigen::MatrixXd D_internal(n - 2, 3);
        D_internal.col(0) = D_internal_x;
        D_internal.col(1) = D_internal_y;
        D_internal.col(2) = D_internal_z;

        //     构建最终控制点矩阵
        D = Eigen::MatrixXd::Zero(n, 3);
        D.row(0) = D0_start;                    // 起始点
        D.block(1, 0, n - 2, 3) = D_internal;   // 内部点
        D.row(n - 1) = D0_end;                  // 终止点

        //    创建新的B样条曲线
        aResultCurve = CreateNewBSplineCurve(m_OriginalCurve, D);

        // 计算新旧曲线的控制点偏差
        Standard_Real aHausdorffDistance = 0;
        if (m_OriginalFitPoints.size() > 0)
        {
            aHausdorffDistance = CurveFair::GetFitPointsCurveHausdorffDistance(m_FitPointParameters, aResultCurve);
            m_FitPointParameters = ReCalculateFitPointParameters(m_OriginalFitPoints, m_OriginalCurve);
        }
        else
        {
            aHausdorffDistance = CurveFair::GetCurveCurveHausdorffDistance(m_OriginalCurve, aResultCurve);
        }

        m_HausdorffDistanceResult = aHausdorffDistance;
        if (aHausdorffDistance <= m_HausdorffDistanceTol || aCompouteTimes >= 100)
        {
            break;
        }
        else
        {
            theAlpha /= 2;
        }
    }
    return aResultCurve;
}

void CurveFair::Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D)
{
    // 分步移动参数
    Standard_Integer n = D0.rows();
    Standard_Real stepSize = 10;
    Standard_Integer aStepCnt = 0;
    while (Standard_True)
    {
        // 计算下一步的位置
        D0 += (D - D0) / 10;
        aStepCnt++;

        // 生成新曲线
        Handle(Geom_BSplineCurve) aCurve = CreateNewBSplineCurve(m_OriginalCurve, D0);
        M = ComputeEnergyMatrix(aCurve, m_OriginalCurve->Degree());
        //C = ComputeContinuityMatrix(aCurve, m_OriginalCurve->Degree());
        SetControlPointWeightMatrix(aCurve, V);
        aCurve = GetTempFairCurveWithTangentConstraint(aCurve, M, V, D0, D, m_Alpha);

        if (aStepCnt == stepSize)
        {
            break;
        }
        else
        {
            m_ResultCurve = aCurve;
        }
    }
}


Handle(Geom_BSplineCurve) CurveFair::SampleAndFitBSpline(
    const Handle(Geom_BSplineCurve)& theOriginalCurve,
    Standard_Integer theSampleNum,
    std::vector<gp_Pnt>& theFitPoints,
    Standard_Integer theMaxDegree,
    GeomAbs_Shape theContinuity,
    Standard_Real theTolerance)
{
    if (theSampleNum < 2 || theOriginalCurve.IsNull())
    {
        std::cerr << "Invalid input: too few points or null curve." << std::endl;
        return nullptr;
    }

    GeomAdaptor_Curve adaptor(theOriginalCurve);
    Standard_Real totalLength = GCPnts_AbscissaPoint::Length(adaptor);

    TColgp_Array1OfPnt sampledPoints(1, theSampleNum);
    for (Standard_Integer i = 0; i < theSampleNum; ++i)
    {
        Standard_Real targetLength = (totalLength * i) / (theSampleNum - 1);
        GCPnts_AbscissaPoint gap(adaptor, targetLength, theOriginalCurve->FirstParameter());
        Standard_Real u = gap.Parameter();
        gp_Pnt pt = theOriginalCurve->Value(u);
        sampledPoints.SetValue(i + 1, pt);
    }
    theFitPoints = OccArrayConvertoVector(sampledPoints);

    Handle(Geom_BSplineCurve) aBSplineCurve;

    GeomAPI_PointsToBSpline fitter(
        sampledPoints,
        theMaxDegree,       // 最小次数 = 最大次数
        theMaxDegree,
        theContinuity,
        theTolerance
    );
    aBSplineCurve = fitter.Curve();


    return aBSplineCurve;
}

std::vector<gp_Pnt> CurveFair::SampleCurveWithArcLength(const Handle(Geom_BSplineCurve)& bsplineCurve, Standard_Integer numSamples)
{
    Handle(GeomAdaptor_Curve) curve = new GeomAdaptor_Curve(bsplineCurve);
    GCPnts_QuasiUniformAbscissa sampler(*curve, numSamples);
    std::vector<gp_Pnt> sampledPointsWithArcLength;

    for (Standard_Integer i = 1; i <= sampler.NbPoints(); ++i)
    {
        Standard_Real param = sampler.Parameter(i);
        gp_Pnt point = curve->Value(param); // 获取点坐标

        // 更新点
        sampledPointsWithArcLength.emplace_back(point);
    }

    return sampledPointsWithArcLength;
}

Handle(Geom_BSplineCurve) CurveFair::SampleAndFitBSpline(
    std::vector<gp_Pnt>& theFitPoints,
    Standard_Real theTolerance,
    Standard_Boolean theResortFlag,
    Standard_Integer theMaxDegree,
    GeomAbs_Shape theContinuity)
{
    if (theFitPoints.size() < 2)
    {
        std::cerr << "Invalid input: too few points." << std::endl;
        return nullptr;
    }


    if (theResortFlag)
    {
        ReorderPointsNearestNeighbor(theFitPoints);
    }

    Handle(Geom_BSplineCurve) aBSplineCurve;

    // 转换为 OpenCASCADE 点数组
    TColgp_Array1OfPnt aSamplePointArray(1, static_cast<Standard_Integer>(theFitPoints.size()));
    for (Standard_Integer i = 0; i < static_cast<Standard_Integer>(theFitPoints.size()); ++i)
    {
        aSamplePointArray.SetValue(i + 1, theFitPoints[i]);
    }
    GeomAPI_PointsToBSpline fitter(
        aSamplePointArray,
        Approx_ChordLength,
        theMaxDegree,
        theMaxDegree,
        theContinuity,
        theTolerance
    );
    aBSplineCurve = fitter.Curve();

    return aBSplineCurve;
}

Handle(Geom_BSplineCurve) CurveFair::RefineCurveByCurvatureAuto(
    const Handle(Geom_BSplineCurve)& theCurve,
    const std::vector<Standard_Real>& myKnotSeq,
    const Standard_Integer baseInsertNum)
{
    if (theCurve.IsNull()) return theCurve;

    Handle(Geom_BSplineCurve) refinedCurve = Handle(Geom_BSplineCurve)::DownCast(theCurve->Copy());
    const Standard_Integer degree = refinedCurve->Degree();
    const Standard_Integer nPoles = refinedCurve->NbPoles();
    const TColStd_Array1OfReal& knotArray = refinedCurve->Knots();

    const Standard_Integer nbKnots = knotArray.Length();
    Standard_Integer lower = knotArray.Lower();
    Standard_Integer upper = knotArray.Upper();
    Standard_Real KnotUpper = knotArray(knotArray.Upper());
    Standard_Real KnotLower = knotArray(knotArray.Lower());
    std::vector<Standard_Real> curvatureRates(nbKnots);
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(refinedCurve);

    // Step 1: 计算每个节点处的曲率指标
    for (Standard_Integer i = lower; i <= upper; ++i)
    {
        Standard_Real u = knotArray(i);
        Standard_Real s = f_inverse(refinedCurve, u);
        Standard_Real f0 = f(s);
        Standard_Real f1 = f(s, 1);
        Standard_Real f2 = f(s, 2);
        Standard_Real f3 = f(s, 3);

        Standard_Real sum = 0.0;
        for (Standard_Integer j = 0; j < nPoles; ++j)
        {
            Standard_Real k = BasisFunctionDerivative(f0, j, degree, 3, myKnotSeq) * std::pow(f1, 3)
                + 3 * BasisFunctionDerivative(f0, j, degree, 2, myKnotSeq) * f1 * f2
                + BasisFunctionDerivative(f0, j, degree, 1, myKnotSeq) * f3;
            sum += k * k;
        }

        curvatureRates[i - lower] = sum;
    }

    // Step 2: 计算平均曲率
    Standard_Real sumCurv = 0.0;
    for (auto c : curvatureRates) sumCurv += c;
    Standard_Real avgCurv = sumCurv / curvatureRates.size();

    // Step 3: 收集所有需要插入的 u 值（曲率大于阈值的中点区间）
    std::vector<Standard_Real> insertParams;

    for (Standard_Integer i = lower; i < upper; ++i)
    {
        if (curvatureRates[i - lower] > avgCurv)
        {
            Standard_Real left = knotArray(i);
            Standard_Real right = knotArray(i + 1);

            if (right - left < (KnotUpper - KnotLower) / 100) continue; // 避免在非常小的区间内插值

            // 自适应插值个数（可按比例加权）
            Standard_Integer insertNum = baseInsertNum;
            if (curvatureRates[i - lower] > 2 * avgCurv)
                insertNum *= 3;

            for (Standard_Integer j = 1; j <= insertNum; ++j)
            {
                Standard_Real u = left + (right - left) * j / (insertNum + 1);
                insertParams.push_back(u);
            }
        }
    }

    // Step 4: 插入所有收集到的 knot 参数值
    if (!insertParams.empty())
    {
        TColStd_Array1OfReal uArray(1, static_cast<Standard_Integer>(insertParams.size()));
        TColStd_Array1OfInteger mults(1, static_cast<Standard_Integer>(insertParams.size()));
        for (Standard_Integer i = 0; i < insertParams.size(); ++i)
        {
            uArray.SetValue(i + 1, insertParams[i]);
            mults.SetValue(i + 1, 1);
        }

        refinedCurve->InsertKnots(uArray, mults, 1e-7, Standard_True);
    }

    return refinedCurve;
}
Handle(Geom_BSplineCurve) CurveFair::RefineCurveByFitPoints(
    const Handle(Geom_BSplineCurve)& theCurve,
    const Standard_Integer baseInsertNum)
{
    if (theCurve.IsNull() || baseInsertNum <= 0)
        return theCurve;          // 兜底

    Handle(Geom_BSplineCurve) refinedCurve =
        Handle(Geom_BSplineCurve)::DownCast(theCurve->Copy());

    const TColStd_Array1OfReal& knotArray = refinedCurve->Knots();
    const Standard_Integer lower = knotArray.Lower();
    const Standard_Integer upper = knotArray.Upper();

    constexpr Standard_Real tol = 1e-12;               // 节点比较公差
    std::map<std::pair<Standard_Integer, Standard_Integer>, int> spanCounts; // 区间计数

    //=== 1) 统计每个区间的拟合点数量 =========================
    for (const auto& fp : m_FitPointParameters)        // (gp_Pnt , uFit)
    {
        Standard_Real uFit = fp.second;

        if (uFit < knotArray(lower) - tol || uFit > knotArray(upper) + tol)
            continue;

        // 线性扫描找 span
        Standard_Integer span = -1;
        for (Standard_Integer i = lower; i < upper; ++i)
        {
            if (uFit + tol < knotArray(i)) break;      // 过早
            if (uFit >= knotArray(i) - tol &&
                uFit < knotArray(i + 1) - tol)
            {
                span = i;
                break;
            }
        }

        if (span == -1) continue;

        Standard_Real uL = knotArray(span);
        Standard_Real uR = knotArray(span + 1);
        if (uR - uL < tol) continue;                   // 极小区间跳过

        spanCounts[{span, span + 1}]++;  // 统计区间出现次数
    }

    //=== 2) 根据区间频率计算插入节点 =========================
    std::vector<Standard_Real> insertParams;
    for (const auto& entry : spanCounts)
    {
        Standard_Integer span = entry.first.first;
        Standard_Real uL = knotArray(span);
        Standard_Real uR = knotArray(span + 1);
        int count = entry.second;

        // 根据区间频率决定插入节点数
        int actualInsertNum = baseInsertNum * count;
        for (Standard_Integer j = 1; j <= actualInsertNum; ++j)
        {
            Standard_Real uNew = uL + (uR - uL) * j / (actualInsertNum + 1);
            insertParams.push_back(uNew);
        }
    }

    if (!insertParams.empty())
    {
        std::sort(insertParams.begin(), insertParams.end());

        // unique - 确保没有重复节点
        insertParams.erase(std::unique(
            insertParams.begin(), insertParams.end(),
            [](Standard_Real a, Standard_Real b) { return std::abs(a - b) < 1e-10; }),
            insertParams.end());

        TColStd_Array1OfReal    uArray(1, static_cast<Standard_Integer>(insertParams.size()));
        TColStd_Array1OfInteger mArray(1, static_cast<Standard_Integer>(insertParams.size()));
        for (Standard_Integer i = 0; i < insertParams.size(); ++i)
        {
            uArray.SetValue(i + 1, insertParams[i]);
            mArray.SetValue(i + 1, 1);   // 单重插 knot 即可
        }

        refinedCurve->InsertKnots(uArray, mArray, 1e-7, Standard_True /* scaled */);
    }

    return refinedCurve;
}

Handle(Geom_BSplineCurve) CurveFair::RefineByTangentAngle(
    const Handle(Geom_BSplineCurve)& theCurve,
    const Standard_Integer theBasicInsertNum)
{
    if (theCurve.IsNull()) return theCurve;

    Handle(Geom_BSplineCurve) refinedCurve = Handle(Geom_BSplineCurve)::DownCast(theCurve->Copy());
    std::vector<Standard_Real> knots = OccArrayConvertoVector(refinedCurve->Knots());

    // 统计每个节点的左右导数夹角
    std::vector<Standard_Real> angles(knots.size(), 0.0);
    for (size_t i = 1; i < knots.size() - 1; ++i)
    {
        gp_Vec leftVec = Dt(refinedCurve, knots[i] - 1e-6, 3);  // 1阶导数
        gp_Vec rightVec = Dt(refinedCurve, knots[i] + 1e-6, 3);
        angles[i] = leftVec.Angle(rightVec) * 180.0 / M_PI;
        std::cout << angles[i] << std::endl;
    }

    Standard_Real avgAngle = 0.0;
    int validCount = 0;
    for (size_t i = 1; i < angles.size() - 1; ++i)
    {
        avgAngle += angles[i];
        ++validCount;
    }

    avgAngle /= std::max(validCount, 1);
    // 统计每个区间的插值数（左右均分）
    std::map<std::pair<Standard_Real, Standard_Real>, int> spanInsertCount;

    for (size_t i = 1; i < knots.size() - 1; ++i)
    {
        int insertNum = theBasicInsertNum + (angles[i] / avgAngle - 1);
        if (insertNum > 0)
        {
            Standard_Real uLeft = knots[i - 1];
            Standard_Real uMid = knots[i];
            Standard_Real uRight = knots[i + 1];

            spanInsertCount[{uLeft, uMid}] += insertNum;
            spanInsertCount[{uMid, uRight}] += insertNum;
        }
    }

    // 统一收集插值参数
    std::vector<Standard_Real> insertParams;
    for (const auto& kv : spanInsertCount)
    {
        Standard_Real u1 = kv.first.first;
        Standard_Real u2 = kv.first.second;
        int count = kv.second;

        for (int i = 1; i <= count; ++i)
        {
            Standard_Real u = u1 + i * (u2 - u1) / (count + 1);
            insertParams.push_back(u);
        }
    }

    if (insertParams.empty())
        return refinedCurve;

    // 去重 + 插值
    std::sort(insertParams.begin(), insertParams.end());
    insertParams.erase(std::unique(insertParams.begin(), insertParams.end(),
        [](Standard_Real a, Standard_Real b) { return std::abs(a - b) < 1e-10; }), insertParams.end());

    // 准备参数数组
    TColStd_Array1OfReal    uArray(1, static_cast<Standard_Integer>(insertParams.size()));
    TColStd_Array1OfInteger mArray(1, static_cast<Standard_Integer>(insertParams.size()));

    for (Standard_Integer i = 0; i < insertParams.size(); ++i)
    {
        uArray.SetValue(i + 1, insertParams[i]);
        mArray.SetValue(i + 1, 1); // multiplicity = 1
    }

    // 插入 knot，使用 scaled 模式
    refinedCurve->InsertKnots(uArray, mArray, 1e-7, Standard_True);

    return refinedCurve;
}

std::vector<Standard_Real> CurveFair::ComputeThirdDerivativeAngles(
    const Handle(Geom_BSplineCurve)& theCurve)
{
    std::vector<Standard_Real> angleResults;

    if (theCurve.IsNull())
        return angleResults;

    // 拷贝一份曲线，避免原始数据被修改
    Handle(Geom_BSplineCurve) curveCopy = Handle(Geom_BSplineCurve)::DownCast(theCurve->Copy());

    // 获取 knots
    std::vector<Standard_Real> knots = OccArrayConvertoVector(curveCopy->Knots());

    if (knots.size() < 3)
        return angleResults; // 没有内部节点

    angleResults.resize(knots.size(), 0.0);

    for (size_t i = 1; i < knots.size() - 1; ++i)
    {
        Standard_Real u = knots[i];
        gp_Vec leftVec = Dt(curveCopy, u - 1e-6, 3);   // 三阶导数（左极限）
        gp_Vec rightVec = Dt(curveCopy, u + 1e-6, 3);  // 三阶导数（右极限）

        if (leftVec.Magnitude() > gp::Resolution() && rightVec.Magnitude() > gp::Resolution())
        {
            angleResults[i] = leftVec.Angle(rightVec) * 180.0 / M_PI; // 返回弧度值
        }
        else
        {
            angleResults[i] = 0.0; // 有一侧为零向量，不计算夹角
        }
    }

    return angleResults;
}

Handle(Geom_BSplineCurve) CurveFair::InsertKnotsBetweenKnotSpan(
    const Handle(Geom_BSplineCurve)& theCurve,
    Standard_Real t1,
    Standard_Real t2,
    Standard_Integer nInsert)
{
    if (theCurve.IsNull() || nInsert <= 0 || t1 >= t2)
    {
        std::cerr << "Invalid input to InsertKnotsBetween." << std::endl;
        return theCurve;
    }

    // 获取原始节点向量
    TColStd_Array1OfReal knots = theCurve->Knots();
    TColStd_Array1OfInteger mults = theCurve->Multiplicities();

    Standard_Integer nbKnots = theCurve->NbKnots();
    Standard_Real knotMin = knots.First();
    Standard_Real knotMax = knots.Last();

    // 检查 t1 和 t2 是否在 knot 范围内
    if (t1 < knotMin || t2 > knotMax)
    {
        std::cerr << "t1 or t2 out of knot vector range." << std::endl;
        return theCurve;
    }

    // 生成 nInsert 个均匀新节点
    std::vector<Standard_Real> newKnots;
    for (int i = 1; i <= nInsert; ++i)
    {
        Standard_Real newKnot = t1 + i * (t2 - t1) / (nInsert + 1);
        newKnots.push_back(newKnot);
    }

    // 复制曲线并插入节点
    Handle(Geom_BSplineCurve) newCurve = Handle(Geom_BSplineCurve)::DownCast(theCurve->Copy());
    for (const auto& u : newKnots)
    {
        newCurve->InsertKnot(u, 1, 1e-7); // multiplicity=1, tolerance可调
    }

    return newCurve;
}

Handle(Geom_BSplineCurve) CurveFair::InsertKnots(
    const Handle(Geom_BSplineCurve)& theCurve,
    const std::vector<Standard_Real>& params,
    Standard_Integer mult)
{
    if (theCurve.IsNull() || params.empty() || mult <= 0)
    {
        std::cerr << "Invalid input to InsertKnotsAtParameters." << std::endl;
        return theCurve;
    }

    Handle(Geom_BSplineCurve) newCurve = Handle(Geom_BSplineCurve)::DownCast(theCurve->Copy());

    for (const auto& u : params)
    {
        newCurve->InsertKnot(u, mult, 1e-7); // 插入 knot，允许重复
    }

    return newCurve;
}

Handle(Geom_BSplineCurve) CurveFair::InsertKnot(
    const Handle(Geom_BSplineCurve)& theCurve,
    Standard_Real u,
    Standard_Integer mult)
{
    if (theCurve.IsNull() || mult <= 0)
    {
        std::cerr << "Invalid input to InsertSingleKnot." << std::endl;
        return theCurve;
    }

    Handle(Geom_BSplineCurve) newCurve = Handle(Geom_BSplineCurve)::DownCast(theCurve->Copy());
    newCurve->InsertKnot(u, mult, 1e-7); // 插入节点，默认容差
    return newCurve;
}


//! @brief 在每一对相邻节点区间内插入 nInsert 个均匀分布的新节点
//! @param [in] theCurve     原始曲线
//! @param [in] nInsert      每个区间插入的新节点数量
//! @return 插值后的新曲线
Handle(Geom_BSplineCurve) CurveFair::InsertUniformKnotsInAllSpans(
    const Handle(Geom_BSplineCurve)& theCurve,
    Standard_Integer nInsert)
{
    if (theCurve.IsNull() || nInsert <= 0)
    {
        std::cerr << "Invalid input to InsertUniformKnotsInAllSpans." << std::endl;
        return theCurve;
    }

    Handle(Geom_BSplineCurve) newCurve = Handle(Geom_BSplineCurve)::DownCast(theCurve->Copy());
    TColStd_Array1OfReal knotsArray = newCurve->Knots();
    Standard_Integer nbKnots = newCurve->NbKnots();
    Standard_Real upperKnot = knotsArray.Upper();
    Standard_Real lowerKnot = knotsArray.Lower();
    std::vector<Standard_Real> newKnots;
    for (Standard_Integer i = 1; i < nbKnots; ++i)
    {
        Standard_Real t1 = knotsArray.Value(i);
        Standard_Real t2 = knotsArray.Value(i + 1);

        Standard_Integer inserKnot = nInsert;
        if (t2 - t1 <= (upperKnot - lowerKnot) / 400)
        {
            inserKnot = 0;
        }
        else if (t2 - t1 <= (upperKnot - lowerKnot) / 200)
        {
            inserKnot = 1;
        }
        else if (t2 - t1 <= (upperKnot - lowerKnot) / 100)
        {
            inserKnot = nInsert / 3;

        }
        else if (t2 - t1 <= (upperKnot - lowerKnot) / 50)
        {
            inserKnot = nInsert / 2;
        }
        for (Standard_Integer j = 1; j <= inserKnot; ++j)
        {
            Standard_Real knot = t1 + j * (t2 - t1) / (inserKnot + 1);
            newKnots.push_back(knot);
        }
    }

    // 插入所有新节点
    for (const auto& u : newKnots)
    {
        newCurve->InsertKnot(u, 1, 1e-7); // multiplicity = 1
    }

    return newCurve;
}


void CurveFair::SetControlPointWeightMatrix(const Handle(Geom_BSplineCurve)& theCurve, Eigen::MatrixXd& V)
{

    V.resize(theCurve->NbPoles(), theCurve->NbPoles());
    V.setZero();
    for (Standard_Integer i = 0; i < theCurve->NbPoles(); i++)
    {
        // 暂定全为 1
        V(i, i) = 1;
    }

    return; 

    V.resize(theCurve->NbPoles(), theCurve->NbPoles());
    V.setZero();
    std::vector<Standard_Real> myKnotSeq = OccArrayConvertoVector(theCurve->KnotSequence());
    std::vector<Standard_Real> aGrevilleAbscissaeArray(theCurve->NbPoles());
    for (Standard_Integer i = 0; i < theCurve->NbPoles(); i++)
    {
        Standard_Real sum = 0.0;
        for (Standard_Integer j = 1; j <= theCurve->Degree(); j++)
        {
            Standard_Real aKnotValue = myKnotSeq[i + j];
            //aKnotValue = f_inverse(theCurve, aKnotValue);
            sum += aKnotValue;
        }
        aGrevilleAbscissaeArray[i] = sum / theCurve->Degree();
    }

    std::vector<Standard_Real> curvatureRates(theCurve->NbPoles());
    for (Standard_Integer i = 0; i < theCurve->NbPoles(); i++)
    {
        Standard_Real s = f_inverse(theCurve, aGrevilleAbscissaeArray[i]);
        Standard_Real f0 = f(s);
        Standard_Real f1 = f(s, 1);
        Standard_Real f2 = f(s, 2);
        Standard_Real f3 = f(s, 3);

        // Ni
        Standard_Real aCurvature = BasisFunctionDerivative(f0, i, theCurve->Degree(), 3, myKnotSeq) * pow(f1, 3)
            + 3 * BasisFunctionDerivative(f0, i, theCurve->Degree(), 2, myKnotSeq) * f1 * f2
            + BasisFunctionDerivative(f0, i, theCurve->Degree(), 1, myKnotSeq) * f3;

        curvatureRates[i] = aCurvature * aCurvature;
    }

    Standard_Real aMinRate = *std::min_element(curvatureRates.begin(), curvatureRates.end());
    Standard_Real aMaxRate = *std::max_element(curvatureRates.begin(), curvatureRates.end());

    for (Standard_Integer i = 0; i < theCurve->NbPoles(); i++)
    {
        // 归一化 curvatureRate 到 [0,1]
        Standard_Real normalized = (curvatureRates[i] - aMinRate) / (aMaxRate - aMinRate + 1e-6);
        // 映射到 [1.0, 0.1]
        Standard_Real maxWeight = 1;
        Standard_Real weight = maxWeight - maxWeight * normalized;

        // 暂定全为 1
        V(i, i) = 1;
    }
}

std::vector<gp_Pnt> CurveFair::ReorderPointsNearestNeighbor(const std::vector<gp_Pnt>& points)
{
    if (points.empty()) return {};

    std::vector<gp_Pnt> ordered;
    std::vector<bool> visited(points.size(), false);
    ordered.reserve(points.size());

    int currentIdx = 0;
    ordered.push_back(points[currentIdx]);
    visited[currentIdx] = true;

    for (std::size_t i = 1; i < points.size(); ++i)
    {
        Standard_Real minDist = std::numeric_limits<Standard_Real>::max();
        int nextIdx = -1;

        for (std::size_t j = 0; j < points.size(); ++j)
        {
            if (!visited[j])
            {
                Standard_Real d = points[currentIdx].Distance(points[j]);
                if (d < minDist)
                {
                    minDist = d;
                    nextIdx = j;
                }
            }
        }

        if (nextIdx >= 0)
        {
            ordered.push_back(points[nextIdx]);
            visited[nextIdx] = true;
            currentIdx = nextIdx;
        }
    }

    return ordered;
}
void CurveFair::Perform()
{
    M = ComputeEnergyMatrix(m_OriginalCurve, m_OriginalCurve->Degree());
    //C = ComputeContinuityMatrix(theCurve, theCurve->Degree());
    SetControlPointWeightMatrix(m_OriginalCurve, V);

    TColgp_Array1OfPnt Poles = m_OriginalCurve->Poles();
    Standard_Integer n = m_OriginalCurve->NbPoles();

    Eigen::MatrixXd D0(n, 3); // 原始控制点
    Eigen::MatrixXd D(n, 3); // 新控制点
    for (Standard_Integer i = 0; i < n; i++)
    {
        D0(i, 0) = Poles.Value(i + 1).X();
        D0(i, 1) = Poles.Value(i + 1).Y();
        D0(i, 2) = Poles.Value(i + 1).Z();
    }

    // 获取优化曲线
    m_ResultCurve = GetTempFairCurveWithTangentConstraint(m_OriginalCurve, M, V, D0, D, m_Alpha);

    Iterator(D0, D);
}