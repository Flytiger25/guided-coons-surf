#include "CurveFair.h"
#include <BSplineAlgo/Fitting/BSCrvFitting.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

// 定义静态成员变量
sggk::BSplineCurve3DPtr CurveFair::m_ArcLengthMappingFunction;

const double M_PI_VAL = 3.14159265358979323846;
const double RESOLUTION = 1e-12;

// Gauss-Legendre 30 点积分表（区间 [-1,1]）
static const double GAUSS_PTS[30] = {
    9.96893484074649505e-01, 9.83668123279747175e-01, 9.60021864968307548e-01,
    9.26200047429274309e-01, 8.82560535792052625e-01, 8.29565762382768357e-01,
    7.67777432104826185e-01, 6.97850494793315734e-01, 6.20526182989242892e-01,
    5.36624148142019863e-01, 4.47033769538089154e-01, 3.52704725530878116e-01,
    2.54636926167889854e-01, 1.53869913608583542e-01, 5.14718425553176914e-02,
    -5.14718425553176984e-02, -1.53869913608583542e-01, -2.54636926167889854e-01,
    -3.52704725530878116e-01, -4.47033769538089210e-01, -5.36624148142019863e-01,
    -6.20526182989242892e-01, -6.97850494793315734e-01, -7.67777432104826185e-01,
    -8.29565762382768357e-01, -8.82560535792052625e-01, -9.26200047429274309e-01,
    -9.60021864968307548e-01, -9.83668123279747175e-01, -9.96893484074649505e-01
};
static const double GAUSS_WTS[30] = {
    7.96819249616664317e-03, 1.84664683110909757e-02, 2.87847078833234002e-02,
    3.87991925696270917e-02, 4.84026728305940665e-02, 5.74931562176190861e-02,
    6.59742298821804907e-02, 7.37559747377051628e-02, 8.07558952294201160e-02,
    8.68997872010829897e-02, 9.21225222377861086e-02, 9.63687371746443644e-02,
    9.95934205867951977e-02, 1.01762389748405471e-01, 1.02852652893558924e-01,
    1.02852652893558924e-01, 1.01762389748405471e-01, 9.95934205867951977e-02,
    9.63687371746443644e-02, 9.21225222377861086e-02, 8.68997872010829897e-02,
    8.07558952294201160e-02, 7.37559747377051628e-02, 6.59742298821804907e-02,
    5.74931562176190861e-02, 4.84026728305940665e-02, 3.87991925696270917e-02,
    2.87847078833234002e-02, 1.84664683110909757e-02, 7.96819249616664317e-03
};

bool isEqual(double x, double y, double epsilon = 1e-10)
{
    return std::fabs(x - y) < epsilon;
}
bool isGreaterThan(double x, double y, double epsilon = 1e-10)
{
    return (x - y) > epsilon;
}
bool isLessThan(double x, double y, double epsilon = 1e-10)
{
    return (y - x) > epsilon;
}
bool isGreaterThanOrEqual(double x, double y, double epsilon = 1e-10)
{
    return (x - y) > -epsilon;
}
bool isLessThanOrEqual(double x, double y, double epsilon = 1e-10)
{
    return (y - x) > -epsilon;
}

// 拷贝一条 B 样条曲线
static sggk::BSplineCurve3DPtr CloneBSpline(const sggk::BSplineCurve3DPtr& c)
{
    if (!c) return nullptr;
    return std::make_shared<sggk::BSplineCurve3D>(*c);
}

// 计算曲线在 [a,b] 区间内的弧长（Gauss-Legendre 30 点数值积分 |D1|）
double CurveFair::ComputeCurveLength(const sggk::BSplineCurve3DPtr& curve, double tol)
{
    if (!curve) return 0.0;
    double a = curve->MinParam();
    double b = curve->MaxParam();
    double len = 0.0;
    for (int i = 0; i < 30; ++i)
    {
        double u = (b - a) * GAUSS_PTS[i] / 2.0 + (a + b) / 2.0;
        sggk::Vector3D d1 = curve->CalcDeriv1(u);
        len += GAUSS_WTS[i] * d1.Length();
    }
    return len * (b - a) / 2.0;
}

void CurveFair::UniformCurve(sggk::BSplineCurve3DPtr& curve)
{
    if (!curve) return;
    double k0 = curve->Knots().front();
    double k1 = curve->Knots().back();
    if (!(k0 == 0 && k1 == 1))
    {
        curve->AdjustKnots(sggk::Interval(0.0, 1.0));
    }
}

// To compute the value of a b-spline basic function value
static double OneBasicFun(
    const double u,
    const int i,
    const int p,
    const std::vector<double>& Knots)
{
    double Nip, uleft, uright, saved, temp;
    int m = (int)Knots.size() - 1;
    std::vector<double>N(p + 1);
    if ((i == 0 && isEqual(u, Knots[0])) || (i == m - p - 1 && isEqual(u, Knots[m])))
    {
        return 1.0;
    }

    if (isLessThan(u, Knots[i]) || isGreaterThanOrEqual(u, Knots[i + p + 1]))
    {
        return 0.0;
    }

    for (size_t j = 0; j <= (size_t)p; j++)
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
    for (size_t k = 1; k <= (size_t)p; k++)
    {
        if (N[0] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = ((u - Knots[i]) * N[0]) / (Knots[i + k] - Knots[i]);
        }
        for (size_t j = 0; j < (size_t)p - k + 1; j++)
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

sggk::BSplineCurve3DPtr CurveFair::GetArclengthParameterMapping(const sggk::BSplineCurve3DPtr& theCurve, const double theTolerance)
{
    // bspline属性
    double aFirstParam = theCurve->MinParam();
    double aLastParam = theCurve->MaxParam();

    // bspline: 计算弧长
    double aBsplineLen = 0;
    double aAvgLen = 0;
    double aParamStep = 0;
    try
    {
        aBsplineLen = ComputeCurveLength(theCurve, theTolerance);
        aAvgLen = aBsplineLen / (m_ArcLengthMappingSampleNum - 1.0);
        aParamStep = (aLastParam - aFirstParam) / (m_ArcLengthMappingSampleNum - 1.0);
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

    // 弧长参数化: 计算新的节点向量与采样点
    std::vector<double> reparamKnots;
    std::vector<sggk::Point3D> reparamPoints;
    reparamKnots.push_back(aFirstParam);
    reparamPoints.push_back(sggk::Point3D(aFirstParam, aFirstParam, aFirstParam));

    double curLen = 0.0;
    double curParam = aFirstParam;
    for (int i = 1; i < (int)m_ArcLengthMappingSampleNum - 1; i++)
    {
        curLen += aAvgLen;
        curParam += aParamStep;

        try
        {
            sggk::Point3D p;
            double tparam = theCurve->CalcParaByLength(curLen, aFirstParam, p, theTolerance);
            reparamKnots.push_back(curParam);
            reparamPoints.push_back(sggk::Point3D(tparam, tparam, tparam));
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

    reparamKnots.push_back(aLastParam);
    reparamPoints.push_back(sggk::Point3D(aLastParam, aLastParam, aLastParam));

    // 用 BSCrvFitting 插值拟合弧长映射函数 f(s) = t
    sggk::BSCrv3DInterpolationOpts opts;
    opts.degree = theCurve->Degree();
    opts.params = reparamKnots;
    sggk::BSCrvFitting3DResult res = sggk::BSCrvFitting::Interpolation3D(reparamPoints, opts);
    return res.curve;
}

double CurveFair::f(
    const double theParameter,
    const int k)
{
    if (k == 0) return m_ArcLengthMappingFunction->CalcPoint(theParameter).X();
    if (k == 1) return m_ArcLengthMappingFunction->CalcDeriv1(theParameter).X();
    if (k == 2) return m_ArcLengthMappingFunction->CalcDeriv2(theParameter).X();
    if (k == 3) return m_ArcLengthMappingFunction->CalcDeriv3(theParameter).X();
    return 0.0;
}

sggk::Vector3D CurveFair::Ds(sggk::BSplineCurve3DPtr& theBSplineCurve, double sParameter, int k)
{
    if (!theBSplineCurve)
        return sggk::Vector3D(0, 0, 0);

    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve);
    double u = f(sParameter, 0);

    double f1 = f(sParameter, 1);
    double f2 = f(sParameter, 2);
    double f3 = f(sParameter, 3);

    sggk::Vector3D D1 = theBSplineCurve->CalcDeriv1(u);
    if (k == 1)
        return D1 * f1;

    sggk::Vector3D D2 = theBSplineCurve->CalcDeriv2(u);
    if (k == 2)
        return D2 * (f1 * f1) + D1 * f2;

    sggk::Vector3D D3 = theBSplineCurve->CalcDeriv3(u);
    if (k == 3)
        return D3 * (f1 * f1 * f1) + D2 * (3.0 * f1 * f2) + D1 * f3;

    if (k == 0)
    {
        sggk::Point3D p = theBSplineCurve->CalcPoint(u);
        return sggk::Vector3D(p.X(), p.Y(), p.Z());
    }

    return sggk::Vector3D(0, 0, 0);
}

sggk::Vector3D CurveFair::Dt(sggk::BSplineCurve3DPtr& theBSplineCurve, double tParameter, int k)
{
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve);
    return Ds(theBSplineCurve, f_inverse(theBSplineCurve, tParameter), k);
}

double CurveFair::BasisFunctionDerivative(
    const double u,
    const int i,
    const int p,
    const int k,
    const std::vector<double>& Knots)
{
    if (k == 0)
    {
        return OneBasicFun(u, i, p, Knots);
    }

    if (k > p || p == 0)
    {
        return 0.0;
    }

    double Term1 = 0.0;
    double Denominator1 = Knots[i + p] - Knots[i];
    if (std::abs(Denominator1) > std::numeric_limits<double>::epsilon())
    {
        Term1 = (p / Denominator1) * BasisFunctionDerivative(u, i, p - 1, k - 1, Knots);
    }

    double Term2 = 0.0;
    double Denominator2 = Knots[i + p + 1] - Knots[i + 1];
    if (std::abs(Denominator2) > std::numeric_limits<double>::epsilon())
    {
        Term2 = (p / Denominator2) * BasisFunctionDerivative(u, i + 1, p - 1, k - 1, Knots);
    }

    return Term1 - Term2;
}

// 自实现 BSplCLib::EvalBsplineBasis（NURBS Book DersBasisFuns 算法）
void CurveFair::EvalBsplineBasis(
    int maxDeriv, int deg, const std::vector<double>& U, double u, int& firstIndex,
    Eigen::MatrixXd& basis)
{
    int m = (int)U.size() - 1;
    int n = m - deg - 1;

    // FindSpan
    int span;
    if (u >= U[n + 1] - 1e-15) span = n;
    else if (u <= U[deg] + 1e-15) span = deg;
    else
    {
        int low = deg, high = n + 1;
        int mid = (low + high) / 2;
        while (u < U[mid] || u >= U[mid + 1])
        {
            if (u < U[mid]) high = mid;
            else low = mid;
            mid = (low + high) / 2;
        }
        span = mid;
    }
    firstIndex = span - deg + 1;  // 1-indexed（与 OCC 一致）

    // DersBasisFuns
    std::vector<std::vector<double>> ndu(deg + 1, std::vector<double>(deg + 1, 0.0));
    ndu[0][0] = 1.0;
    std::vector<double> left(deg + 1), right(deg + 1);
    for (int j = 1; j <= deg; j++)
    {
        left[j] = u - U[span + 1 - j];
        right[j] = U[span + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; r++)
        {
            ndu[j][r] = right[r + 1] + left[j - r];
            double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    std::vector<std::vector<double>> ders(maxDeriv + 1, std::vector<double>(deg + 1, 0.0));
    for (int j = 0; j <= deg; j++) ders[0][j] = ndu[j][deg];

    std::vector<std::vector<double>> a(2, std::vector<double>(deg + 1, 0.0));
    for (int r = 0; r <= deg; r++)
    {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= maxDeriv; k++)
        {
            double d = 0.0;
            int rk = r - k;
            int pk = deg - k;
            if (r >= k)
            {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r - 1 <= pk) ? k - 1 : deg - r;
            for (int j = j1; j <= j2; j++)
            {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk)
            {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::swap(s1, s2);
        }
    }
    int rfact = deg;
    for (int k = 1; k <= maxDeriv; k++)
    {
        for (int j = 0; j <= deg; j++) ders[k][j] *= rfact;
        rfact *= (deg - k);
    }

    basis.resize(maxDeriv + 1, deg + 1);
    basis.setZero();
    for (int k = 0; k <= maxDeriv; k++)
        for (int j = 0; j <= deg; j++)
            basis(k, j) = ders[k][j];
}

sggk::BSplineCurve3DPtr CurveFair::CreateNewBSplineCurve(
    const sggk::BSplineCurve3DPtr& theOriginalCurve,
    const Eigen::MatrixXd& newD)
{
    int degree = theOriginalCurve->Degree();
    sggk::RealArray knots = theOriginalCurve->Knots();
    sggk::UIntArray multiplicities = theOriginalCurve->Mults();

    sggk::Point3DArray poles(newD.rows());
    for (int i = 0; i < newD.rows(); i++)
    {
        poles[i] = sggk::Point3D(newD(i, 0), newD(i, 1), newD(i, 2));
    }

    sggk::BSplineCurve3DPtr newCurve;
    try
    {
        newCurve = std::make_shared<sggk::BSplineCurve3D>(degree, poles, knots, multiplicities);
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

double CurveFair::GetCurveCurveHausdorffDistance(
    const sggk::BSplineCurve3DPtr theOriginalCurve,
    const sggk::BSplineCurve3DPtr theOperateCurve)
{
    double firstParameter = theOriginalCurve->MinParam();
    double lastParameter = theOriginalCurve->MaxParam();
    const int sampleNum = 100;
    const double step = (lastParameter - firstParameter) / sampleNum;
    double hausdorffDistanceResult = INT_MIN;
    for (double t = firstParameter; t <= lastParameter; t += step)
    {
        sggk::Point3D aPntOnOriginalCurve = theOriginalCurve->CalcPoint(t);
        double param;
        sggk::Point3D projected = theOperateCurve->CalcNearestPoint(aPntOnOriginalCurve, param);
        double minDistance = aPntOnOriginalCurve.DistanceTo(projected);
        hausdorffDistanceResult = std::max(hausdorffDistanceResult, minDistance);
    }
    return hausdorffDistanceResult;
}

double CurveFair::GetFitPointsCurveHausdorffDistance(
    const std::vector<std::pair<sggk::Point3D, double>> theFitPointParams,
    const sggk::BSplineCurve3DPtr& theOperateCurve)
{
    sggk::BSplineCurve3DPtr mappingCurve = GetArclengthParameterMapping(theOperateCurve);
    double hausdorffDistance = -1.0;
    std::vector<sggk::Point3D> pointsOnCurve;
    for (auto p : theFitPointParams)
    {
        sggk::Point3D fitPoint = p.first;
        double s = p.second;
        double t = mappingCurve->CalcPoint(s).X();
        sggk::Point3D pointOnCurve = theOperateCurve->CalcPoint(t);
        pointsOnCurve.push_back(pointOnCurve);
        double distance = fitPoint.DistanceTo(pointOnCurve);
        hausdorffDistance = std::max(hausdorffDistance, distance);
    }
    return hausdorffDistance;
}

double CurveFair::GetPointCurveHausdorffDistance(
    const std::vector<sggk::Point3D>& thePoints,
    const sggk::BSplineCurve3DPtr& theCurve)
{
    if (thePoints.empty() || !theCurve)
        return -1.0;

    double hausdorffDistance = -1.0;

    for (const auto& pt : thePoints)
    {
        double param;
        sggk::Point3D projected = theCurve->CalcNearestPoint(pt, param);
        double dist = pt.DistanceTo(projected);
        hausdorffDistance = std::max(hausdorffDistance, dist);
    }

    return hausdorffDistance;
}

double CurveFair::f_inverse(const sggk::BSplineCurve3DPtr& theBSplineCurve, double t)
{
    double leftParam = theBSplineCurve->MinParam();
    double rightParam = theBSplineCurve->MaxParam();
    if (t == leftParam) return leftParam;
    if (t == rightParam) return rightParam;
    while (true)
    {
        double mid = leftParam + (rightParam - leftParam) / 2;
        double t_mid = f(mid);
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

std::vector<std::pair<sggk::Point3D, double>> CurveFair::ReCalculateFitPointParameters(
    const std::vector<sggk::Point3D>& theFitPoints,
    const sggk::BSplineCurve3DPtr& theCurve
)
{
    std::vector<std::pair<sggk::Point3D, double>> result;
    if (theFitPoints.empty())
    {
        return result;
    }

    if (theCurve)
    {
        double first = theCurve->MinParam();
        double last = theCurve->MaxParam();
        double totalLength = ComputeCurveLength(theCurve);

        for (const auto& pt : theFitPoints)
        {
            double u = 0;
            sggk::Point3D projected = theCurve->CalcNearestPoint(pt, u);
            double partialLength = ComputeCurveLengthBetweenParameters(theCurve, first, u);
            double normalized = (partialLength / totalLength) * (last - first) + first;
            result.emplace_back(pt, normalized);
        }
    }
    else
    {
        double totalLength = 0.0;
        std::vector<double> cumulativeLengths(theFitPoints.size(), 0.0);
        for (size_t i = 1; i < theFitPoints.size(); ++i)
        {
            double d = theFitPoints[i].DistanceTo(theFitPoints[i - 1]);
            totalLength += d;
            cumulativeLengths[i] = cumulativeLengths[i - 1] + d;
        }

        for (size_t i = 0; i < theFitPoints.size(); ++i)
        {
            double ratio = (totalLength > 0) ? cumulativeLengths[i] / totalLength : 0.0;
            result.emplace_back(theFitPoints[i], ratio);
        }
    }

    return result;
}

double CurveFair::ComputeCurveLengthBetweenParameters(const sggk::BSplineCurve3DPtr& theCurve, double theParameter1, double theParameter2)
{
    if (!theCurve) return -1.0;
    double len = 0.0;
    for (int i = 0; i < 30; ++i)
    {
        double u = (theParameter2 - theParameter1) * GAUSS_PTS[i] / 2.0 + (theParameter1 + theParameter2) / 2.0;
        sggk::Vector3D d1 = theCurve->CalcDeriv1(u);
        len += GAUSS_WTS[i] * d1.Length();
    }
    return len * (theParameter2 - theParameter1) / 2.0;
}

double CurveFair::ComputeCurveLengthBetweenParameters(const sggk::BSplineCurve3DPtr& theCurve, sggk::Point3D thePnt1, sggk::Point3D thePnt2)
{
    double theParameter1 = GetPntParameterOnCurve(theCurve, thePnt1);
    double theParameter2 = GetPntParameterOnCurve(theCurve, thePnt2);

    if (theParameter1 > theParameter2)
    {
        std::swap(theParameter1, theParameter2);
    }
    return ComputeCurveLengthBetweenParameters(theCurve, theParameter1, theParameter2);
}

double CurveFair::GetPntParameterOnCurve(const sggk::BSplineCurve3DPtr& theCurve, const sggk::Point3D& thePoint)
{
    double param = -1.0;
    theCurve->CalcNearestPoint(thePoint, param);
    return param;
}

std::pair<std::vector<sggk::Point3D>, std::vector<sggk::Point3D>> CurveFair::SampleCurveWithArclengthMapping(const sggk::BSplineCurve3DPtr& theCurve,
    const int nSamples)
{
    std::vector<sggk::Point3D> originPnts;
    std::vector<sggk::Point3D> mappingPnts;

    sggk::BSplineCurve3DPtr mappedCurve = CurveFair::GetArclengthParameterMapping(theCurve);
    if (!mappedCurve)
    {
        std::cerr << "Mapping failed." << std::endl;
    }

    double uStart = theCurve->MinParam();
    double uEnd = theCurve->MaxParam();
    std::vector<double> originLengths;
    std::vector<double> mappingLengths;

    for (int i = 0; i < nSamples; ++i)
    {
        double u = uStart + (uEnd - uStart) * i / (nSamples - 1);

        originLengths.push_back(ComputeCurveLengthBetweenParameters(m_OriginalCurve, uStart, u));
        sggk::Point3D pOriginal = theCurve->CalcPoint(u);

        double uMapped = mappedCurve->CalcPoint(u).X();
        mappingLengths.push_back(ComputeCurveLengthBetweenParameters(m_OriginalCurve, uStart, uMapped));
        sggk::Point3D pMapped = theCurve->CalcPoint(uMapped);

        originPnts.push_back(pOriginal);
        mappingPnts.push_back(pMapped);
    }
    return std::make_pair(originPnts, mappingPnts);
}

sggk::BSplineCurve3DPtr CurveFair::GetTempFairCurve(
    const sggk::BSplineCurve3DPtr& theCurve,
    Eigen::MatrixXd M,
    Eigen::MatrixXd V,
    Eigen::MatrixXd& D0,
    Eigen::MatrixXd& D,
    double theAlpha)
{
    int n = D0.rows();

    Eigen::Vector3d D0_start = D0.row(0);
    Eigen::Vector3d D0_end = D0.row(n - 1);

    Eigen::MatrixXd D0_internal = D0.block(1, 0, n - 2, 3);
    Eigen::VectorXd D0_internal_x = D0.col(0).segment(1, n - 2);
    Eigen::VectorXd D0_internal_y = D0.col(1).segment(1, n - 2);
    Eigen::VectorXd D0_internal_z = D0.col(2).segment(1, n - 2);

    Eigen::MatrixXd M_internal = M.block(1, 1, n - 2, n - 2);
    Eigen::MatrixXd V_internal = V.block(1, 1, n - 2, n - 2);

    int aCompouteTimes = 0;
    sggk::BSplineCurve3DPtr aResultCurve = nullptr;

    while (true)
    {
        aCompouteTimes++;
        Eigen::MatrixXd A = theAlpha * M_internal + V_internal;

        Eigen::VectorXd b_x = Eigen::VectorXd::Zero(n - 2);
        Eigen::VectorXd b_y = Eigen::VectorXd::Zero(n - 2);
        Eigen::VectorXd b_z = Eigen::VectorXd::Zero(n - 2);

        for (int k = 0; k < n - 2; k++)
        {
            double D0_term_x = M(0, k + 1) * D0_start.x();
            double D0_term_y = M(0, k + 1) * D0_start.y();
            double D0_term_z = M(0, k + 1) * D0_start.z();
            double Dn_term_x = M(n - 1, k + 1) * D0_end.x();
            double Dn_term_y = M(n - 1, k + 1) * D0_end.y();
            double Dn_term_z = M(n - 1, k + 1) * D0_end.z();

            double reg_term_x = V_internal(k, k) * D0_internal_x(k);
            double reg_term_y = V_internal(k, k) * D0_internal_y(k);
            double reg_term_z = V_internal(k, k) * D0_internal_z(k);

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

        Eigen::MatrixXd D_internal(n - 2, 3);
        D_internal.col(0) = D_internal_x;
        D_internal.col(1) = D_internal_y;
        D_internal.col(2) = D_internal_z;

        D = Eigen::MatrixXd::Zero(n, 3);
        D.row(0) = D0_start;
        D.block(1, 0, n - 2, 3) = D_internal;
        D.row(n - 1) = D0_end;

        aResultCurve = CreateNewBSplineCurve(m_OriginalCurve, D);

        double aHausdorffDistance = 0;
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
    const sggk::BSplineCurve3DPtr& theBSplineCurve,
    const int p)
{
    const int maxDerivate = 3;
    const int aNum = (int)theBSplineCurve->ControlPoints().size();
    const sggk::RealArray& aKnots = theBSplineCurve->Knots();
    int aDeg = theBSplineCurve->Degree();
    std::vector<double> KnotSeqVector = GetKnotSequence(theBSplineCurve);
    std::vector<double> KnotVector = aKnots;

    std::vector<double> internalKnots(KnotVector.begin() + 1, KnotVector.end() - 1);

    int m = (int)internalKnots.size();
    Eigen::MatrixXd C(m, aNum);
    C.setZero();

    for (int j = 0; j < m; ++j)
    {
        double u_j = internalKnots[j];
        double s_j = f_inverse(theBSplineCurve, u_j);
        double Left_uj = u_j - 1e-7;
        double Right_uj = u_j + 1e-7;
        double f0 = f(s_j, 0);
        double f1 = f(s_j, 1);
        double f2 = f(s_j, 2);
        double f3 = f(s_j, 3);

        Eigen::MatrixXd leftBasis, rightBasis;
        int leftFirstIndex, rightFirstIndex;

        EvalBsplineBasis(maxDerivate, aDeg, KnotSeqVector, Left_uj, leftFirstIndex, leftBasis);
        EvalBsplineBasis(maxDerivate, aDeg, KnotSeqVector, Right_uj, rightFirstIndex, rightBasis);

        for (int localIdx = 0; localIdx < aDeg + 1; ++localIdx)
        {
            int leftGlobalIdx = leftFirstIndex + localIdx - 1;
            int rightGlobalIdx = rightFirstIndex + localIdx - 1;

            if (leftGlobalIdx >= 0 && leftGlobalIdx < aNum)
            {
                double d1Ni_df1 = leftBasis(1, localIdx);
                double d2Ni_df2 = leftBasis(2, localIdx);
                double d3Ni_df3 = leftBasis(3, localIdx);

                double D3Ni_Left = d3Ni_df3 * std::pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;
                C(j, leftGlobalIdx) += D3Ni_Left;
            }
            if (rightGlobalIdx >= 0 && rightGlobalIdx < aNum)
            {
                double d1Ni_df1 = rightBasis(1, localIdx);
                double d2Ni_df2 = rightBasis(2, localIdx);
                double d3Ni_df3 = rightBasis(3, localIdx);

                double D3Ni_Right = d3Ni_df3 * std::pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;
                C(j, rightGlobalIdx) -= D3Ni_Right;
            }
        }
    }
    return C;
}

Eigen::MatrixXd CurveFair::ComputeEnergyMatrix(
    const sggk::BSplineCurve3DPtr& theBSplineCurve,
    const int p,
    const double tol)
{
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve);
    const int maxDerivate = 3;
    const int aGaussNum = 30;
    int aNum = (int)theBSplineCurve->ControlPoints().size();
    const sggk::RealArray& aKnots = theBSplineCurve->Knots();
    int aDeg = theBSplineCurve->Degree();
    std::vector<double> aKnotSeqVector = GetKnotSequence(theBSplineCurve);
    std::vector<double> aKnotVector = aKnots;
    Eigen::MatrixXd M(aNum, aNum);
    M.setZero();

    for (size_t i = 0; i + 1 < aKnots.size(); ++i)
    {
        double sStart = f_inverse(theBSplineCurve, aKnots[i]);
        double sEnd = f_inverse(theBSplineCurve, aKnots[i + 1]);
        for (int GaussIndex = 0; GaussIndex < aGaussNum; ++GaussIndex)
        {
            double s = (sEnd - sStart) * GAUSS_PTS[GaussIndex] / 2.0 + (sStart + sEnd) / 2.0;
            int aFirstIndex;
            Eigen::MatrixXd aBsplineBasis;
            EvalBsplineBasis(maxDerivate, aDeg, aKnotSeqVector, f(s), aFirstIndex, aBsplineBasis);

            double f0 = f(s, 0);
            double f1 = f(s, 1);
            double f2 = f(s, 2);
            double f3 = f(s, 3);
            for (int m = 0; m < aDeg + 1; ++m)
            {
                int globalI = m + aFirstIndex - 1;
                double d1Ni_df1 = aBsplineBasis(1, m);
                double d2Ni_df2 = aBsplineBasis(2, m);
                double d3Ni_df3 = aBsplineBasis(3, m);

                double D3Ni_Dt3 = d3Ni_df3 * std::pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;

                for (int n = 0; n < aDeg + 1; ++n)
                {
                    int globalJ = n + aFirstIndex - 1;
                    double d1Nj_df1 = aBsplineBasis(1, n);
                    double d2Nj_df2 = aBsplineBasis(2, n);
                    double d3Nj_df3 = aBsplineBasis(3, n);

                    double D3Nj_Dt3 = d3Nj_df3 * std::pow(f1, 3) + 3 * d2Nj_df2 * f1 * f2 + d1Nj_df1 * f3;

                    double aElement = D3Ni_Dt3 * D3Nj_Dt3 * GAUSS_WTS[GaussIndex] * (sEnd - sStart) / 2.0;
                    M(globalI, globalJ) += aElement;
                }
            }
        }
    }
    return M;
}

Eigen::SparseMatrix<double> CurveFair::ComputeSparseEnergyMatrix(
    const sggk::BSplineCurve3DPtr& theBSplineCurve,
    const int p,
    const double theTolerance)
{
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(theBSplineCurve);
    const int maxDerivate = 3;
    const int aGaussNum = 30;
    int aNum = (int)theBSplineCurve->ControlPoints().size();
    const sggk::RealArray& aKnots = theBSplineCurve->Knots();
    int aDeg = theBSplineCurve->Degree();
    std::vector<double> aKnotSeqVector = GetKnotSequence(theBSplineCurve);

    std::vector<Eigen::Triplet<double>> triplets;

    for (size_t i = 0; i + 1 < aKnots.size(); ++i)
    {
        double sStart = f_inverse(theBSplineCurve, aKnots[i]);
        double sEnd = f_inverse(theBSplineCurve, aKnots[i + 1]);

        for (int GaussIndex = 0; GaussIndex < aGaussNum; ++GaussIndex)
        {
            double s = (sEnd - sStart) * GAUSS_PTS[GaussIndex] / 2.0 + (sStart + sEnd) / 2.0;

            int aFirstIndex;
            Eigen::MatrixXd aBsplineBasis;
            EvalBsplineBasis(maxDerivate, aDeg, aKnotSeqVector, f(s), aFirstIndex, aBsplineBasis);

            double f1 = f(s, 1);
            double f2 = f(s, 2);
            double f3 = f(s, 3);

            for (int m = 0; m < aDeg + 1; ++m)
            {
                int globalI = m + aFirstIndex - 1;
                double d1Ni_df1 = aBsplineBasis(1, m);
                double d2Ni_df2 = aBsplineBasis(2, m);
                double d3Ni_df3 = aBsplineBasis(3, m);

                double D3Ni_Dt3 = d3Ni_df3 * std::pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;

                for (int n = 0; n < aDeg + 1; ++n)
                {
                    int globalJ = n + aFirstIndex - 1;
                    double d1Nj_df1 = aBsplineBasis(1, n);
                    double d2Nj_df2 = aBsplineBasis(2, n);
                    double d3Nj_df3 = aBsplineBasis(3, n);

                    double D3Nj_Dt3 = d3Nj_df3 * std::pow(f1, 3) + 3 * d2Nj_df2 * f1 * f2 + d1Nj_df1 * f3;

                    double aElement = D3Ni_Dt3 * D3Nj_Dt3 * GAUSS_WTS[GaussIndex] * (sEnd - sStart) / 2.0;

                    if (std::abs(aElement) > theTolerance)
                        triplets.emplace_back(globalI, globalJ, aElement);
                }
            }
        }
    }

    Eigen::SparseMatrix<double> M(aNum, aNum);
    M.setFromTriplets(triplets.begin(), triplets.end());
    return M;
}

double CurveFair::ComputeCurveFairEnergy(const sggk::BSplineCurve3DPtr& theCurve)
{
    const sggk::Point3DArray& Poles = theCurve->ControlPoints();
    const int aNumPoles = (int)Poles.size();

    Eigen::MatrixXd M = ComputeEnergyMatrix(theCurve, 3);
    std::cout << M << std::endl;
    Eigen::MatrixXd P(aNumPoles, 3);
    for (int i = 0; i < aNumPoles; ++i)
    {
        P(i, 0) = Poles[i].X();
        P(i, 1) = Poles[i].Y();
        P(i, 2) = Poles[i].Z();
    }

    return (P.transpose() * M * P).trace();
}

Eigen::MatrixXd CurveFair::ComputeConstraintMatrix(
    const Eigen::MatrixXd& theD0,
    Eigen::VectorXd& theH)
{
    int n = theD0.rows();
    int num_unknowns = n - 2;

    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(4, 3 * num_unknowns);
    theH = Eigen::VectorXd::Zero(4);

    if (num_unknowns <= 0)
    {
        return C;
    }

    Eigen::Vector3d P0 = theD0.row(0);
    Eigen::Vector3d P1 = theD0.row(1);
    Eigen::Vector3d Pn_1 = theD0.row(n - 2);
    Eigen::Vector3d Pn = theD0.row(n - 1);

    Eigen::Vector3d T_start_0 = P1 - P0;
    if (T_start_0.norm() > 1e-9)
    {
        Eigen::Vector3d v_arb = (std::abs(T_start_0.x()) < 0.9) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d u1 = T_start_0.cross(v_arb).normalized();
        Eigen::Vector3d u2 = T_start_0.cross(u1).normalized();

        C(0, 0) = u1.x(); C(0, num_unknowns) = u1.y(); C(0, 2 * num_unknowns) = u1.z();
        C(1, 0) = u2.x(); C(1, num_unknowns) = u2.y(); C(1, 2 * num_unknowns) = u2.z();

        theH(0) = P0.dot(u1);
        theH(1) = P0.dot(u2);
    }

    Eigen::Vector3d T_end_0 = Pn - Pn_1;
    if (T_end_0.norm() > 1e-9)
    {
        Eigen::Vector3d v_arb = (std::abs(T_end_0.x()) < 0.9) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d v1 = T_end_0.cross(v_arb).normalized();
        Eigen::Vector3d v2 = T_end_0.cross(v1).normalized();

        int last_idx = num_unknowns - 1;
        C(2, last_idx) = -v1.x(); C(2, last_idx + num_unknowns) = -v1.y(); C(2, last_idx + 2 * num_unknowns) = -v1.z();
        C(3, last_idx) = -v2.x(); C(3, last_idx + num_unknowns) = -v2.y(); C(3, last_idx + 2 * num_unknowns) = -v2.z();

        theH(2) = -Pn.dot(v1);
        theH(3) = -Pn.dot(v2);
    }

    return C;
}

Eigen::SparseMatrix<double> CurveFair::ComputeSparseConstraintMatrix(
    const Eigen::MatrixXd& theD0,
    Eigen::VectorXd& theH)
{
    int n = theD0.rows();
    int num_unknowns = n - 2;

    Eigen::SparseMatrix<double> C(4, 3 * num_unknowns);
    theH = Eigen::VectorXd::Zero(4);

    if (num_unknowns <= 0)
    {
        return C;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(12);

    Eigen::Vector3d P0 = theD0.row(0);
    Eigen::Vector3d P1 = theD0.row(1);
    Eigen::Vector3d Pn_1 = theD0.row(n - 2);
    Eigen::Vector3d Pn = theD0.row(n - 1);

    Eigen::Vector3d T_start_0 = P1 - P0;
    if (T_start_0.norm() > 1e-9)
    {
        Eigen::Vector3d v_arb = (std::abs(T_start_0.x()) < 0.9) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d u1 = T_start_0.cross(v_arb).normalized();
        Eigen::Vector3d u2 = T_start_0.cross(u1).normalized();

        triplets.emplace_back(0, 0, u1.x());
        triplets.emplace_back(0, num_unknowns, u1.y());
        triplets.emplace_back(0, 2 * num_unknowns, u1.z());

        triplets.emplace_back(1, 0, u2.x());
        triplets.emplace_back(1, num_unknowns, u2.y());
        triplets.emplace_back(1, 2 * num_unknowns, u2.z());

        theH(0) = P0.dot(u1);
        theH(1) = P0.dot(u2);
    }

    Eigen::Vector3d T_end_0 = Pn - Pn_1;
    if (T_end_0.norm() > 1e-9)
    {
        Eigen::Vector3d v_arb = (std::abs(T_end_0.x()) < 0.9) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
        Eigen::Vector3d v1 = T_end_0.cross(v_arb).normalized();
        Eigen::Vector3d v2 = T_end_0.cross(v1).normalized();

        int last_idx = num_unknowns - 1;

        triplets.emplace_back(2, last_idx, -v1.x());
        triplets.emplace_back(2, last_idx + num_unknowns, -v1.y());
        triplets.emplace_back(2, last_idx + 2 * num_unknowns, -v1.z());

        triplets.emplace_back(3, last_idx, -v2.x());
        triplets.emplace_back(3, last_idx + num_unknowns, -v2.y());
        triplets.emplace_back(3, last_idx + 2 * num_unknowns, -v2.z());

        theH(2) = -Pn.dot(v1);
        theH(3) = -Pn.dot(v2);
    }

    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}

sggk::BSplineCurve3DPtr CurveFair::GetTempFairCurveWithTangentConstraint(
    const sggk::BSplineCurve3DPtr& theCurve,
    Eigen::MatrixXd M,
    Eigen::MatrixXd V,
    Eigen::MatrixXd& D0,
    Eigen::MatrixXd& D,
    double theAlpha)
{
    int n = D0.rows();
    if (n <= 2) return nullptr;

    int num_unknowns = n - 2;
    Eigen::Vector3d D0_start = D0.row(0);
    Eigen::Vector3d D0_end = D0.row(n - 1);

    Eigen::VectorXd D0_internal_x = D0.col(0).segment(1, num_unknowns);
    Eigen::VectorXd D0_internal_y = D0.col(1).segment(1, num_unknowns);
    Eigen::VectorXd D0_internal_z = D0.col(2).segment(1, num_unknowns);

    Eigen::MatrixXd M_internal = M.block(1, 1, num_unknowns, num_unknowns);
    Eigen::MatrixXd V_internal = V.block(1, 1, num_unknowns, num_unknowns);

    Eigen::VectorXd h;
    Eigen::MatrixXd C_mat = ComputeConstraintMatrix(D0, h);

    int aCompouteTimes = 0;
    sggk::BSplineCurve3DPtr aResultCurve = nullptr;

    while (true)
    {
        aCompouteTimes++;

        Eigen::MatrixXd A = theAlpha * M_internal + V_internal;

        Eigen::VectorXd b_x(num_unknowns), b_y(num_unknowns), b_z(num_unknowns);
        for (int k = 0; k < num_unknowns; k++)
        {
            double D0_term_x = M(0, k + 1) * D0_start.x();
            double Dn_term_x = M(n - 1, k + 1) * D0_end.x();
            double reg_term_x = V_internal(k, k) * D0_internal_x(k);
            b_x(k) = -theAlpha * (D0_term_x + Dn_term_x) + reg_term_x;

            double D0_term_y = M(0, k + 1) * D0_start.y();
            double Dn_term_y = M(n - 1, k + 1) * D0_end.y();
            double reg_term_y = V_internal(k, k) * D0_internal_y(k);
            b_y(k) = -theAlpha * (D0_term_y + Dn_term_y) + reg_term_y;

            double D0_term_z = M(0, k + 1) * D0_start.z();
            double Dn_term_z = M(n - 1, k + 1) * D0_end.z();
            double reg_term_z = V_internal(k, k) * D0_internal_z(k);
            b_z(k) = -theAlpha * (D0_term_z + Dn_term_z) + reg_term_z;
        }

        int total_vars = 3 * num_unknowns;
        int system_size = total_vars + 4;

        Eigen::SparseMatrix<double> BigLHS_sparse(system_size, system_size);
        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(system_size * 10);

        for (int i = 0; i < num_unknowns; i++)
        {
            for (int j = 0; j < num_unknowns; j++)
            {
                if (std::abs(A(i, j)) > 1e-12)
                {
                    tripletList.push_back(Eigen::Triplet<double>(i, j, A(i, j)));
                    tripletList.push_back(Eigen::Triplet<double>(i + num_unknowns, j + num_unknowns, A(i, j)));
                    tripletList.push_back(Eigen::Triplet<double>(i + 2 * num_unknowns, j + 2 * num_unknowns, A(i, j)));
                }
            }
        }

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < total_vars; j++)
            {
                if (std::abs(C_mat(i, j)) > 1e-12)
                {
                    tripletList.push_back(Eigen::Triplet<double>(total_vars + i, j, C_mat(i, j)));
                    tripletList.push_back(Eigen::Triplet<double>(j, total_vars + i, C_mat(i, j)));
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
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
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

        double aHausdorffDistance = 0;
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

sggk::BSplineCurve3DPtr CurveFair::GetTempContinuniousFairCurve(
    const sggk::BSplineCurve3DPtr& theCurve,
    Eigen::MatrixXd M,
    Eigen::MatrixXd C,
    Eigen::MatrixXd V,
    Eigen::MatrixXd& D0,
    Eigen::MatrixXd& D,
    double theAlpha)
{
    int n = D0.rows();
    int m = C.rows();
    int internal_n = n - 2;
    Eigen::Vector3d D0_start = D0.row(0);
    Eigen::Vector3d D0_end = D0.row(n - 1);

    Eigen::MatrixXd D0_internal = D0.block(1, 0, internal_n, 3);
    Eigen::VectorXd D0_internal_x = D0.col(0).segment(1, internal_n);
    Eigen::VectorXd D0_internal_y = D0.col(1).segment(1, internal_n);
    Eigen::VectorXd D0_internal_z = D0.col(2).segment(1, internal_n);

    Eigen::MatrixXd M_internal = M.block(1, 1, internal_n, internal_n);
    Eigen::MatrixXd V_internal = V.block(1, 1, internal_n, internal_n);
    Eigen::MatrixXd C_internal = C.block(0, 1, m, internal_n);

    Eigen::VectorXd h_x(m);
    Eigen::VectorXd h_y(m);
    Eigen::VectorXd h_z(m);
    for (int j = 0; j < m; ++j)
    {
        h_x(j) = -C(j, 0) * D0_start.x() - C(j, n - 1) * D0_end.x();
        h_y(j) = -C(j, 0) * D0_start.y() - C(j, n - 1) * D0_end.y();
        h_z(j) = -C(j, 0) * D0_start.z() - C(j, n - 1) * D0_end.z();
    }

    Eigen::MatrixXd extendedA(internal_n + m, internal_n + m);
    extendedA.topRightCorner(internal_n, m) = C_internal.transpose();
    extendedA.bottomLeftCorner(m, internal_n) = C_internal;
    extendedA.bottomRightCorner(m, m).setZero();
    int aCompouteTimes = 0;
    sggk::BSplineCurve3DPtr aResultCurve;
    while (++aCompouteTimes <= 100)
    {
        Eigen::MatrixXd A = theAlpha * M_internal + V_internal;
        extendedA.topLeftCorner(internal_n, internal_n) = A;

        Eigen::VectorXd b_x = Eigen::VectorXd::Zero(internal_n);
        Eigen::VectorXd b_y = Eigen::VectorXd::Zero(internal_n);
        Eigen::VectorXd b_z = Eigen::VectorXd::Zero(internal_n);

        for (int k = 0; k < n - 2; k++)
        {
            double D0_term_x = M(0, k + 1) * D0_start.x();
            double D0_term_y = M(0, k + 1) * D0_start.y();
            double D0_term_z = M(0, k + 1) * D0_start.z();
            double Dn_term_x = M(n - 1, k + 1) * D0_end.x();
            double Dn_term_y = M(n - 1, k + 1) * D0_end.y();
            double Dn_term_z = M(n - 1, k + 1) * D0_end.z();

            double reg_term_x = V_internal(k, k) * D0_internal_x(k);
            double reg_term_y = V_internal(k, k) * D0_internal_y(k);
            double reg_term_z = V_internal(k, k) * D0_internal_z(k);

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

        Eigen::MatrixXd D_internal(n - 2, 3);
        D_internal.col(0) = D_internal_x;
        D_internal.col(1) = D_internal_y;
        D_internal.col(2) = D_internal_z;

        D = Eigen::MatrixXd::Zero(n, 3);
        D.row(0) = D0_start;
        D.block(1, 0, n - 2, 3) = D_internal;
        D.row(n - 1) = D0_end;

        aResultCurve = CreateNewBSplineCurve(m_OriginalCurve, D);

        double aHausdorffDistance = 0;
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
    int n = D0.rows();
    double stepSize = 10;
    int aStepCnt = 0;
    while (true)
    {
        D0 += (D - D0) / 10;
        aStepCnt++;

        sggk::BSplineCurve3DPtr aCurve = CreateNewBSplineCurve(m_OriginalCurve, D0);
        M = ComputeEnergyMatrix(aCurve, m_OriginalCurve->Degree());
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

sggk::BSplineCurve3DPtr CurveFair::SampleAndFitBSpline(
    const sggk::BSplineCurve3DPtr& theOriginalCurve,
    int theSampleNum,
    std::vector<sggk::Point3D>& theFitPoints,
    int theMaxDegree,
    int theContinuity,
    double theTolerance)
{
    if (theSampleNum < 2 || !theOriginalCurve)
    {
        std::cerr << "Invalid input: too few points or null curve." << std::endl;
        return nullptr;
    }

    double totalLength = ComputeCurveLength(theOriginalCurve);
    double firstParam = theOriginalCurve->MinParam();

    std::vector<sggk::Point3D> sampledPoints(theSampleNum);
    for (int i = 0; i < theSampleNum; ++i)
    {
        double targetLength = (totalLength * i) / (theSampleNum - 1);
        sggk::Point3D p;
        double u = theOriginalCurve->CalcParaByLength(targetLength, firstParam, p);
        sampledPoints[i] = theOriginalCurve->CalcPoint(u);
    }
    theFitPoints = sampledPoints;

    sggk::BSCrv3DInterpolationOpts opts;
    opts.degree = theMaxDegree;
    sggk::BSCrvFitting3DResult res = sggk::BSCrvFitting::Interpolation3D(sampledPoints, opts);
    return res.curve;
}

std::vector<sggk::Point3D> CurveFair::SampleCurveWithArcLength(const sggk::BSplineCurve3DPtr& bsplineCurve, int numSamples)
{
    std::vector<sggk::Point3D> sampledPointsWithArcLength;
    if (!bsplineCurve || numSamples <= 0) return sampledPointsWithArcLength;

    double totalLength = ComputeCurveLength(bsplineCurve);
    double firstParam = bsplineCurve->MinParam();

    for (int i = 0; i < numSamples; ++i)
    {
        double targetLength = (totalLength * i) / (numSamples - 1);
        sggk::Point3D p;
        double u = bsplineCurve->CalcParaByLength(targetLength, firstParam, p);
        sggk::Point3D point = bsplineCurve->CalcPoint(u);
        sampledPointsWithArcLength.emplace_back(point);
    }

    return sampledPointsWithArcLength;
}

sggk::BSplineCurve3DPtr CurveFair::SampleAndFitBSpline(
    std::vector<sggk::Point3D>& theFitPoints,
    double theTolerance,
    bool theResortFlag,
    int theMaxDegree,
    int theContinuity)
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

    sggk::BSCrv3DInterpolationOpts opts;
    opts.degree = theMaxDegree;
    sggk::BSCrvFitting3DResult res = sggk::BSCrvFitting::Interpolation3D(theFitPoints, opts);
    return res.curve;
}

sggk::BSplineCurve3DPtr CurveFair::RefineCurveByCurvatureAuto(
    const sggk::BSplineCurve3DPtr& theCurve,
    const std::vector<double>& myKnotSeq,
    const int baseInsertNum)
{
    if (!theCurve) return theCurve;

    sggk::BSplineCurve3DPtr refinedCurve = CloneBSpline(theCurve);
    const int degree = refinedCurve->Degree();
    const int nPoles = (int)refinedCurve->ControlPoints().size();
    const sggk::RealArray& knotArray = refinedCurve->Knots();

    const int nbKnots = (int)knotArray.size();
    double KnotUpper = knotArray.back();
    double KnotLower = knotArray.front();
    std::vector<double> curvatureRates(nbKnots);
    m_ArcLengthMappingFunction = GetArclengthParameterMapping(refinedCurve);

    for (int i = 0; i < nbKnots; ++i)
    {
        double u = knotArray[i];
        double s = f_inverse(refinedCurve, u);
        double f0 = f(s);
        double f1 = f(s, 1);
        double f2 = f(s, 2);
        double f3 = f(s, 3);

        double sum = 0.0;
        for (int j = 0; j < nPoles; ++j)
        {
            double k = BasisFunctionDerivative(f0, j, degree, 3, myKnotSeq) * std::pow(f1, 3)
                + 3 * BasisFunctionDerivative(f0, j, degree, 2, myKnotSeq) * f1 * f2
                + BasisFunctionDerivative(f0, j, degree, 1, myKnotSeq) * f3;
            sum += k * k;
        }
        curvatureRates[i] = sum;
    }

    double sumCurv = 0.0;
    for (auto c : curvatureRates) sumCurv += c;
    double avgCurv = sumCurv / curvatureRates.size();

    std::vector<double> insertParams;

    for (int i = 0; i < nbKnots - 1; ++i)
    {
        if (curvatureRates[i] > avgCurv)
        {
            double left = knotArray[i];
            double right = knotArray[i + 1];

            if (right - left < (KnotUpper - KnotLower) / 100) continue;

            int insertNum = baseInsertNum;
            if (curvatureRates[i] > 2 * avgCurv)
                insertNum *= 3;

            for (int j = 1; j <= insertNum; ++j)
            {
                double u = left + (right - left) * j / (insertNum + 1);
                insertParams.push_back(u);
            }
        }
    }

    if (!insertParams.empty())
    {
        for (double u : insertParams)
        {
            refinedCurve->InsertKnots(u, 1);
        }
    }

    return refinedCurve;
}

sggk::BSplineCurve3DPtr CurveFair::RefineCurveByFitPoints(
    const sggk::BSplineCurve3DPtr& theCurve,
    const int baseInsertNum)
{
    if (!theCurve || baseInsertNum <= 0)
        return theCurve;

    sggk::BSplineCurve3DPtr refinedCurve = CloneBSpline(theCurve);

    const sggk::RealArray& knotArray = refinedCurve->Knots();
    const int lower = 0;
    const int upper = (int)knotArray.size() - 1;

    constexpr double tol = 1e-12;
    std::map<std::pair<int, int>, int> spanCounts;

    for (const auto& fp : m_FitPointParameters)
    {
        double uFit = fp.second;

        if (uFit < knotArray[lower] - tol || uFit > knotArray[upper] + tol)
            continue;

        int span = -1;
        for (int i = lower; i < upper; ++i)
        {
            if (uFit + tol < knotArray[i]) break;
            if (uFit >= knotArray[i] - tol && uFit < knotArray[i + 1] - tol)
            {
                span = i;
                break;
            }
        }

        if (span == -1) continue;

        double uL = knotArray[span];
        double uR = knotArray[span + 1];
        if (uR - uL < tol) continue;

        spanCounts[{span, span + 1}]++;
    }

    std::vector<double> insertParams;
    for (const auto& entry : spanCounts)
    {
        int span = entry.first.first;
        double uL = knotArray[span];
        double uR = knotArray[span + 1];
        int count = entry.second;

        int actualInsertNum = baseInsertNum * count;
        for (int j = 1; j <= actualInsertNum; ++j)
        {
            double uNew = uL + (uR - uL) * j / (actualInsertNum + 1);
            insertParams.push_back(uNew);
        }
    }

    if (!insertParams.empty())
    {
        std::sort(insertParams.begin(), insertParams.end());
        insertParams.erase(std::unique(insertParams.begin(), insertParams.end(),
            [](double a, double b) { return std::abs(a - b) < 1e-10; }), insertParams.end());

        for (double u : insertParams)
        {
            refinedCurve->InsertKnots(u, 1);
        }
    }

    return refinedCurve;
}

sggk::BSplineCurve3DPtr CurveFair::RefineByTangentAngle(
    const sggk::BSplineCurve3DPtr& theCurve,
    const int theBasicInsertNum)
{
    if (!theCurve) return theCurve;

    sggk::BSplineCurve3DPtr refinedCurve = CloneBSpline(theCurve);
    std::vector<double> knots = refinedCurve->Knots();

    std::vector<double> angles(knots.size(), 0.0);
    for (size_t i = 1; i < knots.size() - 1; ++i)
    {
        sggk::Vector3D leftVec = Dt(refinedCurve, knots[i] - 1e-6, 3);
        sggk::Vector3D rightVec = Dt(refinedCurve, knots[i] + 1e-6, 3);
        angles[i] = leftVec.CalcAngle(rightVec) * 180.0 / M_PI_VAL;
        std::cout << angles[i] << std::endl;
    }

    double avgAngle = 0.0;
    int validCount = 0;
    for (size_t i = 1; i < angles.size() - 1; ++i)
    {
        avgAngle += angles[i];
        ++validCount;
    }

    avgAngle /= std::max(validCount, 1);
    std::map<std::pair<double, double>, int> spanInsertCount;

    for (size_t i = 1; i < knots.size() - 1; ++i)
    {
        int insertNum = theBasicInsertNum + (int)(angles[i] / avgAngle - 1);
        if (insertNum > 0)
        {
            double uLeft = knots[i - 1];
            double uMid = knots[i];
            double uRight = knots[i + 1];

            spanInsertCount[{uLeft, uMid}] += insertNum;
            spanInsertCount[{uMid, uRight}] += insertNum;
        }
    }

    std::vector<double> insertParams;
    for (const auto& kv : spanInsertCount)
    {
        double u1 = kv.first.first;
        double u2 = kv.first.second;
        int count = kv.second;

        for (int i = 1; i <= count; ++i)
        {
            double u = u1 + i * (u2 - u1) / (count + 1);
            insertParams.push_back(u);
        }
    }

    if (insertParams.empty())
        return refinedCurve;

    std::sort(insertParams.begin(), insertParams.end());
    insertParams.erase(std::unique(insertParams.begin(), insertParams.end(),
        [](double a, double b) { return std::abs(a - b) < 1e-10; }), insertParams.end());

    for (double u : insertParams)
    {
        refinedCurve->InsertKnots(u, 1);
    }

    return refinedCurve;
}

std::vector<double> CurveFair::ComputeThirdDerivativeAngles(
    const sggk::BSplineCurve3DPtr& theCurve)
{
    std::vector<double> angleResults;

    if (!theCurve)
        return angleResults;

    sggk::BSplineCurve3DPtr curveCopy = CloneBSpline(theCurve);

    std::vector<double> knots = curveCopy->Knots();

    if (knots.size() < 3)
        return angleResults;

    angleResults.resize(knots.size(), 0.0);

    for (size_t i = 1; i < knots.size() - 1; ++i)
    {
        double u = knots[i];
        sggk::Vector3D leftVec = Dt(curveCopy, u - 1e-6, 3);
        sggk::Vector3D rightVec = Dt(curveCopy, u + 1e-6, 3);

        if (leftVec.Length() > RESOLUTION && rightVec.Length() > RESOLUTION)
        {
            angleResults[i] = leftVec.CalcAngle(rightVec) * 180.0 / M_PI_VAL;
        }
        else
        {
            angleResults[i] = 0.0;
        }
    }

    return angleResults;
}

sggk::BSplineCurve3DPtr CurveFair::InsertKnotsBetweenKnotSpan(
    const sggk::BSplineCurve3DPtr& theCurve,
    double t1,
    double t2,
    int nInsert)
{
    if (!theCurve || nInsert <= 0 || t1 >= t2)
    {
        std::cerr << "Invalid input to InsertKnotsBetween." << std::endl;
        return theCurve;
    }

    sggk::RealArray knots = theCurve->Knots();
    double knotMin = knots.front();
    double knotMax = knots.back();

    if (t1 < knotMin || t2 > knotMax)
    {
        std::cerr << "t1 or t2 out of knot vector range." << std::endl;
        return theCurve;
    }

    std::vector<double> newKnots;
    for (int i = 1; i <= nInsert; ++i)
    {
        double newKnot = t1 + i * (t2 - t1) / (nInsert + 1);
        newKnots.push_back(newKnot);
    }

    sggk::BSplineCurve3DPtr newCurve = CloneBSpline(theCurve);
    for (const auto& u : newKnots)
    {
        newCurve->InsertKnots(u, 1);
    }

    return newCurve;
}

sggk::BSplineCurve3DPtr CurveFair::InsertKnots(
    const sggk::BSplineCurve3DPtr& theCurve,
    const std::vector<double>& params,
    int mult)
{
    if (!theCurve || params.empty() || mult <= 0)
    {
        std::cerr << "Invalid input to InsertKnotsAtParameters." << std::endl;
        return theCurve;
    }

    sggk::BSplineCurve3DPtr newCurve = CloneBSpline(theCurve);

    for (const auto& u : params)
    {
        newCurve->InsertKnots(u, mult);
    }

    return newCurve;
}

sggk::BSplineCurve3DPtr CurveFair::InsertKnot(
    const sggk::BSplineCurve3DPtr& theCurve,
    double u,
    int mult)
{
    if (!theCurve || mult <= 0)
    {
        std::cerr << "Invalid input to InsertSingleKnot." << std::endl;
        return theCurve;
    }

    sggk::BSplineCurve3DPtr newCurve = CloneBSpline(theCurve);
    newCurve->InsertKnots(u, mult);
    return newCurve;
}

sggk::BSplineCurve3DPtr CurveFair::InsertUniformKnotsInAllSpans(
    const sggk::BSplineCurve3DPtr& theCurve,
    int nInsert)
{
    if (!theCurve || nInsert <= 0)
    {
        std::cerr << "Invalid input to InsertUniformKnotsInAllSpans." << std::endl;
        return theCurve;
    }

    sggk::BSplineCurve3DPtr newCurve = CloneBSpline(theCurve);
    sggk::RealArray knotsArray = newCurve->Knots();
    int nbKnots = (int)knotsArray.size();
    double upperKnot = knotsArray.back();
    double lowerKnot = knotsArray.front();
    std::vector<double> newKnots;
    for (int i = 0; i < nbKnots - 1; ++i)
    {
        double t1 = knotsArray[i];
        double t2 = knotsArray[i + 1];

        int inserKnot = nInsert;
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
        for (int j = 1; j <= inserKnot; ++j)
        {
            double knot = t1 + j * (t2 - t1) / (inserKnot + 1);
            newKnots.push_back(knot);
        }
    }

    for (const auto& u : newKnots)
    {
        newCurve->InsertKnots(u, 1);
    }

    return newCurve;
}

void CurveFair::SetControlPointWeightMatrix(const sggk::BSplineCurve3DPtr& theCurve, Eigen::MatrixXd& V)
{
    int nb = (int)theCurve->ControlPoints().size();
    V.resize(nb, nb);
    V.setZero();
    for (int i = 0; i < nb; i++)
    {
        V(i, i) = 1;
    }
    return;
}

std::vector<sggk::Point3D> CurveFair::ReorderPointsNearestNeighbor(const std::vector<sggk::Point3D>& points)
{
    if (points.empty()) return {};

    std::vector<sggk::Point3D> ordered;
    std::vector<bool> visited(points.size(), false);
    ordered.reserve(points.size());

    int currentIdx = 0;
    ordered.push_back(points[currentIdx]);
    visited[currentIdx] = true;

    for (std::size_t i = 1; i < points.size(); ++i)
    {
        double minDist = std::numeric_limits<double>::max();
        int nextIdx = -1;

        for (std::size_t j = 0; j < points.size(); ++j)
        {
            if (!visited[j])
            {
                double d = points[currentIdx].DistanceTo(points[j]);
                if (d < minDist)
                {
                    minDist = d;
                    nextIdx = (int)j;
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
    SetControlPointWeightMatrix(m_OriginalCurve, V);

    const sggk::Point3DArray& Poles = m_OriginalCurve->ControlPoints();
    int n = (int)Poles.size();

    Eigen::MatrixXd D0(n, 3);
    Eigen::MatrixXd D(n, 3);
    for (int i = 0; i < n; i++)
    {
        D0(i, 0) = Poles[i].X();
        D0(i, 1) = Poles[i].Y();
        D0(i, 2) = Poles[i].Z();
    }

    m_ResultCurve = GetTempFairCurveWithTangentConstraint(m_OriginalCurve, M, V, D0, D, m_Alpha);

    Iterator(D0, D);
}
