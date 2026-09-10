// #ifndef CURVEFAIR_H
// #define CURVEFAIR_H

#include <Geometry/3D/Curve/BSplineCurve3D.h>
#include <GeomBase/Point3D.h>
#include <GeomBase/Vector3D.h>

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <mutex>
#include <map>

// 并行计算矩阵
//#include <omp.h>

const double HAUSDORFFDISTANCETOL = 50;
const double ALPHA = 1;
const double ALPHARATIO = 0.5;
const int PARANUM = 100;
const double FITTOLERANCE = 50;
class CurveFair
{
public:
    //! @brief 光顺算法
    //! @note 输入的 theBsplineCurve 和 theFitPoints 不能同时为空, 如果曲线为空, 则根据 theFitPoints 重新拟合曲线
    //! @param [in]theBSplineCurve 输入的曲线
    //! @param [in]theFitPoints 输入的拟合点
    //! @param [in]theHausdorffDistance 距离容差(默认值为50，若型值点非空，则光顺后曲线与型值点距离需在容差内; 否则，最终曲线与原始曲线的距高在容差内)
    //! @param [in] theKnotsInsertFlag 是否自动插入节点
    //! @param [in] theReFitFlag 是否重新采样曲线拟合新曲线
    CurveFair(
        sggk::BSplineCurve3DPtr& theBSplineCurve,
        std::vector<sggk::Point3D>& theFitPoints,
        double theHausdorffDistanceTol = HAUSDORFFDISTANCETOL)
    {
        if (!theBSplineCurve && theFitPoints.size() == 0)
        {
            m_errorCode.push_back("CONSTRUCT ERROR::theBsplineCurve and theFitPoints cannot be both empty.");
            return;
        }

        if (theBSplineCurve)
        {
            // 检查次数
            if (theBSplineCurve->Degree() != 3)
            {
                m_errorCode.push_back("CONSTRUCT ERROR::The BSplineCurve Degree must be 3!");
                return;
            }

            // 检查节点重复性
            const sggk::RealArray& knots = theBSplineCurve->Knots();
            const sggk::UIntArray& mults = theBSplineCurve->Mults();

            // 前后端点重复次数
            int firstMult = mults.front();
            int lastMult = mults.back();

            if (firstMult < (int)theBSplineCurve->Degree() + 1 || lastMult < (int)theBSplineCurve->Degree() + 1)
            {
                m_errorCode.push_back("CONSTRUCT ERROR::The BSplineCurve must have clamped knot vector with first and last knot multiplicity = degree + 1.");
                return;
            }
        }

        m_ArcLengthMappingSampleNum = PARANUM;
        m_Alpha = ALPHA;
        m_AlphaRatio = ALPHARATIO;
        m_HausdorffDistanceTol = theHausdorffDistanceTol;
        m_FitTolerance = theHausdorffDistanceTol;

        // 输入数据初始化
        m_OriginalCurve = theBSplineCurve;
        m_OriginalFitPoints = theFitPoints;

        if (!m_OriginalCurve && m_OriginalFitPoints.size() > 0)
        {
            m_OriginalCurve = SampleAndFitBSpline(m_OriginalFitPoints, FITTOLERANCE);
        }

        UpdateCurveAttribute();
    }

    //! @brief 默认构造函数
    CurveFair()
    {
        this->m_ArcLengthMappingSampleNum = PARANUM;
        this->m_Alpha = ALPHA;
        this->m_AlphaRatio = ALPHARATIO;
    }

public:

    //! @brief 执行光顺算法
    void Perform();

    //! @brief 根据原始曲线重新采样并拟合 B 样条曲线
    //! @param [in] theFitTolerance 拟合容差，默认值为 50
    //! @param [in] theSampleNum 采样点数，默认值为 100
    inline bool ReFitCurve(double theFitTolerance = 50, int theSampleNum = 100)
    {
        if (!m_OriginalCurve)
        {
            m_errorCode.push_back("ERR_NULL_CURVE");
            return false;
        }
        if (theSampleNum <= 2)
        {
            m_errorCode.push_back("ERR_INVALID_SAMPLE_NUM");
            return false;
        }
        if (theFitTolerance <= 0)
        {
            m_errorCode.push_back("ERR_INVALID_TOLERANCE");
            return false;
        }

        try
        {
            m_OriginalFitPoints = SampleCurveWithArcLength(m_OriginalCurve, theSampleNum);
            if (m_OriginalFitPoints.empty())
            {
                m_errorCode.push_back("ERR_SAMPLE_FAILED");
                return false;
            }

            m_OriginalCurve = SampleAndFitBSpline(m_OriginalFitPoints, theFitTolerance);
            if (!m_OriginalCurve)
            {
                m_errorCode.push_back("ERR_FIT_FAILED");
                return false;
            }

            UpdateCurveAttribute();
        }
        catch (const std::exception& e)
        {
            m_errorCode.push_back(std::string("ERR_EXCEPTION: ") + e.what());
            return false;
        }

        return true;
    }

    //! @brief 按切向角变化插入节点细化曲线
    //! @param [in] theInsertLevel 节点插入数量，默认值为 1
    inline bool InsertKnots(int theInsertLevel = 1)
    {
        if (!m_OriginalCurve)
        {
            m_errorCode.push_back("ERR_NULL_CURVE");
            return false;
        }
        if (theInsertLevel <= 0)
        {
            m_errorCode.push_back("ERR_INVALID_INSERT_LEVEL");
            return false;
        }

        try
        {
            m_OriginalCurve = RefineByTangentAngle(m_OriginalCurve, theInsertLevel);
            if (!m_OriginalCurve)
            {
                m_errorCode.push_back("ERR_REFINE_FAILED");
                return false;
            }

            UpdateCurveAttribute();
        }
        catch (const std::exception& e)
        {
            m_errorCode.push_back(std::string("ERR_EXCEPTION: ") + e.what());
            return false;
        }
        return true;
    }

public:
    //! @brief 获取光顺曲线
    inline sggk::BSplineCurve3DPtr GetResult() { return m_ResultCurve; }

    //! @brief 是否有错误信息
    inline bool HasErrorMessage() const
    {
        return m_errorCode.size() > 0;
    }

    //! @brief 获取错误信息
    inline std::vector<std::string> GetErrorMessage() const
    {
        return m_errorCode;
    }

    //! @brief 获取原始曲线
    inline sggk::BSplineCurve3DPtr& GetOriginalCurve() { return m_OriginalCurve; }

    //! @brief 获取原始型值点 
    inline std::vector<sggk::Point3D>& GetFitPoints() { return m_OriginalFitPoints; }

    //! @brief 获取Hausdorff距离结果
    inline double GetHausdorffDistanceResult() const { return m_HausdorffDistanceResult; }

public:

    //! @brief 计算曲线的光顺能量
    double ComputeCurveFairEnergy(const sggk::BSplineCurve3DPtr& theCurve);

    //! @brief 计算曲线三阶导数对应的角度值
    std::vector<double> ComputeThirdDerivativeAngles(
        const sggk::BSplineCurve3DPtr& theCurve);

    // static
   //! @brief 计算曲线与曲线之间的Hausdorff距离
    static double GetCurveCurveHausdorffDistance(
        const sggk::BSplineCurve3DPtr theOriginalCurve,
        const sggk::BSplineCurve3DPtr theOperateCurve);

    //! @brief 计算一组点到曲线的 Hausdorff 距离
    static double GetPointCurveHausdorffDistance(
        const std::vector<sggk::Point3D>& thePoints,
        const sggk::BSplineCurve3DPtr& theCurve);

    //! @brief 最近邻排序算法
    static std::vector<sggk::Point3D> ReorderPointsNearestNeighbor(const std::vector<sggk::Point3D>& points);

    //! @brief 根据型值点拟合曲线
    static sggk::BSplineCurve3DPtr SampleAndFitBSpline(
        const sggk::BSplineCurve3DPtr& theOriginalCurve,
        int theSampleNum,
        std::vector<sggk::Point3D>& theFitPoints,
        int theMaxDegree = 3,
        int theContinuity = 0,
        double theTolerance = FITTOLERANCE);

    //! @brief 对单条 B 样条曲线进行等参线等距采样
    static std::vector<sggk::Point3D> SampleCurveWithArcLength(const sggk::BSplineCurve3DPtr& bsplineCurve, int numSamples);

    //! @brief 根据型值点拟合曲线
    static sggk::BSplineCurve3DPtr SampleAndFitBSpline(
        std::vector<sggk::Point3D>& theFitPoints,
        double theTolerance = FITTOLERANCE,
        bool theResortFlag = true,
        int theMaxDegree = 3,
        int theContinuity = 0);
    //! @brief 对曲线节点向量[0, 1]化
    void UniformCurve(sggk::BSplineCurve3DPtr& curve);

    //! @brief 更新当前曲线的基本属性
    inline void UpdateCurveAttribute()
    {
        this->m_FitPointParameters = ReCalculateFitPointParameters(m_OriginalFitPoints, m_OriginalCurve);
        this->m_FirstPole = m_OriginalCurve->ControlPoints().front();
        this->m_LastPole = m_OriginalCurve->ControlPoints().back();
        this->m_OriginalCurveKnotSequence = GetKnotSequence(m_OriginalCurve);
        this->m_OriginalCurveKnots = m_OriginalCurve->Knots();
    }

    //! @brief 迭代光顺算法
    void Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D);

    //! @brief 根据曲率自动细分曲线节点向量
    sggk::BSplineCurve3DPtr RefineCurveByCurvatureAuto(
        const sggk::BSplineCurve3DPtr& theCurve,
        const std::vector<double>& myKnotSeq,
        const int baseInsertNum);

    //! @brief 根据型值点对曲线节点向量进行细分
    sggk::BSplineCurve3DPtr RefineCurveByFitPoints(
        const sggk::BSplineCurve3DPtr& theCurve,
        const int baseInsertNum);

    //! @brief 在指定节点区间内插入若干个节点
    sggk::BSplineCurve3DPtr InsertKnotsBetweenKnotSpan(
        const sggk::BSplineCurve3DPtr& theCurve,
        double t1,
        double t2,
        int nInsert);

    //! @brief 根据曲线切向角度自动插入节点细分曲线
    sggk::BSplineCurve3DPtr RefineByTangentAngle(
        const sggk::BSplineCurve3DPtr& theCurve,
        const int theBasicInsertNum = 1);

    //! @brief 在曲线上插入指定参数值的节点
    static sggk::BSplineCurve3DPtr InsertKnots(
        const sggk::BSplineCurve3DPtr& theCurve,
        const std::vector<double>& params,
        int mult);

    //! @brief 在曲线中插入一个节点（参数值）
    sggk::BSplineCurve3DPtr InsertKnot(
        const sggk::BSplineCurve3DPtr& theCurve,
        double u,
        int mult = 1);

    //! @brief 在每一对相邻节点区间内插入 nInsert 个均匀分布的新节点
    sggk::BSplineCurve3DPtr InsertUniformKnotsInAllSpans(
        const sggk::BSplineCurve3DPtr& theCurve,
        int nInsert);

    //! @brief 获取弧长参数映射函数 f(s)
    sggk::BSplineCurve3DPtr GetArclengthParameterMapping(
        const sggk::BSplineCurve3DPtr& theCurve,
        const double theTolerance = 1e-7);

    //! @brief 给定s，计算t = f(s)
    static double f(const double theParameter, const int k = 0);

    //! @brief 计算 B 样条曲线对弧长参数 s 的导向向量 Ds
    sggk::Vector3D Ds(sggk::BSplineCurve3DPtr& theBSplineCurve, double sParameter, int k = 0);

    //! @brief 计算 B 样条曲线对曲线参数 t 的导向向量 Dt
    sggk::Vector3D Dt(sggk::BSplineCurve3DPtr& theBSplineCurve, double tParameter, int k = 0);

    //! @brief 给定t，计算s = f逆(t)
    double f_inverse(const sggk::BSplineCurve3DPtr& theBSplineCurve, double t);

    //! @brief 计算控制点权重矩阵
    static std::vector<std::pair<sggk::Point3D, double>> ReCalculateFitPointParameters(
        const std::vector<sggk::Point3D>& theFitPoints,
        const sggk::BSplineCurve3DPtr& theCurve = nullptr);

    //! @brief 计算 B 样条曲线在两个参数值之间的曲线长度
    static double ComputeCurveLengthBetweenParameters(const sggk::BSplineCurve3DPtr& theCurve,
        double theParameter1, double theParameter2);

    //! @brief 根据给定点找到 B 样条曲线上的最近点参数
    double GetPntParameterOnCurve(const sggk::BSplineCurve3DPtr& theCurve, const sggk::Point3D& thePoint);

    //! @brief 计算 B 样条曲线在两个给定点之间的曲线长度
    double ComputeCurveLengthBetweenParameters(const sggk::BSplineCurve3DPtr& theCurve, sggk::Point3D thePnt1, sggk::Point3D thePnt2);

    //! @brief 对曲线进行弧长参数映射采样
    std::pair<std::vector<sggk::Point3D>, std::vector<sggk::Point3D>> SampleCurveWithArclengthMapping(
        const sggk::BSplineCurve3DPtr& theCurve,
        const int nSamples = 10);

    //! @brief 获取光顺曲线中间结果
    sggk::BSplineCurve3DPtr GetTempFairCurve(const sggk::BSplineCurve3DPtr& theCurve,
        Eigen::MatrixXd M,
        Eigen::MatrixXd V,
        Eigen::MatrixXd& D0,
        Eigen::MatrixXd& D,
        double theAlpha);

    //! @brief 计算光顺曲线带切向约束中间结果
    sggk::BSplineCurve3DPtr GetTempFairCurveWithTangentConstraint(
        const sggk::BSplineCurve3DPtr& theCurve,
        Eigen::MatrixXd M,
        Eigen::MatrixXd V,
        Eigen::MatrixXd& D0,
        Eigen::MatrixXd& D,
        double theAlpha);

    //! @note 利用连续性矩阵计算光顺曲线
    sggk::BSplineCurve3DPtr GetTempContinuniousFairCurve(
        const sggk::BSplineCurve3DPtr& theCurve,
        Eigen::MatrixXd M,
        Eigen::MatrixXd C,
        Eigen::MatrixXd V,
        Eigen::MatrixXd& D0,
        Eigen::MatrixXd& D,
        double alpha);

    //! @brief 计算连续性矩阵
    Eigen::MatrixXd ComputeContinuityMatrix(
        const sggk::BSplineCurve3DPtr& theBSplineCurve,
        const int p);
public:
    //! @brief 计算能量矩阵
    Eigen::MatrixXd ComputeEnergyMatrix(
        const sggk::BSplineCurve3DPtr& theBSplineCurve,
        const int p,
        const double theTolerance = 1e-6);

    //! @brief 计算稀疏能量矩阵
    Eigen::SparseMatrix<double> ComputeSparseEnergyMatrix(
        const sggk::BSplineCurve3DPtr& theBSplineCurve,
        const int p,
        const double theTolerance = 1e-6);

    //! @brief 根据端点切向约束计算约束矩阵 C 和右端项 h
    Eigen::MatrixXd ComputeConstraintMatrix(
        const Eigen::MatrixXd& theD0,
        Eigen::VectorXd& theH);

    //! @brief 计算稀疏约束矩阵 C 和右端向量 H
    Eigen::SparseMatrix<double> ComputeSparseConstraintMatrix(
        const Eigen::MatrixXd& theD0,
        Eigen::VectorXd& theH);

    //! @brief 设置控制点权重矩阵
    void SetControlPointWeightMatrix(
        const sggk::BSplineCurve3DPtr& theCurve,
        Eigen::MatrixXd& V);

    //! @brief 计算基函数导数
    static double BasisFunctionDerivative(
        const double u,
        const int i,
        const int p,
        const int k,
        const std::vector<double>& Knots);

    //! @brief 创建新的B样条曲线
    sggk::BSplineCurve3DPtr CreateNewBSplineCurve(
        const sggk::BSplineCurve3DPtr& theOriginalCurve,
        const Eigen::MatrixXd& D);

    //! @brief 计算曲线与型值点之间的Hausdorff距离
    double GetFitPointsCurveHausdorffDistance(
        const std::vector<std::pair<sggk::Point3D, double>> theFitPointParams,
        const sggk::BSplineCurve3DPtr& theOperateCurve);

    //! @brief 获取弧长参数映射函数
    inline sggk::BSplineCurve3DPtr GetArcLengthMappingFunction() { return m_ArcLengthMappingFunction; }

    //! @brief 展开节点序列（knots + mults → 完整序列）
    inline static std::vector<double> GetKnotSequence(const sggk::BSplineCurve3DPtr& curve)
    {
        std::vector<double> seq;
        if (!curve) return seq;
        const auto& k = curve->Knots();
        const auto& m = curve->Mults();
        for (size_t i = 0; i < k.size(); ++i)
            for (unsigned j = 0; j < m[i]; ++j)
                seq.push_back(k[i]);
        return seq;
    }

    //! @brief 求 B 样条基函数及其导数（自实现，替代 BSplCLib::EvalBsplineBasis）
    //! @return basis 尺寸 (maxDeriv+1) x (deg+1)，basis(d, j) = 第 j 个非零基函数的 (d-1) 阶导数（1-indexed 语义）
    static void EvalBsplineBasis(
        int maxDeriv, int deg, const std::vector<double>& knots, double u, int& firstIndex,
        Eigen::MatrixXd& basis);

    //! @brief 计算曲线总弧长（自实现，数值积分 |D1|）
    static double ComputeCurveLength(const sggk::BSplineCurve3DPtr& curve, double tol = 1e-6);

private:
    // 弧长参数映射采样点数
    double m_ArcLengthMappingSampleNum = PARANUM;

    // 弧长参数映射函数
    static sggk::BSplineCurve3DPtr m_ArcLengthMappingFunction;

    // 原始曲线
    sggk::BSplineCurve3DPtr m_OriginalCurve;
    // 原始型值点
    std::vector<sggk::Point3D> m_OriginalFitPoints;
    // 原始控制点
    sggk::Point3DArray m_OriginalPoles;
    // 原始节点向量
    std::vector<double> m_OriginalCurveKnots;
    // 原始节点向量序列
    std::vector<double> m_OriginalCurveKnotSequence;
    // 原始曲线首控制点
    sggk::Point3D m_FirstPole;
    // 原始曲线尾控制点
    sggk::Point3D m_LastPole;

    // 光顺曲线
    sggk::BSplineCurve3DPtr m_ResultCurve;

    // 如果需要分段光顺的话
    std::vector<sggk::BSplineCurve3DPtr> m_ResultCurveArray;

    // 能量矩阵
    Eigen::MatrixXd M;
    // 连续性矩阵
    Eigen::MatrixXd C;
    // 控制点权重矩阵
    Eigen::MatrixXd V;
    // 光顺能量权重
    double m_Alpha;
    // 距离容差
    double m_HausdorffDistanceTol;
    // 最终Hausdorff距离
    double m_HausdorffDistanceResult;

    double m_FitTolerance;

    // 光顺迭代下降Ratio
    double m_AlphaRatio;
    std::vector<std::pair<sggk::Point3D, double>> m_FitPointParameters;

private:
    std::vector<std::string> m_errorCode;
};
// #endif
