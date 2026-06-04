// #ifndef CURVEFAIR_H
// #define CURVEFAIR_H

#include "Geom_BSplineCurve.hxx"
#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "GeomAPI_Interpolate.hxx"
#include "Geom_BSplineCurve.hxx"
#include "GeomAPI_PointsToBSpline.hxx"
#include "GeomAdaptor_Curve.hxx"
#include "GeomConvert_BSplineCurveToBezierCurve.hxx"
#include "GCPnts_AbscissaPoint.hxx"
#include "GeomConvert_CompCurveToBSplineCurve.hxx"
#include "GCPnts_QuasiUniformAbscissa.hxx"

#include "BRep_Tool.hxx"
#include "BRepTools.hxx"
#include "BRep_Builder.hxx"
#include "BRepPrim_FaceBuilder.hxx"

#include "math_GaussSingleIntegration.hxx"
#include "math_Matrix.hxx"
#include "math.hxx"

#include "TopoDS_Shape.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Edge.hxx"     
#include "TopExp_Explorer.hxx"

#include "Interface_Static.hxx"

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <mutex>
#include <map>

// 并行计算矩阵
//#include <omp.h>

const Standard_Real HAUSDORFFDISTANCETOL = 50;
const Standard_Real ALPHA = 1;
const Standard_Real ALPHARATIO = 0.5;
const Standard_Integer PARANUM = 100;
const Standard_Real FITTOLERANCE = 50;
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
        Handle(Geom_BSplineCurve)& theBSplineCurve,
        std::vector<gp_Pnt>& theFitPoints,
        Standard_Real theHausdorffDistanceTol = HAUSDORFFDISTANCETOL)
    {
        if (theBSplineCurve.IsNull() && theFitPoints.size() == 0)
        {
            m_errorCode.push_back("CONSTRUCT ERROR::theBsplineCurve and theFitPoints cannot be both empty.");
            return;
        }

        if (!theBSplineCurve.IsNull()) 
        {
            // 检查次数
            if (theBSplineCurve->Degree() != 3) 
            {
                m_errorCode.push_back("CONSTRUCT ERROR::The BSplineCurve Degree must be 3!");
                return;
            }

            // 检查节点重复性
            const TColStd_Array1OfReal& knots = theBSplineCurve->Knots();
            const TColStd_Array1OfInteger& mults = theBSplineCurve->Multiplicities();

            // 前后端点重复次数z
            int firstMult = mults.Value(mults.Lower());
            int lastMult = mults.Value(mults.Upper());

            if (firstMult < theBSplineCurve->Degree() + 1 || lastMult < theBSplineCurve->Degree() + 1) 
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

        if (m_OriginalCurve.IsNull() && m_OriginalFitPoints.size() > 0)
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
    //! @note 按弧长对原始曲线进行采样，并用采样点拟合新的 B 样条曲线，
    //!       同时更新拟合点的参数值和曲线属性
    inline bool ReFitCurve(Standard_Real theFitTolerance = 50, Standard_Integer theSampleNum = 100)
    {
        // 基本参数检查
        if (m_OriginalCurve.IsNull())
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
            if (m_OriginalCurve.IsNull())
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
    //! @note 根据切向角变化对曲线进行节点细分，提升局部形状控制精度
    inline bool InsertKnots(Standard_Integer theInsertLevel = 1)
    {
        if (m_OriginalCurve.IsNull())
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
            if (m_OriginalCurve.IsNull())
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
    //! @return 光顺曲线
    inline Handle(Geom_BSplineCurve) GetResult() { return m_ResultCurve; }

    //! @brief 是否有错误信息
    //! @return 错误信息
    inline Standard_Boolean HasErrorMessage() const
    {
        return m_errorCode.size() > 0;
    }

    //! @brief 获取错误信息
    //! @return 错误信息
    inline std::vector<std::string> GetErrorMessage() const
    {
        return m_errorCode;
    }

    //! @brief 获取原始曲线
    //! @return 原始曲线
    inline Handle(Geom_BSplineCurve)& GetOriginalCurve() { return m_OriginalCurve; }

    //! @brief 获取原始型值点 
    //! @return 原始型值点
    inline std::vector<gp_Pnt>& GetFitPoints() { return m_OriginalFitPoints; }

    //! @brief 获取Hausdorff距离结果
    //! @return Hausdorff距离结果
    inline Standard_Real GetHausdorffDistanceResult() const { return m_HausdorffDistanceResult; }

public:

    //! @brief 计算曲线的光顺能量
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @return 曲线的光顺能量值
    Standard_Real ComputeCurveFairEnergy(const Handle(Geom_BSplineCurve)& theCurve);

    //! @brief 计算曲线三阶导数对应的角度值
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @return 每个节点向量的三阶导数角度值数组
    std::vector<Standard_Real> ComputeThirdDerivativeAngles(
        const Handle(Geom_BSplineCurve)& theCurve);

    // static
   //! @brief 计算曲线与曲线之间的Hausdorff距离
   //! @param [in]theOriginalCurve 原始曲线
   //! @param [in]theOperateCurve 操作曲线
   //! @return 曲线与曲线之间的Hausdorff距离
    static Standard_Real GetCurveCurveHausdorffDistance(
        const Handle(Geom_BSplineCurve) theOriginalCurve,
        const Handle(Geom_BSplineCurve) theOperateCurve);

    //! @brief 计算一组点到曲线的 Hausdorff 距离
    //! @param [in] thePoints 输入点集
    //! @param [in] theCurve 输入曲线
    //! @return 点集到曲线的最大距离（Hausdorff 距离）
    static Standard_Real GetPointCurveHausdorffDistance(
        const std::vector<gp_Pnt>& thePoints,
        const Handle(Geom_BSplineCurve)& theCurve);

    //! @brief 最近邻排序算法
    //! @param [in]points 输入的点集
    //! @return 排序后的点集
    static std::vector<gp_Pnt> ReorderPointsNearestNeighbor(const std::vector<gp_Pnt>& points);

    //! @brief 根据型值点拟合曲线
    //! @param [in]theOriginalCurve 输入的曲线
    //! @param [in]theSampleNum 采样点数
    //! @param [in]theFitPoints 采样的拟合点
    //! @param [in]theOccfitFlag 是否利用Occ型值点拟合，若为false，则使用自编算法
    //! @param [in]theMaxDegree 曲线最大度数
    //! @param [in]theContinuity 曲线连续性
    //! @param [in]theTolerance 拟合精度
    //! @return 拟合后FT的曲线
    static Handle(Geom_BSplineCurve) SampleAndFitBSpline(
        const Handle(Geom_BSplineCurve)& theOriginalCurve,
        Standard_Integer theSampleNum,
        std::vector<gp_Pnt>& theFitPoints,
        Standard_Integer theMaxDegree = 3,
        GeomAbs_Shape theContinuity = GeomAbs_G2,
        Standard_Real theTolerance = FITTOLERANCE);

    //! @brief 对单条 B 样条曲线进行等参线等距采样
    //! @note 该函数会在曲线上生成 numSamples 个均匀分布的点，用于拟合操作
    //! @param [in] bsplineCurve 输入的 B 样条曲线
    //! @param [in] numSamples 采样点数量，即希望沿曲线均匀生成的点数
    //! @return std::vector<gp_Pnt> 返回采样得到的点集，每个 gp_Pnt 表示曲线上的一个三维点
    static std::vector<gp_Pnt> SampleCurveWithArcLength(const Handle(Geom_BSplineCurve)& bsplineCurve, Standard_Integer numSamples);

    //! @brief 根据型值点拟合曲线
    //! @param [in]theFitPoints 输入的型值点
    //! @param [in]theOccfitFlag 是否利用Occ型值点拟合，若为false，则使用自编算法
    //! @param [in]theResortFlag 是否重新排序型值点，若为true，则重新排序型值点
    //! @param [in]theMaxDegree 曲线最大度数
    //! @param [in]theContinuity 曲线连续性
    //! @param [in]theTolerance 拟合精度
    //! @return 拟合后的曲线
    static Handle(Geom_BSplineCurve) SampleAndFitBSpline(
        std::vector<gp_Pnt>& theFitPoints,
        Standard_Real theTolerance = FITTOLERANCE,
        Standard_Boolean theResortFlag = Standard_True,
        Standard_Integer theMaxDegree = 3,
        GeomAbs_Shape theContinuity = GeomAbs_C3);
    //! @brief 对曲线节点向量[0, 1]化
    //! @param [in,out] curve 待处理曲线
    void UniformCurve(Handle(Geom_BSplineCurve)& curve);

    //! @brief 更新当前曲线的基本属性
    //! @note 提取首末控制点、节点序列和主节点数组等信息，用于后续计算或显示
    inline void UpdateCurveAttribute()
    {
        this->m_FitPointParameters = ReCalculateFitPointParameters(m_OriginalFitPoints, m_OriginalCurve);
        this->m_FirstPole = m_OriginalCurve->Pole(1);
        this->m_LastPole = m_OriginalCurve->Pole(m_OriginalCurve->NbPoles());
        this->m_OriginalCurveKnotSequence = OccArrayConvertoVector(m_OriginalCurve->KnotSequence());
        this->m_OriginalCurveKnots = OccArrayConvertoVector(m_OriginalCurve->Knots());
    }

    //! @brief 迭代光顺算法
    //! @param [in] D0 原始控制点
    //! @param [in] D 新生成控制点
    void Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D);

    //! @brief 根据曲率自动细分曲线节点向量
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] myKnotSeq 已知的节点序列
    //! @param [in] baseInsertNum 每个区间基础插入节点数
    //! @return 细分后的新 B 样条曲线
    Handle(Geom_BSplineCurve) RefineCurveByCurvatureAuto(
        const Handle(Geom_BSplineCurve)& theCurve,
        const std::vector<Standard_Real>& myKnotSeq,
        const Standard_Integer baseInsertNum);

    //! @brief 根据型值点对曲线节点向量进行细分
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] baseInsertNum 每个区间基础插入节点数
    //! @return 细分后的新 B 样条曲线
    Handle(Geom_BSplineCurve) RefineCurveByFitPoints(
        const Handle(Geom_BSplineCurve)& theCurve,
        const Standard_Integer baseInsertNum);

    //! @brief 在指定节点区间内插入若干个节点
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] t1 区间起始参数值
    //! @param [in] t2 区间结束参数值
    //! @param [in] nInsert 插入节点数量
    //! @return 插入节点后的新曲线
    Handle(Geom_BSplineCurve) InsertKnotsBetweenKnotSpan(
        const Handle(Geom_BSplineCurve)& theCurve,
        Standard_Real t1,
        Standard_Real t2,
        Standard_Integer nInsert);

    //! @brief 根据曲线切向角度自动插入节点细分曲线
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] theBasicInsertNum 每个区间基础插入节点数，默认值为 1
    //! @return 细分后的新曲线
    Handle(Geom_BSplineCurve) RefineByTangentAngle(
        const Handle(Geom_BSplineCurve)& theCurve,
        const Standard_Integer theBasicInsertNum = 1);

    //! @brief 在曲线上插入指定参数值的节点
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] params 待插入的参数值集合
    //! @param [in] mult 插入节点的重数
    //! @return 插入节点后的新曲线
    static Handle(Geom_BSplineCurve) InsertKnots(
        const Handle(Geom_BSplineCurve)& theCurve,
        const std::vector<Standard_Real>& params,
        Standard_Integer mult);

    //! @brief 在曲线中插入一个节点（参数值）
    //! @param [in] theCurve  原始 B 样条曲线
    //! @param [in] u         插入的参数值（节点）
    //! @param [in] mult      插入的重数（默认 1）
    //! @return 插入节点后的新曲线
    Handle(Geom_BSplineCurve) InsertKnot(
        const Handle(Geom_BSplineCurve)& theCurve,
        Standard_Real u,
        Standard_Integer mult = 1);

    //! @brief 在每一对相邻节点区间内插入 nInsert 个均匀分布的新节点
    //! @param [in] theCurve     原始曲线
    //! @param [in] nInsert      每个区间插入的新节点数量
    //! @return 插值后的新曲线
    Handle(Geom_BSplineCurve) InsertUniformKnotsInAllSpans(
        const Handle(Geom_BSplineCurve)& theCurve,
        Standard_Integer nInsert);

    //! @brief 获取弧长参数映射函数 f(s)
    //! @note 原始曲线C(t)不满足弧长参数化,计算弧长映射函数 使得C(f(s)) 参数 s 满足弧长参数化
    //! @param [in]theCurve 输入的曲线
    //! @param [in]theTolerance 拟合精度
    //! @return 弧长参数映射函数
    Handle(Geom_BSplineCurve) GetArclengthParameterMapping(
        const Handle(Geom_BSplineCurve)& theCurve,
        const Standard_Real theTolerance = 1e-7);

    //! @brief 给定s，计算t = f(s)
    //! @param [in]theParameter 参数
    //! @param [in]k 阶数
    //! @return t = f(s)
    static Standard_Real f(const Standard_Real theParameter, const Standard_Integer k = 0);

    //! @brief 计算 B 样条曲线对弧长参数 s 的导向向量 Ds
    //! @param [in] theBSplineCurve 输入的 B 样条曲线
    //! @param [in] sParameter 弧长参数 s
    //! @param [in] k 导数阶数，默认值为 0
    //! @return 曲线在 s 参数处的导向向量 Ds
    gp_Vec Ds(Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real sParameter, Standard_Integer k = 0);

    //! @brief 计算 B 样条曲线对曲线参数 t 的导向向量 Dt
    //! @param [in] theBSplineCurve 输入的 B 样条曲线
    //! @param [in] tParameter 曲线参数 t
    //! @param [in] k 导数阶数，默认值为 0
    //! @return 曲线在 t 参数处的导向向量 Dt
    gp_Vec Dt(Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real tParameter, Standard_Integer k = 0);

    //! @brief 给定t，计算s = f逆(t)
    //! @param [in]theBSplineCurve 输入的曲线
    //! @param [in]t 参数
    //! @return t参数对应的弧长参数s
    Standard_Real f_inverse(const Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real t);

    //! @brief 计算控制点权重矩阵
    //! @param [in]theCurve 输入的曲线
    //! @param [in]theFirPoints 输入的型值点
    //! @return 型值点对应在曲线上的弧长参数值
    static std::vector<std::pair<gp_Pnt, Standard_Real>> ReCalculateFitPointParameters(
        const std::vector<gp_Pnt>& theFitPoints,
        const Handle(Geom_BSplineCurve)& theCurve = nullptr);

    //! @brief 计算 B 样条曲线在两个参数值之间的曲线长度
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] theParameter1 起始参数值
    //! @param [in] theParameter2 终止参数值
    //! @return 曲线在指定参数区间内的长度
    Standard_Real ComputeCurveLengthBetweenParameters(const Handle(Geom_BSplineCurve)& theCurve, 
        Standard_Real theParameter1, Standard_Real theParameter2);

    //! @brief 根据给定点 gp_Pnt 找到 B 样条曲线上的最近点参数
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] thePoint 指定点
    //! @return 曲线上距离 thePoint 最近的参数值
    Standard_Real GetPntParameterOnCurve(const Handle(Geom_BSplineCurve)& theCurve, const gp_Pnt& thePoint);

    //! @brief 计算 B 样条曲线在两个给定点之间的曲线长度
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] thePnt1 起始点
    //! @param [in] thePnt2 终止点
    //! @return 曲线在指定点之间的长度
    Standard_Real ComputeCurveLengthBetweenParameters(const Handle(Geom_BSplineCurve)& theCurve, gp_Pnt thePnt1, gp_Pnt thePnt2);

    //! @brief 对曲线进行弧长参数映射采样
    //! @param [in] theCurve 输入的 B 样条曲线
    //! @param [in] nSamples 采样点数量，默认值为 10
    //! @return pair，第一个为采样点列表，第二个为对应的参数值列表
    std::pair<std::vector<gp_Pnt>, std::vector <gp_Pnt>> SampleCurveWithArclengthMapping(
        const Handle(Geom_BSplineCurve)& theCurve,
        const Standard_Integer nSamples = 10);

    //! @brief 获取光顺曲线中间结果
    //! @param [in]theCurve 输入的曲线
    //! @param [in]M 能量矩阵
    //! @param [in]V 控制点权重矩阵
    //! @param [in]D0 原始控制点
    //! @param [out]D 新生成控制点
    //! @param [in]theAlpha 能量权重
    //! @return 光顺曲线中间结果
    Handle(Geom_BSplineCurve) GetTempFairCurve(const Handle(Geom_BSplineCurve)& theCurve,
        Eigen::MatrixXd M,
        Eigen::MatrixXd V,
        Eigen::MatrixXd& D0,
        Eigen::MatrixXd& D,
        Standard_Real theAlpha);

    //! @brief 计算光顺曲线带切向约束中间结果
    //! @param [in]theCurve 输入的曲线
    //! @param [in]M 能量矩阵
    //! @param [in]V 控制点权重矩阵
    //! @param [in]D0 原始控制点
    //! @param [out]D 新生成控制点
    //! @param [in]theAlpha 能量权重
    //! @return 光顺曲线中间结果
    Handle(Geom_BSplineCurve) GetTempFairCurveWithTangentConstraint(
        const Handle(Geom_BSplineCurve)& theCurve,
        Eigen::MatrixXd M,
        Eigen::MatrixXd V,
        Eigen::MatrixXd& D0,
        Eigen::MatrixXd& D,
        Standard_Real theAlpha);

    //! @note 利用连续性矩阵计算光顺曲线
    //! @param [in]theCurve 输入的曲线
    //! @param [in]M 能量矩阵
    //! @param [in]C 连续性矩阵
    //! @param [in]V 控制点权重矩阵
    //! @param [in]D0 原始控制点
    //! @param [out]D 新生成控制点
    //! @param [in]theAlpha 能量权重
    //! @return 光顺曲线中间结果
    Handle(Geom_BSplineCurve) GetTempContinuniousFairCurve(
        const Handle(Geom_BSplineCurve)& theCurve,
        Eigen::MatrixXd M,          // n×n 矩阵
        Eigen::MatrixXd C,          // m×n 矩阵（包含首尾控制点）
        Eigen::MatrixXd V,          // n×n 矩阵
        Eigen::MatrixXd& D0,        // n×3 初始控制点
        Eigen::MatrixXd& D,         // n×3 输出控制点
        Standard_Real alpha);

    //! @brief 计算连续性矩阵
    //! @param [in]theBSplineCurve 输入的曲线
    //! @param [in]p 曲线阶数
    //! @return 连续性矩阵
    Eigen::MatrixXd ComputeContinuityMatrix(
        const Handle(Geom_BSplineCurve)& theBSplineCurve,
        const Standard_Integer p);
public:
    //! @brief 计算能量矩阵
    //! @param [in]theBSplineCurve 输入的曲线
    //! @param [in]p 曲线阶数
    //! @param [in]tol 拟合精度
    //! @return 能量矩阵
    Eigen::MatrixXd ComputeEnergyMatrix(
        const Handle(Geom_BSplineCurve)& theBSplineCurve,
        const Standard_Integer p,
        const Standard_Real theTolerance = 1e-6);

    //! @brief 计算稀疏能量矩阵
    //! @param [in] theBSplineCurve 输入的 B 样条曲线
    //! @param [in] p 曲线阶数
    //! @param [in] theTolerance 积分或计算的容差，默认值 1e-6
    //! @return 稀疏能量矩阵
    Eigen::SparseMatrix<Standard_Real> ComputeSparseEnergyMatrix(
        const Handle(Geom_BSplineCurve)& theBSplineCurve,
        const Standard_Integer p,
        const Standard_Real theTolerance = 1e-6);

    //! @brief 根据端点切向约束计算约束矩阵 C 和右端项 h
    //! @param [in]theD0 原始控制点矩阵 (n x 3)
    //! @param [out]theH  输出参数，约束方程的右端项 h (4 x 1)
    //! @return Eigen::MatrixXd 约束矩阵 C (4 x 3*(n-2))
    Eigen::MatrixXd ComputeConstraintMatrix(
        const Eigen::MatrixXd& theD0,
        Eigen::VectorXd& theH);

    //! @brief 计算稀疏约束矩阵 C 和右端向量 H
    //! @param [in] theD0 原始控制点矩阵 (n x 3)
    //! @param [out] theH 输出约束向量 (4 x 1)
    //! @return 稀疏约束矩阵 C (4 x 3*(n-2))
    Eigen::SparseMatrix<Standard_Real> ComputeSparseConstraintMatrix(
        const Eigen::MatrixXd& theD0,
        Eigen::VectorXd& theH);

    //! @brief 设置控制点权重矩阵
    //! @note 利用 GrevilleAbscissae 计算控制点权重
    //! @param [in]theCurve 输入的曲线
    //! @param [in]V 控制点权重矩阵
    void SetControlPointWeightMatrix(
        const Handle(Geom_BSplineCurve)& theCurve,
        Eigen::MatrixXd& V);

    //! @brief 计算基函数导数
    //! @note  Nurbs Book 公式(2 - 9)
    //! @param [in]u 参数
    //! @param [in]i 控制点索引
    //! @param [in]p 曲线阶数
    //! @param [in]k 导数阶数
    //! @param [in]Knots 节点向量
    static Standard_Real BasisFunctionDerivative(
        const Standard_Real u,
        const Standard_Integer i,
        const Standard_Integer p,
        const Standard_Integer k,
        const std::vector<Standard_Real>& Knots);

    //! @brief 创建新的B样条曲线
    //! @param [in]theOriginalCurve 原始曲线
    //! @param [in]D 控制点
    //! @return 新的B样条曲线
    Handle(Geom_BSplineCurve) CreateNewBSplineCurve(
        const Handle(Geom_BSplineCurve)& theOriginalCurve,
        const Eigen::MatrixXd& D);



    //! @brief 计算曲线与型值点之间的Hausdorff距离
    //! @param [in]theFitPointParams 型值点及其对应弧长参数
    //! @param [in]theOperateCurve 操作曲线
    //! @return 曲线与型值点之间的Hausdorff距离
    Standard_Real GetFitPointsCurveHausdorffDistance(
        const std::vector<std::pair<gp_Pnt, Standard_Real>> theFitPointParams,
        const Handle(Geom_BSplineCurve)& theOperateCurve);

    //! @brief 获取弧长参数映射函数
    //! @return 弧长参数映射函数
    inline Handle(Geom_BSplineCurve) GetArcLengthMappingFunction() { return m_ArcLengthMappingFunction; }

    //! @brief 将Occ数组转换为标准向量
    //! @param [in]theOccArray Occ数组
    inline static std::vector<Standard_Real> OccArrayConvertoVector(TColStd_Array1OfReal theOccArray)
    {
        std::vector<Standard_Real> VectorArray;
        for (Standard_Integer i = theOccArray.Lower(); i <= theOccArray.Upper(); i++)
        {
            VectorArray.push_back(theOccArray.Value(i));
        }
        return VectorArray;
    }

    //! @brief 将Occ数组转换为标准向量
    //! @param [in]theOccArray Occ数组
    inline static std::vector<gp_Pnt> OccArrayConvertoVector(TColgp_Array1OfPnt theOccArray)
    {
        std::vector<gp_Pnt> PntArray;
        for (int i = 1; i <= theOccArray.Size(); i++)
        {
            PntArray.push_back(theOccArray.Value(i));
        }
        return PntArray;
    }

private:
    // 弧长参数映射采样点数
    Standard_Real m_ArcLengthMappingSampleNum = PARANUM;

    // 弧长参数映射函数
    static Handle(Geom_BSplineCurve) m_ArcLengthMappingFunction;

    // 原始曲线
    Handle(Geom_BSplineCurve) m_OriginalCurve;
    // 原始型值点
    std::vector<gp_Pnt> m_OriginalFitPoints;
    // 原始控制点
    TColgp_Array1OfPnt m_OriginalPoles;
    // 原始节点向量
    std::vector<Standard_Real> m_OriginalCurveKnots;
    // 原始节点向量序列
    std::vector<Standard_Real> m_OriginalCurveKnotSequence;
    // 原始曲线首控制点
    gp_Pnt m_FirstPole;
    // 原始曲线尾控制点
    gp_Pnt m_LastPole;

    // 光顺曲线
    Handle(Geom_BSplineCurve) m_ResultCurve;

    // 如果需要分段光顺的话
    std::vector<Handle(Geom_BSplineCurve)> m_ResultCurveArray;

    // 能量矩阵
    Eigen::MatrixXd M;
    // 连续性矩阵
    Eigen::MatrixXd C;
    // 控制点权重矩阵
    Eigen::MatrixXd V;
    // 光顺能量权重
    Standard_Real m_Alpha;
    // 距离容差
    Standard_Real m_HausdorffDistanceTol;
    // 最终Hausdorff距离
    Standard_Real m_HausdorffDistanceResult;

    Standard_Real m_FitTolerance;

    // 光顺迭代下降Ratio
    Standard_Real m_AlphaRatio;
    std::vector<std::pair<gp_Pnt, Standard_Real>> m_FitPointParameters;

private:
    std::vector<std::string> m_errorCode;
};
// #endif
