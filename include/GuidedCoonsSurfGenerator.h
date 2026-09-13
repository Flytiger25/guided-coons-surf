#pragma once
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <STEPControl_Writer.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <GeomConvert.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <TopoDS_Builder.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <GCPnts_AbscissaPoint.hxx>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#include <map>
#include <chrono>
#include <string>

#include "KnotUpdate.h"
#include "CurveFair.h"
#include "LobattoQuadrature.h"

#define TRIM_TOLERANCE 10.0          // 裁剪容差
#define COONS_TOLERANCE 10.0         // Coos边界线容差
#define APPROXIMATE_SAMPLING_NUM 50  // 边界线拟合采样点数量
#define APPROXIMATE_KNOTS_NUM 12     // 边界线拟合初始节点数量

#define MAX_OFFSETDISTANCE 500       // 最大优化偏移距离
#define MAX_ITERATIONS 5             // 最大迭代次数
#define STEP_LENGTH 1                // 步长
#define FAIRNESS_WEIGHT 0.1        // 光顺能量权重
//#define GUIDE_SAMPLING_NUM 20        // 采样点数量
#define GUIDE_SAMPLING_INTERVAL 500  // 引导线采样间隔（暂定每500mm采一个点）
#define KNOT_INTERVAL_THRESHOLD 0.1  // 节点插入间隔区间

//! @class GuidedCoonsSurfGenerator
//! @brief 生成带引导线的Coons曲面
//! @author 胡新宇
//! @date 2025.3.21
class GuidedCoonsSurfGenerator
{
public:
    //! @brief 获取指定阶数 Lobatto 求积数据（数据表独立维护于 LobattoQuadrature.h）
    inline const Lobatto::QuadratureData* getQuadratureData(int n) {
        return Lobatto::GetQuadratureData(n);
    }

    // 检查是否支持指定阶数
    inline bool isSupported(int n) {
        return Lobatto::HasQuadratureData(n);
    }


	//! @brief 构造函数
	//! @param [In] theBoundaryCurves 边界线
	//! @param [In] theGuideCurves 引导线
	//! @param [In] theTol 容差（默认为5.0）
	//! @return 
	GuidedCoonsSurfGenerator(const std::vector<Handle(Geom_BSplineCurve)>& boundaryCurves, const std::vector<Handle(Geom_BSplineCurve)>& guideCurves, Standard_Real theTol = 5.0);

	//! @brief 执行迭代生成曲面算法
	//! @return void
	//! 改返回值errorcode
	void Perform();

	//! @brief get
	//! @return 返回生成的带引导线的曲面
	inline Handle(Geom_BSplineSurface) GuidedSurf() const 
	{
		return m_guidedSurf;
	}

	//! @brief get
	//! @return 返回生成的带引导线的曲面
	inline std::vector<Handle(Geom_BSplineCurve)> GuideCurves() const
	{
		return m_guideCurves;
	}

    //! @brief get
    //! @return 返回生成的初始Coons曲面
    inline Handle(Geom_BSplineSurface) GetCoons() const
    {
        return m_coonsSurf;
    }       

	//! @brief 返回迭代是否结束
	//! @return 结束返回true，未结束返回false
    inline Standard_Boolean IsDone() const
	{
		return m_isDone;
	}

	void SetOriginalSurf(Handle(Geom_BSplineSurface) originalSurf)
	{
		m_originalSurf = originalSurf;
	}

	//! @brief 设置光顺能量方式
	//! @param [In] isCurveFair true=双向等参线光顺，false=格雷维尔坐标Laplace光顺
	inline void SetIsCurveFair(Standard_Boolean isCurveFair)
	{
		m_isCurveFair = isCurveFair;
	}

	//! @brief 设置中间结果（Coons曲面/每轮迭代曲面）STEP 导出目录
	//! @param [In] coonsOutDir 导出目录，为空字符串则不导出
	inline void SetCoonsOutDir(const std::string& coonsOutDir)
	{
		m_coonsOutDir = coonsOutDir;
	}

private:

	//------------------------构造Coons曲面-------------------------------

	//! @brief 构造Coons曲面
	//! @return void
	void ConstructCoonsSurf();

	//! @brief 对于三边情况构造退化边
	//! @param [In] theBoundaryCurves 边界线
	//! @return void
	void AddDegenerateCurve(std::vector<Handle(Geom_BSplineCurve)>& boundaryCurves);

	//! @brief 对给定曲线数组进行Coons曲面G0连续性处理和排列
	//! @param [Out] curveArray 输入的曲线数组，处理后可能会被修改
	//! @param [In] bslpineCurve1 第一条边界曲线
	//! @param [In] bslpineCurve2 第二条边界曲线
	//! @param [In] bslpineCurve3 第三条边界曲线
	//! @param [In] bslpineCurve4 第四条边界曲线
	//! @param [In] tol 容差值
	//! @param [In] isModify 
	//! @return Standard_Integer 返回处理结果状态码
	Standard_Integer Arrange_Coons_G0(std::vector<Handle(Geom_BSplineCurve)>& curveArray, Handle(Geom_BSplineCurve)& bslpineCurve1,
		Handle(Geom_BSplineCurve)& bslpineCurve2, Handle(Geom_BSplineCurve)& bslpineCurve3, Handle(Geom_BSplineCurve)& bslpineCurve4, Standard_Real tol = COONS_TOLERANCE, Standard_Integer isModify = true);

	//! @brief 基于四条边界曲线构造满足G0连续性的Coons曲面
	//! @param [In] curve1 第一条边界曲线
	//! @param [In] curve2 第二条边界曲线
	//! @param [In] curve3 第三条边界曲线
	//! @param [In] curve4 第四条边界曲线
	//! @param [Out] mySurface_coons 构造完成的Coons曲面
	//! @return void
	void Coons_G0(Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4, Handle(Geom_BSplineSurface)& mySurface_coons);

	//! @brief 根据边界曲线对内部曲线进行裁剪处理(730新算法)
	//! @param [Out] guideBSplineCurves 内部B样条曲线数组，处理后保留裁剪结果
	//! @param [In] boundaryCurveArray 边界曲线数组，定义裁剪范围
	//! @param [In] toleranceDistance 距离容差，用于判断曲线是否在边界内
	//! @return void
    void TrimInternalCurves(
		std::vector<Handle(Geom_BSplineCurve)>& theInternalBSplineCurves,
		const std::vector<Handle(Geom_BSplineCurve)>& theBoundaryCurveArray,
		Standard_Real theToleranceDistance = 10);

	Standard_Boolean IsCurveInsideBoundaries(
		const Handle(Geom_BSplineCurve)& theCurve,
		std::vector<Handle(Geom_BSplineCurve)>& theBoundaryCurveArray,
		Standard_Real theToleranceDistance = 10);

	// 辅助函数：判断一个点是否在由一系列二维点构成的多边形内部（使用射线法）。
	Standard_Boolean IsPointInPolygon2D(
		const gp_Pnt2d& theTestPoint,
		const std::vector<gp_Pnt2d>& thePolygon2d,
		Standard_Real theTolerance = 10);

	std::vector<gp_Pnt> DiscretizeBSplineCurve(
		const Handle(Geom_BSplineCurve)& theCurve,
		Standard_Integer numSegments,
		Standard_Boolean theBoundaryFlag = Standard_True);

	Standard_Boolean CurvesConnectedLoop(
		std::vector<Handle(Geom_BSplineCurve)>& theCurves,
		Standard_Real theTolerance = 10);

	Standard_Boolean IsCurveInsideSurface(
		const Handle(Geom_BSplineCurve)& theCurve,
		const Handle(Geom_Surface)& theSurface,
		const Standard_Real theTolerance);

	//! @brief 根据边界曲线对内部曲线进行裁剪处理
	//! @param [Out] guideBSplineCurves 内部B样条曲线数组，处理后保留裁剪结果
	//! @param [In] boundaryCurveArray 边界曲线数组，定义裁剪范围
	//! @param [In] toleranceDistance 距离容差，用于判断曲线是否在边界内
	//! @return void
	void TrimGuideCurves(std::vector<Handle(Geom_BSplineCurve)>& guideBSplineCurves, const std::vector<Handle(Geom_BSplineCurve)>& boundaryCurveArray, Standard_Real toleranceDistance = TRIM_TOLERANCE);

	//! @brief 计算两条B样条曲线之间的最小距离
	//! @param [In] guideCurve 第一条B样条曲线
	//! @param [In] boundaryCurve 第二条B样条曲线
	//! @return Standard_Real 返回两条曲线之间的最小距离
	Standard_Real ComputeCurveCurveDistance(const Handle(Geom_BSplineCurve)& guideCurve, const Handle(Geom_BSplineCurve)& boundaryCurve);

	//! @brief 对B样条曲线数组进行逼近处理以优化曲线表示
	//! @param [Out] curves B样条曲线数组，处理后曲线将被优化
	//! @param [In] samplingNum 采样点数量，用于曲线逼近(默认为50)
	//! @return void
	void ApproximateBoundaryCurves(std::vector<Handle(Geom_BSplineCurve)>& curves, Standard_Integer samplingNum = APPROXIMATE_SAMPLING_NUM);

	//! @brief 根据参数点生成B样条曲线节点矢量
	//! @param [In] params 参数点数组，用于生成节点矢量
	//! @param [In] n 控制点数
	//! @param [In] p 曲线次数
	//! @return std::vector<Standard_Real> 返回生成的节点矢量
	std::vector<Standard_Real> KnotGernerationByParams(const std::vector<Standard_Real>& params, Standard_Integer n, Standard_Integer p);

	//! @brief 通过迭代方法逼近点集生成B样条曲线
	//! @param [In] insertKnots 插入节点数组，用于迭代优化
	//! @param [In] pnts 数据点集
	//! @param [In] pntsParams 数据点对应的参数值
	//! @param [In] initKnots 初始节点矢量
	//! @param [In] degree 曲线次数
	//! @param [In] maxIterNum 最大迭代次数
	//! @param [In] toler 逼近容差
	//! @return Handle(Geom_BSplineCurve) 返回逼近生成的B样条曲线
	Handle(Geom_BSplineCurve) IterateApproximate(std::vector<Standard_Real>& insertKnots, const std::vector<gp_Pnt>& pnts, std::vector<Standard_Real>& pntsParams,
		std::vector<Standard_Real>& initKnots, Standard_Integer degree, Standard_Integer maxIterNum = 10, Standard_Real toler = 1);

	//! @brief 使用给定参数点和节点矢量逼近点集生成B样条曲线
	//! @param [In] pnts 数据点集
	//! @param [In] params 数据点对应的参数值
	//! @param [In] fKnots 最终使用的节点矢量
	//! @param [In] degree 曲线次数
	//! @return Handle(Geom_BSplineCurve) 返回逼近生成的B样条曲线
	Handle(Geom_BSplineCurve) ApproximateCurve(const std::vector<gp_Pnt>& pnts, std::vector<Standard_Real>& params, std::vector<Standard_Real>& fKnots, Standard_Integer degree);

	//! @brief 使用更详细参数设置逼近点集生成B样条曲线
	//! @param [In] pnts 数据点集
	//! @param [In] pntsParams 数据点对应的参数值
	//! @param [In] knots B样条节点数组
	//! @param [In] mutis 节点重复度数组
	//! @param [In] fKnots 最终使用的节点矢量
	//! @param [In] degree 曲线次数
	//! @return Handle(Geom_BSplineCurve) 返回逼近生成的B样条曲线
	Handle(Geom_BSplineCurve) ApproximateCurve(const std::vector<gp_Pnt>& pnts, std::vector<Standard_Real>& pntsParams, TColStd_Array1OfReal& knots,
		TColStd_Array1OfInteger& mutis, std::vector<Standard_Real>& fKnots, Standard_Integer degree);

	//! @brief 计算B样条曲线在指定参数处的残差向量
	//! @param [In] k 参数索引
	//! @param [In] dataPoints 数据点集
	//! @param [In] parameters 参数值数组
	//! @param [In] p 曲线次数
	//! @param [In] knots 节点矢量
	//! @param [In] ctrlPntNum 控制点数
	//! @return gp_Vec 返回计算得到的残差向量
	gp_Vec CalResPnt(Standard_Integer k, const std::vector<gp_Pnt>& dataPoints, const std::vector<Standard_Real>& parameters, Standard_Integer p, std::vector<Standard_Real>& knots, Standard_Integer ctrlPntNum);

	//! @brief 将参数序列转换为B样条节点矢量和重复度
	//! @param [In] sequence 输入的参数序列
	//! @param [Out] knots 输出的节点矢量
	//! @param [Out] multiplicities 输出的节点重复度数组
	//! @return void
	void SequenceToKnots(const std::vector<Standard_Real>& sequence, std::vector<Standard_Real>& knots, std::vector<Standard_Integer>& multiplicities);

	//------------------------带引导线的Coons-------------------------------

	//! @brief 生成的带引导线的Coons曲面
	//! @return void
	void ConstructSurfWithGuideCrvs();

	//! @brief 在引导线上采样得到采样点
	//! @return void
	void GetGuideSamples();

	//! @brief 获取采样点和对应的投影点之间偏移量和投影点参数
	//! @param [Out] thePntParams 投影点参数
	//! @param [Out] theOffsets 偏移量
	//! @param [In] isOriginal 是否为初始曲面（默认为true）
	//! @return void
	void GetSamplesOffset(std::vector<gp_Pnt2d>& thePntParams, std::vector<gp_Pnt>& theOffsets, Standard_Boolean isOriginal = Standard_True);

	//! @brief 对曲线均匀采样
	//! @param [In] theCurve 曲线
	//! @param [In] theSamplesNum 采样点数量
	//! @return 采样点数组
	std::vector<gp_Pnt> SampleGuideCurve(const Handle(Geom_BSplineCurve)& theCurve, Standard_Real startParam, Standard_Real endParam, Standard_Integer theSamplesNum);

	//! @brief 获取点到曲面的投影点和参数，获取偏移量
	//! @param [In] thePoints 点集
	//! @param [In] theSurface 曲面
	//! @param [Out] thePntParams 投影点参数
	//! @return 偏移量数组
	std::vector<gp_Pnt> ProjectPntsToSurf(const std::vector<gp_Pnt>& thePoints, std::vector<gp_Pnt>& theProjectionPoints, const Handle(Geom_BSplineSurface)& theSurface, std::vector<gp_Pnt2d>& thePntParams);

	//! @brief 计算偏移量
	//! @param [In] theSamples 采样点数组
	//! @param [In] theProjections 投影点数组
	//! @return 偏移量数组
	std::vector<gp_Pnt> CalOffsets(const std::vector<gp_Pnt>& theSamples, const std::vector<gp_Pnt>& theProjections);

	//------------------------拟合偏移曲面-------------------------------

	//! @brief 拟合偏移曲面计算控制点
	//! @param [In] theSamplePntOffsets 采样点偏移量向量
	//! @param [In] thePntParamsU u方向参数
	//! @param [In] thePntParamsV v方向参数
	//! @param [In] theUKnots u方向节点
	//! @param [In] theVKnots v方向节点
	//! @param [In] theDegU u方向degree
	//! @param [In] theDegV v方向degree
	//! @param [Out] theCtrlPoints 偏移曲面控制点
	//! @return void 
	void FitOffsetSurface(const std::vector<Eigen::Vector3d>& theSamplePntOffsets, const std::vector<Standard_Real>& thePntParamsU, const std::vector<Standard_Real>& thePntParamsV,
		const std::vector<Standard_Real>& theUKnots, const std::vector<Standard_Real>& theVKnots, Standard_Integer theDegU, Standard_Integer theDegV, std::vector<Eigen::Vector3d>& theCtrlPoints
	);
 
	//! @brief 构建非约束项系数矩阵
	//! @param [In] thePntParamsU u方向参数
	//! @param [In] thePntParamsV v方向参数
	//! @param [In] theUKnots u方向节点
	//! @param [In] theVKnots v方向节点
	//! @param [In] theDegU u方向degree
	//! @param [In] theDegV v方向degree
	//! @param [In] theCtrlPtsUNum u方向控制点数
	//! @param [In] theCtrlPtsVNum v方向控制点数
	//! @param [Out] theMatrixN 非约束项系数矩阵
	//! @return void
	void BuildMatrixUnconstraint(const std::vector<Standard_Real>& thePntParamsU, const std::vector<Standard_Real>& thePntParamsV, const std::vector<Standard_Real>& theUKnots, 
		const std::vector<Standard_Real>& theVKnots,Standard_Integer theDegU, Standard_Integer theDegV, Standard_Integer theCtrlPtsUNum, Standard_Integer theCtrlPtsVNum, Eigen::MatrixXd& theMatrixN
	);

	//! @brief 构建约束项系数矩阵
	//! @param [In] theCtrlPtsUNum u方向控制点数
	//! @param [In] theCtrlPtsVNum v方向控制点数
	//! @param [Out] theMatrixM 约束项系数矩阵
	//! @return void
	void BuildMatrixConstraint(Standard_Integer theCtrlPtsUNum, Standard_Integer theCtrlPtsVNum, Eigen::MatrixXd& theMatrixM);
 
	//! @brief 构建权重系数矩阵
	//! @param [In] thePntParamsU u方向参数
	//! @param [In] thePntParamsV v方向参数
	//! @param [In] theCtrlPtsUNum u方向控制点数
	//! @param [In] theCtrlPtsVNum v方向控制点数
	//! @param [In] alpha 光顺能量权重系数
	//! @param [Out] theMatrixW 权重系数矩阵
	//! @return void
	void BuildMatrixWeight(Standard_Integer thePntParamsSize, Standard_Integer theCtrlPtsUNum, Standard_Integer theCtrlPtsVNum, Standard_Real alpha, Eigen::MatrixXd& theMatrixW);

	//! @brief 计算基函数值
	//! @param [In] param 参数
	//! @param [In] index 下标
	//! @param [In] deg 次数
	//! @param [In] knots 节点向量
	//! @return 基函数值
	Standard_Real CalBasicFunction(Standard_Real param, Standard_Integer index, Standard_Integer deg, const std::vector<Standard_Real>& knots);

	//! @brief 计算Kronecker积
	//! @param [In] theMatA 第一个矩阵
	//! @param [In] theMatB 第二个矩阵
	//! @return 返回两个矩阵的Kronecker积
	Eigen::MatrixXd CalKroneckerProduct(const Eigen::MatrixXd& theMatA, const Eigen::MatrixXd& theMatB);

	//------------------------光顺能量矩阵-------------------------------

	/**
	 * @brief 获取控制点在展开向量中的索引
	 * @param i u方向索引
	 * @param j v方向索引
	 * @param n_v v方向控制点数量减1
	 * @return 在P_vec中的线性索引
	 */
	inline int GetControlPointIndex(int i, int j, int n_v) {
		return i * (n_v + 1) + j;
	}

	/**
	 * @brief 构造单条曲线的光顺能量矩阵（稠密版本）
	 * @param theBSplineCurve 输入的B样条曲线
	 * @param derivative_order 导数阶数（通常为2，表示二阶导数能量）
	 * @param tolerance 计算精度（默认为1e-6）
	 * @return 曲线光顺能量矩阵 M_v 或 M_u
	 *
	 * @note 这是单条曲线的能量矩阵，维度为 (m+1)×(m+1)
	 *       能量形式为: E_curve = Q^T * M * Q
	 *       调用 CurveFair::ComputeEnergyMatrix 实现
	 */
	Eigen::MatrixXd ConstructCurveSmoothingMatrix(
		const Handle(Geom_BSplineCurve)& theBSplineCurve,
		int derivative_order = 2,
		double tolerance = 1e-6
	);

	/**
	 * @brief 构造双向曲面光顺能量矩阵（稠密版本）
	 * @param theBSplineSurface 输入的B样条曲面
	 * @param u_params u-等参线的参数值列表
	 * @param v_params v-等参线的参数值列表
	 * @param derivative_order 导数阶数（默认为2）
	 * @param tolerance 计算精度（默认为1e-6）
	 * @return 总的双向光顺能量矩阵 M_total
	 *
	 * @note 返回的矩阵维度为 ((n_u+1)(n_v+1))×((n_u+1)(n_v+1))
	 *       总能量形式为: E_surface = P_vec^T * M_total * P_vec
	 *       其中 M_total = Σ(C_uk^T * M_v * C_uk) + Σ(D_vl^T * M_u * D_vl)
	 *       所有曲面参数（次数、节点向量等）自动从曲面对象中提取
	 */
	Eigen::MatrixXd ConstructBidirectionalSmoothingMatrix(
		const Handle(Geom_BSplineSurface)& theBSplineSurface,
		const std::vector<double>& u_params,
		const std::vector<double>& v_params,
		int derivative_order = 2,
		double tolerance = 1e-6
	);

	/**
	 * @brief 构造双向光顺能量矩阵的内部实现（模板函数）
	 * @tparam MatrixType 矩阵类型（Eigen::MatrixXd 或 Eigen::SparseMatrix<double>）
	 */
	Eigen::MatrixXd ConstructBidirectionalSmoothingMatrixImpl(
		const Handle(Geom_BSplineSurface)& theBSplineSurface,
		int n_u, int n_v,
		int p_u, int p_v,
		const std::vector<double>& knots_u,
		const std::vector<double>& knots_v,
		const std::vector<double>& u_params,
		const std::vector<double>& v_params,
		int derivative_order,
		double tolerance
	);

	/**
	 * @brief 构造等参线提取矩阵（稠密版本）
	 * @param n_u u方向控制点数量减1
	 * @param n_v v方向控制点数量减1
	 * @param p_u u方向B样条次数
	 * @param p_v v方向B样条次数
	 * @param knots_u u方向节点向量
	 * @param knots_v v方向节点向量
	 * @param direction 等参线方向
	 * @param param_value 固定参数的值（u0或v0）
	 * @return 提取矩阵C（稠密矩阵）
	 */
	Eigen::MatrixXd ConstructDenseMatrix(
		int n_u, int n_v,
		int p_u, int p_v,
		const std::vector<double>& knots_u,
		const std::vector<double>& knots_v,
		int direction,
		double param_value
	);

	/**
	 * @brief 构造u-等参线提取矩阵（内部函数）
	 */
	void ConstructUIsoparam(
		int n_u, int n_v,
		int p_u,
		const std::vector<double>& knots_u,
		double u_value,
		std::vector<Eigen::Triplet<double>>& triplets,
		Eigen::MatrixXd* dense_matrix = nullptr
	);

	/**
	 * @brief 构造v-等参线提取矩阵（内部函数）
	 */
	void ConstructVIsoparam(
		int n_u, int n_v,
		int p_v,
		const std::vector<double>& knots_v,
		double v_value,
		std::vector<Eigen::Triplet<double>>& triplets,
		Eigen::MatrixXd* dense_matrix = nullptr
	);


	Standard_Integer SetSameDistribution(Handle(Geom_BSplineCurve)& C1, Handle(Geom_BSplineCurve)& C2);

	//! @brief 在容差意义下比较 x 是否等于 y
	//! @param [In] x 第一个数
	//! @param [In] y 第二个数
	//! @return x 等于 y 则返回true， 否则返回false
	inline Standard_Boolean IsEqual(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return std::fabs(x - y) < tol;
	}

	//! @brief 在容差意义下比较 x 是否大于 y
	//! @param [In] x 第一个数
	//! @param [In] y 第二个数
	//! @return x 大于 y 则返回true， 否则返回false
	inline Standard_Boolean IsGreater(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return (x - y) > tol;
	}

	//! @brief 在容差意义下比较 x 是否小于 y
	//! @param [In] x 第一个数
	//! @param [In] y 第二个数
	//! @return x 小于 y 则返回true， 否则返回false
	inline Standard_Boolean IsLess(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return (y - x) > tol;
	}

	//! @brief 在容差意义下比较 x 是否大于等于 y
	//! @param [In] x 第一个数
	//! @param [In] y 第二个数
	//! @return x 大于等于 y 则返回true， 否则返回false
	inline Standard_Boolean IsGreaterOrEqual(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return (x - y) > -tol;
	}

	//! @brief 在容差意义下比较 x 是否小于等于 y
	//! @param [In] x 第一个数
	//! @param [In] y 第二个数
	//! @return x 小于等于 y 则返回true， 否则返回false
	inline Standard_Boolean IsLessOrEqual(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return (y - x) > -tol;
	}


	Handle(Geom_BSplineSurface) m_originalSurf; // 初始曲面
	Handle(Geom_BSplineSurface) m_coonsSurf; // Coons曲面
	Handle(Geom_BSplineSurface) m_guidedSurf; // 引导后的曲面
	std::vector<Handle(Geom_BSplineCurve)> m_boundaryCurves; // 输入的边界线
	std::vector<Handle(Geom_BSplineCurve)> m_guideCurves; // 输入的内部引导线
	std::vector<std::vector<std::pair<Standard_Real, Standard_Real>>> m_guideCurvesTrimIntervals; // 引导线裁剪区间
	std::vector<gp_Pnt> m_samples; // 所有引导线采样点

	Standard_Real m_tol; // 逼近的容差精度
	Standard_Integer m_iterateCount; // 迭代次数
	Standard_Boolean m_isDone; // 迭代完成的标志
	Standard_Boolean m_isCurveFair = Standard_True; // 光顺能量方式：true=双向等参线光顺，false=格雷维尔坐标Laplace光顺
	std::string m_coonsOutDir; // 中间结果 STEP 导出目录（空则不导出）
};

