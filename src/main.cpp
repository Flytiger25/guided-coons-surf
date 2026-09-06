#include <iostream>
#include <TopoDS.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Trsf.hxx>
#include <Geom_Line.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Type.hxx>
#include <IFSelect_ReturnStatus.hxx>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "GuidedCoonsSurfGenerator.h"

void Export_step_OCC(TopoDS_Shape shape, std::string filePath)
{
	STEPControl_Writer stepWriter;
	stepWriter.Transfer(shape, STEPControl_AsIs);
	IFSelect_ReturnStatus status = stepWriter.Write(filePath.c_str());
}

void LoadBSplineCurves(const std::string& filePath, std::vector<Handle(Geom_BSplineCurve)>& curveArray)
{
	// 获取文件后缀
	std::string extension = filePath.substr(filePath.find_last_of('.') + 1);

	TopoDS_Shape boundary;
	if (extension == "brep")
	{
		// 初始化边界Shape
		BRep_Builder B1;
		// 从文件读取BRep数据
		BRepTools::Read(boundary, filePath.c_str(), B1);
	}
	else if (extension == "step" || extension == "stp")
	{
		// 创建 STEP 文件读取器
		STEPControl_Reader reader;
		IFSelect_ReturnStatus status = reader.ReadFile(filePath.c_str());

		if (status == IFSelect_ReturnStatus::IFSelect_RetDone)
		{
			// 传输读取的数据
			reader.TransferRoots();
			boundary = reader.OneShape();
		}
	}

	// 遍历Shape中的边
	TopExp_Explorer explorer(boundary, TopAbs_EDGE);
	for (; explorer.More(); explorer.Next())
	{
		TopoDS_Edge edge = TopoDS::Edge(explorer.Current());

		// 获取边的几何表示
		TopLoc_Location loc;
		Standard_Real first, last;
		Handle(Geom_Curve) gcurve = BRep_Tool::Curve(edge, loc, first, last);
		gcurve = Handle(Geom_Curve)::DownCast(gcurve->Copy());

		// 检查曲线类型
		if (gcurve->DynamicType() == STANDARD_TYPE(Geom_Line))
		{
			// 如果是直线，转换为BSpline
			Handle(Geom_TrimmedCurve) aTrimmedLine = new Geom_TrimmedCurve(gcurve, first, last);
			Handle(Geom_BSplineCurve) aGeom_BSplineCurve = GeomConvert::CurveToBSplineCurve(aTrimmedLine);
			if (!aGeom_BSplineCurve.IsNull() && aGeom_BSplineCurve->IsKind(STANDARD_TYPE(Geom_BSplineCurve)))
			{
				curveArray.push_back(aGeom_BSplineCurve);
			}
		}
		else if (gcurve->DynamicType() == STANDARD_TYPE(Geom_BSplineCurve))
		{
			// 如果已经是BSpline，直接处理
			Handle(Geom_BSplineCurve) aGeom_BSplineCurve = Handle(Geom_BSplineCurve)::DownCast(gcurve);
			if (!aGeom_BSplineCurve.IsNull()) {
				aGeom_BSplineCurve->Segment(first, last);
				if (!aGeom_BSplineCurve.IsNull() && aGeom_BSplineCurve->IsKind(STANDARD_TYPE(Geom_BSplineCurve)))
				{
					curveArray.push_back(aGeom_BSplineCurve);
				}
			}
		}
	}
}

int main() 
{
	// 获取guideCurves
	std::vector<Handle(Geom_BSplineCurve)> guideCurves;
	std::string brepName = std::string(GUIDED_COONS_DATA_DIR) + "/input/1_internal.brep";
	// brepName += std::to_string(i);
	// brepName += "_internal.brep";
    LoadBSplineCurves(brepName, guideCurves);

	// 获取boundary
	std::vector<Handle(Geom_BSplineCurve)> boundary;
	brepName = std::string(GUIDED_COONS_DATA_DIR) + "/input/1_boundary.brep";
	// brepName += std::to_string(i);
	// brepName += "_internal.brep";
    LoadBSplineCurves(brepName, boundary);
    
	// try guide
	Handle(Geom_BSplineSurface) guidedSurf;
	GuidedCoonsSurfGenerator msg = GuidedCoonsSurfGenerator(boundary, guideCurves);
	msg.Perform();
	guidedSurf = msg.GuidedSurf();
	if (!msg.IsDone()) 
    {
		std::cout << "Case 1 " << " Guided Failing!" << std::endl;
	}

    TopoDS_Face guidedFace = BRepBuilderAPI_MakeFace(guidedSurf, Precision::Confusion());
    brepName = std::string(GUIDED_COONS_DATA_DIR) + "/output/";
	brepName += "1_guidedCoonsSurf.step";
	Export_step_OCC(guidedFace, brepName);

    return 0;
}