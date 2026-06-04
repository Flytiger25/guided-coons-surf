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
#include <map>

enum KONT_UPDATE_TYPE
{
	MID_KNOT_BY_SINGLE_ERROR,
	MID_KNOT_BY_INTERVAL_ERROR,
	PARAM_BASED_BY_INTERVAL_ERROR,
	ADJUST_KNOT,
	UNIFORM_UP_DATE
};

//this class use to update knot in iterate approximate
class KnotUpdate
{
public:

	Standard_Boolean IsEqual(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return std::fabs(x - y) < tol;
	}

	Standard_Boolean IsGreater(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return (x - y) > tol;
	}

	Standard_Boolean IsLess(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return (y - x) > tol;
	}

	Standard_Boolean IsGreaterOrEqual(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return (x - y) > -tol;
	}

	Standard_Boolean IsLessOrEqual(Standard_Real x, Standard_Real y, Standard_Real tol = Precision::Angular())
	{
		return (y - x) > -tol;
	}

	KnotUpdate(Handle(Geom_BSplineCurve)& Bspline, const std::vector<Standard_Real>& Sequences, const std::vector<gp_Pnt>& Pnts, const std::vector<Standard_Real>& Params);

	Standard_Real SelfSingleUpdate(KONT_UPDATE_TYPE type);

	Standard_Real getMaxError()
	{
		return maxError;
	}

	std::vector<Standard_Real> getSequences()
	{
		return myCurrentSequences;
	}

	~KnotUpdate() {}

private:
	Standard_Real adjustKnots();

	Standard_Real selfUpdateUniform();

	Standard_Real selfUpdateForLspia();

	Standard_Real selfUpdateForMidKnot(Standard_Boolean isSingle = true);

	Standard_Real selfUpdateForMidKnot_IntervalError();
private:

	Standard_Real error(Standard_Real u, const gp_Pnt& P);

	void updateKnotsAndMutis();

	void updateSequences();

	void updateSequences(Standard_Real newKnot);

	Standard_Integer checkNewKnot(Standard_Real knot);

private:
	Handle(Geom_BSplineCurve)& bspline;

	std::vector<Standard_Real> myCurrentSequences;
	std::vector<gp_Pnt> myPnts;
	std::vector<Standard_Real> myParams;

	std::vector<Standard_Real> myCurrentKnots;
	std::vector<Standard_Integer> myCurrentMutis;

	Standard_Real maxError;
};

