#pragma once
#include <Geometry/3D/Curve/BSplineCurve3D.h>
#include <GeomBase/Point3D.h>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>

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

	bool IsEqual(double x, double y, double tol = 1e-12)
	{
		return std::fabs(x - y) < tol;
	}

	bool IsGreater(double x, double y, double tol = 1e-12)
	{
		return (x - y) > tol;
	}

	bool IsLess(double x, double y, double tol = 1e-12)
	{
		return (y - x) > tol;
	}

	bool IsGreaterOrEqual(double x, double y, double tol = 1e-12)
	{
		return (x - y) > -tol;
	}

	bool IsLessOrEqual(double x, double y, double tol = 1e-12)
	{
		return (y - x) > -tol;
	}

	KnotUpdate(const sggk::BSplineCurve3DPtr& Bspline, const std::vector<double>& Sequences, const std::vector<sggk::Point3D>& Pnts, const std::vector<double>& Params);

	double SelfSingleUpdate(KONT_UPDATE_TYPE type);

	double getMaxError()
	{
		return maxError;
	}

	std::vector<double> getSequences()
	{
		return myCurrentSequences;
	}

	~KnotUpdate() {}

private:
	double adjustKnots();

	double selfUpdateUniform();

	double selfUpdateForLspia();

	double selfUpdateForMidKnot(bool isSingle = true);

	double selfUpdateForMidKnot_IntervalError();
private:

	double error(double u, const sggk::Point3D& P);

	void updateKnotsAndMutis();

	void updateSequences();

	void updateSequences(double newKnot);

	int checkNewKnot(double knot);

private:
	const sggk::BSplineCurve3DPtr& bspline;

	std::vector<double> myCurrentSequences;
	std::vector<sggk::Point3D> myPnts;
	std::vector<double> myParams;

	std::vector<double> myCurrentKnots;
	std::vector<int> myCurrentMutis;

	double maxError;
};
