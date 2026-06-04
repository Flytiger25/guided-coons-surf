#include "KnotUpdate.h"

KnotUpdate::KnotUpdate(Handle(Geom_BSplineCurve)& Bspline, const std::vector<Standard_Real>& Sequences, const std::vector<gp_Pnt>& Pnts, const std::vector<Standard_Real>& Params)
	:myParams(Params), myPnts(Pnts), bspline(Bspline), maxError(10000)
{
	std::for_each(Sequences.begin(), Sequences.end(), [&](Standard_Real REALT) {myCurrentSequences.push_back(REALT); });
	updateKnotsAndMutis();
}

Standard_Real KnotUpdate::SelfSingleUpdate(KONT_UPDATE_TYPE type)
{
	switch (type)
	{
	case MID_KNOT_BY_SINGLE_ERROR:
	{
		return selfUpdateForMidKnot();
		break;
	}
	case MID_KNOT_BY_INTERVAL_ERROR:
	{
		return selfUpdateForMidKnot_IntervalError();
		break;
	}
	case PARAM_BASED_BY_INTERVAL_ERROR:
	{
		return selfUpdateForLspia();
		break;
	}
	case UNIFORM_UP_DATE:
	{
		return selfUpdateUniform();
		break;
	}
	case ADJUST_KNOT:
	{
		return adjustKnots();
		break;
	}
	default:
		break;
	}
	return 0.0;
}

Standard_Real KnotUpdate::adjustKnots() {
	// Step 1: Compute differences and find param_max
	Standard_Real max_dist = -1.0;
	Standard_Real param_max = 0.0;
	for (size_t i = 0; i < myParams.size(); ++i) {
		gp_Pnt p_curve = bspline->Value(myParams[i]);
		Standard_Real dist = p_curve.Distance(myPnts[i]);

		if (dist > max_dist) {
			max_dist = dist;
			param_max = myParams[i];
		}
	}
	maxError = max_dist;

	// Step 2: Find the two knots closest to param_max
	// Assume myCurrentKnots is sorted in increasing order
	auto it = std::lower_bound(myCurrentKnots.begin(), myCurrentKnots.end(), param_max);

	size_t index_upper = it - myCurrentKnots.begin(); // Index of the first knot >= param_max
	size_t index_lower = 0;

	if (it == myCurrentKnots.begin()) {
		// param_max is less than the first knot (do not adjust endpoints)
		index_lower = std::numeric_limits<size_t>::max();
		index_upper = (myCurrentKnots.size() > 1) ? 1 : std::numeric_limits<size_t>::max();
	}
	else if (it == myCurrentKnots.end()) {
		// param_max is greater than the last knot (do not adjust endpoints)
		index_upper = std::numeric_limits<size_t>::max();
		index_lower = (myCurrentKnots.size() > 1) ? myCurrentKnots.size() - 2 : std::numeric_limits<size_t>::max();
	}
	else {
		// param_max is between knots
		index_upper = it - myCurrentKnots.begin();
		index_lower = index_upper - 1;
	}

	size_t first_knot_index = 0;
	size_t last_knot_index = myCurrentKnots.size() - 1;

	// Do not adjust endpoint knots
	if (index_lower == first_knot_index) {
		index_lower = std::numeric_limits<size_t>::max();
	}
	if (index_upper == last_knot_index) {
		index_upper = std::numeric_limits<size_t>::max();
	}

	// Step 3: Adjust the knots towards param_max
	if (index_lower != std::numeric_limits<size_t>::max()) {
		Standard_Real knot = myCurrentKnots[index_lower];
		myCurrentKnots[index_lower] = knot + 0.5 * (param_max - knot);
	}
	if (index_upper != std::numeric_limits<size_t>::max()) {
		Standard_Real knot = myCurrentKnots[index_upper];
		myCurrentKnots[index_upper] = knot + 0.5 * (param_max - knot);
	}
	updateSequences();

	return 0.0;
}

Standard_Real KnotUpdate::selfUpdateUniform() {

	return 0.0;
}

Standard_Real KnotUpdate::selfUpdateForLspia()
{
	Standard_Real maxSingleParamError = 0, maxIntervalError = 0, leftKnot = 0, rightKnot = 0, newKnot = 0;
	maxSingleParamError = maxIntervalError = error(myParams[0], myPnts[0]);
	Standard_Integer maxParamIntervalLeftIndex = 0;
	Standard_Integer maxParamIntervalRightIndex = 0;
	Standard_Integer knotIndex;
	Standard_Integer paramIndex = 0;
	//开始遍历节点和参数点，要求传入参数点一定包含于节点。只遍历到最后一个下标的前一个
	for (knotIndex = 0; knotIndex < myCurrentKnots.size() - 1; knotIndex++)
	{
		//获取当前节点区间 [leftKnot,rightKnot)
		leftKnot = myCurrentKnots[knotIndex];
		rightKnot = myCurrentKnots[knotIndex + 1];
		//获取位于当前节点区间内的参数点，并且同步更新最大单点误差，和最大区间误差
		Standard_Real intervalError = 0; //记录当前节点区间误差
		Standard_Integer count = 0;//记录当前区间多少个参数
		while (paramIndex < myParams.size())//确保当前参数下标有效
		{
			if (IsGreater(myParams[paramIndex], rightKnot))//当前参数超过当前区间右端点
			{
				break;
			}
			//当前参数位于区间内
			Standard_Real singleError = error(myParams[paramIndex], myPnts[paramIndex]);
			maxSingleParamError = (maxSingleParamError > singleError) ? maxSingleParamError : singleError;
			intervalError += singleError;
			paramIndex++;
			count++;
		}
		//若参数下标无效（超出索引）

		if (IsGreater(intervalError, maxIntervalError))//该区间为admissible 区间，且最大区间变化，需要更新数据
		{
			maxParamIntervalRightIndex = paramIndex - 1;
			maxParamIntervalLeftIndex = paramIndex - count;
			maxIntervalError = intervalError;
		}
	}

	//全部检查完毕，计算新节点
	for (Standard_Integer i = maxParamIntervalLeftIndex; i <= maxParamIntervalRightIndex; i++)
	{
		newKnot += myParams[i];
	}
	newKnot /= (maxParamIntervalRightIndex - maxParamIntervalLeftIndex + 1);
	//将newKnot插入到当前节点，注意比较是不是和现有节点相等，如果相等需要检测加入后重复度是否超过次数
	auto index = checkNewKnot(newKnot);
	if (index != -1)
	{
		newKnot = myCurrentKnots[index] + (myCurrentKnots[index + 1] - myCurrentKnots[index]) / 20;
	}
	updateSequences(newKnot);
	updateKnotsAndMutis();
	maxError = maxSingleParamError;
	return newKnot;
}

Standard_Real KnotUpdate::selfUpdateForMidKnot(Standard_Boolean isSingle)
{
	Standard_Real maxSingleParamError = 0, leftKnot = 0, rightKnot = 0, newKnot = 0, maxParam = 0;
	newKnot = 0.5;
	maxSingleParamError = error(myParams[0], myPnts[0]);
	Standard_Integer knotIndex = 0, maxKnotIndex = 0;
	Standard_Integer paramIndex = 0;
	//开始遍历节点和参数点，要求传入参数点一定包含于节点。只遍历到最后一个下标的前一个
	for (knotIndex = 0; knotIndex < myCurrentKnots.size() - 1; knotIndex++)
	{
		//获取当前节点区间 [leftKnot,rightKnot)
		leftKnot = myCurrentKnots[knotIndex];
		rightKnot = myCurrentKnots[knotIndex + 1];
		//获取位于当前节点区间内的参数点，并且同步更新最大单点误差，和最大区间误差
		while (paramIndex < myParams.size())//确保当前参数下标有效
		{
			if (IsGreater(myParams[paramIndex], rightKnot))//当前参数超过当前区间右端点
			{
				break;
			}
			//当前参数位于区间内
			Standard_Real singleError = error(myParams[paramIndex], myPnts[paramIndex]);
			if (IsGreater(singleError, maxSingleParamError))
			{
				maxSingleParamError = singleError;
				maxKnotIndex = knotIndex;
				maxParam = myParams[paramIndex];
			}
			paramIndex++;
		}
	}

	if (isSingle)
	{
		if (IsEqual(maxParam, 0))
		{
			newKnot = (myCurrentKnots[0] + myCurrentKnots[1]) / 2;
		}
		else if (IsEqual(maxParam, 1))
		{
			newKnot = (myCurrentKnots[myCurrentKnots.size() - 1] + myCurrentKnots[myCurrentKnots.size() - 2]) / 2;
		}
		else
		{
			newKnot = (myCurrentKnots[maxKnotIndex] + myCurrentKnots[maxKnotIndex + 1]) / 2;
		}
	}
	updateSequences(newKnot);
	updateKnotsAndMutis();
	maxError = maxSingleParamError;
	return newKnot;
}

Standard_Real KnotUpdate::selfUpdateForMidKnot_IntervalError()
{
	Standard_Real maxSingleParamError = 0, maxIntervalError = 0, leftKnot = 0, rightKnot = 0, newKnot = 0;
	maxSingleParamError = maxIntervalError = error(myParams[0], myPnts[0]);
	Standard_Integer knotIndex = 0, maxKnotIndex = 0;
	Standard_Integer paramIndex = 0;
	//开始遍历节点和参数点，要求传入参数点一定包含于节点。只遍历到最后一个下标的前一个
	for (knotIndex = 0; knotIndex < myCurrentKnots.size() - 1; knotIndex++)
	{
		//获取当前节点区间 [leftKnot,rightKnot)
		leftKnot = myCurrentKnots[knotIndex];
		rightKnot = myCurrentKnots[knotIndex + 1];
		//获取位于当前节点区间内的参数点，并且同步更新最大单点误差，和最大区间误差
		Standard_Real intervalError = 0; //记录当前节点区间误差
		while (paramIndex < myParams.size())//确保当前参数下标有效
		{
			if (IsGreater(myParams[paramIndex], rightKnot))//当前参数超过当前区间右端点
			{
				break;
			}
			//当前参数位于区间内
			Standard_Real singleError = error(myParams[paramIndex], myPnts[paramIndex]);
			maxSingleParamError = (maxSingleParamError > singleError) ? maxSingleParamError : singleError;
			intervalError += singleError;
			paramIndex++;
		}
		//若参数下标无效（超出索引）

		if (IsGreater(intervalError, maxIntervalError))//该区间为admissible 区间，且最大区间变化，需要更新数据
		{
			maxIntervalError = intervalError;
			maxKnotIndex = knotIndex;
		}
	}

	//全部检查完毕，计算新节点
	if (IsEqual(maxKnotIndex, myCurrentKnots.size() - 1))
	{
		newKnot = (myCurrentKnots[myCurrentKnots.size() - 1] + myCurrentKnots[myCurrentKnots.size() - 2]) / 2;
	}
	else
	{
		newKnot = (myCurrentKnots[maxKnotIndex] + myCurrentKnots[maxKnotIndex + 1]) / 2;
	}
	updateSequences(newKnot);
	updateKnotsAndMutis();
	maxError = maxSingleParamError;
	return newKnot;
}

Standard_Real KnotUpdate::error(Standard_Real u, const gp_Pnt& P)
{
	return P.Distance(bspline->Value(u));
}

void KnotUpdate::updateKnotsAndMutis()
{
	if (myCurrentSequences.empty()) return;

	std::map<Standard_Real, Standard_Integer> knotMap;

	// 使用map来统计每个节点的重复次数
	for (Standard_Real value : this->myCurrentSequences) {
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
	myCurrentKnots.clear();
	myCurrentMutis.clear();
	// 将map的内容转移到knots和multiplicities向量
	for (const auto& knot : knotMap) {
		myCurrentKnots.push_back(knot.first);
		myCurrentMutis.push_back(knot.second);
	}
}

void KnotUpdate::updateSequences()
{
	myCurrentSequences.clear();
	for (size_t i = 0; i < myCurrentKnots.size(); i++)
	{
		for (Standard_Integer j = 0; j < myCurrentMutis[i]; j++)
		{
			myCurrentSequences.push_back(myCurrentKnots[i]);
		}
	}
}

void KnotUpdate::updateSequences(Standard_Real newKnot)
{
	myCurrentSequences.clear();
	bool notPush = true;
	for (size_t i = 0; i < myCurrentKnots.size(); i++)
	{

		for (Standard_Integer j = 0; j < myCurrentMutis[i]; j++)
		{
			myCurrentSequences.push_back(myCurrentKnots[i]);
		}
		if (notPush && IsGreater(newKnot, myCurrentKnots[i]) && IsLess(newKnot, myCurrentKnots[i + 1]))
		{
			myCurrentSequences.push_back(newKnot);
			notPush = false;
		}
	}
}

Standard_Integer KnotUpdate::checkNewKnot(Standard_Real knot)
{
	for (size_t i = 0; i < myCurrentKnots.size(); i++)
	{
		if (IsEqual(knot, myCurrentKnots[i]))
		{
			return i;
		}
	}
	return -1;
}