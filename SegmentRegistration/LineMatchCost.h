#pragma once
#include "ecmlr.h"
class LineMatchCost
	//: public ceres::CostFunction
{
	public:
		LineMatchCost(ecmlr* pEcmlr) {m_pEcmlr = pEcmlr;};
	ecmlr* m_pEcmlr;
	virtual ~LineMatchCost() {}
	virtual bool Evaluate(double const* const* parameters,
	double* residuals,
	double** jacobians) const;
};