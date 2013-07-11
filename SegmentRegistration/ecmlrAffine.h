#pragma once

#include "ecmlr.h"

class ecmlrAffine : public ecmlr
{
public:
	enum AdjustParamIndex
	{
		a0_param = 0,
		a1_param,
		a2_param,
		b0_param,
		b1_param,
		b2_param,
		NUM_ADJUSTABLE_PARAMS // not an index
	};

	ecmlrAffine( const vector<cvline_polar>& source_lines_polar, const vector<cvline_polar>& reference_lines_polar)
		:ecmlr(source_lines_polar, reference_lines_polar)
	{
		initAdjustableParameters();
	};
	~ecmlrAffine(){};

	//virtual vector<double> getParameterScale() {return m_parameter_scale;};

	virtual cvline_polar forward(const cvline_polar& l, Eigen::VectorXd parameters)const;
	virtual cvline_polar forward(const cvline_polar& l)const;

	virtual int get_PARAMETER_NUM() const{return NUM_ADJUSTABLE_PARAMS;};
	virtual void initAdjustableParameters();
	virtual bool parameters_optimization();
	virtual void updateDeltaSquare();
};