#pragma once

#include "ecmlr.h"

class ecmlrPerspective : public ecmlr
{
public:
	enum AdjustParamIndex
	{
		M11_param = 0,
		M12_param,
		M13_param,
		M21_param,
		M22_param,
		M23_param,
		M31_param,
		M32_param,
		NUM_ADJUSTABLE_PARAMS // not an index
	};

	ecmlrPerspective( const vector<cvline_polar>& reference_lines_polar, const vector<cvline_polar>& source_lines_polar)
		:ecmlr(reference_lines_polar, source_lines_polar)
	{
	};
	~ecmlrPerspective(){};

	//virtual vector<double> getParameterScale() {return m_parameter_scale;};

	virtual cvline_polar forward(const cvline_polar& l, Eigen::VectorXd parameters)const;
	virtual cvline_polar forward(const cvline_polar& l)const;

	virtual int get_PARAMETER_NUM() const{return NUM_ADJUSTABLE_PARAMS;};
	virtual void initAdjustableParameters();
	virtual bool parameters_optimization();
	virtual void updateDeltaSquare();

	bool eigen_levmar_optimization();

	friend struct lmder_functor_perspective;
};

struct lmder_functor_perspective : Functor<double>
{
public:
	lmder_functor_perspective(int np, int nobs): Functor<double>(np,nobs) {}
	int operator()(const VectorXd &x, VectorXd &fvec) const;
	int df(const VectorXd &x, MatrixXd &fjac) const;
	ecmlrPerspective *pData;
};