#pragma once
#include "lineConstant.h"
#include <vector>
using namespace std;

#include <levmar.h>

#include <armadillo>
using namespace arma;
#include "RANSAC/ransac/ransac.h"
#include "RANSAC/estimators/Solver.h"
#include "RANSAC/estimators/affineSolver.h"
#include "RANSAC/estimators/affineError.h"
#include "RANSAC/estimators/lineAffineSolver.h"
#include "RANSAC/estimators/lineAffineError.h"
using namespace groupsac;

#include "func.h"
#include "draw.h"
#include "ecmlr.h"

#include <iostream>
#include <unsupported/Eigen/NonLinearOptimization>
//using namespace Eigen;
class ecmlr;
struct lmder_functor;

class AffineLineRegistration
{
public :
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
	enum AdjustParamIndex2
	{
		c1_param = 0,
		c2_param,
		c3_param,
		c4_param,
		c5_param,
		c6_param,
	};

	AffineLineRegistration( const vector<cvline_polar>& source_lines_polar, const vector<cvline_polar>& reference_lines_polar)
	{
		m_source_lines_polar = source_lines_polar;
		m_reference_lines_polar = reference_lines_polar;
		m_nX = (int)m_source_lines_polar.size();
		m_nY = (int)m_reference_lines_polar.size();

		m_srcData.resize(m_nX, 2);
		m_refData.resize(m_nY	, 2);
		for (int i = 0;i < m_nX;++i)
		{
			double rho = m_source_lines_polar[i].rho;
			double theta = m_source_lines_polar[i].theta;
			double x = -sin(theta) / (rho + DBL_EPSILON);
			double y = cos(theta) / (rho + DBL_EPSILON);
			m_srcData(i,0) = x;
			m_srcData(i,1) = y;
		}
		for (int i = 0;i < m_nY;++i)
		{
			double rho = m_reference_lines_polar[i].rho;
			double theta = m_reference_lines_polar[i].theta;
			double x = -sin(theta) / (rho + DBL_EPSILON);
			double y = cos(theta) / (rho + DBL_EPSILON);
			m_refData(i,0) = x;
			m_refData(i,1) = y;
		}
		m_ndim = 2;
	};

	virtual Eigen::VectorXd getParameters()const;

	//virtual Eigen::Vector2d forward(const Eigen::Vector2d& pt)const
	//{
	//	Eigen::VectorXd parameters = getParameters();
	//	Eigen::Vector2d newPt;
	//	double m = m_parameters2[c5_param]*pt(0)+m_parameters2[c6_param]*pt(1)+1.0;
	//	newPt(0) = (m_parameters2[c1_param]*pt(0) + m_parameters2[c2_param]*pt(1));
	//	newPt(1) = (m_parameters2[c3_param]*pt(0) + m_parameters2[c4_param]*pt(1));
	//	return newPt;
	//};

	//virtual cvline_polar forward(const cvline_polar& l)const
	//{
	//	Eigen::VectorXd parameters = getParameters();
	//	return forward(l, parameters);
	//};

	//virtual cvline_polar forward(const cvline_polar& l, Eigen::VectorXd parameters)const;

	virtual bool solve(Eigen::VectorXd & parameters)
	{
		bool rst = solve();
		parameters = getParameters();
		return rst;
	};
	virtual bool solve(vector<double> & parameters)
	{
		bool rst = solve();
		parameters.clear();
		Eigen::VectorXd parm = getParameters();
		for(int i = 0;i < (int)parm.size();++i)
		{
			parameters.push_back(parm(i));
		}
		return rst;
	};

	virtual int get_PARAMETER_NUM()const {return NUM_ADJUSTABLE_PARAMS;};
	vector<int> m_z;

	string m_referenceImageFile;
	string m_sourceImageFile;

	virtual void setOutFileStream(fstream *fs){m_fs = fs;};
	virtual void updateDeltaSquare();
	virtual void updateTransLines();
	virtual void updateParameters();
	virtual double getDeltaSquare(){return m_delta_2;};
	fstream *m_fs;
	
	vector<cvline_polar> m_matched_src_ls;
	vector<cvline_polar> m_matched_ref_ls;
	Eigen::MatrixX2d m_srcData;
	Eigen::MatrixX2d m_refData;
	Eigen::MatrixX2d m_transData;
	//Eigen::Matrix<double, -1, 2, 0, -1, -1>
	vector<double> m_parameters;
	int get_nX(){return m_nX;};
	int get_nY(){return m_nY;};
	int get_nDim(){return m_ndim;};

protected:
	//Eigen::VectorXd m_parameters;	// registration parameters
	Eigen::MatrixXd m_alpha;	// the posterior probabilities
	Eigen::MatrixXd m_w;		// virtual observations
	vector<double> m_lambda;	// weights of virtual observations
	vector<cvline_polar> m_source_lines_polar;
	vector<cvline_polar> m_reference_lines_polar;
	vector<cvline_polar> m_trans_lines_polar;

	vector<LinePolarPair> m_LinePairList;
	vector<int> m_inliers;
	vector<double> m_parameters2;

	//Eigen::MatrixXd m_modelData;
	//Eigen::MatrixXd m_transData;
	//Eigen::MatrixXd m_observationData;
	double m_delta_2;
	double m_lammda_sum;

	int m_nX;
	int m_nY;
	int m_ndim;

	virtual bool solve();
	virtual bool initialization();
	virtual void initAdjustableParameters();

	virtual bool E_Step();
	virtual bool parameters_optimization();
	virtual bool classification();

	//bool eigen_levmar_optimization();
	bool linear_optimization();

	correspondence_struct assignment_multi();
	correspondence_struct assignment_one();
	
	//friend struct lmder_functor_affine;
};

//struct lmder_functor_affine : Functor<double>
//{
//public:
//	lmder_functor_affine(int np, int nobs): Functor<double>(np,nobs) {}
//	int operator()(const VectorXd &x, VectorXd &fvec) const;
//	int df(const VectorXd &x, MatrixXd &fjac) const;
//	AffineLineRegistration *pData;
//};