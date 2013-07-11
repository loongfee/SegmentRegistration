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

struct _intersection_struct_{
	double theta;
	fPoint pt;
	std::vector<double> ratios;
};

typedef struct _line_struct_
{
	cvline_polar line;
	vector<_intersection_struct_> intersections;
} LineStruct;

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

	AffineLineRegistration( const vector<cvline_polar>& source_lines_polar, const vector<cvline_polar>& reference_lines_polar)
	{
		m_source_lines_polar = source_lines_polar;
		m_reference_lines_polar = reference_lines_polar;
	};

	double m_srcMaxSize;
	double m_refMaxSize;
	vector<fPoint> m_srcIntersections;
	vector<fPoint> m_refIntersections;

	virtual Eigen::VectorXd getParameters()const;

	//virtual Eigen::Vector2d forward(const cvline_polar& l)const;
	cvline_polar forward(const cvline_polar& l)const;
	fPoint forward(const fPoint& pt)const;

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
	virtual double getDeltaSquare(){return m_delta_2;};
	fstream *m_fs;
	
	vector<cvline_polar> m_matched_src_ls;
	vector<cvline_polar> m_matched_ref_ls;
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

	bool intersection(const cvline_polar &l1, const cvline_polar &l2, vector<fPoint>& outPointList, double maxSize = 2000);
	bool intersection(const cvline_polar &l, const vector<cvline_polar> &lineList, vector<LineStruct>& outLineList, double threshold_dist = 2.0);

	static bool intersection_comparison(const _intersection_struct_& node1, const _intersection_struct_& node2);
	//bool sortIntersections(vector<_intersection_struct_>& pts);

	bool ratio_forward_match(const vector<double>& ratioList1, const vector<double>& ratioList2);
	bool ratio_inverse_match(const vector<double>& ratioList1, const vector<double>& ratioList2);
	bool pairwise_distance(const vector<LineStruct>& lineList1, const vector<LineStruct>& lineList2, Eigen::MatrixXd& distMatrix);

	vector<LinePolarPair> m_LinePairList;
	vector<int> m_inliers;

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