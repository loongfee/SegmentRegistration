#pragma once
#include "lineConstant.h"
#include <vector>
using namespace std;

#include <levmar.h>

#include <armadillo>
using namespace arma;
//#include <nlopt.hpp>
//#include "alglib\optimization.h"
//using namespace alglib;
#include "RANSAC/ransac/ransac.h"
#include "RANSAC/estimators/Solver.h"
#include "RANSAC/estimators/affineSolver.h"
#include "RANSAC/estimators/affineError.h"
#include "RANSAC/estimators/lineAffineSolver.h"
#include "RANSAC/estimators/lineAffineError.h"
using namespace groupsac;

#include "func.h"
#include "draw.h"

#include <iostream>
#include <unsupported/Eigen/NonLinearOptimization>
//using namespace Eigen;
class ecmlr;
struct lmder_functor;

//#define GLOG_NO_ABBREVIATED_SEVERITIES
//#include "ceres/ceres.h"
//#include "gflags/gflags.h"
//#include "glog/logging.h"
//using namespace ceres;
//#include <splm.h>
// General Solver class.
// -- Compute a list of model that fit the candidates data.
#define scale_angle 50.0
#define lammda_rho	10.0
class ecmlr
{
public :
	ecmlr( const vector<cvline_polar>& source_lines_polar, const vector<cvline_polar>& reference_lines_polar);

	enum{PARAMETER_NUMBER = 0};

	virtual Eigen::VectorXd getParameters()const;

	virtual cvline_polar forward(const cvline_polar& l)const
	{
		Eigen::VectorXd parameters = getParameters();
		return forward(l, parameters);
	};

	virtual cvline_polar forward(const cvline_polar& l, Eigen::VectorXd parameters)const = 0;
	//virtual vector<double> getParameterScale() = 0;

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

	virtual int get_PARAMETER_NUM() const =0;
	vector<int> m_z;		//classification

	string m_referenceImageFile;
	string m_sourceImageFile;

	virtual void setOutFileStream(fstream *fs){m_fs = fs;};
	virtual void updateDeltaSquare();
	virtual void updateTransLines();
	virtual double getDeltaSquare(){return m_delta_2;};
	virtual bool setMatched();
	//virtual double distance(const cvline_polar& l1, const cvline_polar& l2, Eigen::VectorXd parameters) = 0;
	virtual double distance(const cvline_polar& l1, const cvline_polar& l2);
	virtual fPoint distance2(const cvline_polar& l1, const cvline_polar& l2);

	virtual double point2polarline(const fPoint& p, const cvline_polar& l);
	virtual fPoint segment2polarline(const cvline_polar& segment, const cvline_polar& line);
	fstream *m_fs;
	
	vector<cvline_polar> m_matched_src_ls;
	vector<cvline_polar> m_matched_ref_ls;
	vector<double> m_parameters;
	int get_nX(){return m_nX;};
	int get_nY(){return m_nY;};
	int get_nDim(){return m_ndim;};
	virtual void setParameters(const vector<double> & parameters);

protected:
	//Eigen::VectorXd m_parameters;	// registration parameters
	Eigen::MatrixXd m_alpha;	// the posterior probabilities
	vector<cvline_polar> m_w;		// virtual observations
	vector<double> m_lambda;	// weights of virtual observations
	vector<Eigen::MatrixXd> m_covariances;	// the covariance matrices
	vector<cvline_polar> m_source_lines_polar;
	vector<cvline_polar> m_reference_lines_polar;
	vector<cvline_polar> m_trans_lines_polar;

	vector<LinePolarPair> m_LinePairList;
	vector<unsigned int> m_inliers;

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
	virtual void initAdjustableParameters() = 0;

	virtual bool E_Step();
	virtual bool parameters_optimization();
	virtual bool classification();
	virtual bool classification1();

	static void funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata);
	static void jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata);

	static void levmar_function_fvec(double *param, double *hx, int nparameter, int nequation, void *adata);
	static void levmar_function_jac(double *param, double *jac, int nparameter, int nequation, void *adata);

	//static void  alglib_function_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr);
	//static void  alglib_function_jac(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr);

	//static void sparseLM_func(double *p, double *hx, int m, int n, void *adata);
	//static void sparseLM_jac(double *p, struct splm_crsm *jac, int m, int n, void *adata);

	//static void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr);
	static double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);

	//bool alglib_optimization();
	//bool spaseLM_optimization();
	bool levmar_optimization();
	bool eigen_levmar_optimization();
	//bool ceres_optimization();

	static double sumVec(const Eigen::VectorXd& v)
	{
		double sum = 0.0;
		for (int i = 0;i < (int)v.size();++i)
		{
			sum += v(i);
		}
		return sum;
	};

	virtual bool MetropolisHastings();
	virtual correspondence_struct getValidCorrespondence();
	virtual correspondence_struct MarkovNext(const correspondence_struct& correspond);
	virtual double getCorrepondenceProbility(const correspondence_struct& correspond);
	virtual double acceptance_ratio(const correspondence_struct& c, const correspondence_struct& c_prime);

	correspondence_struct assignment_multi();
	correspondence_struct assignment_one();
	
	friend class LineMatchCost;
	friend struct lmder_functor;
};

// Generic functor
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct Functor
{
	typedef _Scalar Scalar;
	enum {
		InputsAtCompileTime = NX,
		ValuesAtCompileTime = NY
	};
	typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
	typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
	typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

	const int m_inputs, m_values;

	Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
	Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

	int inputs() const { return m_inputs; }
	int values() const { return m_values; }
	
	// you should define that in the subclass :
	//  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};

struct lmder_functor : Functor<double>
{
public:
	lmder_functor(int np, int nobs): Functor<double>(np,nobs) {}
	int operator()(const VectorXd &x, VectorXd &fvec) const;
	int df(const VectorXd &x, MatrixXd &fjac) const;
	ecmlr *pData;
};