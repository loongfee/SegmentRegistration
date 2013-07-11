#pragma once
#include "lineConstant.h"
#include <vector>
using namespace std;

#include <levmar.h>

#include <armadillo>
using namespace arma;
// General Solver class.
// -- Compute a list of model that fit the candidates data.
class mcempr
{
public :
	mcempr(const arma::mat & modelData, const arma::mat & observationData)
	{
		m_modelData = modelData;
		m_observationData = observationData;
	};

	enum{PARAMETER_NUMBER = 0};


	virtual arma::vec forward(const arma::vec& original)const
	{
		return forward(original, m_parameters);
	};

	virtual arma::vec forward(const arma::vec&, arma::vec parameters)const = 0;
	virtual bool initialParameters() = 0;

	virtual bool solve(arma::vec & parameters)
	{
		bool rst = solve();
		parameters = m_parameters;
		return rst;
	};
	virtual bool solve(vector<double> & parameters)
	{
		bool rst = solve();
		parameters.clear();
		for(int i = 0;i < (int)m_parameters.size();++i)
		{
			parameters.push_back(m_parameters(i));
		}
		return rst;
	};

	virtual int get_PARAMETER_NUM() const =0;
	vector<int> m_z;		//classification

protected:
	arma::vec m_parameters;	// registration parameters
	arma::mat m_alpha;	// the posterior probabilities
	arma::mat m_w;		// virtual observations
	arma::vec m_lambda;	// weights of virtual observations
	vector<arma::mat> m_covariances;	// the covariance matrices
	arma::mat m_modelData;
	arma::mat m_transData;
	arma::mat m_observationData;
	double m_delta_2;

	virtual bool solve();
	virtual bool initialization();

	virtual bool E_Step();
	virtual bool parameters_optimization();
	virtual bool classification();
	virtual bool parameters_refine();
	virtual bool updateCovariances();
	virtual bool MetropolisHastings();
	virtual correspondence_struct getValidCorrespondence();
	virtual correspondence_struct MarkovNext(const correspondence_struct& correspond);
	virtual double getCorrepondenceProbility(const correspondence_struct& correspond);
	virtual double acceptance_ratio(const correspondence_struct& c, const correspondence_struct& c_prime);
	static void funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata);
	static void jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata);

	static double sumVec(const arma::vec& v)
	{
		double sum = 0.0;
		for (int i = 0;i < (int)v.size();++i)
		{
			sum += v(i);
		}
		return sum;
	};

	static double norm2_Vec(const arma::vec& v)
	{
		double sum = 0.0;
		for (int i = 0;i < (int)v.size();++i)
		{
			sum += v(i)*v(i);
		}
		return sqrt(sum);
	}

	static double MahalanobisDistance(const arma::vec& a, const arma::vec& b, const arma::mat& covariance)
	{
		if(a.is_col && b.is_col)
		{
			arma::vec delta = a - b;
			//if(fabs(det(covariance)) < FLOAT_EPSILON)
			//{
				return det(delta.t() * delta);
			//}
			//return det(delta.t() * inv(covariance) * delta);
		}
		else if(a.is_row && b.is_row)
		{
			arma::rowvec delta = a - b;
			if(fabs(det(covariance)) < DBL_EPSILON)
			{
				return det(delta.t() * delta);
			}
			return det(delta * inv(covariance) * delta.t());
		}
		else if(a.is_col && b.is_row)
		{
			return MahalanobisDistance(a, b.t(), covariance);
		}
		else
		{
			return MahalanobisDistance(a.t(), b, covariance);
		}
	}

};