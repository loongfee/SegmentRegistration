#pragma once
#include <cassert>
#include <iostream>
#include <vector>
#include <armadillo>
#include "affineError.h"
#include "Solver.h"
using namespace std;
using namespace arma;

namespace groupsac  {
namespace estimators  {

/// Fit a 2D line with to a set of points
/// Specifically, find a and b in the model y = ax + b
///
/// Input data must be typed as follow :
/// X0 Y0
/// X1 Y1
/// X... Y ...
/// Internal extractor function allow to extract sampled data.
template<typename T = mat, typename Model = vec>
class lineAffineSolver : public Solver<T,Model>
{
  /// At least two point are necessary to solve the line equation y = ax+b
  enum { MINIMUM_SAMPLES = 3 };
public :

  int get_MINIMUM_SAMPLES() const {return MINIMUM_SAMPLES;}
  /// See groupsac::estimators::Solver
  bool solve(const T & candidates, vector<Model> & model) const
  {
      // Build matrices to solve Ax = b problem:

	  //vec b(candidates.n_rows * 2);
      //mat A(candidates.n_rows * 2, 6);
      //b = candidates.col(6);
	VectorXd b(6);
	MatrixXd A(6, 6);
	b.fill(0.0);
	A.fill(0.0);
	for(int i = 0;i < (int)candidates.n_rows;++i)
	{
		double rho = candidates(i, 0);
		double theta = candidates(i, 1);
		double x1 = candidates(i, 2);
		double y1 = candidates(i, 3);
		double x2 = candidates(i, 4);
		double y2 = candidates(i, 5);
		double sin_theta = sin(theta);
		double cos_theta = cos(theta);

		
		vector<double> coeff(6);
		coeff[0] = -sin_theta;
		coeff[1] = -sin_theta * x1;
		coeff[2] = -sin_theta * y1;
		coeff[3] = cos_theta;
		coeff[4] = cos_theta * x1;
		coeff[5] = cos_theta * y1;

		for(int p1=0;p1<6;++p1)
		{        
			b[p1] += coeff[p1] * rho;
			for(int p2=0;p2<6;++p2)
			{
				A(p1,p2) += coeff[p1] * coeff[p2];
			}
		}

		coeff[0] = -sin_theta;
		coeff[1] = -sin_theta * x2;
		coeff[2] = -sin_theta * y2;
		coeff[3] = cos_theta;
		coeff[4] = cos_theta * x2;
		coeff[5] = cos_theta * y2;
		for(int p1=0;p1<6;++p1)
		{        
			b[p1] += coeff[p1] * rho;
			for(int p2=0;p2<6;++p2)
			{
				A(p1,p2) += coeff[p1] * coeff[p2];
			}
		}
	}
	
	//Eigen::MatrixXd damper = MatrixXd::Identity(6, 6);
	//A = A + damper * 1.0E-6;
	double tmp = A.determinant();
	if ( fabs(tmp) > DBL_EPSILON )
	{
		Eigen::VectorXd parameter = A.inverse() * b;
		vec X(6);
		for(int i=0;i<6;++i) X[i] = parameter[i];
        model.push_back(X);
		return true;
	}
	else
	{
		cout << "ERROR : cannot find a solution" << endl;
		return false;
	}
  }

  /**
  * Return the candidates that are estimated as inliers to the best model
  *
  * \param[in] model (The model(s) that fit the data).
  * \param[in] candidates (The input data).
  * \param[in] threshold (To evaluate if the candidates is an Inlier or Outlier)
  *
  * \return The list of point that are considered as inliers
  */
  static vector<int> defaultEvaluator(vector<Model> & model,
                                      const T & candidates,
                                      double threshold)
  {
	  if(model.size() < 1)
	  {
		  vector<int> inlier;
		  return inlier;
	  }
    assert(model.size() > 0);
    vector< vector<int> > inliers(model.size());
    int bestIndex = 0;
    // For each model compute the number of inliers and the index of the inliers
    // Return the longest inliers vector.
    // must use the pointToLineDist.h file.
    for (size_t i = 0; i < model.size(); ++i)
    {
      const Model & modelToTest = model[i];
      for (size_t j = 0; j < candidates.n_rows; ++j)
      {
        double dist = lineAffineError( modelToTest, trans(candidates.row(j)) );

        if ( abs(dist) < threshold)
          inliers[i].push_back(j);
      }
      if ( i > 0 && inliers[bestIndex].size() < inliers[i].size())
      {
        bestIndex = i;
      }
    }
    return inliers[bestIndex];
  }

  /**
    * Extract the sampled indices from the data container.
    *
    * \param[in] data (The input data).
    * \param[in] samples (The indices of data to extract (line or row)).
    *
    * \return The sampled data.
    */
  static T extractor(const T & data, const vector<int> & sampled)
  {
    mat test;
    test.zeros(sampled.size(), data.n_cols);
    for(size_t i=0; i < sampled.size(); ++i)
      test.row(i) = data.row( sampled[i] );
    return test;
  }
};

}; // namespace estimators
}; // namespace groupsac

