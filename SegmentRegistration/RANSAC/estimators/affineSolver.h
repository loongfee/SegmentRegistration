#ifndef GROUPSAC_ESTIMATORS_AFFINESOLVER_H
#define GROUPSAC_ESTIMATORS_AFFINESOLVER_H

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
class affineSolver : public Solver<T,Model>
{
  /// At least two point are necessary to solve the line equation y = ax+b
  enum { MINIMUM_SAMPLES = 3 };
public :

  int get_MINIMUM_SAMPLES() const {return MINIMUM_SAMPLES;}
  /// See groupsac::estimators::Solver
  bool solve(const T & candidates, vector<Model> & model) const
  {
      // Build matrices to solve Ax = b problem:
      vec b(candidates.n_rows * 2);
      mat A(candidates.n_rows * 2, 6);
      //b = candidates.col(6);
	  for(int i = 0;i < (int)candidates.n_rows;++i)
	  {
		  A(2*i, 0) = 1.0;
		  A(2*i, 1) = candidates(i, 0);
		  A(2*i, 2) = candidates(i, 1);
		  A(2*i, 3) = 0.0;
		  A(2*i, 4) = 0.0;
		  A(2*i, 5) = 0.0;
		  b(2*i) = candidates(i, 2);

		  A(2*i+1, 0) = 0.0;
		  A(2*i+1, 1) = 0.0;
		  A(2*i+1, 2) = 0.0;
		  A(2*i+1, 3) = 1.0;
		  A(2*i+1, 4) = candidates(i, 0);
		  A(2*i+1, 5) = candidates(i, 1);
		  b(2*i+1) = candidates(i, 3);
	  }
	  //A.col(0) = candidates.col(0);
	  //A.col(1) = candidates.col(1);
	  //A.col(2) = candidates.col(2);
	  //A.col(3) = candidates.col(3);
	  //A.col(4) = candidates.col(4);
	  //A.col(5) = candidates.col(5);

      // Compute least-squares solution:
      vec X;

	  //arma::mat AA = trans(A)*A;
	  //A.save("debug_A.txt", arma::arma_ascii);
	  //AA.save("debug_AA.txt", arma::arma_ascii);
	  //b.save("debug_b.txt", arma::arma_ascii);
	  //double detA = det(AA);
	  //if(fabs(detA) < FLOAT_EPSILON)
	  //{
		 // cout << "ERROR : cannot find a solution" << endl;
		 // return false;
	  //}
	  //arma::mat bb = trans(A)*b;
	  //X = inv(AA) * bb;
	  //X.save("debug_X.txt", arma::arma_ascii);
	  //return true;

      if ( ::solve(X, A, b) )
      {
		//  X.save("debug_X.txt", arma::arma_ascii);
        model.push_back(X);
        return true;
      }
      else
      {
		  cout << "ERROR : cannot find a solution" << endl;
		  candidates.save("debug_candidates.txt", arma::arma_ascii);
		  A.save("debug_A.txt", arma::arma_ascii);
		  arma::mat AA = trans(A)*A;
		  AA.save("debug_AA.txt", arma::arma_ascii);
        return false;
      }
      return false;
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
        double dist = affineError( modelToTest, trans(candidates.row(j)) );

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

#endif //GROUPSAC_ESTIMATORS_LINEFITTINGSOLVER_H
