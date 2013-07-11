#ifndef GROUPSAC_ESTIMATORS_LINEFITTINGSOLVER_H
#define GROUPSAC_ESTIMATORS_LINEFITTINGSOLVER_H

#include <cassert>
#include <iostream>
#include <vector>
#include "armadillo/include/armadillo"
#include "pointToLineDist.h"
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
class lineFittingSolver : public Solver<T,Model>
{
  /// At least two point are necessary to solve the line equation y = ax+b
  enum { MINIMUM_SAMPLES = 2 };
public :

  int get_MINIMUM_SAMPLES() const {return MINIMUM_SAMPLES;}
  /// See groupsac::estimators::Solver
  bool solve(const T & candidates, vector<Model> & model) const
  {
      // Build matrices to solve Ax = b problem:
      vec b(candidates.n_rows);
      mat A(candidates.n_rows, 2);
      b = candidates.col(1);
      A.col(0) = candidates.col(0);
      A.col(1).fill(1.0);

      // Compute least-squares solution:
      vec X;
      if ( ::solve(X, A, b) )
      {
        model.push_back(X);
        return true;
      }
      else
      {
        cout << "ERROR : cannot find a solution" << endl;
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
        double dist = pt2LineDist( modelToTest, trans(candidates.row(j)) );
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
