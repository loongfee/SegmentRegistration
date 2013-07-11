#ifndef GROUPSAC_ESTIMATORS_FUNDAMENTAL_7PT_H
#define GROUPSAC_ESTIMATORS_FUNDAMENTAL_7PT_H

#include <cassert>
#include <iostream>
#include <set>
#include <vector>
#include "armadillo/include/armadillo"
#include "pointToLineDist.h"
#include "Solver.h"
using namespace std;
using namespace arma;
#define _USE_MATH_DEFINES
#include <math.h>

// This Fundamental matrix solver is based on LIBMV open Source
// computer vision library.
// http://code.google.com/p/libmv/

namespace groupsac  {
namespace estimators  {

// Solve the cubic polynomial
//
//   x^3 + a*x^2 + b*x + c = 0
//
// The number of roots (from zero to three) is returned. If the number of roots
// is less than three, then higher numbered x's are not changed. For example,
// if there are 2 roots, only x0 and x1 are set.
//
// The GSL cubic solver was used as a reference for this routine.
// TODO(pmoulon) move to poly.h
template<typename Real>
int SolveCubicPolynomial(Real a, Real b, Real c,
                         Real *x0, Real *x1, Real *x2) {
  Real q = a * a - 3 * b;
  Real r = 2 * a * a * a - 9 * a * b + 27 * c;

  Real Q = q / 9;
  Real R = r / 54;

  Real Q3 = Q * Q * Q;
  Real R2 = R * R;

  Real CR2 = 729 * r * r;
  Real CQ3 = 2916 * q * q * q;

  if (R == 0 && Q == 0) {
    // Tripple root in one place.
    *x0 = *x1 = *x2 = -a / 3 ;
    return 3;

  } else if (CR2 == CQ3) {
    // This test is actually R2 == Q3, written in a form suitable for exact
    // computation with integers.
    //
    // Due to finite precision some double roots may be missed, and considered
    // to be a pair of complex roots z = x +/- epsilon i close to the real
    // axis.
    Real sqrtQ = sqrt (Q);
    if (R > 0) {
      *x0 = -2 * sqrtQ - a / 3;
      *x1 =      sqrtQ - a / 3;
      *x2 =      sqrtQ - a / 3;
    } else {
      *x0 =     -sqrtQ - a / 3;
      *x1 =     -sqrtQ - a / 3;
      *x2 =  2 * sqrtQ - a / 3;
    }
    return 3;

  } else if (CR2 < CQ3) {
    // This case is equivalent to R2 < Q3.
    Real sqrtQ = sqrt (Q);
    Real sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
    Real theta = acos (R / sqrtQ3);
    Real norm = -2 * sqrtQ;
    *x0 = norm * cos (theta / 3) - a / 3;
    *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
    *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;

    // Put the roots in ascending order.
    if (*x0 > *x1) {
      std::swap(*x0, *x1);
    }
    if (*x1 > *x2) {
      std::swap(*x1, *x2);
      if (*x0 > *x1) {
        std::swap(*x0, *x1);
      }
    }
    return 3;
  }
  Real sgnR = (R >= 0 ? 1 : -1);
  Real A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
  Real B = Q / A ;
  *x0 = A + B - a / 3;
  return 1;
}

// The coefficients are in ascending powers, i.e. coeffs[N]*x^N.
template<typename Real>
int SolveCubicPolynomial(const Real *coeffs, Real *solutions) {
  if (coeffs[0] == 0.0) {
    // This is a quadratic not a cubic. Implement a quadratic
    // solver!
    return 0;
  }
  Real a = coeffs[2] / coeffs[3];
  Real b = coeffs[1] / coeffs[3];
  Real c = coeffs[0] / coeffs[3];
  return SolveCubicPolynomial(a, b, c,
                              solutions + 0,
                              solutions + 1,
                              solutions + 2);
}

// Return the square of a number.
template<typename T>
inline T Square(T x) {
  return x * x;
}

struct SampsonError {
  static double Error(const mat &F, const vec &x1, const vec &x2) {
    vec x(3),y(3);
    x(0) = x1(0); x(1) = x1(1); x(2) = 1.0;
    y(0) = x2(0); y(1) = x2(1); y(2) = 1.0;
    swap(x,y);
    // See page 287 equation (11.9) of HZ.
    vec F_x = F * x;
    vec Ft_y = trans(F) * y;
    return Square(dot(y,F_x)) /
      (Square(F_x(0)) + Square(F_x(1))
      + Square(Ft_y(0)) + Square(Ft_y(1)));
  }
};

struct SymmetricEpipolarDistanceError {
  static double Error(const mat &F, const vec &x1, const vec &x2) {
  vec x(3),y(3);
    x(0) = x1(0); x(1) = x1(1); x(2) = 1.0;
    y(0) = x2(0); y(1) = x2(1); y(2) = 1.0;
    swap(x,y);
    // See page 288 equation (11.10) of HZ.
    vec F_x = F * x;
    vec Ft_y = trans(F) * y;
    return Square(dot(y,F_x)) *
      (1.0/(Square(F_x(0)) + Square(F_x(1)))
      + 1.0/(Square(Ft_y(0)) + Square(Ft_y(1))))
      / 4.0;// The divide by 4 is to make this match the sampson distance.
  }
};

struct EpipolarDistanceError {
  static double Error(const mat &F, const vec &x1, const vec &x2) {
	vec x(3),y(3);
    x(0) = x1(0); x(1) = x1(1); x(2) = 1.0;
    y(0) = x2(0); y(1) = x2(1); y(2) = 1.0;
    // See page 287 equation (11.9) of HZ.
    return dot( x, F * y);
  }
};

/// Compute the fundamental matrix from 7 points.
/// It's the minimal case.
/// It could compute one to three solution.
///
/// Input data must be typed as follow (Image a and Image B):
/// Xa0 Xa1
/// Ya0 Ya1
/// Xb0 Xb1
/// Yb0 Yb1
/// Internal extractor function allow to extract sampled data.
template<typename T = mat, typename Model = mat>
class Fundamental7ptSolver : public Solver<T,Model>
{
  /// At least seven points are necessary to solve F.
  enum { MINIMUM_SAMPLES = 7 };
public :

  int get_MINIMUM_SAMPLES() const {return MINIMUM_SAMPLES;}
  /// See groupsac::estimators::Solver
  bool solve(const T & candidates, vector<Model> & model) const
  {
      // Set up the homogeneous system Af = 0 from the equations x'T*F*x = 0.
      mat A(candidates.n_cols, 9);
      for (int i = 0; i < candidates.n_cols; ++i) {
        A(i, 0) = candidates(2, i) * candidates(0, i);  // 0 represents x coords,
        A(i, 1) = candidates(2, i) * candidates(1, i);  // 1 represents y coords.
        A(i, 2) = candidates(2, i);
        A(i, 3) = candidates(3, i) * candidates(0, i);
        A(i, 4) = candidates(3, i) * candidates(1, i);
        A(i, 5) = candidates(3, i);
        A(i, 6) = candidates(0, i);
        A(i, 7) = candidates(1, i);
        A(i, 8) = 1.0;
      }

      // Find the two F matrices in the nullspace of A.
      mat U,V;
      colvec s;
      if (svd(U,s,V,A))
      {
        mat F1 = reshape(V.col(A.n_cols - 1),3,3);
        mat F2 = reshape(V.col(A.n_cols - 2),3,3);

        // Then, use the condition det(F) = 0 to determine F. In other words, solve
        // det(F1 + a*F2) = 0 for a.
        double a = F1(0, 0), j = F2(0, 0),
               b = F1(0, 1), k = F2(0, 1),
               c = F1(0, 2), l = F2(0, 2),
               d = F1(1, 0), m = F2(1, 0),
               e = F1(1, 1), n = F2(1, 1),
               f = F1(1, 2), o = F2(1, 2),
               g = F1(2, 0), p = F2(2, 0),
               h = F1(2, 1), q = F2(2, 1),
               i = F1(2, 2), r = F2(2, 2);

        // Run fundamental_7point_coeffs.py to get the below coefficients.
        // The coefficients are in ascending powers of alpha, i.e. P[N]*x^N.
        double P[4] = {
          a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g,
          a*e*r + a*i*n + b*f*p + b*g*o + c*d*q + c*h*m + d*h*l + e*i*j + f*g*k -
          a*f*q - a*h*o - b*d*r - b*i*m - c*e*p - c*g*n - d*i*k - e*g*l - f*h*j,
          a*n*r + b*o*p + c*m*q + d*l*q + e*j*r + f*k*p + g*k*o + h*l*m + i*j*n -
          a*o*q - b*m*r - c*n*p - d*k*r - e*l*p - f*j*q - g*l*n - h*j*o - i*k*m,
          j*n*r + k*o*p + l*m*q - j*o*q - k*m*r - l*n*p,
        };

        // Solve for the roots of P[3]*x^3 + P[2]*x^2 + P[1]*x + P[0] = 0.
        double roots[3];
        int num_roots = SolveCubicPolynomial(P, roots);

        // Build the fundamental matrix for each solution.
        for (int kk = 0; kk < num_roots; ++kk)  {
          model.push_back(F1 + roots[kk] * F2);
          model[kk] /= model[kk](2,2);
        }
        return true;
      }
      else
      {
        cout << "ERROR : cannot find a solution" << endl;
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
    vector< vector<int> > vec_inliers(model.size());
    int bestIndex = 0;
    // For each model compute the number of inliers and the index of the inliers
    // Return the longest inliers vector.
    for (size_t i = 0; i < model.size(); ++i)
    {
      const Model & modelToTest = model[i];
      for (size_t j = 0; j < candidates.n_cols; ++j)
      {
        vec x1 = candidates.submat(0,j,1,j);
        vec x2 = candidates.submat(2,j,3,j);
        double dist = SampsonError::Error(modelToTest,x1,x2);
        if ( abs(dist) < threshold)
          vec_inliers[i].push_back(j);
      }
      if ( i > 0 && vec_inliers[bestIndex].size() < vec_inliers[i].size())
        bestIndex = i;
    }
    return vec_inliers[bestIndex];
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
    mat test = zeros(data.n_rows,sampled.size());
    for(size_t i = 0; i < sampled.size(); ++i)
    {
      test.col(i) = data.col( sampled[i] );
    }
    return test;
  }
};

}; // namespace estimators
}; // namespace groupsac

#endif //GROUPSAC_ESTIMATORS_FUNDAMENTAL_7PT_H
