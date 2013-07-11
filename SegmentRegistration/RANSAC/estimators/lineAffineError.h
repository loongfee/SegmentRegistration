#pragma once

#include <armadillo>
using namespace arma;

namespace groupsac  {
namespace estimators  {

/**
* The distance between the given point and the line y=ax+b, i.e. ax-y+b=0
*
* \param[in] ab (ax+b line equation).
* \param[in] pt (The point on which distance will be computed).  
*
* \return The distance from the point to the line.
* http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
*/
inline double lineAffineError(const vec & model,const vec & obs)  {
	//double dist = 0.0;
	//int k = 0;
	//for( k;k < (int)model.n_rows;k++ )
	//{
	//	dist += obs(k) * model(k);
	//}
	
	double rho = obs(0);
	double theta = obs(1);
	double x1 = obs(2);
	double y1 = obs(3);
	double x2 = obs(4);
	double y2 = obs(5);
	double sin_theta = sin(theta);
	double cos_theta = cos(theta);

	double d1 = -x1*sin_theta + y1*cos_theta - rho;
	double d2 = -x2*sin_theta + y2*cos_theta - rho;


	//dist += obs(k);
	return std::sqrt(d1*d1 + d2*d2);
  //return -1.0f; ///Todo(pmoulon) wait the matrix framework
}

}; // namespace estimators
}; // namespace groupsac

