#ifndef GROUPSAC_ESTIMATORS_AFFINE_ERROR_H
#define GROUPSAC_ESTIMATORS_AFFINE_ERROR_H

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
inline double affineError(const vec & model,const vec & obs)  {
	//double dist = 0.0;
	//int k = 0;
	//for( k;k < (int)model.n_rows;k++ )
	//{
	//	dist += obs(k) * model(k);
	//}
	double d1 = 1.0*model(0) + obs(0)*model(1) + obs(1)*model(2) - obs(2);
	double d2 = 1.0*model(3) + obs(0)*model(4) + obs(1)*model(5) - obs(3);

	//dist += obs(k);
	return std::sqrt(d1*d1 + d2*d2);
  //return -1.0f; ///Todo(pmoulon) wait the matrix framework
}

}; // namespace estimators
}; // namespace groupsac

#endif //GROUPSAC_ESTIMATORS_POINTTOLINEDIST_H
