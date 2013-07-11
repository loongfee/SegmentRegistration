#ifndef GROUPSAC_ESTIMATORS_POINTTOLINEDIST_H
#define GROUPSAC_ESTIMATORS_POINTTOLINEDIST_H

#include "armadillo/include/armadillo"
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
double pt2LineDist(const vec & ab,const vec & pt)  {
  return abs(ab(0)*pt(0) - pt(1) + ab(1)) / sqrt(ab(0) * ab(0) + 1);
  //return -1.0f; ///Todo(pmoulon) wait the matrix framework
}

}; // namespace estimators
}; // namespace groupsac

#endif //GROUPSAC_ESTIMATORS_POINTTOLINEDIST_H
