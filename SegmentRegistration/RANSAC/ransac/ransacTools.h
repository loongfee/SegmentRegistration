#pragma once
#include <cassert>

namespace groupsac  {
namespace ransac  {

  
//compute the threshold used in RANSAC
inline double ransac_threshold(int codimension,float sigma)
{
  double threshold = 0.0;
  double sigma2 = sigma * sigma;
  switch (codimension)
  {
    case 1:
      threshold = 3.84 * sigma2;
      break;
    case 2:
      threshold = 5.99 * sigma2;
      break;
    case 3:
      threshold = 7.81 * sigma2;
      break;
    case 4:
      threshold = 9.49 * sigma2;
      break;
    default:
      assert(codimension <= 4);
  }
  return threshold;
}

}; // namespace ransac
}; // namespace groupsac
