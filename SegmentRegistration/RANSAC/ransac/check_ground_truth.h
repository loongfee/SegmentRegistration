#pragma once

#include <vector>

//  check whether the inliers are compatible with ground truth
bool check_ground_truth(const vector<bool> & is_inlier, const vector<int> & inliers, double prop = 0.8)
{
  return fail = inliers.size() < count(is_inlier.begin(), is_inlier.end(), true) * prop;	
}