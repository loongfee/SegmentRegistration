#pragma once

#include "lineConstant.h"

#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <vector>
using namespace cv;
using namespace std;

bool surf_extract(string srcImage, vector<fPoint> &pointList, CvSURFParams params = cvSURFParams(500, 1));
bool surf_match(const char* srcImage, const char* refImage, vector<PointPair>& pointPairList, CvSURFParams params = cvSURFParams(500, 1));

bool harris_extract(const char* srcImage, vector<fPoint> &pointList, int threshold=10, bool nonmaxSuppression=true);

bool orb_extract(const char* srcImage, vector<fPoint> &pointList);
bool orb_match(const char* srcImage, const char* refImage, vector<PointPair>& pointPairList);