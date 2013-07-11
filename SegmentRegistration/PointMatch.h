#pragma once

#ifndef  _POINT_MATCH_H
#define  _POINT_MATCH_H

#include <algorithm>
#include <iostream>
//#include <intsafe.h>
#include <math.h>

#include "lineConstant.h"

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <vector>
using namespace std;

class PointMatch
{
public:
	struct DataStruct{
		PointMatch* pThis;
		vector<PointPair> pointPairList;
		DataStruct()
		{
			pThis = NULL;
		}
	};

	PointMatch();
	~PointMatch(){};

	vector<double> modelParameters;

	fPoint transModel(double x, double y);
	fPoint transModel(fPoint pt);

	fPoint getForwardDeriv( int p , fPoint pt, double pstep_scale);

	static void funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata);
	static void jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata);

	static double huber(double v, double c);

	PointPair findNearestPointPair(const fPoint &basePoint, const vector<fPoint> &pointList);
	void findBestModel(std::vector<PointPair> pointPairList);
	fPoint pointPairDistance(PointPair pointPair);
	void registration(const vector<fPoint>& src_points, const vector<fPoint>& ref_points);
	void registration(std::vector<PointPair> pointPairList);
	void drawResult(const std::vector<PointPair>& pointPairList, const vector<int>& inliers, const char* szFilename);

};

#endif
