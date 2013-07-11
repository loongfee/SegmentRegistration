#pragma once

#ifndef  _LINE_MATCH_H
#define  _LINE_MATCH_H

#include <algorithm>
#include <iostream>
#include <intsafe.h>
#include <math.h>

#include "lineConstant.h"

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <vector>
using namespace std;

class LineMatch
{
public:
	struct DataStruct{
		LineMatch* pThis;
		vector<LinePolarPair> linePairList;
		DataStruct()
		{
			pThis = NULL;
		}
	};
	LineMatch(const vector<cvline_polar>& src_lines, const vector<cvline_polar>& ref_lines);
	~LineMatch(){};
	vector<cvline_polar> m_src_lines;
	vector<cvline_polar> m_ref_lines;

	vector<double> modelParameters;

	fPoint transModel(double x, double y);

	fPoint getForwardDeriv( int p , fPoint pt, double pstep_scale);

	static void funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata);
	static void jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata);

	static double huber(double v, double c);

	LinePolarPair findNearestLinePair(const cvline_polar &baseLine, const vector<cvline_polar> &lineList);
	void findBestModel(std::vector<LinePolarPair> linePairList);
	fPoint linePairDistance(LinePolarPair linePair);
	void registration();

};

#endif
