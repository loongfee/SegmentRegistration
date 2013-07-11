#ifndef LINE_CONSTANT_H
#define LINE_CONSTANT_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <intsafe.h>
#include <math.h>
#include <stdio.h>
#include <tchar.h>
#include <float.h>
using namespace std;

//#include "targetver.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include "surf_match.h"

//#include "RANSAC/ransac/ransac.h"
//#include "RANSAC/estimators/Solver.h"
//#include "RANSAC/estimators/affineSolver.h"
//#include "RANSAC/estimators/affineError.h"
//#include <armadillo>
//using namespace groupsac;

#include <Eigen/Dense>
using namespace Eigen;
using Eigen::MatrixXd;

//#include <Eigen/Dense>
//
////using Eigen::MatrixXd;
//using namespace Eigen;
//using namespace Eigen::internal;
//using namespace Eigen::Architecture;

#define PI 3.14159265f
#define MAX_LINE_NUM 1000



/*
 Definitions for unit type.
*/
enum UnitType
{
   UNIT_UNKNOWN    = 0,
   METERS          = 1,
   FEET            = 2,
   US_SURVEY_FEET  = 3,
   DEGREES         = 4,
   RADIANS         = 5,
   NAUTICAL_MILES  = 6,
   SECONDS         = 7,
   MINUTES         = 8,
   PIXEL           = 9,
   MILES           = 10,
   MILLIMETERS     = 11,
   MICRONS         = 12,
   CENTIMETERS     = 13,
   YARDS           = 14,
   INCHES          = 15,
   KILOMETERS      = 16
};


typedef struct POINT_FLOAT
{
	double x;
	double y;
	POINT_FLOAT()
	{
		x = 0.0;
		y = 0.0;
	};
	POINT_FLOAT(double xx, double yy)
	{
		x = xx;
		y = yy;
	};
} fPoint;


typedef struct CORRESPONDENCE_STRUCT
{
	vector<int> forward_map;
	vector<int> inverse_map;

	CORRESPONDENCE_STRUCT()
	{
	};

	CORRESPONDENCE_STRUCT(const CORRESPONDENCE_STRUCT& o)
	{
		this->forward_map = o.forward_map;
		this->inverse_map = o.inverse_map;
	};

	CORRESPONDENCE_STRUCT& operator=(const CORRESPONDENCE_STRUCT& o)
	{
		this->forward_map = o.forward_map;
		this->inverse_map = o.inverse_map;
		return *this;
	};
} correspondence_struct;

/************************************************************************/
/* 
	rho = -x * sin(theta) + y * cos(theta)
*/
/************************************************************************/
typedef struct LINE_POLAR
{
	fPoint pt1;
	fPoint pt2;
	double rho;
	double theta;
} cvline_polar;

typedef struct LINE_POLAR_PAIR
{
	// -¦Ð~+¦Ð
	//¦È
	//¦È'
	cvline_polar line;
	int index;
	cvline_polar line_prime;
	int index_prime;
} LinePolarPair;

typedef struct POINT_PAIR
{
	fPoint point;
	fPoint point_prime;
} PointPair;

#endif