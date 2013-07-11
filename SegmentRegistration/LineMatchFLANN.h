#pragma once

#include "lineFeature.h"
#include "cv.h"
#include <stdio.h>
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "lbd/LineDescriptor.hh"
using namespace cv;

//void LineMatch(IplImage* ref_seg,IplImage* src_seg,vector<vector<double>> &ref_LineFeature,vector<vector<double>> &src_LineFeature,
//	vector<int>&ref_LineNum,vector<int>&src_LineNum,vector<cvline_polar>&ref_lines_list,vector<cvline_polar>&src_lines_list);

void LineMatch(vector<vector<double>> &ref_LineFeature,vector<vector<double>> &src_LineFeature,
	vector<int>&ref_LineNum,vector<int>&src_LineNum,vector<cvline_polar>&ref_lines_list,vector<cvline_polar>&src_lines_list,
	vector<cvline_polar>& matched_ref_lines_list, vector<cvline_polar>& matched_src_lines_list);

void LineMatch(const ScaleLines &ref_LineFeature, const ScaleLines &src_LineFeature,
			   vector<cvline_polar>& matched_ref_lines_list, vector<cvline_polar>& matched_src_lines_list);