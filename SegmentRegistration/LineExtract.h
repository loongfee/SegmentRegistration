#pragma once

#include "lineConstant.h"
#include "lp\LineMatcher.h"

IplImage* get_lines(IplImage* img,vector<cvline_polar>& vec_lines, double scale = 0.8, int seletion_num = 200);
IplImage* get_lines(const char* strImage, vector<cvline_polar>& vec_lines, double scale = 0.8, int seletion_num = 200);

bool isBoundary(IplImage* img, cvline_polar line, cv::Scalar bgColor = CV_RGB(0,0,0));
bool isBoundary(cv::Mat img, cvline_polar line, cv::Scalar bgColor = CV_RGB(0,0,0));


bool GetImGray(IplImage* im, double x, double y, double &gray);
bool lp_LineDetect(const char* strImage, std::vector<Line>& lineList, double scale = 0.8, int seletion_num = 200);

IplImage* get_lines_lswms(const char* strImage, vector<cvline_polar>& vec_lines, int numMaxLSegs = 200);