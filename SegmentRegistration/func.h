#pragma once

#include "lineConstant.h"
#include <Eigen/Eigen>

static bool approximate(double a, double b, double threshold = 1e-6)
{
	if (fabs(a-b) < threshold)
	{
		return true;
	}
	return false;
};

bool fexists(const char *filename);

bool fexists(string filename);

void normalize_segment(cvline_polar& segment, double len = 50.0);

double segment_length(const cvline_polar& l);

cvline_polar segment_move_vertical(const cvline_polar& l, double dist);

fPoint line2point(const cvline_polar& l);

double distance_(const fPoint& pt1, const fPoint& pt2);

cv::Scalar random_color(CvRNG* rng);

/************************************************************************/
/* 
	rho = -x * sin(theta) + y * cos(theta)
*/
/************************************************************************/
cvline_polar line2polar(CvPoint* line);

cvline_polar line2polar(fPoint* line);

cvline_polar line2polar(fPoint pt1, fPoint pt2);

cvline_polar mergeline(cvline_polar l1, cvline_polar l2);

bool pt_comparison(const fPoint& node1, const fPoint& node2);
bool sortPoints(vector<fPoint>& pts);

void drawMatchedSegments(string foldName, string sourceName, string referenceName, string linesName);
void drawMatchedSegments2(string foldName, string sourceName, string referenceName, string linesName);

void drawResult(const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const char* szFilename);

//void drawResultLine(IplImage* pic, const cvline_polar& line, CvScalar color, const char* strId, int offset = 100, int thickness = 1);
void drawResultLine(cv::Mat& pic, const cvline_polar& line, cv::Scalar color, const char* strId, int offset = 0, int thickness = 1);

void drawResult(Eigen::MatrixX2d ptList1, Eigen::MatrixX2d ptList2, const char* szFilename);

void drawLastResult(const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const vector<int>& matcher, const char* szFilename1, const char* szFilename2);
void drawLastResult(string strImage1, string strImage2, vector<double> param,
					const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const vector<int>& matcher, const char* szFilename1, const char* szFilename2);
void drawLastResult(string strImage1, string strImage2,
					const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const vector<int>& matcher, const char* szFilename1, const char* szFilename2);


void drawPolarResult(const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const char* szFilename);

cv::Size getNewSize(int iRows, int iCols, cv::Mat transMat);
cv::Size getPerspectiveNewSize(int iRows, int iCols, cv::Mat transMat);

Eigen::VectorXd affine_ransac(vector<cvline_polar> matched_ref_lines_list,
	vector<cvline_polar> matched_src_lines_list, vector<int> &inliers);

template<class T>
static void affine_warp_Image(string image_in, string image_out, T transParameter)
{
	cv::Mat img = cv::imread(image_in);
	int iRows = img.rows;
	int iCols = img.cols;

	cv::Mat affineTrans(2, 3, CV_64F);

	//for(int i = 0;i < 2;i++)
	//{
	//	for(int j = 0;j < 3;j++)
	//	{
	//		affineTrans.at<double>(i, (j-1)%3) = transParameter[i * 3 + j];
	//	}
	//}
	affineTrans.at<double>(0, 0) = transParameter[1];
	affineTrans.at<double>(0, 1) = transParameter[2];
	affineTrans.at<double>(0, 2) = transParameter[0];
	affineTrans.at<double>(1, 0) = transParameter[4];
	affineTrans.at<double>(1, 1) = transParameter[5];
	affineTrans.at<double>(1, 2) = transParameter[3];


	//double mm = transParameter[0] * transParameter[4] - transParameter[1] * transParameter[3];
	//affineTrans.at<double>(0, 0) = transParameter[4] / (mm + FLOAT_EPSILON);
	//affineTrans.at<double>(0, 1) = -transParameter[1] / (mm + FLOAT_EPSILON);
	//affineTrans.at<double>(0, 2) = (transParameter[1]*transParameter[5] - transParameter[2]*transParameter[4]) / (mm + FLOAT_EPSILON);
	//affineTrans.at<double>(1, 0) = transParameter[3] / (-mm + FLOAT_EPSILON);
	//affineTrans.at<double>(1, 1) = -transParameter[0] / (-mm + FLOAT_EPSILON);
	//affineTrans.at<double>(1, 2) = (transParameter[0]*transParameter[5] - transParameter[2]*transParameter[3]) / (-mm + FLOAT_EPSILON);

	cv::Size newSize = getNewSize(iRows, iCols, affineTrans);

	cv::Mat newImg;
	warpAffine(img, newImg, affineTrans , newSize, cv::INTER_CUBIC);

	//ÏÔÊ¾Í¼Ïñ
	//imshow("image", newImg);

	imwrite(image_out, newImg);
}

template<class T>
static void invert_affine_warp_Image(string image_in, string image_out, T transParameter)
{
	cv::Mat img = cv::imread(image_in);
	int iRows = img.rows;
	int iCols = img.cols;

	cv::Mat affineTrans(2, 3, CV_64F);

	//for(int i = 0;i < 2;i++)
	//{
	//	for(int j = 0;j < 3;j++)
	//	{
	//		affineTrans.at<double>(i, (j-1)%3) = transParameter[i * 3 + j];
	//	}
	//}
	affineTrans.at<double>(0, 0) = transParameter[1];
	affineTrans.at<double>(0, 1) = transParameter[2];
	affineTrans.at<double>(0, 2) = transParameter[0];
	affineTrans.at<double>(1, 0) = transParameter[4];
	affineTrans.at<double>(1, 1) = transParameter[5];
	affineTrans.at<double>(1, 2) = transParameter[3];

	cv::Mat invertAffineTrans;
	invertAffineTransform(affineTrans, invertAffineTrans);


	//double mm = transParameter[0] * transParameter[4] - transParameter[1] * transParameter[3];
	//affineTrans.at<double>(0, 0) = transParameter[4] / (mm + FLOAT_EPSILON);
	//affineTrans.at<double>(0, 1) = -transParameter[1] / (mm + FLOAT_EPSILON);
	//affineTrans.at<double>(0, 2) = (transParameter[1]*transParameter[5] - transParameter[2]*transParameter[4]) / (mm + FLOAT_EPSILON);
	//affineTrans.at<double>(1, 0) = transParameter[3] / (-mm + FLOAT_EPSILON);
	//affineTrans.at<double>(1, 1) = -transParameter[0] / (-mm + FLOAT_EPSILON);
	//affineTrans.at<double>(1, 2) = (transParameter[0]*transParameter[5] - transParameter[2]*transParameter[3]) / (-mm + FLOAT_EPSILON);

	cv::Size newSize = getNewSize(iRows, iCols, invertAffineTrans);

	cv::Mat newImg;
	warpAffine(img, newImg, invertAffineTrans , newSize, cv::INTER_CUBIC);

	//ÏÔÊ¾Í¼Ïñ
	//imshow("image", newImg);

	imwrite(image_out, newImg);
}

template<class T>
static cv::Mat affine_warp_Image(string image_in, T transParameter)
{
	cv::Mat img = cv::imread(image_in);
	int iRows = img.rows;
	int iCols = img.cols;

	cv::Mat affineTrans(2, 3, CV_64F);

	//for(int i = 0;i < 2;i++)
	//{
	//	for(int j = 0;j < 3;j++)
	//	{
	//		affineTrans.at<double>(i, (j-1)%3) = transParameter[i * 3 + j];
	//	}
	//}
	affineTrans.at<double>(0, 0) = transParameter[1];
	affineTrans.at<double>(0, 1) = transParameter[2];
	affineTrans.at<double>(0, 2) = transParameter[0];
	affineTrans.at<double>(1, 0) = transParameter[4];
	affineTrans.at<double>(1, 1) = transParameter[5];
	affineTrans.at<double>(1, 2) = transParameter[3];

	cv::Size newSize = getNewSize(iRows, iCols, affineTrans);

	cv::Mat newImg;
	warpAffine(img, newImg, affineTrans , newSize, cv::INTER_CUBIC);

	return newImg;
	//cvSaveImage(image_out, newImg);
}


template<class T>
static cv::Mat perspective_warp_Image(string image_in, T transParameter)
{
	cv::Mat img = cv::imread(image_in);
	int iRows = img.rows;
	int iCols = img.cols;

	cv::Mat affineTrans(2, 3, CV_64F);

	//for(int i = 0;i < 2;i++)
	//{
	//	for(int j = 0;j < 3;j++)
	//	{
	//		affineTrans.at<double>(i, (j-1)%3) = transParameter[i * 3 + j];
	//	}
	//}
	affineTrans.at<double>(0, 0) = transParameter[1];
	affineTrans.at<double>(0, 1) = transParameter[2];
	affineTrans.at<double>(0, 2) = transParameter[0];
	affineTrans.at<double>(1, 0) = transParameter[4];
	affineTrans.at<double>(1, 1) = transParameter[5];
	affineTrans.at<double>(1, 2) = transParameter[3];

	cv::Size newSize = getNewSize(iRows, iCols, affineTrans);

	cv::Mat newImg;
	warpPerspective(img, newImg, affineTrans , newSize, cv::INTER_CUBIC);

	return newImg;
	//cvSaveImage(image_out, newImg);
}