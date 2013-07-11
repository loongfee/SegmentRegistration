//// hough.cpp : Defines the entry point for the console application.
////
#include "LineExtract.h"
#include "lineConstant.h"
#include "func.h"
//#include "LineMatch.h"
//#include "PointExtrat.h"
//#include "PointMatch.h"
////
//#include "surf_match.h"
//#include "ecmprAffine.h"
#include "ecmlrAffine.h"
#include "ecmlrPerspective.h"
#include "ecmlr.h"
//#include "mcemprAffine.h"
#include "lineFeature.h"
#include "LineMatchFLANN.h"
#include "lbd/EDLineDetector.hh"
#include "lbd/LineDescriptor.hh"
#include "lbd/PairwiseLineMatching.hh"
#include "lp/LineMatcher.h"
#include "lp/lp.h"
#include <iostream>
#include <Eigen/Eigen>
#include <Windows.h>
#include <direct.h>
using namespace std;
#pragma comment(lib,"LineMatcher.lib")
vector<double> parameters(6);
//
////#define GLOG_NO_ABBREVIATED_SEVERITIES
////#include "ceres/ceres.h"
//
void warp(std::string sourceImage, std::string resultImage)
{
	cv::Mat img = cv::imread(sourceImage.c_str());
	int iRows = img.rows;
	int iCols = img.cols;

	cv::Mat affineTrans(2, 3, CV_64F);

	//affineTrans.at<double>(0, 2) = 100.0;
	//affineTrans.at<double>(0, 0) = 1.2;
	//affineTrans.at<double>(0, 1) = 0.25;
	//affineTrans.at<double>(1, 2) = 162.0;
	//affineTrans.at<double>(1, 0) = 0.2;
	//affineTrans.at<double>(1, 1) = 1.2;
	affineTrans.at<double>(0, 2) = -0.0;
	affineTrans.at<double>(0, 0) = 0.8;
	affineTrans.at<double>(0, 1) = 0.2;
	affineTrans.at<double>(1, 2) = -0.0;
	affineTrans.at<double>(1, 0) = 0.1;
	affineTrans.at<double>(1, 1) = 0.9;

	cv::Mat invert_affineTrans;
	invertAffineTransform(affineTrans, invert_affineTrans);

	cv::Size newSize = getNewSize(iRows, iCols, invert_affineTrans);

	cv::Mat newImg;
	warpAffine(img, newImg, invert_affineTrans , newSize, cv::INTER_CUBIC);


	imwrite(resultImage.c_str(), newImg);
}

cvline_polar affineTransLine(cv::Mat affineTrans, const cvline_polar& inLine)
{
	double a0 = affineTrans.at<double>(0, 2);
	double a1 = affineTrans.at<double>(0, 0);
	double a2 = affineTrans.at<double>(0, 1);
	double b0 = affineTrans.at<double>(1, 2);
	double b1 = affineTrans.at<double>(1, 0);
	double b2 = affineTrans.at<double>(1, 1);
	cvline_polar outLine;
	fPoint line[2];
	line[0].x = a0 + a1 * inLine.pt1.x + a2 * inLine.pt1.y;
	line[0].y = b0 + b1 * inLine.pt1.x + b2 * inLine.pt1.y;
	line[1].x = a0 + a1 * inLine.pt2.x + a2 * inLine.pt2.y;
	line[1].y = b0 + b1 * inLine.pt2.x + b2 * inLine.pt2.y;
	outLine = line2polar(line);
	return outLine;
}

void perspective_warp(std::string sourceImage, std::string resultImage)
{
	cv::Mat img = cv::imread(sourceImage.c_str());
	int iRows = img.rows;
	int iCols = img.cols;

	cv::Mat perspectiveTrans(3, 3, CV_64F);

	perspectiveTrans.at<double>(0, 0) = 1.3;
	perspectiveTrans.at<double>(0, 1) = 0.2;
	perspectiveTrans.at<double>(0, 2) = 10.0;
	perspectiveTrans.at<double>(1, 0) = 0.1;
	perspectiveTrans.at<double>(1, 1) = 0.8;
	perspectiveTrans.at<double>(1, 2) = 20.2;
	perspectiveTrans.at<double>(2, 0) = 0.0;
	perspectiveTrans.at<double>(2, 1) = 0.0;
	perspectiveTrans.at<double>(2, 2) = 0.7;

	//cv::Mat invert_affineTrans;
	//invertAffineTransform(affineTrans, invert_affineTrans);

	//cv::Size newSize = getNewSize(iRows, iCols, perspectiveTrans);

	cv::Mat invertPerspectiveTrans;
	cv::invert(perspectiveTrans, invertPerspectiveTrans);
	cv::Size newSize = getPerspectiveNewSize(iRows, iCols, invertPerspectiveTrans);

	cv::Mat newImg;
	warpPerspective(img, newImg, invertPerspectiveTrans , newSize, cv::INTER_CUBIC);
	//warpAffine(img, newImg, invert_affineTrans , newSize, cv::INTER_CUBIC);


	imwrite(resultImage.c_str(), newImg);
}

//void plotLineAngle(string outName, const vector<cvline_polar>& keylines)
//{
//	int nTotal = (int)keylines.size();
//	int hist[360] = {0};
//	for (int i = 0;i < nTotal;++i)
//	{
//		double dDegree = keylines[i].theta * 180 / PI;
//		int iDegree = (int)round(dDegree);
//		hist[iDegree]++;
//	}
//	int maxCount = 0;
//	for (int i = 0;i < 360;++i)
//	{
//		if (maxCount < hist[i])
//		{
//			maxCount = hist[i];
//		}
//	}
//}

void outLineAngle(string outName, const vector<cvline_polar>& reference_lines, const vector<cvline_polar>& source_lines)
{
	const int nHistogram = 180;
	int nTotal_ref = (int)reference_lines.size();
	int hist_ref[nHistogram] = {0};
	for (int i = 0;i < nTotal_ref;++i)
	{
		double dDegree = reference_lines[i].theta * 180 / PI;
		int iDegree = (int)round(dDegree)+90;
		iDegree = (iDegree < 0)?0:iDegree;
		iDegree = (iDegree > nHistogram-1)?(nHistogram-1):iDegree;
		hist_ref[iDegree]++;
	}

	int nTotal_src = (int)source_lines.size();
	int hist_src[nHistogram] = {0};
	for (int i = 0;i < nTotal_src;++i)
	{
		double dDegree = source_lines[i].theta * 180 / PI;
		int iDegree = (int)round(dDegree)+90;
		iDegree = (iDegree < 0)?0:iDegree;
		iDegree = (iDegree > (nHistogram-1))?(nHistogram-1):iDegree;
		hist_src[iDegree]++;
	}

	double accumulate_ref[nHistogram] = {0.0};
	double accumulate_src[nHistogram] = {0.0};
	accumulate_ref[0] = hist_ref[0];
	accumulate_src[0] = hist_src[0];
	for (int i = 1;i < nHistogram;++i)
	{
		accumulate_ref[i] = accumulate_ref[i-1]+hist_ref[i];
		accumulate_src[i] = accumulate_src[i-1]+hist_src[i];
	}
	for (int i = 1;i < nHistogram;++i)
	{
		accumulate_ref[i] = accumulate_ref[i] / (accumulate_ref[nHistogram-1] + DBL_EPSILON);
		accumulate_src[i] = accumulate_src[i] / (accumulate_src[nHistogram-1] + DBL_EPSILON);
	}

	int angleMap[nHistogram] = {0};
	for (int i = 0;i < nHistogram;++i)
	{
		double minDist = DBL_MAX;
		int _map = 0;
		for (int j = 0;j < nHistogram;++j)
		{
			double dist = fabs(accumulate_src[i] - accumulate_ref[j]);
			if (dist < minDist)
			{
				minDist = dist;
				_map = j;
			}
		}
		angleMap[i] = _map;
	}

	fstream fs;
	fs.open(outName.c_str(), std::ios_base::out);
	for (int i = 0;i < nHistogram;++i)
	{
		fs<<i-90<<"\t"<<hist_ref[i]<<"\t"<<hist_src[i]<<"\t"<<hist_src[angleMap[i]]<<endl;
	}
	fs.close();
	//int maxCount = 0;
	//for (int i = 0;i < 360;++i)
	//{
	//	if (maxCount < hist[i])
	//	{
	//		maxCount = hist[i];
	//	}
	//}
}

void drawOutliers(string foldName, string sourceName, string referenceName, string linesName)
{
	string SourceImage = foldName + "\\" + sourceName;//待配准影像
	string ReferImage = foldName + "\\" + referenceName;//参考影像
	string linesText = foldName + "\\" + linesName;//参考影像
	string resultImage = foldName + "\\matched.png";//结果影像
	string outlierImage = foldName + "\\outlier.png";//结果影像
	string outlierLog = foldName + "\\outlier.txt";
	fstream fs;
	fs.open(outlierLog.c_str(), std::ios_base::out);

	std::vector<cvline_polar> source_ls;
	std::vector<cvline_polar> reference_ls;

	//FILE* pf;
	ifstream ifs;//(linesName);
	ifs.open(linesText, std::ios_base::in);
	int id_pt1;
	int id_pt2;
	double pt1_source_x;
	double pt1_source_y;
	double pt2_source_x;
	double pt2_source_y;
	double pt1_reference_x;
	double pt1_reference_y;
	double pt2_reference_x;
	double pt2_reference_y;
	double temp;
	while (ifs>>id_pt1>>pt1_source_x>>pt1_source_y>>pt1_reference_x>>pt1_reference_y>>temp
		&& ifs>>id_pt2>>pt2_source_x>>pt2_source_y>>pt2_reference_x>>pt2_reference_y>>temp)
	{
		if(id_pt1 != id_pt2)
		{
			continue;
		}
		fPoint point_source[2];
		fPoint point_reference[2];
		point_source[0] = fPoint(pt1_source_x, pt1_source_y);
		point_reference[0] = fPoint(pt1_reference_x, pt1_reference_y);
		point_source[1] = fPoint(pt2_source_x, pt2_source_y);
		point_reference[1] = fPoint(pt2_reference_x, pt2_reference_y);
		cvline_polar line_source = line2polar(point_source);
		cvline_polar line_reference = line2polar(point_reference);
		source_ls.push_back(line_source);
		reference_ls.push_back(line_reference);
	}
	ifs.close();

	int nMatched = (int)source_ls.size();
	ecmlrAffine* ecmlr = new ecmlrAffine(reference_ls, source_ls);
	ecmlr->m_referenceImageFile = ReferImage;
	ecmlr->m_sourceImageFile = SourceImage;
	ecmlr->setOutFileStream(&fs);
	ecmlr->setMatched();
	ecmlr->parameters_optimization();
	ecmlr->updateTransLines();
	ecmlr->updateDeltaSquare();
	double delta_square = ecmlr->getDeltaSquare();
	Eigen::VectorXd modelParameters = ecmlr->getParameters();
	int nOutlier = 0;
	double outlier_threshold = 2.0;
	double outlier_threshold_squre = outlier_threshold*outlier_threshold;
	vector<cvline_polar> outlier_ref_lines_list;
	vector<cvline_polar> outlier_src_lines_list;
	for (int i = 0;i < nMatched;++i)
	{
		cvline_polar outPt = ecmlr->forward(reference_ls[i]);
		fPoint dist = ecmlr->segment2polarline(source_ls[i], outPt);
		double dist2_squre = dist.x*dist.x + dist.y*dist.y;

		fPoint center_src((source_ls[i].pt1.x+source_ls[i].pt2.x)*0.5, (source_ls[i].pt1.y+source_ls[i].pt2.y)*0.5);
		fPoint center_ref((outPt.pt1.x+outPt.pt2.x)*0.5, (outPt.pt1.y+outPt.pt2.y)*0.5);
		double cen_dis = (center_src.x-center_ref.x)*(center_src.x-center_ref.x) + (center_src.y-center_ref.y)*(center_src.y-center_ref.y);
		//if (cen_dis > 1000.0)
		//{
		//	outlier_ref_lines_list.push_back(matched_ref_lines_list[i]);
		//	outlier_src_lines_list.push_back(matched_src_lines_list[i]);
		//	nOutlier++;
		//	continue;
		//}
		if (dist2_squre > outlier_threshold_squre*delta_square)
		{	
			outlier_ref_lines_list.push_back(reference_ls[i]);
			outlier_src_lines_list.push_back(source_ls[i]);
			nOutlier++;
			continue;
		}
	}


	fs<<"all matched:"<<nMatched<<endl;
	fs<<"outliers:"<<nOutlier<<endl;
	fs<<endl;

	for (int i = 0;i < (int)modelParameters.size();++i)
	{
		cout<<modelParameters[i]<<endl;
	}
	cout<<"finished!"<<endl;

	fs.close();


	cv::Mat srcImage = cv::imread(SourceImage);
	cv::Mat refImage = cv::imread(ReferImage);
	cv::Mat img_matches;
	cv::Mat img_outlier;

	int thickness = 2;
	

	drawLineMatches(srcImage, source_ls, refImage, reference_ls,
		img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);

	drawLineMatches(srcImage, outlier_src_lines_list, refImage, outlier_ref_lines_list,
		img_outlier, cv::Scalar::all(-1), cv::Scalar::all(-1),
		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	
	cv::imwrite(resultImage, img_matches);
	cv::imwrite(outlierImage, img_outlier);
}

//
//void point_match()
//{
//	string foldName = "E:\\Dropbox\\Programs\\data\\temple";
//	//const char*  SourceImage = "D:\\workspace\\testdata\\features\\base.jpg";//待配准影像
//	string SourceImage = foldName + "\\2.png";//待配准影像
//	//const char*  SourceImage = "D:\\workspace\\testdata\\Landsat\\std\\p122r025_7t20000914.tif";//待配准影像
//	//const char*  ReferImage = "D:\\workspace\\testdata\\features\\warp.jpg";//参考影像
//	string ReferImage = foldName + "\\1.png";//参考影像
//	string resultImage = foldName + "\\registrated_surf.png";//参考影像
//
//	cv::initModule_nonfree();
//	IplImage* source = cvLoadImage( SourceImage.c_str(), CV_LOAD_IMAGE_GRAYSCALE );
//	IplImage* refer = cvLoadImage( ReferImage.c_str(), CV_LOAD_IMAGE_GRAYSCALE );
//	if( !source || !refer)
//	{
//		return ;
//	}
//
//	CvMemStorage* storage = cvCreateMemStorage(0);
//	CvSURFParams params = cvSURFParams(500, 1);
//
//	CvSeq* sourceKeypoints = 0, *sourceDescriptors = 0;
//	CvSeq* referKeypoints = 0, *referDescriptors = 0;
//
//	cvExtractSURF( source, 0, &sourceKeypoints, &sourceDescriptors, storage, params );
//	cvExtractSURF( refer, 0, &referKeypoints, &referDescriptors, storage, params );
//
//	//vector<fPoint> src_points;
//	//vector<fPoint> ref_points;
//
//	//for( int i = 0; i < referKeypoints->total; i++ )
//	//{
//	//	CvPoint* point = (CvPoint*)cvGetSeqElem(referKeypoints, i);
//	//	//cvLine( hough_result, line[0], line[1], CV_RGB(255,0,0), 1, CV_AA, 0);
//	//	ref_points.push_back(fPoint(point->x, point->y));
//	//}
//
//	//for( int i = 0; i < sourceKeypoints->total; i++ )
//	//{
//	//	CvPoint* point = (CvPoint*)cvGetSeqElem(sourceKeypoints, i);
//	//	//cvLine( hough_result, line[0], line[1], CV_RGB(255,0,0), 1, CV_AA, 0);
//	//	src_points.push_back(fPoint(point->x, point->y));
//	//}
//
//	IplImage* source_color = cvCreateImage(cvGetSize(source), 8, 3);
//	cvCvtColor( source, source_color, CV_GRAY2BGR );
//
//	static CvScalar colors[] =
//	{
//		{{0,0,255}},
//		{{0,128,255}},
//		{{0,255,255}},
//		{{0,255,0}},
//		{{255,128,0}},
//		{{255,255,0}},
//		{{255,0,0}},
//		{{255,0,255}},
//		{{255,255,255}}
//	};
//
//#ifdef USE_FLANN
//	printf("Using approximate nearest neighbor search\n");
//#endif
//
//	vector<int> ptpairs;
//#ifdef USE_FLANN
//	flannFindPairs( sourceKeypoints, sourceDescriptors, referKeypoints, referDescriptors, ptpairs );
//#else
//	findPairs( sourceKeypoints, sourceDescriptors, referKeypoints, referDescriptors, ptpairs );
//#endif
//	
//// RANSAC detect outliers
//	auto_ptr< estimators::Solver<mat,vec> > ptrSolver(
//		new estimators::affineSolver<mat,vec>);
//
//	//-- Create input data
//	mat dataPoints(ptpairs.size()/2, 4);// = "0 0; 1 1; 2 2; 3 3";
//	fstream model_file;
//	model_file.open("model.txt", ios_base::out);
//	fstream scene_file;
//	scene_file.open("scene.txt", ios_base::out);
//
//	for(int i = 0;i < (int)ptpairs.size()/2;++i)
//	{
//		CvSURFPoint* r1 = (CvSURFPoint*)cvGetSeqElem( sourceKeypoints, ptpairs[2*i] );
//		CvSURFPoint* r2 = (CvSURFPoint*)cvGetSeqElem( referKeypoints, ptpairs[2*i+1] );
//		dataPoints(i, 0) = r1->pt.x;
//		dataPoints(i, 1) = r1->pt.y;
//		dataPoints(i, 2) = r2->pt.x;
//		dataPoints(i, 3) = r2->pt.y;
//
//		model_file<<r1->pt.x<<"  "<<r1->pt.y<<endl;
//		scene_file<<r2->pt.x<<"  "<<r2->pt.y<<endl;
//	}
//	model_file.close();
//	scene_file.close();
//
//	vector<int> inliers;
//	vector<vec> models;
//
//	ransac::Ransac_Handler ransac_fun_Handler;
//	bool result = ransac::Ransac_RobustEstimator
//		(
//		dataPoints, // the input data
//		estimators::affineSolver<mat,vec>::extractor, // How select sampled point from indices
//		dataPoints.n_rows,  // the number of putatives data
//		*(ptrSolver.get()),  // compute the underlying model given a sample set
//		estimators::affineSolver<mat,vec>::defaultEvaluator,  // the function to evaluate a given model
//		//Ransac Object that contain function:
//		// CandidatesSelector, Sampler and TerminationFunction
//		ransac_fun_Handler, // the basic ransac object
//		1000,  // the maximum rounds for RANSAC routine
//		inliers, // inliers to the final solution
//		models, // models array that fit input data
//		0.95 // the confidence want to achieve at the end
//		);
//
//	for (int i = 0;i < (int)models[0].n_rows;++i)
//	{
//		cout<<models[0](i)<<endl;
//	}
//	cout<<endl;
//
//// calculate the affine model
//	vector<PointPair> pointPairList;
//	for(int i = 0;i < (int)inliers.size();++i)
//	{
//		PointPair pointPair;
//		pointPair.point = fPoint(dataPoints(inliers[i], 0), dataPoints(inliers[i], 1));
//		pointPair.point_prime = fPoint(dataPoints(inliers[i], 2), dataPoints(inliers[i], 3));
//		pointPairList.push_back(pointPair);
//	}
//
//	PointMatch pointMatch;
//	pointMatch.findBestModel(pointPairList);
//	for (int i = 0;i < (int)pointMatch.modelParameters.size();++i)
//	{
//		cout<<pointMatch.modelParameters[i]<<endl;
//	}
//
//	affine_warp_Image< vector<double> >(SourceImage.c_str(), resultImage.c_str(), pointMatch.modelParameters);
//}
//
//void surf_match()
//{
//	string foldName = "E:\\Dropbox\\Programs\\data\\ding";
//	//string foldName = "D:\\workspace";
//	//const char*  SourceImage = "D:\\workspace\\testdata\\features\\base.jpg";//待配准影像
//	string SourceImage = foldName + "\\1.png";//待配准影像
//	//const char*  SourceImage = "D:\\workspace\\testdata\\Landsat\\std\\p122r025_7t20000914.tif";//待配准影像
//	//const char*  ReferImage = "D:\\workspace\\testdata\\features\\warp.jpg";//参考影像
//	string ReferImage = foldName + "\\2.png";//参考影像
//	string resultImage = foldName + "\\registrated_surf.png";//参考影像
//
//	cv::initModule_nonfree();
//	vector<PointPair> pointPairList;
//	if(!surf_match(SourceImage.c_str(), ReferImage.c_str(), pointPairList))
//	{
//		cout<<"surf extract failed for \""<<SourceImage<<"\"!"<<endl;
//		system("pause");
//		return;
//	}
//
//	// RANSAC detect outliers
//	auto_ptr< estimators::Solver<mat,vec> > ptrSolver(
//		new estimators::affineSolver<mat,vec>);
//
//	//-- Create input data
//	mat dataPoints((int)pointPairList.size(), 4);// = "0 0; 1 1; 2 2; 3 3";
//	for(int i = 0;i < (int)pointPairList.size();++i)
//	{
//		dataPoints(i, 0) = pointPairList[i].point.x;
//		dataPoints(i, 1) = pointPairList[i].point.y;
//		dataPoints(i, 2) = pointPairList[i].point_prime.x;
//		dataPoints(i, 3) = pointPairList[i].point_prime.y;
//	}
//
//	vector<int> inliers;
//	vector<vec> models;
//
//	ransac::Ransac_Handler ransac_fun_Handler;
//	bool result = ransac::Ransac_RobustEstimator
//		(
//		dataPoints, // the input data
//		estimators::affineSolver<mat,vec>::extractor, // How select sampled point from indices
//		dataPoints.n_rows,  // the number of putatives data
//		*(ptrSolver.get()),  // compute the underlying model given a sample set
//		estimators::affineSolver<mat,vec>::defaultEvaluator,  // the function to evaluate a given model
//		//Ransac Object that contain function:
//		// CandidatesSelector, Sampler and TerminationFunction
//		ransac_fun_Handler, // the basic ransac object
//		1000,  // the maximum rounds for RANSAC routine
//		inliers, // inliers to the final solution
//		models, // models array that fit input data
//		0.95 // the confidence want to achieve at the end
//		);
//
//	for (int i = 0;i < (int)models[0].n_rows;++i)
//	{
//		cout<<models[0](i)<<endl;
//	}
//	cout<<endl;
//
//	// calculate the affine model
//	pointPairList.clear();
//	for(int i = 0;i < (int)inliers.size();++i)
//	{
//		PointPair pointPair;
//		pointPair.point = fPoint(dataPoints(inliers[i], 0), dataPoints(inliers[i], 1));
//		pointPair.point_prime = fPoint(dataPoints(inliers[i], 2), dataPoints(inliers[i], 3));
//		pointPairList.push_back(pointPair);
//	}
//
//	PointMatch pointMatch;
//	pointMatch.findBestModel(pointPairList);
//	for (int i = 0;i < (int)pointMatch.modelParameters.size();++i)
//	{
//		cout<<pointMatch.modelParameters[i]<<endl;
//	}
//
//	affine_warp_Image< vector<double> >(SourceImage.c_str(), resultImage.c_str(), pointMatch.modelParameters);
//}
//
//void orb_match()
//{
//
//	string foldName = "E:\\Dropbox\\Programs\\data\\beijing";
//	//const char*  SourceImage = "D:\\workspace\\testdata\\features\\base.jpg";//待配准影像
//	string SourceImage = foldName + "\\2.png";//待配准影像
//	//const char*  SourceImage = "D:\\workspace\\testdata\\Landsat\\std\\p122r025_7t20000914.tif";//待配准影像
//	//const char*  ReferImage = "D:\\workspace\\testdata\\features\\warp.jpg";//参考影像
//	string ReferImage = foldName + "\\1.png";//参考影像
//	string resultImage = foldName + "\\registrated_orb.png";//参考影像
//
//	vector<PointPair> pointPairList;
//	orb_match(SourceImage.c_str(), ReferImage.c_str(), pointPairList);
//
//	// RANSAC detect outliers
//	auto_ptr< estimators::Solver<mat,vec> > ptrSolver(
//		new estimators::affineSolver<mat,vec>);
//
//	//-- Create input data
//	mat dataPoints((int)pointPairList.size(), 4);// = "0 0; 1 1; 2 2; 3 3";
//	for(int i = 0;i < (int)pointPairList.size();++i)
//	{
//		dataPoints(i, 0) = pointPairList[i].point.x;
//		dataPoints(i, 1) = pointPairList[i].point.y;
//		dataPoints(i, 2) = pointPairList[i].point_prime.x;
//		dataPoints(i, 3) = pointPairList[i].point_prime.y;
//	}
//
//	vector<int> inliers;
//	vector<vec> models;
//
//	ransac::Ransac_Handler ransac_fun_Handler;
//	bool result = ransac::Ransac_RobustEstimator
//		(
//		dataPoints, // the input data
//		estimators::affineSolver<mat,vec>::extractor, // How select sampled point from indices
//		dataPoints.n_rows,  // the number of putatives data
//		*(ptrSolver.get()),  // compute the underlying model given a sample set
//		estimators::affineSolver<mat,vec>::defaultEvaluator,  // the function to evaluate a given model
//		//Ransac Object that contain function:
//		// CandidatesSelector, Sampler and TerminationFunction
//		ransac_fun_Handler, // the basic ransac object
//		1000,  // the maximum rounds for RANSAC routine
//		inliers, // inliers to the final solution
//		models, // models array that fit input data
//		0.95 // the confidence want to achieve at the end
//		);
//
//	for (int i = 0;i < (int)models[0].n_rows;++i)
//	{
//		cout<<models[0](i)<<endl;
//	}
//	cout<<endl;
//
//	// calculate the affine model
//	pointPairList.clear();
//	for(int i = 0;i < (int)inliers.size();++i)
//	{
//		PointPair pointPair;
//		pointPair.point = fPoint(dataPoints(inliers[i], 0), dataPoints(inliers[i], 1));
//		pointPair.point_prime = fPoint(dataPoints(inliers[i], 2), dataPoints(inliers[i], 3));
//		pointPairList.push_back(pointPair);
//	}
//
//	PointMatch pointMatch;
//	pointMatch.findBestModel(pointPairList);
//	for (int i = 0;i < (int)pointMatch.modelParameters.size();++i)
//	{
//		cout<<pointMatch.modelParameters[i]<<endl;
//	}
//
//	affine_warp_Image< vector<double> >(SourceImage.c_str(), resultImage.c_str(), pointMatch.modelParameters);
//}
//void point_match1()
//{
//	const char*  SourceImage = "D:\\workspace\\testdata\\features\\base.jpg";//待配准影像
//	//const char*  SourceImage = "D:\\workspace\\testdata\\Landsat\\std\\p122r025_7t20000914.tif";//待配准影像
//	const char*  ReferImage = "D:\\workspace\\testdata\\features\\warp.jpg";//参考影像
//	const char*  resultImage = "D:\\workspace\\testdata\\features\\result.jpg";//参考影像
//
//	cv::initModule_nonfree();
//	IplImage* source = cvLoadImage( SourceImage, CV_LOAD_IMAGE_GRAYSCALE );
//	IplImage* refer = cvLoadImage( ReferImage, CV_LOAD_IMAGE_GRAYSCALE );
//	if( !source || !refer)
//	{
//		return ;
//	}
//
//	CvMemStorage* storage = cvCreateMemStorage(0);
//	CvSURFParams params = cvSURFParams(500, 1);
//
//	CvSeq* sourceKeypoints = 0, *sourceDescriptors = 0;
//	CvSeq* referKeypoints = 0, *referDescriptors = 0;
//
//	cvExtractSURF( source, 0, &sourceKeypoints, &sourceDescriptors, storage, params );
//	cvExtractSURF( refer, 0, &referKeypoints, &referDescriptors, storage, params );
//
//	vector<fPoint> src_points;
//	vector<fPoint> ref_points;
//
//	for( int i = 0; i < referKeypoints->total; i++ )
//	{
//		CvSURFPoint* point = (CvSURFPoint*)cvGetSeqElem( referKeypoints, i );
//		//cvLine( hough_result, line[0], line[1], CV_RGB(255,0,0), 1, CV_AA, 0);
//		ref_points.push_back(fPoint(point->pt.x, point->pt.y));
//	}
//
//	for( int i = 0; i < sourceKeypoints->total; i++ )
//	{
//		CvSURFPoint* point = (CvSURFPoint*)cvGetSeqElem( sourceKeypoints, i );
//		//cvLine( hough_result, line[0], line[1], CV_RGB(255,0,0), 1, CV_AA, 0);
//		src_points.push_back(fPoint(point->pt.x, point->pt.y));
//	}
//
//	PointMatch pointMatch;
//	pointMatch.registration(src_points, ref_points);
//	for (int i = 0;i < (int)pointMatch.modelParameters.size();++i)
//	{
//		cout<<pointMatch.modelParameters[i]<<endl;
//	}
//	cout<<"finished!"<<endl;
//}
//
//void point_match2()
//{
//	string foldName = "E:\\Dropbox\\Programs\\data\\beijing";
//	//const char*  SourceImage = "D:\\workspace\\testdata\\features\\base.jpg";//待配准影像
//	string SourceImage = foldName + "\\1.png";//待配准影像
//	//const char*  SourceImage = "D:\\workspace\\testdata\\Landsat\\std\\p122r025_7t20000914.tif";//待配准影像
//	//const char*  ReferImage = "D:\\workspace\\testdata\\features\\warp.jpg";//参考影像
//	string ReferImage = foldName + "\\2.png";//参考影像
//	string resultImage = foldName + "\\registrated_ecm.png";//参考影像
//
//	CvSURFParams params = cvSURFParams(10500, 1);
//
//	vector<fPoint> src_points;
//	vector<fPoint> ref_points;
//
//	//surf
//	if(!surf_extract(SourceImage, src_points, params))
//	{
//		cout<<"surf extract failed for \""<<SourceImage<<"\"!"<<endl;
//		system("pause");
//		return;
//	}
//	if(!surf_extract(ReferImage, ref_points, params))
//	{
//		cout<<"surf extract failed for \""<<ReferImage<<"\"!"<<endl;
//		system("pause");
//		return;
//	}
//
//	//// orb
//	//if(!orb_extract(SourceImage.c_str(), src_points))
//	//{
//	//	cout<<"surf extract failed for \""<<SourceImage<<"\"!"<<endl;
//	//	system("pause");
//	//	return;
//	//}
//	//if(!orb_extract(ReferImage.c_str(), ref_points))
//	//{
//	//	cout<<"surf extract failed for \""<<ReferImage<<"\"!"<<endl;
//	//	system("pause");
//	//	return;
//	//}
//
//	//// Harris
//	//if(!harris_extract(SourceImage, src_points, 40, true))
//	//{
//	//	cout<<"harris extract failed for \""<<SourceImage<<"\"!"<<endl;
//	//	system("pause");
//	//}
//	//if(!harris_extract(ReferImage, ref_points, 40, true))
//	//{
//	//	cout<<"harris extract failed for \""<<ReferImage<<"\"!"<<endl;
//	//	system("pause");
//	//}
//
//	Eigen::MatrixXd modelData((int)ref_points.size(), 2);
//	for( int i = 0; i < (int)ref_points.size(); i++ )
//	{
//		modelData(i, 0) = ref_points[i].x;
//		modelData(i, 1) = ref_points[i].y;
//	}
//
//	Eigen::MatrixXd observationData((int)src_points.size(), 2);
//	for( int i = 0; i < (int)src_points.size(); i++ )
//	{
//		observationData(i, 0) = src_points[i].x;
//		observationData(i, 1) = src_points[i].y;
//	}
//
//	//fstream ofs;
//	//ofs.open("model.txt", ios_base::out);
//	//for (int i = 0;i < (int)ref_points.size();++i)
//	//{
//	//	ofs<<ref_points[i].x<<"\t"<<ref_points[i].y<<endl;
//	//}
//	//ofs.close();
//	//ofs.open("scene.txt", ios_base::out);
//	//for (int i = 0;i < (int)src_points.size();++i)
//	//{
//	//	ofs<<src_points[i].x<<"\t"<<src_points[i].y<<endl;
//	//}
//	//ofs.close();
//
//	ecmprAffine* mcem = new ecmprAffine(modelData, observationData);
//	mcem->m_sourceImageFile = ReferImage;
//	vector<double> modelParameters;
//	mcem->solve(modelParameters);
//
//	//mcemprAffine* mcem = new mcemprAffine(modelData, observationData);
//	//vector<double> modelParameters;
//	//mcem->solve(modelParameters);
//
//	// calculate the affine model
//	vector<PointPair> pointPairList;
//	for(int i = 0;i < (int)mcem->m_z.size();++i)
//	{
//		if(mcem->m_z[i] == ref_points.size()) continue;
//		PointPair pointPair;
//		pointPair.point = ref_points[mcem->m_z[i]];
//		pointPair.point_prime = src_points[i];
//		pointPairList.push_back(pointPair);
//	}
//
//	PointMatch pointMatch;
//	pointMatch.findBestModel(pointPairList);
//	for (int i = 0;i < (int)pointMatch.modelParameters.size();++i)
//	{
//		cout<<pointMatch.modelParameters[i]<<endl;
//	}
//
//	affine_warp_Image< vector<double> >(ReferImage.c_str(), resultImage.c_str(), pointMatch.modelParameters);
//	//affine_warp_Image< vector<double> >(SourceImage, resultImage, pointMatch.modelParameters);
//
//	//icp matcher(ref_points, src_points);
//	//matcher.solve();
//}
//

//void lbd_line_match(string foldName, string sourceName, string referenceName)
//{
//	string SourceImage = foldName + "\\" + sourceName;//待配准影像
//	string ReferImage = foldName + "\\" + referenceName;//参考影像
//	string registratedImage = foldName + "\\lbd_registrated.png";//结果影像
//	string referImage_line = foldName + "\\lined1.png";//参考影像
//	string sourceImage_line = foldName + "\\lined2.png";//待配准影像
//	string matchedImage = foldName + "\\lbd_segment_matched.png";//结果影像
//	string outlierImage = foldName + "\\lbd_segment_outlier.png";//结果影像
//	string strLog = foldName + "\\lbd_log.txt";
//	fstream fs;
//	fs.open(strLog.c_str(), std::ios_base::out);
//
//	//用clock()来计时  毫秒
//	clock_t  clockBegin, clockEnd;
//	clockBegin = clock();
//
//	//extract lines, compute their descriptors and match lines
//	LineDescriptor lineDesc;
//	PairwiseLineMatching lineMatch;
//
//	ScaleLines   linesInReference;
//	ScaleLines   linesInSource;
//	std::vector<unsigned int> matchResult;
//
//	cv::Mat refImage = cv::imread(ReferImage.c_str());
//	cv::Mat srcImage = cv::imread(SourceImage.c_str());
//	cvtColor(refImage, refImage, CV_RGB2GRAY );
//	cvtColor(srcImage, srcImage, CV_RGB2GRAY );
//
//	//cv::Mat debugImage;
//	//cv::Sobel(srcImage, debugImage, CV_8U, 1, 0, 3, 1.0, 0.0, cv::BORDER_CONSTANT);
//	//cv::imwrite("E:\\debug.png", debugImage);
//
//	lineDesc.GetLineDescriptor(refImage, linesInReference);
//	lineDesc.GetLineDescriptor(srcImage, linesInSource);
//
//	lineMatch.LineMatching(linesInReference,linesInSource,matchResult);
//	
//	vector<cvline_polar> matched_ref_lines_list;
//	vector<cvline_polar> matched_src_lines_list;
//
//
//	int nMatched = (int)matchResult.size()/2;
//	for(int i = 0; i < nMatched; i++)
//	{
//		cvline_polar refLine;
//		cvline_polar srcLine;
//
//		fPoint line[2];
//		//line[0].x = linesInReference[matchRefIndex[i]][0].startPointX;
//		//line[0].y = linesInReference[matchRefIndex[i]][0].startPointY;
//		//line[1].x = linesInReference[matchRefIndex[i]][0].endPointX;
//		//line[1].y = linesInReference[matchRefIndex[i]][0].endPointY;
//		line[0].x = linesInReference[matchResult[2*i]][0].sPointInOctaveX;
//		line[0].y = linesInReference[matchResult[2*i]][0].sPointInOctaveY;
//		line[1].x = linesInReference[matchResult[2*i]][0].ePointInOctaveX;
//		line[1].y = linesInReference[matchResult[2*i]][0].ePointInOctaveY;
//		refLine = line2polar(line);
//
//		//line[0].x = linesInSource[matchSrcIndex[i]][0].startPointX;
//		//line[0].y = linesInSource[matchSrcIndex[i]][0].startPointY;
//		//line[1].x = linesInSource[matchSrcIndex[i]][0].endPointX;
//		//line[1].y = linesInSource[matchSrcIndex[i]][0].endPointY;
//		line[0].x = linesInSource[matchResult[2*i+1]][0].sPointInOctaveX;
//		line[0].y = linesInSource[matchResult[2*i+1]][0].sPointInOctaveY;
//		line[1].x = linesInSource[matchResult[2*i+1]][0].ePointInOctaveX;
//		line[1].y = linesInSource[matchResult[2*i+1]][0].ePointInOctaveY;
//		srcLine = line2polar(line);
//
//		matched_ref_lines_list.push_back(refLine);
//		matched_src_lines_list.push_back(srcLine);
//	}
//
//	//std::vector<short> matchRefIndex;
//	//std::vector<short> matchSrcIndex;
//	//lineDesc.MatchLineByDescriptor(linesInReference, 	linesInSource,
//	//	matchRefIndex, matchSrcIndex,
//	//	LineDescriptor::NearestNeighbor);
//
//	//int nMatched = (int)matchRefIndex.size();
//	//for(int i = 0; i < nMatched; i++)
//	//{
//	//	cvline_polar refLine;
//	//	cvline_polar srcLine;
//
//	//	fPoint line[2];
//	//	//line[0].x = linesInReference[matchRefIndex[i]][0].startPointX;
//	//	//line[0].y = linesInReference[matchRefIndex[i]][0].startPointY;
//	//	//line[1].x = linesInReference[matchRefIndex[i]][0].endPointX;
//	//	//line[1].y = linesInReference[matchRefIndex[i]][0].endPointY;
//	//	line[0].x = linesInReference[matchRefIndex[i]][0].sPointInOctaveX;
//	//	line[0].y = linesInReference[matchRefIndex[i]][0].sPointInOctaveY;
//	//	line[1].x = linesInReference[matchRefIndex[i]][0].ePointInOctaveX;
//	//	line[1].y = linesInReference[matchRefIndex[i]][0].ePointInOctaveY;
//	//	refLine = line2polar(line);
//
//	//	//line[0].x = linesInSource[matchSrcIndex[i]][0].startPointX;
//	//	//line[0].y = linesInSource[matchSrcIndex[i]][0].startPointY;
//	//	//line[1].x = linesInSource[matchSrcIndex[i]][0].endPointX;
//	//	//line[1].y = linesInSource[matchSrcIndex[i]][0].endPointY;
//	//	line[0].x = linesInSource[matchSrcIndex[i]][0].sPointInOctaveX;
//	//	line[0].y = linesInSource[matchSrcIndex[i]][0].sPointInOctaveY;
//	//	line[1].x = linesInSource[matchSrcIndex[i]][0].ePointInOctaveX;
//	//	line[1].y = linesInSource[matchSrcIndex[i]][0].ePointInOctaveY;
//	//	srcLine = line2polar(line);
//
//	//	matched_ref_lines_list.push_back(refLine);
//	//	matched_src_lines_list.push_back(srcLine);
//	//}
//
//	//LineMatch(linesInReference, linesInSource, matched_ref_lines_list, matched_src_lines_list);
//
//	fs<<"source image:"<<SourceImage<<endl;
//
//	fs<<"reference image:"<<ReferImage<<endl;
//
//	ecmlrAffine* ecmlr = new ecmlrAffine(matched_ref_lines_list, matched_src_lines_list);
//	ecmlr->m_referenceImageFile = ReferImage;
//	ecmlr->m_sourceImageFile = SourceImage;
//	ecmlr->setOutFileStream(&fs);
//	ecmlr->setMatched();
//	ecmlr->parameters_optimization();
//	ecmlr->updateTransLines();
//	ecmlr->updateDeltaSquare();
//	double delta_square = ecmlr->getDeltaSquare();
//	Eigen::VectorXd modelParameters = ecmlr->getParameters();
//	int nOutlier = 0;
//	double outlier_threshold = 2.0;
//	double outlier_threshold_squre = outlier_threshold*outlier_threshold;
//	vector<cvline_polar> outlier_ref_lines_list;
//	vector<cvline_polar> outlier_src_lines_list;
//	for (int i = 0;i < nMatched;++i)
//	{
//		cvline_polar outPt = ecmlr->forward(matched_ref_lines_list[i]);
//		fPoint dist = ecmlr->segment2polarline(matched_src_lines_list[i], outPt);
//		double dist2_squre = dist.x*dist.x + dist.y*dist.y;
//
//		fPoint center_src((matched_src_lines_list[i].pt1.x+matched_src_lines_list[i].pt2.x)*0.5, (matched_src_lines_list[i].pt1.y+matched_src_lines_list[i].pt2.y)*0.5);
//		fPoint center_ref((outPt.pt1.x+outPt.pt2.x)*0.5, (outPt.pt1.y+outPt.pt2.y)*0.5);
//		double cen_dis = (center_src.x-center_ref.x)*(center_src.x-center_ref.x) + (center_src.y-center_ref.y)*(center_src.y-center_ref.y);
//		//if (cen_dis > 1000.0)
//		//{
//		//	outlier_ref_lines_list.push_back(matched_ref_lines_list[i]);
//		//	outlier_src_lines_list.push_back(matched_src_lines_list[i]);
//		//	nOutlier++;
//		//	continue;
//		//}
//		if (dist2_squre > outlier_threshold_squre*delta_square)
//		{	
//			outlier_ref_lines_list.push_back(matched_ref_lines_list[i]);
//			outlier_src_lines_list.push_back(matched_src_lines_list[i]);
//			nOutlier++;
//			continue;
//		}
//	}
//
//
//	fs<<"all matched:"<<nMatched<<endl;
//	fs<<"outliers:"<<nOutlier<<endl;
//	fs<<endl;
//
//	clockEnd = clock();
//	fs.setf(ios::fixed, ios::floatfield);
//	fs.precision(3);
//	fs<<"time consuming:"<<(clockEnd - clockBegin)*1e-3<<"second"<<endl;
//
//	//invert_affine_warp_Image< Eigen::VectorXd >(SourceImage.c_str(), registratedImage.c_str(), modelParameters);
//
//
//	for (int i = 0;i < (int)modelParameters.size();++i)
//	{
//		cout<<modelParameters[i]<<endl;
//	}
//	cout<<"finished!"<<endl;
//
//	fs.close();
//
//
//	srcImage = cv::imread(SourceImage);
//	refImage = cv::imread(ReferImage);
//	cv::Mat img_matches;
//	cv::Mat img_outlier;
//
//	int thickness = 2;
//	drawLineMatches(srcImage, matched_src_lines_list, refImage, matched_ref_lines_list,
//		img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
//		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
//
//	drawLineMatches(srcImage, outlier_src_lines_list, refImage, outlier_ref_lines_list,
//		img_outlier, cv::Scalar::all(-1), cv::Scalar::all(-1),
//		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
//
//	cv::imwrite(matchedImage, img_matches);
//	cv::imwrite(outlierImage, img_outlier);
//}

void lp_line_match(string foldName, string sourceName, string referenceName)
{
	cv::initModule_nonfree();
	string SourceImage = foldName + "\\" + sourceName;//待配准影像
	string ReferImage = foldName + "\\" + referenceName;//参考影像
	string registratedImage = foldName + "\\lp_registrated.png";//结果影像
	string matchedImage = foldName + "\\lp_segment_matched.png";//结果影像
	string outlierImage = foldName + "\\lp_segment_outlier.png";//结果影像
	string strLog = foldName + "\\lp_log.txt";
	fstream fs;
	fs.open(strLog.c_str(), std::ios_base::out);
	double sourceimage_p1 = 0.7;
	double sourceimage_p2 = 5000;
	double referimage_p1 = 0.7;
	double referimage_p2 = 5000;

	//用clock()来计时  毫秒
	clock_t  clockBegin, clockEnd;
	clockBegin = clock();

	int nLine1, nLine2;
	std::vector<Line> LineReferenceList;
	std::vector<Line> LineSourceList;
	lp_LineDetect(ReferImage.c_str(), LineReferenceList, referimage_p1, referimage_p2);
	lp_LineDetect(SourceImage.c_str(), LineSourceList, sourceimage_p1, sourceimage_p2);
	nLine1 = (int)LineReferenceList.size();
	nLine2 = (int)LineSourceList.size();
	cout << "extracted lines : " << nLine1 << " , " << nLine2 << endl;

	vector<MatchPt> PtMatches;
	int nPtMatch = PtMatch(ReferImage.c_str(),SourceImage.c_str(),PtMatches);
	cout << nPtMatch << " matched points" << endl;

	vector<MatchLine> outLM;
	//int n;
	int nMatched = MatchingLines(&LineReferenceList[0],&LineSourceList[0],nLine1,nLine2,0,PtMatches,outLM,false);

	vector<cvline_polar> matched_ref_lines_list;
	vector<cvline_polar> matched_src_lines_list;

	for(int i = 0; i < nMatched; i++)
	{
		int id1 = outLM[i].ID1;
		int id2 = outLM[i].ID2;
		cvline_polar refLine;
		cvline_polar srcLine;

		fPoint line[2];
		line[0].x = LineReferenceList[id1].StartPt.x;
		line[0].y = LineReferenceList[id1].StartPt.y;
		line[1].x = LineReferenceList[id1].EndPt.x;
		line[1].y = LineReferenceList[id1].EndPt.y;
		refLine = line2polar(line);

		line[0].x = LineSourceList[id2].StartPt.x;
		line[0].y = LineSourceList[id2].StartPt.y;
		line[1].x = LineSourceList[id2].EndPt.x;
		line[1].y = LineSourceList[id2].EndPt.y;
		srcLine = line2polar(line);

		matched_ref_lines_list.push_back(refLine);
		matched_src_lines_list.push_back(srcLine);
	}


	fs<<"source image:"<<SourceImage<<endl;

	fs<<"reference image:"<<ReferImage<<endl;

	//std::vector<int> inliers;
	//Eigen::VectorXd tmpParameters = affine_ransac(matched_ref_lines_list, matched_src_lines_list, inliers);
	//
	//vector<cvline_polar> inlier_ref_lines_list;
	//vector<cvline_polar> inlier_src_lines_list;
	//for (int i = 0; i < (int)inliers.size(); i++)
	//{
	//	inlier_ref_lines_list.push_back(matched_ref_lines_list[inliers[i]]);
	//	inlier_src_lines_list.push_back(matched_src_lines_list[inliers[i]]);
	//}
	ecmlrAffine* ecmlr = new ecmlrAffine(matched_ref_lines_list, matched_src_lines_list);
	ecmlr->m_referenceImageFile = ReferImage;
	ecmlr->m_sourceImageFile = SourceImage;
	ecmlr->setOutFileStream(&fs);
	ecmlr->setMatched();
	ecmlr->parameters_optimization();
	ecmlr->updateTransLines();
	ecmlr->updateDeltaSquare();
	//vector<double> tempParameters;
	//ecmlr->solve(tempParameters);
	double delta_square = ecmlr->getDeltaSquare();
	Eigen::VectorXd modelParameters = ecmlr->getParameters();
	int nOutlier = 0;
	double outlier_threshold = 2.0;
	double outlier_threshold_squre = outlier_threshold*outlier_threshold;
	vector<cvline_polar> outlier_ref_lines_list;
	vector<cvline_polar> outlier_src_lines_list;
	for (int i = 0;i < nMatched;++i)
	{
		cvline_polar outPt = ecmlr->forward(matched_ref_lines_list[i]);
		fPoint dist = ecmlr->segment2polarline(matched_src_lines_list[i], outPt);
		double dist2_squre = dist.x*dist.x + dist.y*dist.y;

		fPoint center_src((matched_src_lines_list[i].pt1.x+matched_src_lines_list[i].pt2.x)*0.5, (matched_src_lines_list[i].pt1.y+matched_src_lines_list[i].pt2.y)*0.5);
		fPoint center_ref((outPt.pt1.x+outPt.pt2.x)*0.5, (outPt.pt1.y+outPt.pt2.y)*0.5);
		double cen_dis = (center_src.x-center_ref.x)*(center_src.x-center_ref.x) + (center_src.y-center_ref.y)*(center_src.y-center_ref.y);
		//if (cen_dis > 1000.0)
		//{
		//	outlier_ref_lines_list.push_back(matched_ref_lines_list[i]);
		//	outlier_src_lines_list.push_back(matched_src_lines_list[i]);
		//	nOutlier++;
		//	continue;
		//}
		if (dist2_squre > outlier_threshold_squre*delta_square)
		{	
			outlier_ref_lines_list.push_back(matched_ref_lines_list[i]);
			outlier_src_lines_list.push_back(matched_src_lines_list[i]);
			nOutlier++;
			continue;
		}
	}


	fs<<"all matched:"<<nMatched<<endl;
	fs<<"outliers:"<<nOutlier<<endl;
	fs<<endl;

	clockEnd = clock();
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(3);
	fs<<"time consuming:"<<(clockEnd - clockBegin)*1e-3<<"second"<<endl;

	invert_affine_warp_Image< Eigen::VectorXd >(SourceImage.c_str(), registratedImage.c_str(), modelParameters);


	for (int i = 0;i < (int)modelParameters.size();++i)
	{
		cout<<modelParameters[i]<<endl;
	}
	cout<<"finished!"<<endl;

	fs.close();


	cv::Mat refImage = cv::imread(ReferImage.c_str());
	cv::Mat srcImage = cv::imread(SourceImage.c_str());
	cv::Mat img_matches;
	cv::Mat img_outlier;

	int thickness = 2;
	drawLineMatches(srcImage, matched_src_lines_list, refImage, matched_ref_lines_list,
		img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);

	drawLineMatches(srcImage, outlier_src_lines_list, refImage, outlier_ref_lines_list,
		img_outlier, cv::Scalar::all(-1), cv::Scalar::all(-1),
		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);

	cv::imwrite(matchedImage, img_matches);
	cv::imwrite(outlierImage, img_outlier);
}

void msld_line_match(string foldName, string sourceName, string referenceName)
{
	string SourceImage = foldName + "\\" + sourceName;//待配准影像
	string ReferImage = foldName + "\\" + referenceName;//参考影像
	string resultImage = foldName + "\\0.png";//结果影像
	string referImage_line = foldName + "\\lined1.png";//参考影像
	string sourceImage_line = foldName + "\\lined2.png";//待配准影像

	IplImage* ref = cvLoadImage(ReferImage.c_str(), 0 );
	IplImage* src = cvLoadImage(SourceImage.c_str(), 0 );

	if (!ref||!src)
	{
		cout<<"open image failed !"<<endl;
		system("pause");
		return;
	}

	vector<cvline_polar> ref_lines_polar_List;
	vector<cvline_polar> src_lines_polar_List;
	double sourceimage_p1 = 0.4;
	double sourceimage_p2 = 500;
	double referimage_p1 = 0.9;
	double referimage_p2 = 500;

	cvSaveImage(referImage_line.c_str(), get_lines(ref, ref_lines_polar_List, referimage_p1, referimage_p2), 0);
	cvSaveImage(sourceImage_line.c_str(), get_lines(src, src_lines_polar_List, sourceimage_p1, sourceimage_p2), 0);//提取出直线特征

	vector<vector<double>> ref_LineFeature;
	vector<vector<double>> src_LineFeature;
	vector<int>ref_LineNum;
	vector<int>src_LineNum;//直线下标（靠近图像R边缘的直线未求其平行域）
	LineFeature(ref,ref_lines_polar_List,ref_LineFeature,ref_LineNum);
	LineFeature(src,src_lines_polar_List,src_LineFeature,src_LineNum);//得到直线描述子

	//LineMatch(ref,src,ref_LineFeature,src_LineFeature,ref_LineNum,src_LineNum,ref_lines_polar_List,src_lines_polar_List);
	vector<cvline_polar> matched_ref_lines_list;
	vector<cvline_polar> matched_src_lines_list;
	LineMatch(ref_LineFeature,src_LineFeature,ref_LineNum,src_LineNum,ref_lines_polar_List,src_lines_polar_List, matched_ref_lines_list, matched_src_lines_list);

	cvReleaseImage(&ref);
	cvReleaseImage(&src);

	cv::Mat img_source = cv::imread(SourceImage);
	cv::Mat img_reference = cv::imread(ReferImage);
	cv::Mat img_matches;

	int thickness = 2;
	drawLineMatches(img_source, matched_src_lines_list, img_reference, matched_ref_lines_list,
		img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);

	resultImage = foldName + "\\segment_matched.png";//结果影像
	cv::imwrite(resultImage, img_matches);
}

void line_match(string foldName, string sourceName, string referenceName)
{
	string SourceImage = foldName + "\\" + sourceName;//待配准影像
	string ReferImage = foldName + "\\" + referenceName;//参考影像
	string resultImage = foldName + "\\0.png";//结果影像
	string sourceMatch = foldName + "\\line_" + sourceName;
	string referenceMatch = foldName + "\\line_" + referenceName;
	string matchedImage = foldName + "\\matched.png";//结果影像
	//const char*  SourceImage = "D:\\workspace\\testdata\\Landsat\\std\\p122r025_7t20000914.tif";//待配准影像
	//const char*  ReferImage = "D:\\workspace\\testdata\\features\\warp.jpg";//参考影像
	string referImage_line = foldName + "\\lined1.png";
	string sourceImage_line = foldName + "\\lined2.png";
	string line_angle = foldName + "\\angles.txt";
	string strLog = foldName + "\\processinglog.txt";
	fstream fs;
	fs.open(strLog.c_str(), std::ios_base::out);

	//warp(ReferImage, SourceImage);
	//perspective_warp(ReferImage, SourceImage);

	vector<cvline_polar> src_lines_polar;
	vector<cvline_polar> ref_lines_polar;

	double sourceimage_p1 = 0.8;
	double sourceimage_p2 = 500;
	double referimage_p1 = 0.9;
	double referimage_p2 = 5000;

	//cvSaveImage(sourceImage_line.c_str(), get_lines(src, src_lines_polar, 0.2), 0);
	//cvSaveImage(referImage_line.c_str(), get_lines(ref, ref_lines_polar, 0.8), 0);
	IplImage* src_line = get_lines(SourceImage.c_str(), src_lines_polar, sourceimage_p1, sourceimage_p2);
	IplImage* ref_line = get_lines(ReferImage.c_str(), ref_lines_polar, referimage_p1, referimage_p2);
	//IplImage* src_line = get_lines_lswms(SourceImage.c_str(), src_lines_polar, sourceimage_p2);
	//IplImage* ref_line = get_lines_lswms(ReferImage.c_str(), ref_lines_polar, referimage_p2);
	
	if (!src_line || !ref_line)
	{
		return;
	}
	//outLineAngle(line_angle, ref_lines_polar, src_lines_polar);
	cvSaveImage(sourceImage_line.c_str(), src_line, 0);
	cvSaveImage(referImage_line.c_str(), ref_line, 0);

	cvReleaseImage(&src_line);
	cvReleaseImage(&ref_line);

	//IplImage* src = cvLoadImage(SourceImage.c_str(), 0 );
	//if (!src)
	//{
	//	cout<<"open image failed for \""<<SourceImage<<"\"!"<<endl;
	//	system("pause");
	//	return;
	//}

	//IplImage* ref = cvLoadImage(ReferImage.c_str(), 0 );
	//if (!ref)
	//{
	//	cout<<"open image failed for \""<<ReferImage<<"\"!"<<endl;
	//	system("pause");
	//	return;
	//}

	//cvSaveImage(sourceImage_line.c_str(), get_lines(src, src_lines_polar, 0.8), 0);
	//cvSaveImage(referImage_line.c_str(), get_lines(ref, ref_lines_polar, 0.8), 0);


	fs<<"source image:"<<SourceImage<<endl;
	fs<<"scale parameter:"<<sourceimage_p1<<"\t number parameter:"<<sourceimage_p2<<endl;
	fs<<"detected lines:"<<src_lines_polar.size()<<endl;
	fs<<endl;

	fs<<"reference image:"<<ReferImage<<endl;
	fs<<"scale parameter:"<<referimage_p1<<"\t number parameter:"<<referimage_p2<<endl;
	fs<<"detected lines:"<<ref_lines_polar.size()<<endl;
	fs<<endl;

	cout<<src_lines_polar.size()<<" segments are detected in source image."<<endl;
	cout<<ref_lines_polar.size()<<" segments are detected in reference image."<<endl;
	ecmlrAffine* ecmlr = new ecmlrAffine(src_lines_polar, ref_lines_polar);
	ecmlr->setParameters(parameters);
	//ecmlrPerspective* ecmlr = new ecmlrPerspective(ref_lines_polar, src_lines_polar);
	ecmlr->m_referenceImageFile = SourceImage;
	ecmlr->m_sourceImageFile = ReferImage;
	ecmlr->setOutFileStream(&fs);
	vector<double> modelParameters;

	//用clock()来计时  毫秒
	clock_t  clockBegin, clockEnd;
	clockBegin = clock();
	ecmlr->solve(modelParameters);
	clockEnd = clock();
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(3);
	fs<<"time consuming:"<<(clockEnd - clockBegin)*1e-3<<"second"<<endl;

	//invert_affine_warp_Image< vector<double> >(ReferImage.c_str(), resultImage.c_str(), modelParameters);
	affine_warp_Image< vector<double> >(SourceImage.c_str(), resultImage.c_str(), modelParameters);


	for (int i = 0;i < (int)modelParameters.size();++i)
	{
		cout<<modelParameters[i]<<endl;
	}
	cout<<"finished!"<<endl;
	fs.close();
	
	int thickness = 2;
	//cv::Mat img_source = cv::imread(SourceImage);
	//cv::Mat img_reference = cv::imread(ReferImage);
	//drawLineMatches2(img_source, ecmlr->m_matched_src_ls, img_reference, ecmlr->m_matched_ref_ls,
	//	cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
	//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//cv::imwrite(sourceMatch, img_source);
	//cv::imwrite(referenceMatch, img_reference);

	//cv::Mat img_source = cv::imread(SourceImage);
	//cv::Mat img_reference = cv::imread(ReferImage);
	//cv::Mat img_outImage;
	//drawLineMatches(img_source, ecmlr->m_matched_src_ls, img_reference, ecmlr->m_matched_ref_ls,
	//	img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
	//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//cv::imwrite(matchedImage, img_outImage);
}

#include "AffineLineRegistration.h"
void line_match2(string foldName, string sourceName, string referenceName)
{
	string SourceImage = foldName + "\\" + sourceName;//待配准影像
	string ReferImage = foldName + "\\" + referenceName;//参考影像
	string resultImage = foldName + "\\0.png";//结果影像
	string sourceMatch = foldName + "\\line_" + sourceName;
	string referenceMatch = foldName + "\\line_" + referenceName;
	string referImage_line = foldName + "\\lined1.png";
	string sourceImage_line = foldName + "\\lined2.png";
	string line_angle = foldName + "\\angles.txt";
	string strLog = foldName + "\\processinglog.txt";
	fstream fs;
	fs.open(strLog.c_str(), std::ios_base::out);

	//warp(ReferImage, SourceImage);
	//perspective_warp(ReferImage, SourceImage);

	vector<cvline_polar> src_lines_polar;
	vector<cvline_polar> ref_lines_polar;

	double sourceimage_p1 = 0.8;
	double sourceimage_p2 = 500;
	double referimage_p1 = 0.8;
	double referimage_p2 = 500;

	//cvSaveImage(sourceImage_line.c_str(), get_lines(src, src_lines_polar, 0.2), 0);
	//cvSaveImage(referImage_line.c_str(), get_lines(ref, ref_lines_polar, 0.8), 0);
	IplImage* src_line = get_lines(SourceImage.c_str(), src_lines_polar, sourceimage_p1, sourceimage_p2);
	IplImage* ref_line = get_lines(ReferImage.c_str(), ref_lines_polar, referimage_p1, referimage_p2);
	if (!src_line || !ref_line)
	{
		return;
	}
	//outLineAngle(line_angle, ref_lines_polar, src_lines_polar);
	cvSaveImage(sourceImage_line.c_str(), src_line, 0);
	cvSaveImage(referImage_line.c_str(), ref_line, 0);

	cvReleaseImage(&src_line);
	cvReleaseImage(&ref_line);


	fs<<"source image:"<<SourceImage<<endl;
	fs<<"scale parameter:"<<sourceimage_p1<<"\t number parameter:"<<sourceimage_p2<<endl;
	fs<<"detected lines:"<<src_lines_polar.size()<<endl;
	fs<<endl;

	fs<<"reference image:"<<ReferImage<<endl;
	fs<<"scale parameter:"<<referimage_p1<<"\t number parameter:"<<referimage_p2<<endl;
	fs<<"detected lines:"<<ref_lines_polar.size()<<endl;
	fs<<endl;

	cout<<src_lines_polar.size()<<" segments are detected in source image."<<endl;
	cout<<ref_lines_polar.size()<<" segments are detected in reference image."<<endl;
	AffineLineRegistration* ecmlr = new AffineLineRegistration(src_lines_polar, ref_lines_polar);
	//ecmlrPerspective* ecmlr = new ecmlrPerspective(ref_lines_polar, src_lines_polar);
	ecmlr->m_referenceImageFile = SourceImage;
	ecmlr->m_sourceImageFile = ReferImage;
	ecmlr->setOutFileStream(&fs);
	vector<double> modelParameters;

	//用clock()来计时  毫秒
	clock_t  clockBegin, clockEnd;
	clockBegin = clock();
	ecmlr->solve(modelParameters);
	clockEnd = clock();
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(3);
	fs<<"time consuming:"<<(clockEnd - clockBegin)*1e-3<<"second"<<endl;

	invert_affine_warp_Image< vector<double> >(SourceImage.c_str(), resultImage.c_str(), modelParameters);


	for (int i = 0;i < (int)modelParameters.size();++i)
	{
		cout<<modelParameters[i]<<endl;
	}
	cout<<"finished!"<<endl;
	fs.close();

	int thickness = 2;
	//cv::Mat img_source = cv::imread(SourceImage);
	//cv::Mat img_reference = cv::imread(ReferImage);
	//drawLineMatches2(img_source, ecmlr->m_matched_src_ls, img_reference, ecmlr->m_matched_ref_ls,
	//	cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
	//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//cv::imwrite(sourceMatch, img_source);
	//cv::imwrite(referenceMatch, img_reference);

	cv::Mat img_source = cv::imread(SourceImage);
	cv::Mat img_reference = cv::imread(ReferImage);
	cv::Mat img_outImage;
	drawLineMatches(img_source, ecmlr->m_matched_src_ls, img_reference, ecmlr->m_matched_ref_ls,
		img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	cv::imwrite(resultImage, img_outImage);
}

//
//void line_match2()
//{
//	//string SourImage = "D:\\MyPaper\\MyTest\\data\\my_test\\MIHT\\qb\\qb_s.tif";//待配准影像
//	//string ReferImage = "D:\\MyPaper\\MyTest\\data\\my_test\\MIHT\\qb\\qb_r.tif";//参考影像
//	//string SourImage = "D:\\workspace\\testdata\\straight-line\\sitongqiao.jpg";//待配准影像
//	//string SourImage = "D:\\workspace\\testdata\\straight-line\\alos2.jpg";//待配准影像
//	//string ReferImage = "D:\\workspace\\testdata\\straight-line\\sitongqiao4.jpg";//参考影像
//	//string SourImage = "D:\\Loong\\Programs\\CV\\testdata\\polygons_de.jpg";//待配准影像
//	//string ReferImage = "D:\\Loong\\Programs\\CV\\testdata\\polygons.jpg";//参考影像
//	//string SourImage = "D:\\workspace\\testdata\\features\\1.jpg";//待配准影像
//	//string ReferImage = "D:\\workspace\\testdata\\features\\2.jpg";//参考影像
//	//string SourImage = "E:\\testdata\\straight-line\\sitongqiao.jpg";//待配准影像
//	//string ReferImage = "E:\\testdata\\straight-line\\sitongqiao4.jpg";//参考影像
//
//	//string SourImage = "D:\\MyPaper\\MyTest\\data\\my_test\\MIHT\\qb25.tif";//待配准影像
//	//string ReferImage = "D:\\MyPaper\\MyTest\\data\\my_test\\MIHT\\qb2m.tif";//参考影像
//	string foldName = "F:\\Dropbox\\Programs\\data\\ding";
//	string SourceImage = foldName + "\\1.png";//参考影像
//	string ReferImage = foldName + "\\3.png";//待配准影像
//	string resultImage = foldName + "\\0.png";//结果影像
//	//const char*  SourceImage = "D:\\workspace\\testdata\\Landsat\\std\\p122r025_7t20000914.tif";//待配准影像
//	//const char*  ReferImage = "D:\\workspace\\testdata\\features\\warp.jpg";//参考影像
//	string referImage_line = foldName + "\\lined1.png";//参考影像
//	string sourceImage_line = foldName + "\\lined2.png";//参考影像
//
//
//	//string SourOpen =Get_cannyOpen(SourImage,40);   //边缘提取
//	//string ReferOpen =Get_cannyOpen(ReferImage,40);
//
//	vector<cvline_polar> src_lines_polar;
//	vector<cvline_polar> ref_lines_polar;
//
//
//	IplImage* src = cvLoadImage(SourceImage.c_str(), 0 );
//	if (!src)
//	{
//		cout<<"open image failed for \""<<SourceImage<<"\"!"<<endl;
//		system("pause");
//		return;
//	}
//
//	IplImage* ref = cvLoadImage(ReferImage.c_str(), 0 );
//	if (!ref)
//	{
//		cout<<"open image failed for \""<<ReferImage<<"\"!"<<endl;
//		system("pause");
//		return;
//	}
//
//	cvSaveImage(sourceImage_line.c_str(), get_lines(src, src_lines_polar), 0);
//	cvSaveImage(referImage_line.c_str(), get_lines(ref, ref_lines_polar), 0);
//
//	//hough_single(SourceImage, src_lines_polar);
//	//hough_single(ReferImage, ref_lines_polar);
//
//	//Get_hough(SourImage,SourOpen,ReferImage,ReferOpen, src_lines_polar, ref_lines_polar);   //hough检测直线段，并存放入对应的数组中
//
//
//	Eigen::MatrixXd modelData((int)ref_lines_polar.size(), 2);
//	for( int i = 0; i < (int)ref_lines_polar.size(); i++ )
//	{
//		fPoint pt = line2point(ref_lines_polar[i]);
//		modelData(i, 0) = pt.x;
//		modelData(i, 1) = pt.y;
//	}
//
//	Eigen::MatrixXd observationData((int)src_lines_polar.size(), 2);
//	for( int i = 0; i < (int)src_lines_polar.size(); i++ )
//	{
//		fPoint pt = line2point(src_lines_polar[i]);
//		observationData(i, 0) = pt.x;
//		observationData(i, 1) = pt.y;
//	}
//
//	//fstream ofs;
//	//ofs.open("model.txt", ios_base::out);
//	//for (int i = 0;i < (int)ref_points.size();++i)
//	//{
//	//	ofs<<ref_points[i].x<<"\t"<<ref_points[i].y<<endl;
//	//}
//	//ofs.close();
//	//ofs.open("scene.txt", ios_base::out);
//	//for (int i = 0;i < (int)src_points.size();++i)
//	//{
//	//	ofs<<src_points[i].x<<"\t"<<src_points[i].y<<endl;
//	//}
//	//ofs.close();
//
//	ecmprAffine* mcem = new ecmprAffine(modelData, observationData);
//	mcem->m_sourceImageFile = ReferImage;
//	vector<double> modelParameters;
//	mcem->solve(modelParameters);
//
//	//ecmlrAffine* ecmlr = new ecmlrAffine(src_lines_polar, ref_lines_polar);
//	//ecmlr->m_sourceImageFile = SourceImage;
//	//vector<double> modelParameters;
//	//ecmlr->solve(modelParameters);
//
//	affine_warp_Image< vector<double> >(ReferImage.c_str(), resultImage.c_str(), modelParameters);
//}

void line_match(string foldName, string sourceName, string referenceName, 
				 const vector<cvline_polar>& src_lines_polar, const vector<cvline_polar>& ref_lines_polar)
{
	string SourceImage = foldName + "\\" + sourceName;//待配准影像
	string ReferImage = foldName + "\\" + referenceName;//参考影像
	string resultImage = foldName + "\\0.png";//结果影像
	string sourceMatch = foldName + "\\line_" + sourceName;
	string referenceMatch = foldName + "\\line_" + referenceName;
	string matchedImage = foldName + "\\matched.png";//结果影像
	string strLog = foldName + "\\processinglog.txt";
	fstream fs;
	fs.open(strLog.c_str(), std::ios_base::out);

	ecmlrAffine* ecmlr = new ecmlrAffine(src_lines_polar, ref_lines_polar);
	ecmlr->setParameters(parameters);
	//ecmlrPerspective* ecmlr = new ecmlrPerspective(ref_lines_polar, src_lines_polar);
	ecmlr->m_referenceImageFile = SourceImage;
	ecmlr->m_sourceImageFile = ReferImage;
	ecmlr->setOutFileStream(&fs);
	vector<double> modelParameters;

	//用clock()来计时  毫秒
	clock_t  clockBegin, clockEnd;
	clockBegin = clock();
	ecmlr->solve(modelParameters);
	clockEnd = clock();
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(3);
	fs<<"time consuming:"<<(clockEnd - clockBegin)*1e-3<<"second"<<endl;

	affine_warp_Image< vector<double> >(SourceImage.c_str(), resultImage.c_str(), modelParameters);


	for (int i = 0;i < (int)modelParameters.size();++i)
	{
		cout<<modelParameters[i]<<endl;
	}
	cout<<"finished!"<<endl;
	fs.close();

	int thickness = 2;

	cv::Mat img_source = cv::imread(SourceImage);
	cv::Mat img_reference = cv::imread(ReferImage);
	cv::Mat img_outImage;
	drawLineMatches(img_source, ecmlr->m_matched_src_ls, img_reference, ecmlr->m_matched_ref_ls,
		img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	cv::imwrite(matchedImage, img_outImage);
}


void rotationMatch(string foldName, double angle, string referenceName,
				   const vector<cvline_polar>& src_lines_polar, const vector<cvline_polar>& ref_lines_polar)
{
	if (angle > 360.0)
	{
		angle = angle - 360.0;
	}
	std::ostringstream s;
	s << angle;
	string strAngle = s.str();
	string subFoldName = foldName + "\\" + strAngle;
	if (!fexists(subFoldName))
	{
		_mkdir(subFoldName.c_str());
	}

	cv::Mat img = cv::imread(foldName + "\\" +referenceName);
	int iRows = img.rows;
	int iCols = img.cols;
	int y0 = iRows / 2;
	int x0 = iCols / 2;

	cv::Mat affineTrans(2, 3, CV_64F);
	double angle_ = angle * 3.14159265 / 180.0;
	double cos_theta = cos(angle_);
	double sin_theta = sin(angle_);

	double a0 = x0*cos_theta+y0*sin_theta;
	double b0 = -x0*sin_theta+y0*cos_theta;
	affineTrans.at<double>(0, 2) = x0-a0;
	affineTrans.at<double>(0, 0) = cos_theta;
	affineTrans.at<double>(0, 1) = sin_theta;
	affineTrans.at<double>(1, 2) = y0-b0;
	affineTrans.at<double>(1, 0) = -sin_theta;
	affineTrans.at<double>(1, 1) = cos_theta;
	//parameters[0] = affineTrans.at<double>(0, 2) = x0-a0;
	//parameters[1] = affineTrans.at<double>(0, 0) = cos_theta;
	//parameters[2] = affineTrans.at<double>(0, 1) = sin_theta;
	//parameters[3] = affineTrans.at<double>(1, 2) = y0-b0;
	//parameters[4] = affineTrans.at<double>(1, 0) = -sin_theta;
	//parameters[5] = affineTrans.at<double>(1, 1) = cos_theta;

	//parameters[0] = 0.0;
	//parameters[3] = 0.0;

	cv::Mat invert_affineTrans;
	invertAffineTransform(affineTrans, invert_affineTrans);

	cv::Size newSize = getNewSize(iRows, iCols, invert_affineTrans);
	std::vector<cvline_polar> new_src_line_polar;
	for (int i = 0;i < (int)src_lines_polar.size();++i)
	{
		new_src_line_polar.push_back(affineTransLine(invert_affineTrans, src_lines_polar[i]));
	}

	cv::Mat newImg;
	warpAffine(img, newImg, invert_affineTrans , cv::Size(iRows, iCols), cv::INTER_CUBIC);


	imwrite(subFoldName + "\\2.png", newImg);
	CopyFile((foldName + "\\" +referenceName).c_str(),(subFoldName + "\\" +referenceName).c_str(), false);

	line_match(subFoldName, "2.png", referenceName, new_src_line_polar, ref_lines_polar);
}

void rotationMatch(string foldName, double angle, string referenceName)
{
	if (angle > 360.0)
	{
		angle = angle - 360.0;
	}
	std::ostringstream s;
	s << angle;
	string strAngle = s.str();
	string subFoldName = foldName + "\\" + strAngle;
	if (!fexists(subFoldName))
	{
		_mkdir(subFoldName.c_str());
	}

	cv::Mat img = cv::imread(foldName + "\\" +referenceName);
	int iRows = img.rows;
	int iCols = img.cols;
	int y0 = iRows / 2;
	int x0 = iCols / 2;

	cv::Mat affineTrans(2, 3, CV_64F);
	double angle_ = angle * 3.14159265 / 180.0;
	double cos_theta = cos(angle_);
	double sin_theta = sin(angle_);

	double a0 = x0*cos_theta+y0*sin_theta;
	double b0 = -x0*sin_theta+y0*cos_theta;
	affineTrans.at<double>(0, 2) = x0-a0;
	affineTrans.at<double>(0, 0) = cos_theta;
	affineTrans.at<double>(0, 1) = sin_theta;
	affineTrans.at<double>(1, 2) = y0-b0;
	affineTrans.at<double>(1, 0) = -sin_theta;
	affineTrans.at<double>(1, 1) = cos_theta;
	//parameters[0] = affineTrans.at<double>(0, 2) = x0-a0;
	//parameters[1] = affineTrans.at<double>(0, 0) = cos_theta;
	//parameters[2] = affineTrans.at<double>(0, 1) = sin_theta;
	//parameters[3] = affineTrans.at<double>(1, 2) = y0-b0;
	//parameters[4] = affineTrans.at<double>(1, 0) = -sin_theta;
	//parameters[5] = affineTrans.at<double>(1, 1) = cos_theta;

	//parameters[0] = 0.0;
	//parameters[3] = 0.0;

	cv::Mat invert_affineTrans;
	invertAffineTransform(affineTrans, invert_affineTrans);

	cv::Size newSize = getNewSize(iRows, iCols, invert_affineTrans);

	cv::Mat newImg;
	warpAffine(img, newImg, invert_affineTrans , cv::Size(iRows, iCols), cv::INTER_CUBIC);


	imwrite(subFoldName + "\\2.png", newImg);
	CopyFile((foldName + "\\" +referenceName).c_str(),(subFoldName + "\\" +referenceName).c_str(), false);

	line_match(subFoldName, "2.png", referenceName);
}

void get_line_candidate(double angle, string foldName, string referenceName,
						std::vector<cvline_polar>& src_lines_polar, std::vector<cvline_polar>& ref_lines_polar)
{
	cv::Mat img = cv::imread(foldName + "\\" +referenceName);
	int iRows = img.rows;
	int iCols = img.cols;
	int y0 = iRows / 2;
	int x0 = iCols / 2;

	cv::Mat affineTrans(2, 3, CV_64F);
	double angle_ = angle * 3.14159265 / 180.0;
	double cos_theta = cos(angle_);
	double sin_theta = sin(angle_);

	double a0 = x0*cos_theta+y0*sin_theta;
	double b0 = -x0*sin_theta+y0*cos_theta;
	affineTrans.at<double>(0, 2) = x0-a0;
	affineTrans.at<double>(0, 0) = cos_theta;
	affineTrans.at<double>(0, 1) = sin_theta;
	affineTrans.at<double>(1, 2) = y0-b0;
	affineTrans.at<double>(1, 0) = -sin_theta;
	affineTrans.at<double>(1, 1) = cos_theta;

	cv::Mat invert_affineTrans;
	invertAffineTransform(affineTrans, invert_affineTrans);

	cv::Size newSize = getNewSize(iRows, iCols, invert_affineTrans);

	cv::Mat newImg;
	warpAffine(img, newImg, invert_affineTrans , cv::Size(iRows, iCols), cv::INTER_CUBIC);

	string SourceImage = foldName + "\\2.png";
	string ReferImage = foldName + "\\" +referenceName;
	imwrite(SourceImage, newImg);

	double sourceimage_p1 = 0.7;
	double sourceimage_p2 = 500;
	double referimage_p1 = 0.7;
	double referimage_p2 = 500;
	vector<cvline_polar> temp_src_lines_polar;

	src_lines_polar.clear();
	ref_lines_polar.clear();
	get_lines(SourceImage.c_str(), temp_src_lines_polar, sourceimage_p1, sourceimage_p2);
	get_lines(ReferImage.c_str(), ref_lines_polar, referimage_p1, referimage_p2);
	//get_lines_lswms(SourceImage.c_str(), temp_src_lines_polar, sourceimage_p2);
	//get_lines_lswms(ReferImage.c_str(), ref_lines_polar, referimage_p2);

	for (int i = 0;i < (int)temp_src_lines_polar.size();++i)
	{
		src_lines_polar.push_back(affineTransLine(affineTrans, temp_src_lines_polar[i]));
	}
}

void rotationTest(double initial_angle, string foldName, string referenceName)
{
	std::ostringstream s;
	s << initial_angle<<"_initial";
	string strAngle = s.str();
	string subFoldName = foldName + "\\" + strAngle;
	if (!fexists(subFoldName))
	{
		_mkdir(subFoldName.c_str());
	}
	CopyFile((foldName + "\\" +referenceName).c_str(),(subFoldName + "\\" +referenceName).c_str(), false);

	cv::Mat img = cv::imread(foldName + "\\" +referenceName);
	int iRows = img.rows;
	int iCols = img.cols;
	int y0 = iRows / 2;
	int x0 = iCols / 2;
	cv::Mat affineTrans(2, 3, CV_64F);
	double angle_ = initial_angle * 3.14159265 / 180.0;
	double cos_theta = cos(angle_);
	double sin_theta = sin(angle_);

	double a0 = x0*cos_theta+y0*sin_theta;
	double b0 = -x0*sin_theta+y0*cos_theta;
	parameters[0] = x0-a0;
	parameters[1] = cos_theta;
	parameters[2] = sin_theta;
	parameters[3] = y0-b0;
	parameters[4] = -sin_theta;
	parameters[5] = cos_theta;

	for (int i = 0;i < 7;++i)
	{
		//rotationMatch(subFoldName, initial_angle+i*10.0, referenceName);
	}
	for (int i = 35;i > 29;--i)
	{
		rotationMatch(subFoldName, initial_angle+i*10.0, referenceName);
	}

	//vector<cvline_polar> src_lines_polar;
	//vector<cvline_polar> ref_lines_polar;
	//get_line_candidate(20.0, subFoldName, referenceName,
	//	 src_lines_polar, ref_lines_polar);
	//for (int i = 0;i < 6;++i)
	//{
	//	rotationMatch(subFoldName, initial_angle+i*10.0, referenceName, src_lines_polar, ref_lines_polar);
	//}
	//get_line_candidate(-30.0, subFoldName, referenceName,
	//	src_lines_polar, ref_lines_polar);
	//for (int i = 35;i > 30;--i)
	//{
	//	rotationMatch(subFoldName, initial_angle+i*10.0, referenceName, src_lines_polar, ref_lines_polar);
	//}
}

void rotationTest(string foldName, string referenceName)
{
	// 0
	rotationTest(0.0, foldName, referenceName);

	// 90
	//rotationTest(90.0, foldName, referenceName);


	//180
	//rotationTest(180.0, foldName, referenceName);


	//270
	//rotationTest(270.0, foldName, referenceName);

}

void salt(cv::Mat& image, double percentage)  
{  
	int nCols = image.cols;
	int nRows = image.rows;
	int n = percentage * nCols * nRows;
	for(int k=0; k<n; k++)  
	{  
		int i = rand()%image.cols;  
		int j = rand()%image.rows;  

		if(image.channels() == 1)  
		{  
			image.at<uchar>(j,i) = 255;  
		}  
		else  
		{  
			image.at<Vec3b>(j,i)[0] = 255;  
			image.at<Vec3b>(j,i)[1] = 255;  
			image.at<Vec3b>(j,i)[2] = 255;  
		}  
	}  
} 

void noiseMatch(string foldName, double dev, string referenceName)
{
	std::ostringstream s;
	s << dev;
	string strAngle = s.str();
	string subFoldName = foldName + "\\" + strAngle;
	if (!fexists(subFoldName))
	{
		_mkdir(subFoldName.c_str());
	}

	cv::Mat img = cv::imread(foldName + "\\" + referenceName);
	cv::Mat gaussian_noise = img.clone();
	randn(gaussian_noise, 0, dev);
	img += gaussian_noise;
	//salt(img, dev*0.01) ;

	// rotate
	int iRows = img.rows;
	int iCols = img.cols;
	int y0 = iRows / 2;
	int x0 = iCols / 2;

	cv::Mat affineTrans(2, 3, CV_64F);
	double angle_ = 20.0 * 3.14159265 / 180.0;
	double cos_theta = cos(angle_);
	double sin_theta = sin(angle_);

	double a0 = x0*cos_theta+y0*sin_theta;
	double b0 = -x0*sin_theta+y0*cos_theta;
	affineTrans.at<double>(0, 2) = x0-a0;
	affineTrans.at<double>(0, 0) = cos_theta;
	affineTrans.at<double>(0, 1) = sin_theta;
	affineTrans.at<double>(1, 2) = y0-b0;
	affineTrans.at<double>(1, 0) = -sin_theta;
	affineTrans.at<double>(1, 1) = cos_theta;

	cv::Mat invert_affineTrans;
	invertAffineTransform(affineTrans, invert_affineTrans);

	cv::Size newSize = getNewSize(iRows, iCols, invert_affineTrans);

	cv::Mat newImage;
	warpAffine(img, newImage, invert_affineTrans , cv::Size(iRows, iCols), cv::INTER_CUBIC);
	
	imwrite(subFoldName + "\\2.png", newImage);

	CopyFile((foldName + "\\" +referenceName).c_str(),(subFoldName + "\\" +referenceName).c_str(), false);

	line_match(subFoldName, "2.png", referenceName);
}

void noiseTest(string foldName, string referenceName)
{
	//noiseMatch(foldName, 14.0, referenceName);
	noiseMatch(foldName, 46.0, referenceName);
	for (int i = 55;i < 61;++i)
	{
		noiseMatch(foldName, i*2.0, referenceName);
	}
}

void blurMatch(string foldName, int ksize, string referenceName)
{
	std::ostringstream s;
	s << ksize;
	string strAngle = s.str();
	string subFoldName = foldName + "\\" + strAngle;
	if (!fexists(subFoldName))
	{
		_mkdir(subFoldName.c_str());
	}

	cv::Mat img = cv::imread(foldName + "\\2.png");
	cv::GaussianBlur(img, img, cv::Size(ksize, ksize), 0.0);
	imwrite(subFoldName + "\\2.png", img);
	CopyFile((foldName + "\\" +referenceName).c_str(),(subFoldName + "\\" +referenceName).c_str(), false);

	line_match(subFoldName, "2.png", referenceName);
}

void blurTest(string foldName, string referenceName)
{
	blurMatch(foldName, 9, referenceName);
	for (int i = 0;i < 20;++i)
	{
		blurMatch(foldName, i*2+1, referenceName);
	}
}

void brightMatch(string foldName, double alpha, double beta, string referenceName)
{
	std::ostringstream s;
	s << alpha<<"_"<<beta;
	string strAngle = s.str();
	string subFoldName = foldName + "\\" + strAngle;
	if (!fexists(subFoldName))
	{
		_mkdir(subFoldName.c_str());
	}

	cv::Mat img = cv::imread(foldName + "\\"+referenceName);
	/// Do the operation new_image(i,j) = alpha*image(i,j) + beta
	for( int y = 0; y < img.rows; y++ )
	{ 
		for( int x = 0; x < img.cols; x++ )
		{ 
			for( int c = 0; c < 3; c++ )
			{
				img.at<Vec3b>(y,x)[c] =
				saturate_cast<uchar>( alpha*( img.at<Vec3b>(y,x)[c] ) + beta );
			}
		}
	}
	// rotate
	int iRows = img.rows;
	int iCols = img.cols;
	int y0 = iRows / 2;
	int x0 = iCols / 2;

	cv::Mat affineTrans(2, 3, CV_64F);
	double angle_ = 20.0 * 3.14159265 / 180.0;
	double cos_theta = cos(angle_);
	double sin_theta = sin(angle_);

	double a0 = x0*cos_theta+y0*sin_theta;
	double b0 = -x0*sin_theta+y0*cos_theta;
	affineTrans.at<double>(0, 2) = x0-a0;
	affineTrans.at<double>(0, 0) = cos_theta;
	affineTrans.at<double>(0, 1) = sin_theta;
	affineTrans.at<double>(1, 2) = y0-b0;
	affineTrans.at<double>(1, 0) = -sin_theta;
	affineTrans.at<double>(1, 1) = cos_theta;

	cv::Mat invert_affineTrans;
	invertAffineTransform(affineTrans, invert_affineTrans);

	cv::Size newSize = getNewSize(iRows, iCols, invert_affineTrans);

	cv::Mat newImage;
	warpAffine(img, newImage, invert_affineTrans , cv::Size(iRows, iCols), cv::INTER_CUBIC);
	
	imwrite(subFoldName + "\\2.png", newImage);
	CopyFile((foldName + "\\" +referenceName).c_str(),(subFoldName + "\\" +referenceName).c_str(), false);

	line_match(subFoldName, "2.png", referenceName);
}

void brightTest(string foldName, string referenceName)
{
	for (int i = -11;i < 12;++i)
	{
		brightMatch(foldName, 0.25, i*20.0, referenceName);
	}
	//for (int i = -10;i < 12;++i)
	//{
	//	if(i>=-10 && i <10) continue;
	//	brightMatch(foldName, 0.3, i*20.0, referenceName);
	//}
	//for (int i = -10;i < 12;++i)
	//{
	//	if(i>=-10 && i <10) continue;
	//	brightMatch(foldName, 0.5, i*20.0, referenceName);
	//}
	//for (int i = -10;i < 12;++i)
	//{
	//	if(i>=-10 && i <10) continue;
	//	brightMatch(foldName, 1.0, i*20.0, referenceName);
	//}
	//for (int i = -5;i < 6;++i)
	//{
	//	brightMatch(foldName, 2.0, i*20.0, referenceName);
	//}
	//for (int i = -5;i < 6;++i)
	//{
	//	brightMatch(foldName, 3.0, i*20.0, referenceName);
	//}
}

void main(int argc, char** argv)
{
	string dropbox = getenv("DROPBOX");
	//string foldName = dropbox + "\\Programs\\data\\nx\\compare\\em\\noRansac";
	string foldName = dropbox + "\\Programs\\data\\envisat";
	//string foldName = dropbox + "\\Programs\\data\\invariant\\noise";
	//string foldName = dropbox + "\\Programs\\data\\msld\\bright";

	//string foldName = dropbox + "\\Programs\\data\\uva";
	//string foldName = dropbox + "\\Document\\Notes\\Research\\Straight-line_Match\\JSTARS\\pic\\nx\\lsr";
	string sourceName = "2.tif";
	string referenceName = "1.tif";
	parameters[0] =0.0;
	parameters[1] =1.0;
	parameters[2] =0.0;
	parameters[3] =0.0;
	parameters[4] =0.0;
	parameters[5] =1.0;
	//fft();
	//orb_match();
	//surf_match();
	//point_match2();
	//foldName = dropbox + "\\Programs\\matching\\compare\\TestLMDll\\Release";
	//sourceName = "L503_A.jpg";
	//referenceName = "L503_B.jpg";
	line_match(foldName, sourceName, referenceName);
	//rotationMatch(foldName, 90.0, referenceName);
	//line_match2(foldName, sourceName, referenceName);
	//lbd_line_match(foldName, sourceName, referenceName);
	//lp_line_match(foldName, sourceName, referenceName);
	//msld_line_match(foldName, sourceName, referenceName);
	//line_match2();
	//cout<<"Hello~"<<endl;
	//rotationTest(foldName, referenceName);
	//brightTest(foldName, referenceName);
	//blurTest(foldName, referenceName);
	//noiseTest(foldName, referenceName);

	////string resultFold = "E:\\Dropbox\\Document\\Notes\\Research\\Straight-line_Match\\JSTARS\\pic\\tj";
	//string resultFold =  dropbox + "\\Programs\\line\\SegmentRegistration\\SegmentRegistration\\pic\\aerial";
	//string resultFold =  dropbox + "\\Programs\\data\\nx\\ls";
	//drawMatchedSegments(foldName, sourceName, referenceName, "lp_lines.txt");
	//drawMatchedSegments2(foldName, sourceName, referenceName, "lines.txt");
	drawOutliers(foldName, sourceName, referenceName, "lines.txt");
	system("pause");
}
