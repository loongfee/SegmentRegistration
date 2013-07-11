#include "lineFeature.h"
#include "cv.h"
#include <stdio.h>
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "func.h"
#include "LineMatchFLANN.h"
using namespace cv;

//void LineMatch(IplImage* ref_seg,IplImage* src_seg,vector<vector<double>> &ref_LineFeature,vector<vector<double>> &src_LineFeature,
//	vector<int>&ref_LineNum,vector<int>&src_LineNum,vector<cvline_polar>&ref_lines_list,vector<cvline_polar>&src_lines_list,
//	vector<cvline_polar>& matched_ref_lines_list, vector<cvline_polar>& matched_src_lines_list)
void LineMatch(vector<vector<double>> &ref_LineFeature,vector<vector<double>> &src_LineFeature,
	vector<int>&ref_LineNum,vector<int>&src_LineNum,vector<cvline_polar>&ref_lines_list,vector<cvline_polar>&src_lines_list,
	vector<cvline_polar>& matched_ref_lines_list, vector<cvline_polar>& matched_src_lines_list)
{
	//string RefImageLine = "data\\r1_s.png";//提取出直线的参考影像
	//string SrcImageLine = "data\\r2_s.png";//提取出直线的源影像
	//ref_seg = cvLoadImage(RefImageLine.c_str(), 1 );
	//src_seg = cvLoadImage(SrcImageLine.c_str(), 1 );

	//if (!ref_seg||!src_seg)
	//{
	//	cout<<"open image after segmentation failed !"<<endl;
	//	system("pause");
	//	return;
	//}
	
	int N_ref=ref_LineFeature.size();//参考影像直线个数
	int N_src=src_LineFeature.size();//源影像直线个数
	int M=ref_LineFeature[0].size();//描述子维数
	int i,j;
   FlannBasedMatcher matcher;
   std::vector< DMatch > matches;
  CvMat* descriptors_ref=cvCreateMat(N_ref,M,CV_32FC1);
  CvMat* descriptors_src=cvCreateMat(N_src,M,CV_32FC1);
  for(i=0;i!=N_ref;i++)
  {
    for(j=0;j!=M;j++)
	{
		cvmSet(descriptors_ref,i,j,ref_LineFeature[i][j]);
	}
  }
  
  for(i=0;i!=N_src;i++)
  {
    for(j=0;j!=M;j++)
	{
		cvmSet(descriptors_src,i,j,src_LineFeature[i][j]);
	}
  }
  matcher.match( descriptors_ref, descriptors_src, matches );

  double max_dist = 0; double min_dist = 100;

  //-- Quick calculation of max and min distances between keypoints
  for(i = 0; i < N_ref; i++ )
  { double dist = matches[i].distance;
    if( dist < min_dist ) min_dist = dist;
    if( dist > max_dist ) max_dist = dist;
  }

  printf("-- Max dist : %f \n", max_dist );
  printf("-- Min dist : %f \n", min_dist );

  //-- Draw only "good" matches (i.e. whose distance is less than 2*min_dist )
  //-- PS.- radiusMatch can also be used here.
  std::vector< DMatch > good_matches;
  for( i = 0; i < N_ref; i++ )
  { 
	  if( matches[i].distance <1.3*min_dist )
    { good_matches.push_back( matches[i]); }
  }

  //-- Draw only "good" matches
  for( i = 0; i < good_matches.size(); i++ )
  { printf( "-- Good Match [%d] ref_seg: %d  -- src_seg: %d  \n", i, good_matches[i].queryIdx, good_matches[i].trainIdx ); }
  CvFont font;
  cvInitFont(&font,CV_FONT_HERSHEY_SIMPLEX,0.5,0.5,0,1,8);

  matched_ref_lines_list.clear();
  matched_src_lines_list.clear();
  for(i = 0; i < good_matches.size(); i++)
  {
	  matched_ref_lines_list.push_back(ref_lines_list[ref_LineNum[good_matches[i].queryIdx]]);
	  matched_src_lines_list.push_back(src_lines_list[src_LineNum[good_matches[i].trainIdx]]);
	 //double rx=( ref_lines_list[ref_LineNum[good_matches[i].queryIdx]].pt1.x+ref_lines_list[ref_LineNum[good_matches[i].queryIdx]].pt2.x)/2;
	 //double ry=( ref_lines_list[ref_LineNum[good_matches[i].queryIdx]].pt1.y+ref_lines_list[ref_LineNum[good_matches[i].queryIdx]].pt2.y)/2;
	 //double sx=( src_lines_list[src_LineNum[good_matches[i].trainIdx]].pt1.x+src_lines_list[src_LineNum[good_matches[i].trainIdx]].pt2.x)/2;
	 //double sy=( src_lines_list[src_LineNum[good_matches[i].trainIdx]].pt1.y+src_lines_list[src_LineNum[good_matches[i].trainIdx]].pt2.y)/2;
	 //char num[20];
  //   itoa(i,num,10);
	 //cvPutText(ref_seg,num,cvPoint(rx,ry),&font,CV_RGB(255,0,0));
	 //cvPutText(src_seg,num,cvPoint(sx,sy),&font,CV_RGB(255,0,0));
	 //cvSaveImage("data\\r1_match.png",ref_seg);
	 //cvSaveImage("data\\r2_match.png",src_seg);
  }
  //   cvNamedWindow("ref",0);
	 //cvNamedWindow("src",0);
	 //cvShowImage("ref",ref_seg);
	 //cvShowImage("src",src_seg);
	 //cvWaitKey(0);
	 //system("pause");
	 //cvReleaseImage(&ref_seg);
	 //cvReleaseImage(&src_seg);
	 cvDestroyAllWindows();
    
}


void LineMatch(const ScaleLines &ref_LineFeature, const ScaleLines &src_LineFeature,
	vector<cvline_polar>& matched_ref_lines_list, vector<cvline_polar>& matched_src_lines_list)
{

	int N_ref = (int)ref_LineFeature.size();			//参考影像直线个数
	int N_src = (int)src_LineFeature.size();			//源影像直线个数
	int nOctave = (int)ref_LineFeature[0].size();	//尺度个数
	int M=0;																//描述子维数
	for (int i = 0;i < nOctave;++i)
	{
		M += (int)ref_LineFeature[0][i].descriptor.size();
	}

	FlannBasedMatcher matcher;
	std::vector< DMatch > matches;
	CvMat* descriptors_ref=cvCreateMat(N_ref, M, CV_32FC1);
	CvMat* descriptors_src=cvCreateMat(N_src, M, CV_32FC1);
	for(int i=0;i<N_ref;i++)
	{
		int nCount = 0;
		for (int j = 0;j < nOctave;++j)
		{
			int nDes = (int)ref_LineFeature[i][j].descriptor.size();
			for (int k = 0; k < nDes; k++)
			{
				cvmSet(descriptors_ref, i, nCount++, ref_LineFeature[i][j].descriptor[k]);
			}
		}
	}

	for(int i=0;i<N_src;i++)
	{
		int nCount = 0;
		for (int j = 0;j < nOctave;++j)
		{
			int nDes = (int)src_LineFeature[i][j].descriptor.size();
			for (int k = 0; k < nDes; k++)
			{
				cvmSet(descriptors_src, i, nCount++, src_LineFeature[i][j].descriptor[k]);
			}
		}
	}

	matcher.match( descriptors_ref, descriptors_src, matches );

	double max_dist = 0; double min_dist = 100;

	//-- Quick calculation of max and min distances between keypoints
	for(int i = 0; i < N_ref; i++ )
	{ double dist = matches[i].distance;
	if( dist < min_dist ) min_dist = dist;
	if( dist > max_dist ) max_dist = dist;
	}

	printf("-- Max dist : %f \n", max_dist );
	printf("-- Min dist : %f \n", min_dist );

	//-- Draw only "good" matches (i.e. whose distance is less than 2*min_dist )
	//-- PS.- radiusMatch can also be used here.
	std::vector< DMatch > good_matches;
	for(int i = 0; i < N_ref; i++ )
	{ 
		if( matches[i].distance <1.7*min_dist )
		{ good_matches.push_back( matches[i]); }
	}
	good_matches = matches;

	//-- Draw only "good" matches
	for(int i = 0; i < good_matches.size(); i++ )
	{ printf( "-- Good Match [%d] ref_seg: %d  -- src_seg: %d  \n", i, good_matches[i].queryIdx, good_matches[i].trainIdx ); }
	CvFont font;
	cvInitFont(&font,CV_FONT_HERSHEY_SIMPLEX,0.5,0.5,0,1,8);

	matched_ref_lines_list.clear();
	matched_src_lines_list.clear();
	for(int i = 0; i < good_matches.size(); i++)
	{
		cvline_polar refLine;
		cvline_polar srcLine;

		fPoint line[2];
		line[0].x = ref_LineFeature[good_matches[i].queryIdx][0].startPointX;
		line[0].y = ref_LineFeature[good_matches[i].queryIdx][0].startPointY;
		line[1].x = ref_LineFeature[good_matches[i].queryIdx][0].endPointX;
		line[1].y = ref_LineFeature[good_matches[i].queryIdx][0].endPointY;
		refLine = line2polar(line);

		line[0].x = src_LineFeature[good_matches[i].trainIdx][0].startPointX;
		line[0].y = src_LineFeature[good_matches[i].trainIdx][0].startPointY;
		line[1].x = src_LineFeature[good_matches[i].trainIdx][0].endPointX;
		line[1].y = src_LineFeature[good_matches[i].trainIdx][0].endPointY;
		srcLine = line2polar(line);

		matched_ref_lines_list.push_back(refLine);
		matched_src_lines_list.push_back(srcLine);
	}

}