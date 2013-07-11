/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

#include "draw.h"

using namespace std;

//const int draw_shift_bits = 4;
//const int draw_multiplier = 1 << draw_shift_bits;

using namespace cv;

///*
// * Functions to draw keypoints and matches.
// */
//void _drawKeyline( cv::Mat& img, const cvline_polar& l, const cv::Scalar& color, int flags, int thickness/*=1*/)
//{
//	cv::line(img, cvPoint(cvRound(l.pt1.x), cvRound(l.pt1.y)),
//		cvPoint(cvRound(l.pt2.x), cvRound(l.pt2.y)), color, thickness, CV_AA);
//}
//
//void drawKeylines( const cv::Mat& image, const vector<cvline_polar>& keylines, cv::Mat& outImg,
//                    const Scalar& _color, int flags, int thickness/*=1*/)
//{
//    if( !(flags & DrawMatchesFlags::DRAW_OVER_OUTIMG) )
//        cvtColor( image, outImg, CV_GRAY2BGR );
//
//    RNG& rng=theRNG();
//    bool isRandColor = _color == Scalar::all(-1);
//
//    for( vector<cvline_polar>::const_iterator i = keylines.begin(), ie = keylines.end(); i != ie; ++i )
//    {
//        Scalar color = isRandColor ? Scalar(rng(256), rng(256), rng(256)) : _color;
//        _drawKeyline( outImg, *i, color, flags, thickness);
//    }
//}
//
//
//void _prepareImgAndDrawKeylines( const cv::Mat& img1, const vector<cvline_polar>& keylines1,
//									   const cv::Mat& img2, const vector<cvline_polar>& keylines2,
//									   cv::Mat& outImg, cv::Mat& outImg1, cv::Mat& outImg2,
//									   const Scalar& singlePointColor, int flags )
//{
//	Size size( img1.cols + img2.cols, MAX(img1.rows, img2.rows) );
//	if( flags & DrawMatchesFlags::DRAW_OVER_OUTIMG )
//	{
//		if( size.width > outImg.cols || size.height > outImg.rows )
//			CV_Error( CV_StsBadSize, "outImg has size less than need to draw img1 and img2 together" );
//		outImg1 = outImg( Rect(0, 0, img1.cols, img1.rows) );
//		outImg2 = outImg( Rect(img1.cols, 0, img2.cols, img2.rows) );
//	}
//	else
//	{
//		outImg.create( size, CV_MAKETYPE(img1.depth(), 3) );
//		outImg1 = outImg( Rect(0, 0, img1.cols, img1.rows) );
//		outImg2 = outImg( Rect(img1.cols, 0, img2.cols, img2.rows) );
//
//		if( img1.type() == CV_8U )
//			cvtColor( img1, outImg1, CV_GRAY2BGR );
//		else
//			img1.copyTo( outImg1 );
//
//		if( img2.type() == CV_8U )
//			cvtColor( img2, outImg2, CV_GRAY2BGR );
//		else
//			img2.copyTo( outImg2 );
//	}
//
//	// draw keypoints
//	if( !(flags & DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS) )
//	{
//		cv::Mat outImg1 = outImg( Rect(0, 0, img1.cols, img1.rows) );
//		drawKeylines( outImg1, keylines1, outImg1, singlePointColor, flags + DrawMatchesFlags::DRAW_OVER_OUTIMG );
//
//		cv::Mat outImg2 = outImg( Rect(img1.cols, 0, img2.cols, img2.rows) );
//		drawKeylines( outImg2, keylines2, outImg2, singlePointColor, flags + DrawMatchesFlags::DRAW_OVER_OUTIMG );
//	}
//}
//
//void _drawLineMatch( cv::Mat& outImg, cv::Mat& outImg1, cv::Mat& outImg2 ,
//							  const cvline_polar& l1, const cvline_polar& l2, const Scalar& matchColor, int flags, int thickness/*=1*/)
//{
//	Scalar color(0, 0, 255);
//	//RNG& rng = theRNG();
//	//bool isRandMatchColor = matchColor == Scalar::all(-1);
//	//color = isRandMatchColor ? Scalar( rng(256), rng(256), rng(256) ) : matchColor;
//
//	_drawKeyline( outImg1, l1, color, flags );
//	_drawKeyline( outImg2, l2, color, flags );
//
//	cv::Point pt1((l1.pt1.x + l1.pt2.x)*0.5, (l1.pt1.y + l1.pt2.y)*0.5);
//	cv::Point pt2( std::min((l2.pt1.x + l2.pt2.x)*0.5+outImg1.cols, double(outImg.cols-1)), (l2.pt1.y + l2.pt2.y)*0.5 );
//	//cv::line(outImg, pt1, pt2, color, thickness, CV_AA);
//}
//
//
//void _drawLabel( cv::Mat& outImg, const cvline_polar& line,const cv::Scalar& matchColor, int flags, string strName/*=""*/)
//{
//	Point2f pt1((line.pt1.x + line.pt2.x)*0.5, (line.pt1.y + line.pt2.y)*0.5);
//
//	int baseline=0;
//	int fontFace = CV_FONT_HERSHEY_SIMPLEX;
//	double fontScale = 0.4;
//	int thickness = 1;
//	Size textSize = getTextSize(strName, fontFace,
//		fontScale, thickness, &baseline);
//	baseline += thickness;
//	baseline = 2;
//
//	rectangle(outImg, cv::Point(pt1.x, pt1.y) + Point(0, baseline),
//		cv::Point(pt1.x, pt1.y) + Point(textSize.width, -textSize.height) - Point(0, baseline),
//		Scalar(255,255,255), CV_FILLED);
//
//	putText(outImg, strName, cv::Point(pt1.x, pt1.y), fontFace, fontScale,
//        Scalar::all(0), thickness, CV_AA);
//}
//
//void _drawMatchLabel( cv::Mat& outImg, cv::Mat& outImg1, cv::Mat& outImg2 ,
//							  const cvline_polar& l1, const cvline_polar& l2, const Scalar& matchColor, int flags, string strName/*=""*/)
//{
//	Point2f pt1((l1.pt1.x + l1.pt2.x)*0.5, (l1.pt1.y + l1.pt2.y)*0.5);
//	Point2f pt2( std::min((l2.pt1.x + l2.pt2.x)*0.5+outImg1.cols, double(outImg.cols-1)), (l2.pt1.y + l2.pt2.y)*0.5 );
//
//	int baseline=0;
//	int fontFace = CV_FONT_HERSHEY_SIMPLEX;
//	double fontScale = 0.4;
//	int thickness = 1;
//	Size textSize = getTextSize(strName, fontFace,
//		fontScale, thickness, &baseline);
//	baseline += thickness;
//	baseline = 2;
//	//line(outImg, cv::Point(pt1.x, pt1.y) + Point(0, thickness),
// //    cv::Point(pt1.x, pt1.y) + Point(textSize.width, thickness),
// //    Scalar(0, 0, 255));
//	rectangle(outImg, cv::Point(pt1.x, pt1.y) + Point(0, baseline),
//		cv::Point(pt1.x, pt1.y) + Point(textSize.width, -textSize.height) - Point(0, baseline),
//		Scalar(255,255,255), CV_FILLED);
//	//line(outImg, cv::Point(pt2.x, pt2.y) + Point(0, thickness),
// //    cv::Point(pt2.x, pt2.y) + Point(textSize.width, thickness),
// //    Scalar(0, 0, 255));
//	rectangle(outImg, cv::Point(pt2.x, pt2.y) + Point(0, baseline),
//		cv::Point(pt2.x, pt2.y) + Point(textSize.width, -textSize.height) - Point(0, baseline),
//		Scalar(255,255,255), CV_FILLED);
//
//	//cv::putText(outImg, strName, cv::Point(pt1.x, pt1.y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, CV_AA );
//	//cv::putText(outImg, strName, cv::Point(dpt2.x, dpt2.y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, CV_AA );
//	putText(outImg, strName, cv::Point(pt1.x, pt1.y), fontFace, fontScale,
//        Scalar::all(0), thickness, CV_AA);
//	putText(outImg, strName, cv::Point(pt2.x, pt2.y), fontFace, fontScale,
//        Scalar::all(0), thickness, CV_AA);
//}

void drawLineMatches( const cv::Mat& img1, const vector<cvline_polar>& keylines1,
					 const cv::Mat& img2, const vector<cvline_polar>& keylines2,
					 cv::Mat& outImg,
					 const Scalar& matchColor, const Scalar& singlePointColor,
					 const vector< vector< char > >& matchesMask, int flags, int thickness/*=1*/)
{
	cv::Mat outImg1, outImg2;
	_prepareImgAndDrawKeylines( img1, keylines1, img2, keylines2,
		outImg, outImg1, outImg2, singlePointColor, flags );

	// draw matches
	//for( size_t i = 0; i < matches1to2.size(); i++ )
	//{
	//	for( size_t j = 0; j < matches1to2[i].size(); j++ )
	//	{
	//		int i1 = matches1to2[i][j].queryIdx;
	//		int i2 = matches1to2[i][j].trainIdx;
	//		if( matchesMask.empty() || matchesMask[i][j] )
	//		{
	//			const cvline_polar &l1 = keylines1[i1], &l2 = keylines2[i2];
	//			_drawLineMatch( outImg, outImg1, outImg2, l1, l2, matchColor, flags );
	//		}
	//	}
	//}

	int nCount = keylines1.size();
	//nCount = 20;
	for( size_t i = 0; i < nCount; i++ )
	{
		const cvline_polar &l1 = keylines1[i], &l2 = keylines2[i];
		char temp[64];
		sprintf_s(temp, "%d", i+1);
		string strName(temp);
		_drawLineMatch( outImg, outImg1, outImg2, l1, l2, matchColor, flags);
	}
	for( size_t i = 0; i < nCount; i++ )
	{
		const cvline_polar &l1 = keylines1[i], &l2 = keylines2[i];
		char temp[64];
		sprintf_s(temp, "%d", i+1);
		string strName(temp);
		_drawMatchLabel( outImg, outImg1, outImg2, l1, l2, matchColor, flags, strName);
	}
}

void drawLineMatches2( cv::Mat& img1, const vector<cvline_polar>& keylines1,
					 cv::Mat& img2, const vector<cvline_polar>& keylines2,
					 const Scalar& matchColor, const Scalar& singlePointColor,
					 const vector< vector< char > >& matchesMask, int flags, int thickness/*=1*/)
{
	size_t nCount = keylines1.size();
	for( size_t i = 0; i < nCount; i++ )
	{
		_drawKeyline( img1, keylines1[i], matchColor, cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
		_drawKeyline( img2, keylines2[i], matchColor, cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	}
	
	for( size_t i = 0; i < nCount; i++ )
	{
		char temp[64];
		sprintf_s(temp, "%d", i+1);
		string strName(temp);
		_drawLabel( img1, keylines1[i],matchColor, cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, strName);
		_drawLabel( img2, keylines2[i],matchColor, cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, strName);
	}
}