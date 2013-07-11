#ifndef _DRAW_H_
#define _DRAW_H_

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <vector>
#include "lineConstant.h"
using namespace std;

const int draw_shift_bits = 4;
const int draw_multiplier = 1 << draw_shift_bits;

using namespace cv;
/*
 * Functions to draw keypoints and matches.
 */
static void _drawKeyline( cv::Mat& img, const cvline_polar& l, const cv::Scalar& color, int flags, int thickness=1)
{
	cv::line(img, cvPoint(cvRound(l.pt1.x), cvRound(l.pt1.y)),
		cvPoint(cvRound(l.pt2.x), cvRound(l.pt2.y)), color, thickness, CV_AA);
}

static void drawKeylines( const cv::Mat& image, const vector<cvline_polar>& keylines, cv::Mat& outImg,
                    const cv::Scalar& _color, int flags, int thickness=1)
{
    if( !(flags & DrawMatchesFlags::DRAW_OVER_OUTIMG) )
        cvtColor( image, outImg, CV_GRAY2BGR );

    RNG& rng=theRNG();
    bool isRandColor = _color == Scalar::all(-1);

    for( vector<cvline_polar>::const_iterator i = keylines.begin(), ie = keylines.end(); i != ie; ++i )
    {
        Scalar color = isRandColor ? Scalar(rng(256), rng(256), rng(256)) : _color;
        _drawKeyline( outImg, *i, color, flags, thickness);
    }
}

static void _prepareImgAndDrawKeylines( const cv::Mat& img1, const vector<cvline_polar>& keylines1,
									   const cv::Mat& img2, const vector<cvline_polar>& keylines2,
									   cv::Mat& outImg, cv::Mat& outImg1, cv::Mat& outImg2,
									   const cv::Scalar& singlePointColor, int flags )
{
	Size size( img1.cols + img2.cols, MAX(img1.rows, img2.rows) );
	if( flags & DrawMatchesFlags::DRAW_OVER_OUTIMG )
	{
		if( size.width > outImg.cols || size.height > outImg.rows )
			CV_Error( CV_StsBadSize, "outImg has size less than need to draw img1 and img2 together" );
		outImg1 = outImg( Rect(0, 0, img1.cols, img1.rows) );
		outImg2 = outImg( Rect(img1.cols, 0, img2.cols, img2.rows) );
	}
	else
	{
		outImg.create( size, CV_MAKETYPE(img1.depth(), 3) );
		outImg1 = outImg( Rect(0, 0, img1.cols, img1.rows) );
		outImg2 = outImg( Rect(img1.cols, 0, img2.cols, img2.rows) );

		if( img1.type() == CV_8U )
			cvtColor( img1, outImg1, CV_GRAY2BGR );
		else
			img1.copyTo( outImg1 );

		if( img2.type() == CV_8U )
			cvtColor( img2, outImg2, CV_GRAY2BGR );
		else
			img2.copyTo( outImg2 );
	}

	// draw keypoints
	if( !(flags & DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS) )
	{
		cv::Mat outImg1 = outImg( Rect(0, 0, img1.cols, img1.rows) );
		drawKeylines( outImg1, keylines1, outImg1, singlePointColor, flags + DrawMatchesFlags::DRAW_OVER_OUTIMG );

		cv::Mat outImg2 = outImg( Rect(img1.cols, 0, img2.cols, img2.rows) );
		drawKeylines( outImg2, keylines2, outImg2, singlePointColor, flags + DrawMatchesFlags::DRAW_OVER_OUTIMG );
	}
}

static void _drawLineMatch( cv::Mat& outImg, cv::Mat& outImg1, cv::Mat& outImg2 ,
							  const cvline_polar& l1, const cvline_polar& l2, const cv::Scalar& matchColor, int flags, int thickness=1)
{
	Scalar color(0, 0, 255);
	//RNG& rng = theRNG();
	//bool isRandMatchColor = matchColor == Scalar::all(-1);
	//color = isRandMatchColor ? Scalar( rng(256), rng(256), rng(256) ) : matchColor;

	_drawKeyline( outImg1, l1, color, flags );
	_drawKeyline( outImg2, l2, color, flags );

	cv::Point pt1((l1.pt1.x + l1.pt2.x)*0.5, (l1.pt1.y + l1.pt2.y)*0.5);
	cv::Point pt2( std::min((l2.pt1.x + l2.pt2.x)*0.5+outImg1.cols, double(outImg.cols-1)), (l2.pt1.y + l2.pt2.y)*0.5 );
	//cv::line(outImg, pt1, pt2, color, thickness, CV_AA);
}


static inline void _drawLabel( cv::Mat& outImg, const cvline_polar& line,const cv::Scalar& matchColor, int flags, string strName="", int thickness=1)
{
	Point2f pt1((line.pt1.x + line.pt2.x)*0.5, (line.pt1.y + line.pt2.y)*0.5);

	int baseline=0;
	int fontFace = CV_FONT_HERSHEY_SIMPLEX;
	double fontScale = 0.35;
	fontScale *= thickness;
	//int thickness = 1;
	Size textSize = getTextSize(strName, fontFace,
		fontScale, thickness, &baseline);
	baseline += thickness;
	baseline = 2;

	rectangle(outImg, cv::Point(pt1.x, pt1.y) + Point(0, baseline),
		cv::Point(pt1.x, pt1.y) + Point(textSize.width, -textSize.height) - Point(0, baseline),
		Scalar(255,255,255), CV_FILLED);

	putText(outImg, strName, cv::Point(pt1.x, pt1.y), fontFace, fontScale,
        Scalar::all(0), thickness, CV_AA);
}

static void _drawMatchLabel( cv::Mat& outImg, cv::Mat& outImg1, cv::Mat& outImg2 ,
							  const cvline_polar& l1, const cvline_polar& l2, const cv::Scalar& matchColor, int flags, string strName="")
{
	Point2f pt1((l1.pt1.x + l1.pt2.x)*0.5, (l1.pt1.y + l1.pt2.y)*0.5);
	Point2f pt2( std::min((l2.pt1.x + l2.pt2.x)*0.5+outImg1.cols, double(outImg.cols-1)), (l2.pt1.y + l2.pt2.y)*0.5 );

	int baseline=0;
	int fontFace = CV_FONT_HERSHEY_SIMPLEX;
	double fontScale = 0.4;
	int thickness = 1;
	Size textSize = getTextSize(strName, fontFace,
		fontScale, thickness, &baseline);
	baseline += thickness;
	baseline = 2;
	//line(outImg, cv::Point(pt1.x, pt1.y) + Point(0, thickness),
 //    cv::Point(pt1.x, pt1.y) + Point(textSize.width, thickness),
 //    Scalar(0, 0, 255));
	rectangle(outImg, cv::Point(pt1.x, pt1.y) + Point(0, baseline),
		cv::Point(pt1.x, pt1.y) + Point(textSize.width, -textSize.height) - Point(0, baseline),
		Scalar(255,255,255), CV_FILLED);
	//line(outImg, cv::Point(pt2.x, pt2.y) + Point(0, thickness),
 //    cv::Point(pt2.x, pt2.y) + Point(textSize.width, thickness),
 //    Scalar(0, 0, 255));
	rectangle(outImg, cv::Point(pt2.x, pt2.y) + Point(0, baseline),
		cv::Point(pt2.x, pt2.y) + Point(textSize.width, -textSize.height) - Point(0, baseline),
		Scalar(255,255,255), CV_FILLED);

	//cv::putText(outImg, strName, cv::Point(pt1.x, pt1.y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, CV_AA );
	//cv::putText(outImg, strName, cv::Point(dpt2.x, dpt2.y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, CV_AA );
	putText(outImg, strName, cv::Point(pt1.x, pt1.y), fontFace, fontScale,
        Scalar::all(0), thickness, CV_AA);
	putText(outImg, strName, cv::Point(pt2.x, pt2.y), fontFace, fontScale,
        Scalar::all(0), thickness, CV_AA);
}

void drawLineMatches( const cv::Mat& img1, const vector<cvline_polar>& keylines1,
					 const cv::Mat& img2, const vector<cvline_polar>& keylines2,
					 cv::Mat& outImg,
					 const cv::Scalar& matchColor, const cv::Scalar& singlePointColor,
					 const vector< vector< char > >& matchesMask, int flags, int thickness=1);
void drawLineMatches2( cv::Mat& img1, const vector<cvline_polar>& keylines1,
					 cv::Mat& img2, const vector<cvline_polar>& keylines2,
					 const Scalar& matchColor, const Scalar& singlePointColor,
					 const vector< vector< char > >& matchesMask, int flags, int thickness=1);

#endif