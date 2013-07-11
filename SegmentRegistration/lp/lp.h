/* 
Author: Bin Fan (NLPR), bfan@nlpr.ia.ac.cn

This is a demo about how to use the 'LineMatcher' library, which conducts line matching by affine invariants. 
Please refer to the following paper for more details about the algorithm.

Bin Fan, Fuchao Wu and Zhanyi Hu. Line Matching Leveraged by Point Correspondences, In CVPR 2010, pp 390-391.

The line segmentation detector (LSD) used in the demo is from the website:
http://www.ipol.im/pub/algo/gjmr_line_segment_detector/

Rafael Grompone von Gioi, J¨¦r¨¦mie Jakubowicz, Jean-Michel Morel, Gregory Randall, LSD: A Fast Line Segment Detector with a False Detection Control, IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 32, no. 4, pp. 722-732, April 2010.

The matched points in the demo are obtained based on the opencv 2.2

If you have questions, please contact: bfan@nlpr.ia.ac.cn
*/

#ifndef LP_LINE_MATCH
#define LP_LINE_MATCH
//extern "C"
//{
//#include "..\lsd.h"
//};
#include <stdio.h>
#include <highgui.h>
#include <cv.h>
#include "LineMatcher.h"
#include <math.h>
#include <iostream>
using namespace std;
using namespace cv;


#ifndef BYTE
typedef unsigned char BYTE;
#endif // !BYTE


//Line* LineDetect(const char* imfile, int &nLines);

void crossCheckMatching( Ptr<cv::DescriptorMatcher>& descriptorMatcher,
						const cv::Mat& descriptors1, const cv::Mat& descriptors2,
						vector<DMatch>& filteredMatches12, int knn=1 );


int PtMatch(const char* imfile1, const char* imfile2, vector<MatchPt> &match);


#endif