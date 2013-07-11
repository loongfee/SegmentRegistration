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
#include "lp.h"

//Line* LineDetect(const char* imfile, int &nLines)
//{
//	IplImage* im = cvLoadImage(imfile,CV_LOAD_IMAGE_GRAYSCALE);
//	image_double image = new_image_double(im->width, im->height);
//	unsigned char* im_src = (unsigned char*)im->imageData;
//	int xsize = image->xsize;
//	int ysize = image->ysize;
//	int y,x;
//	for (y = 0;y < ysize;y++)
//	{
//		for (x = 0;x < xsize;x++)
//		{
//			image->data[y * xsize + x] = im_src[y * im->widthStep + x];
//		}
//	}
//	ntuple_list detected_lines = lsd(image);
//	free_image_double(image);
//
//	nLines = detected_lines->size;
//
//	IplImage* smoothed_im = cvCreateImage(cvSize(im->width,im->height),IPL_DEPTH_8U,1);
//	cvSmooth(im,smoothed_im);
//	cvReleaseImage(&im);
//
//	Line* pLines = new Line[nLines];
//	int nCount = 0;
//	int i,j;
//	int dim = detected_lines->dim;
//	for (i = 0;i < nLines;i++)
//	{
//		double x1 = detected_lines->values[i*dim+0];
//		double y1 = detected_lines->values[i*dim+1];
//		double x2 = detected_lines->values[i*dim+2];
//		double y2 = detected_lines->values[i*dim+3];
//		double len = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
//		len = sqrt(len);
//		if(len <= 20) continue;
//		double a = y1 - y2;
//		double b = -(x1 - x2);
//		double c = x1 * y2 - x2 * y1;
//		vector<CvPoint2D32f> Pts;
//		int nPts = 0;
//		CvPoint2D32f p;
//		p.x = x1;
//		p.y = y1;
//		Pts.push_back(p);
//		nPts++;
//		if ( abs((a+0.0000000000001) / (b+0.0000000000001)) < 1 )
//		{
//			double iterval = (x2 - x1) / len;
//			while (1)
//			{
//				p.x += iterval;
//				if( (iterval > 0 && p.x > x2) || (iterval < 0 && p.x < x2) ) break;
//				p.y = -(a * p.x + c) / b;
//				Pts.push_back(p);
//				nPts++;
//			}
//		}
//		else
//		{
//			double iterval = (y2 - y1) / len;
//			while (1)
//			{
//				p.y += iterval;
//				if( (iterval > 0 && p.y > y2) || (iterval < 0 && p.y < y2) ) break;
//				p.x = -(b * p.y + c) / a;
//				Pts.push_back(p);
//				nPts++;
//			}
//		}
//		double line_grad_x, line_grad_y;
//		int nValid = 0;
//		line_grad_x = line_grad_y = 0;
//		for (j = 0;j < nPts;j++)
//		{
//			double gray2;
//			double gray1;
//			if( !GetImGray(smoothed_im, Pts[j].x + 1, Pts[j].y, gray2) ) continue;
//			if( !GetImGray(smoothed_im, Pts[j].x - 1, Pts[j].y, gray1) ) continue;
//			line_grad_x += (gray2 - gray1) / 2;
//			if( !GetImGray(smoothed_im, Pts[j].x, Pts[j].y + 1, gray2) ) continue;
//			if( !GetImGray(smoothed_im, Pts[j].x, Pts[j].y - 1, gray1) ) continue;
//			line_grad_y += (gray2 - gray1) / 2;
//			nValid++;
//		}
//		if( nValid == 0) continue;
//		line_grad_x /= nValid;
//		line_grad_y /= nValid;
//		double ExProd = (Pts[nPts-1].x - Pts[0].x) * line_grad_y - line_grad_x * (Pts[nPts-1].y - Pts[0].y);
//		if (ExProd > 0)
//		{
//			pLines[nCount].StartPt = Pts[0];
//			pLines[nCount].EndPt = Pts[nPts-1];
//		}
//		else
//		{
//			pLines[nCount].StartPt = Pts[nPts-1];
//			pLines[nCount].EndPt = Pts[0];
//		}
//		pLines[nCount].length = len;
//		pLines[nCount].Center.x = (x1 + x2) / 2;
//		pLines[nCount].Center.y = (y1 + y2) / 2;
//		pLines[nCount].para_a = a;
//		pLines[nCount].para_b = b;
//		pLines[nCount].para_c = c;
//		nCount++;
//	}
//	cvReleaseImage(&smoothed_im);
//
//	free_ntuple_list(detected_lines);
//
//	nLines = nCount;
//	return pLines;
//}

void crossCheckMatching( Ptr<DescriptorMatcher>& descriptorMatcher,
						const Mat& descriptors1, const Mat& descriptors2,
						vector<DMatch>& filteredMatches12, int knn/*=1*/ )
{
	filteredMatches12.clear();
	vector<vector<DMatch> > matches12, matches21;
	descriptorMatcher->knnMatch( descriptors1, descriptors2, matches12, knn );
	descriptorMatcher->knnMatch( descriptors2, descriptors1, matches21, knn );
	for( size_t m = 0; m < matches12.size(); m++ )
	{
		bool findCrossCheck = false;
		for( size_t fk = 0; fk < matches12[m].size(); fk++ )
		{
			DMatch forward = matches12[m][fk];

			for( size_t bk = 0; bk < matches21[forward.trainIdx].size(); bk++ )
			{
				DMatch backward = matches21[forward.trainIdx][bk];
				if( backward.trainIdx == forward.queryIdx )
				{
					filteredMatches12.push_back(forward);
					findCrossCheck = true;
					break;
				}
			}
			if( findCrossCheck ) break;
		}
	}
}


int PtMatch(const char* imfile1, const char* imfile2, vector<MatchPt> &match)
{
	Ptr<FeatureDetector> detector = FeatureDetector::create( "SIFT" );
	Ptr<DescriptorExtractor> descriptorExtractor = DescriptorExtractor::create( "SIFT" );
	Ptr<DescriptorMatcher> descriptorMatcher = DescriptorMatcher::create( "BruteForce" );
	if( detector.empty() || descriptorExtractor.empty() || descriptorMatcher.empty()  )
	{
		cout << "Can not create detector or descriptor exstractor or descriptor matcher of given types" << endl;
		return -1;
	}

	Mat img1 = imread( imfile1 );
	Mat img2 = imread( imfile2 );
	if( img1.empty() || img2.empty() )
	{
		cout << "Can not read images" << endl;
		return -1;
	}

	cout << endl << "Extracting keypoints from first image..." << endl;
	vector<KeyPoint> keypoints1;
	detector->detect( img1, keypoints1 );
	cout << keypoints1.size() << " points" << endl;

	cout << "Computing descriptors for keypoints from first image..." << endl;
	Mat descriptors1;
	descriptorExtractor->compute( img1, keypoints1, descriptors1 );

	cout << endl << "Extracting keypoints from second image..." << endl;
	vector<KeyPoint> keypoints2;
	detector->detect( img2, keypoints2 );
	cout << keypoints2.size() << " points" << endl;

	cout << "Computing descriptors for keypoints from second image..." << endl;
	Mat descriptors2;
	descriptorExtractor->compute( img2, keypoints2, descriptors2 );

	cout << "Matching descriptors..." << endl;
	vector<DMatch> filteredMatches;
	//FlannBasedMatcher matcher;
	//matcher.match( descriptors1, descriptors2, filteredMatches );
	crossCheckMatching( descriptorMatcher, descriptors1, descriptors2, filteredMatches );

	vector<int> queryIdxs( filteredMatches.size() ), trainIdxs( filteredMatches.size() );
	for( size_t i = 0; i < filteredMatches.size(); i++ )
	{
		queryIdxs[i] = filteredMatches[i].queryIdx;
		trainIdxs[i] = filteredMatches[i].trainIdx;
	}

	vector<Point2f> points1; KeyPoint::convert(keypoints1, points1, queryIdxs);
	vector<Point2f> points2; KeyPoint::convert(keypoints2, points2, trainIdxs);

	int nPtMatches = points1.size();
	match.clear();
	for(int i = 0;i < nPtMatches;i++){
		MatchPt tmp;
		tmp.P1 = points1[i];
		tmp.P2 = points2[i];
		match.push_back(tmp);
	}

	filteredMatches.clear();
	return nPtMatches;
}