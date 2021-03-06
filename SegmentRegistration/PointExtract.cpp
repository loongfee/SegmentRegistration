#include "PointExtrat.h"

bool orb_extract(const char* srcImage, vector<fPoint> &pointList)
{
	cv::Mat img = cv::imread(srcImage);
	if (!img.data)
	{
		cout << "error reading images " << endl;
		return false;
	}

	ORB orb;
	vector<KeyPoint> keyPoints;
	cv::Mat descriptors;

	orb(img, cv::Mat(), keyPoints, descriptors);


	pointList.clear();
	for( int i = 0; i < (int)keyPoints.size(); i++ )
	{
		pointList.push_back(fPoint(keyPoints[i].pt.x, keyPoints[i].pt.y));
	}

	//Mat img_matches;
	//drawMatches(img_1, keyPoints_1, img_2, keyPoints_2,
	//	good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
	//	vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
	//imshow( "Match", img_matches);
	//cvWaitKey();
	return true;
}

bool orb_match(const char* srcImage, const char* refImage, vector<PointPair>& pointPairList)
{
	cv::Mat img_1 = imread(srcImage);
	cv::Mat img_2 = imread(refImage);
	if (!img_1.data || !img_2.data)
	{
		cout << "error reading images " << endl;
		return false;
	}

	cv::ORB orb;
	vector<cv::KeyPoint> keyPoints_1, keyPoints_2;
	cv::Mat descriptors_1, descriptors_2;

	orb(img_1, cv::Mat(), keyPoints_1, descriptors_1);
	orb(img_2, cv::Mat(), keyPoints_2, descriptors_2);

	cv::BruteForceMatcher<cv::HammingLUT> matcher;
	vector<cv::DMatch> matches;
	matcher.match(descriptors_1, descriptors_2, matches);

	double max_dist = 0; double min_dist = 100;
	//-- Quick calculation of max and min distances between keypoints
	for( int i = 0; i < descriptors_1.rows; i++ )
	{ 
		double dist = matches[i].distance;
		if( dist < min_dist ) min_dist = dist;
		if( dist > max_dist ) max_dist = dist;
	}
	printf("-- Max dist : %f \n", max_dist );
	printf("-- Min dist : %f \n", min_dist );
	//-- Draw only "good" matches (i.e. whose distance is less than 0.6*max_dist )
	//-- PS.- radiusMatch can also be used here.
	std::vector< cv::DMatch > good_matches;
	for( int i = 0; i < descriptors_1.rows; i++ )
	{ 
		if( matches[i].distance < 0.6*max_dist )
		{ 
			good_matches.push_back( matches[i]); 
		}
	}
	//
	pointPairList.clear();
	for(int i = 0;i < (int)good_matches.size();++i)
	{
		PointPair ptPair;
		ptPair.point = fPoint(keyPoints_1[good_matches[i].queryIdx].pt.x, keyPoints_1[good_matches[i].queryIdx].pt.y);
		ptPair.point_prime = fPoint(keyPoints_2[good_matches[i].trainIdx].pt.x, keyPoints_2[good_matches[i].trainIdx].pt.y);
		pointPairList.push_back(ptPair);
	}
	//

	cv::Mat img_matches;
	drawMatches(img_1, keyPoints_1, img_2, keyPoints_2,
		good_matches, img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
		vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
	imshow( "Match", img_matches);
	cvWaitKey();
	return true;
}

bool surf_extract(string srcImage, vector<fPoint> &pointList, CvSURFParams params/* = cvSURFParams(500, 1)*/)
{
	cv::initModule_nonfree();

	IplImage* source = cvLoadImage( srcImage.c_str(), CV_LOAD_IMAGE_GRAYSCALE );
	if( !source)
	{
		return false;
	}

	CvMemStorage* storage = cvCreateMemStorage(0);

	CvSeq* sourceKeypoints = 0, *sourceDescriptors = 0;

	cvExtractSURF( source, 0, &sourceKeypoints, &sourceDescriptors, storage, params );

	pointList.clear();
	for( int i = 0; i < sourceKeypoints->total; i++ )
	{
		CvSURFPoint* point = (CvSURFPoint*)cvGetSeqElem( sourceKeypoints, i );
		pointList.push_back(fPoint(point->pt.x, point->pt.y));
	}

	return true;
}


bool surf_match(const char* srcImage, const char* refImage, vector<PointPair>& pointPairList, CvSURFParams params/* = cvSURFParams(500, 1)*/)
{
	cv::initModule_nonfree();
	IplImage* source = cvLoadImage( srcImage, CV_LOAD_IMAGE_GRAYSCALE );
	IplImage* refer = cvLoadImage( refImage, CV_LOAD_IMAGE_GRAYSCALE );
	if( !source || !refer)
	{
		return false;
	}

	CvMemStorage* storage = cvCreateMemStorage(0);

	CvSeq* sourceKeypoints = 0, *sourceDescriptors = 0;
	CvSeq* referKeypoints = 0, *referDescriptors = 0;

	cvExtractSURF( source, 0, &sourceKeypoints, &sourceDescriptors, storage, params );
	cvExtractSURF( refer, 0, &referKeypoints, &referDescriptors, storage, params );

	IplImage* source_color = cvCreateImage(cvGetSize(source), 8, 3);
	cvCvtColor( source, source_color, CV_GRAY2BGR );

	static CvScalar colors[] =
	{
		{{0,0,255}},
		{{0,128,255}},
		{{0,255,255}},
		{{0,255,0}},
		{{255,128,0}},
		{{255,255,0}},
		{{255,0,0}},
		{{255,0,255}},
		{{255,255,255}}
	};

#ifdef USE_FLANN
	printf("Using approximate nearest neighbor search\n");
#endif

	vector<int> ptpairs;
#ifdef USE_FLANN
	flannFindPairs( sourceKeypoints, sourceDescriptors, referKeypoints, referDescriptors, ptpairs );
#else
	findPairs( sourceKeypoints, sourceDescriptors, referKeypoints, referDescriptors, ptpairs );
#endif


	pointPairList.clear();
	for(int i = 0;i < (int)ptpairs.size()/2;++i)
	{
		CvSURFPoint* r1 = (CvSURFPoint*)cvGetSeqElem( sourceKeypoints, ptpairs[2*i] );
		CvSURFPoint* r2 = (CvSURFPoint*)cvGetSeqElem( referKeypoints, ptpairs[2*i+1] );

		PointPair ptPair;
		ptPair.point = fPoint(r1->pt.x, r1->pt.y);
		ptPair.point_prime = fPoint(r2->pt.x, r2->pt.y);
		pointPairList.push_back(ptPair);
	}


	//cv::Mat img_matches;
	//drawMatches(source, sourceKeypoints, refer, referKeypoints,
	//	good_matches, img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
	//	vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
	//imshow( "Match", img_matches);
	cvWaitKey();
	return true;
}

bool harris_extract(const char* srcImage, vector<fPoint> &pointList, int threshold/*=10*/, bool nonmaxSuppression/*=true*/)
{
	cv::Mat image, image1 = cv::imread (srcImage);
	
	if( image1.empty())
	{
		return false;
	}
	cv::cvtColor (image1, image, CV_BGR2GRAY);
	//快速角点检测
	std::vector<cv::KeyPoint> keypoints;
	cv::FastFeatureDetector fast(threshold, nonmaxSuppression);
	fast.detect(image, keypoints);
	cv::drawKeypoints (image, keypoints,image,cv::Scalar::all(255),cv::DrawMatchesFlags::DRAW_OVER_OUTIMG);

	pointList.clear();
	
	for( int i = 0; i < (int)keypoints.size(); i++ )
	{
		pointList.push_back(fPoint(keypoints[i].pt.x, keypoints[i].pt.y));
	}
	return true;
}