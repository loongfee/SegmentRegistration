#include "func.h"
#include "draw.h"

#include "RANSAC/ransac/ransac.h"
#include "RANSAC/estimators/Solver.h"
#include "RANSAC/estimators/affineSolver.h"
#include "RANSAC/estimators/affineError.h"
#include "RANSAC/estimators/lineAffineSolver.h"

#include <armadillo>
using namespace groupsac;
using namespace arma;

bool fexists(const char *filename)
{
	ifstream ifile(filename);
	if (ifile) {
		return true;
	}
	return false;
}

bool fexists(string filename)
{
	return fexists(filename.c_str());
}

void normalize_segment(cvline_polar& segment, double len/* = 400.0*/)
{
	fPoint centerPt((segment.pt1.x+segment.pt2.x)*0.5, (segment.pt1.y+segment.pt2.y)*0.5);

	double delta_x = cos(segment.theta) * len * 0.5;
	double delta_y = sin(segment.theta) * len * 0.5;

	segment.pt1.x = centerPt.x - delta_x;
	segment.pt1.y = centerPt.y - delta_y;
	segment.pt2.x = centerPt.x + delta_x;
	segment.pt2.y = centerPt.y + delta_y;

	//  for check
	fPoint line[2];
	line[0] = segment.pt1;
	line[1] = segment.pt2;
	cvline_polar l = line2polar(line);
	double delta_theta = l.theta - segment.theta;
	double delta_rho = l.rho - segment.rho;
}

double segment_length(const cvline_polar& l)
{
	double delta_x = l.pt1.x - l.pt2.x;
	double delta_y = l.pt1.y - l.pt2.y;
	return sqrt(delta_x*delta_x + delta_y*delta_y);
}

cvline_polar segment_move_vertical(const cvline_polar& l, double dist)
{
	fPoint pt1 = l.pt1;
	fPoint pt2 = l.pt2;
	double sin_theta = sin(l.theta);
	double cos_theta = cos(l.theta);
	double rho = l.rho;
	double x1 = pt1.x*cos_theta*cos_theta + pt1.y*sin_theta*cos_theta - sin_theta*(rho - dist);
	double y1 = pt1.x*sin_theta*cos_theta + pt1.y*sin_theta*sin_theta + cos_theta*(rho - dist);
	double x2 = pt2.x*cos_theta*cos_theta + pt2.y*sin_theta*cos_theta - sin_theta*(rho - dist);
	double y2 = pt2.x*sin_theta*cos_theta + pt2.y*sin_theta*sin_theta + cos_theta*(rho - dist);

	return line2polar(fPoint(x1, y1), fPoint(x2, y2));
}

fPoint line2point(const cvline_polar& l)
{
	fPoint pt;
	pt.x = (l.pt2.x - l.pt1.x);
	pt.y = (l.pt2.y - l.pt1.y);
	double len = sqrt(pt.x*pt.x + pt.y*pt.y) / 200.0;
	pt.x /= len;
	pt.y /= len;
	return pt;
}


cv::Scalar random_color(CvRNG* rng)
{
	//int color = cvRandInt(rng);
	return CV_RGB(rand()&255, rand()&255, rand()&255 );
	//return CV_RGB(color&255, (color>>8)&255, (color>>16)&255);
}

/************************************************************************/
/* 
	rho = -x * sin(theta) + y * cos(theta)
*/
/************************************************************************/
cvline_polar line2polar(CvPoint* line)
{
	cvline_polar line_polar;

	line_polar.pt1 = fPoint(line[0].x, line[0].y);
	line_polar.pt2 = fPoint(line[1].x, line[1].y);

	double rho = 0;
	//角度为对应直线与X轴的夹角
	double theta = 0;
	if( line[0].x == line[1].x)
	{
		// 如果与Y轴平行
		rho = (double)line[0].x;
		theta = PI / 2.0;
	} 
	else
	{
		theta = atan((line[1].y - line[0].y) / (double)(line[1].x - line[0].x));
		//float t0 = (line[1].x * line[0].y - line[0].x * line[1].y) * cos(theta);
		//float t1 = -line[0].x * sin(theta) + line[0].y * cos(theta);
		//float t2 = -line[1].x * sin(theta) + line[1].y * cos(theta);
		rho = -line[0].x * sin(theta) + line[0].y * cos(theta);
	}
	line_polar.rho = rho; 
	line_polar.theta = theta;

	return line_polar;
}


cvline_polar line2polar(fPoint* line)
{
	cvline_polar line_polar;

	line_polar.pt1 = line[0];
	line_polar.pt2 = line[1];
	double rho = 0;
	//角度为对应直线与X轴的夹角
	double theta = 0;
	if( line[0].x == line[1].x)
	{
		// 如果与Y轴平行
		rho = (double)line[0].x;
		theta = PI / 2.0f;
	} 
	else
	{
		theta = atan((line[1].y - line[0].y) / (double)(line[1].x - line[0].x));
		//float t0 = (line[1].x * line[0].y - line[0].x * line[1].y) * cos(theta);
		//float t1 = -line[0].x * sin(theta) + line[0].y * cos(theta);
		//float t2 = -line[1].x * sin(theta) + line[1].y * cos(theta);
		rho = -line[0].x * sin(theta) + line[0].y * cos(theta);
	}
	line_polar.rho = rho; 
	line_polar.theta = theta;

	return line_polar;
}

cvline_polar line2polar(fPoint pt1, fPoint pt2)
{
	fPoint line[2];

	line[0] = pt1;
	line[1] = pt2;
	return line2polar(line);
}

cvline_polar mergeline(cvline_polar l1, cvline_polar l2)
{
	//vertical
	if(fabs(l1.pt1.x-l1.pt2.x) < 2.0 && fabs(l2.pt1.x-l2.pt2.x) < 2.0 )
	{
		fPoint up_pt = l1.pt1;
		fPoint low_pt = l1.pt1;
		
		if(up_pt.y > l1.pt2.y)
		{
			up_pt = l1.pt2;
		}
		if(low_pt.y < l1.pt2.y)
		{
			low_pt = l1.pt2;
		}


		if(up_pt.y > l2.pt1.y)
		{
			up_pt = l2.pt1;
		}
		if(low_pt.y < l2.pt1.y)
		{
			low_pt = l2.pt1;
		}

		if(up_pt.y > l2.pt2.y)
		{
			up_pt = l2.pt2;
		}
		if(low_pt.y < l2.pt2.y)
		{
			low_pt = l2.pt2;
		}

		return line2polar(up_pt, low_pt);
	}
	else
	{
		// other
		fPoint left_pt = l1.pt1;
		fPoint right_pt = l1.pt1;


		if(left_pt.x > l1.pt2.x)
		{
			left_pt = l1.pt2;
		}
		if(right_pt.x < l1.pt2.x)
		{
			right_pt = l1.pt2;
		}


		if(left_pt.x > l2.pt1.x)
		{
			left_pt = l2.pt1;
		}
		if(right_pt.x < l2.pt1.x)
		{
			right_pt = l2.pt1;
		}


		if(left_pt.x > l2.pt2.x)
		{
			left_pt = l2.pt2;
		}
		if(right_pt.x < l2.pt2.x)
		{
			right_pt = l2.pt2;
		}

		return line2polar(left_pt, right_pt);
	}
}

double distance_(const fPoint& pt1, const fPoint& pt2)
{
	return sqrt((pt1.x-pt2.x)*(pt1.x-pt2.x) + (pt1.y-pt2.y)*(pt1.y-pt2.y));
}

bool pt_comparison(const fPoint& node1, const fPoint& node2)
{
	if (node1.x < node2.x) return true;
	if (node1.x == node2.y) return node1.y < node2.y;

	return false;
}

bool sortPoints(vector<fPoint>& pts)
{
	std::sort(pts.begin(), pts.end(), pt_comparison);
	return true;
}

void drawMatchedSegments(string foldName, string sourceName, string referenceName, string linesName)
{
	string SourceImage = foldName + "\\" + sourceName;//待配准影像
	string ReferImage = foldName + "\\" + referenceName;//参考影像
	string linesText = foldName + "\\" + linesName;//参考影像
	string resultImage = foldName + "\\segment_matched.png";//结果影像

	cv::Mat img_source = cv::imread(SourceImage);
	cv::Mat img_reference = cv::imread(ReferImage);
	cv::Mat img_matches;
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
	int thickness = 3;
	drawLineMatches2( img_source, source_ls, img_reference, reference_ls,
					 cv::Scalar::all(-1), cv::Scalar::all(-1),
		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//drawLineMatches(img_source, source_ls, img_reference, reference_ls,
	//	img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
	//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);

	cv::imwrite(resultImage, img_matches);
}

void drawMatchedSegments2(string foldName, string sourceName, string referenceName, string linesName)
{
	string SourceImage = foldName + "\\" + sourceName;//待配准影像
	string ReferImage = foldName + "\\" + referenceName;//参考影像
	string linesText = foldName + "\\" + linesName;//参考影像
	string resultImage = foldName + "\\segment_matched.png";//结果影像
	string sourceMatch = foldName + "\\line_" + sourceName;
	string referenceMatch = foldName + "\\line_" + referenceName;

	cv::Mat img_source = cv::imread(SourceImage);
	cv::Mat img_reference = cv::imread(ReferImage);
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
	int thickness = 2;

	size_t nCount = source_ls.size();
	cv::Scalar color(0, 0, 255);
	//RNG& rng = theRNG();
	//bool isRandMatchColor = matchColor == Scalar::all(-1);
	//color = isRandMatchColor ? Scalar( rng(256), rng(256), rng(256) ) : matchColor;
	
	for( size_t i = 0; i < nCount; i++ )
	{
		thickness = 9;
		_drawKeyline( img_source, source_ls[i], color, cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
		thickness = 5;
		_drawKeyline( img_reference, reference_ls[i], color, cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	}
	
	for( size_t i = 0; i < nCount; i++ )
	{
		char temp[64];
		sprintf_s(temp, "%d", i+1);
		string strName(temp);
		_drawLabel( img_source, source_ls[i],color, cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, strName, 3);
		_drawLabel( img_reference, reference_ls[i],color, cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, strName, 2);
	}

	cv::imwrite(sourceMatch, img_source);
	cv::imwrite(referenceMatch, img_reference);
}

void drawResult(const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const char* szFilename)
{
	// plot
	// for debug
	double minX = DBL_MAX;
	double minY = DBL_MAX;
	double maxX = 0;
	double maxY = 0;

	for (int i = 0;i < (int)lineList1.size();++i)
	{
		minX = (minX > lineList1[i].pt1.x) ? lineList1[i].pt1.x : minX;
		minY = (minY > lineList1[i].pt1.y) ? lineList1[i].pt1.y : minY;
		minX = (minX > lineList1[i].pt2.x) ? lineList1[i].pt2.x : minX;
		minY = (minY > lineList1[i].pt2.y) ? lineList1[i].pt2.y : minY;
		maxX = (maxX < lineList1[i].pt1.x) ? lineList1[i].pt1.x : maxX;
		maxY = (maxY < lineList1[i].pt1.y) ? lineList1[i].pt1.y : maxY;
		maxX = (maxX < lineList1[i].pt2.x) ? lineList1[i].pt2.x : maxX;
		maxY = (maxY < lineList1[i].pt2.y) ? lineList1[i].pt2.y : maxY;
	}
	for (int i = 0;i < (int)lineList2.size();++i)
	{
		minX = (minX > lineList2[i].pt1.x) ? lineList2[i].pt1.x : minX;
		minY = (minY > lineList2[i].pt1.y) ? lineList2[i].pt1.y : minY;
		minX = (minX > lineList2[i].pt2.x) ? lineList2[i].pt2.x : minX;
		minY = (minY > lineList2[i].pt2.y) ? lineList2[i].pt2.y : minY;
		maxX = (maxX < lineList2[i].pt1.x) ? lineList2[i].pt1.x : maxX;
		maxY = (maxY < lineList2[i].pt1.y) ? lineList2[i].pt1.y : maxY;
		maxX = (maxX < lineList2[i].pt2.x) ? lineList2[i].pt2.x : maxX;
		maxY = (maxY < lineList2[i].pt2.y) ? lineList2[i].pt2.y : maxY;
	}

	//int offset = 100;
	int offset = 0;
	int offset2 = 2.0*offset;
	int nwidth = (int)(maxX - minX + 0.5) + offset2;
	int nheight = (int)(maxY - minY + 0.5) + offset2;
	nwidth = (nwidth < offset2) ? offset2 : nwidth;
	nheight = (nheight < offset2) ? offset2 : nheight;
	// create a pic
	IplImage* pic = cvCreateImage( cvSize(nwidth, nheight), 8, 3 );
	cvFloodFill(pic, cvPoint(nwidth/2+offset, nheight/2+offset), CV_RGB(255, 255, 255));

	int nCount = 0;
	for (int i = 0;i < (int)lineList1.size();++i)
	{
		cvLine(pic, cvPoint(cvRound(lineList1[i].pt1.x) + offset, cvRound(lineList1[i].pt1.y) + offset),
			cvPoint(cvRound(lineList1[i].pt2.x) + offset, cvRound(lineList1[i].pt2.y) + offset), CV_RGB(255, 0, 0), 2);
	}
	for (int i = 0;i < (int)lineList2.size();++i)
	{
		cvLine(pic, cvPoint(cvRound(lineList2[i].pt1.x) + offset, cvRound(lineList2[i].pt1.y) + offset),
			cvPoint(cvRound(lineList2[i].pt2.x) + offset, cvRound(lineList2[i].pt2.y) + offset), CV_RGB(0, 128, 255), 1);
	}
	cvSaveImage(szFilename, pic);
	cvReleaseImage(&pic);
}

void drawResult(Eigen::MatrixX2d ptList1, Eigen::MatrixX2d ptList2, const char* szFilename)
{
	// plot
	// for debug
	double minX = DBL_MAX;
	double minY = DBL_MAX;
	double maxX = 0;
	double maxY = 0;

	int n1 = ptList1.rows();
	int n2 = ptList2.rows();
	for (int i = 0;i < n1;++i)
	{
		minX = (minX > ptList1(i,0)) ? ptList1(i,0) : minX;
		minY = (minY > ptList1(i,1)) ? ptList1(i,1) : minY;
		maxX = (maxX < ptList1(i,0)) ? ptList1(i,0) : maxX;
		maxY = (maxY < ptList1(i,1)) ? ptList1(i,1) : maxY;
	}
	for (int i = 0;i < n2;++i)
	{
		minX = (minX > ptList2(i,0)) ? ptList2(i,0) : minX;
		minY = (minY > ptList2(i,1)) ? ptList2(i,1) : minY;
		maxX = (maxX < ptList2(i,0)) ? ptList2(i,0) : maxX;
		maxY = (maxY < ptList2(i,1)) ? ptList2(i,1) : maxY;
	}

	//int offset = 100;
	int offset = 0;
	int offset2 = 2.0*offset;
	int nwidth = (int)(maxX - minX + 0.5) + offset2;
	int nheight = (int)(maxY - minY + 0.5) + offset2;
	nwidth = (nwidth < offset2) ? offset2 : nwidth;
	nheight = (nheight < offset2) ? offset2 : nheight;
	// create a pic
	IplImage* pic = cvCreateImage( cvSize(nwidth, nheight), 8, 3 );
	cvFloodFill(pic, cvPoint(nwidth/2+offset, nheight/2+offset), CV_RGB(255, 255, 255));

	int nCount = 0;
	for (int i = 0;i < n1;++i)
	{
		cvCircle(pic, cvPoint(cvRound(ptList1(i,0)) + offset, cvRound(ptList1(i,1)) + offset), 3, CV_RGB(255, 0, 0));
	}
	for (int i = 0;i < n2;++i)
	{
		cvCircle(pic, cvPoint(cvRound(ptList2(i,0)) + offset, cvRound(ptList2(i,1)) + offset), 3, CV_RGB(0, 128, 255));
	}
	cvSaveImage(szFilename, pic);
	cvReleaseImage(&pic);
}

void drawResultLine(cv::Mat& pic, const cvline_polar& line, cv::Scalar color, const char* strId, int offset/* = 100*/, int thickness/* = 1*/)
{
	cv::line(pic, cvPoint(cvRound(line.pt1.x) + offset, cvRound(line.pt1.y) + offset),
		cvPoint(cvRound(line.pt2.x) + offset, cvRound(line.pt2.y) + offset), color, thickness, CV_AA);
	//cvLine(pic, cvPoint(cvRound(line.pt1.x) + offset, cvRound(line.pt1.y) + offset),
	//	cvPoint(cvRound(line.pt2.x) + offset, cvRound(line.pt2.y) + offset), color, thickness);

	double sin_theta = fabs(sin(line.theta));
	double cos_theta = fabs(cos(line.theta));
	double dis = 5.0;
	int centerx = cvRound((line.pt1.x + line.pt2.x)*0.5 + dis*sin_theta + offset);
	int centery = cvRound((line.pt1.y + line.pt2.y)*0.5 - dis*cos_theta + offset);
	//cvPutText(pic, strId, cvPoint(centerx, centery), &font, color);
	cv::putText(pic, strId, cvPoint(centerx, centery), CV_FONT_HERSHEY_SIMPLEX, 0.4, color);
}

void drawLastResult(string strImage1, string strImage2, vector<double> param,
					const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const vector<int>& matcher, const char* szFilename1, const char* szFilename2)
{
	// plot
	// for debug
	cv::Mat image1 = cv::imread(strImage1);
	int nwidth1 = image1.cols;
	int nheight1 = image1.rows;

	cv::Mat image2= affine_warp_Image(strImage2.c_str(), param);
	int nwidth2 = image2.cols;
	int nheight2 = image2.rows;


	// create a pic
	//IplImage* pic1 = cvCreateImage( cvSize(nwidth1, nheight1), 8, 3 );
	//IplImage* pic2 = cvCreateImage( cvSize(nwidth2, nheight2), 8, 3 );

	int n = (int)matcher.size() / 2;
	CvRNG rng;
	for (int i = 0;i < n;++i)
	{
		//cv::Scalar color = random_color(&rng);
		//cv::Scalar color = CV_RGB(255, 0, 0);
		cv::Scalar color = CV_RGB(0, 255, 0);
		char strId[256];
		sprintf_s(strId, "%d", i+1);
		drawResultLine(image1, lineList1[matcher[2*i]], color, strId, 0, 1);
		drawResultLine(image2, lineList2[matcher[2*i+1]], color, strId, 0, 1);
	}
	cv::imwrite(szFilename1, image1);
	cv::imwrite(szFilename2, image2);
	//cvSaveImage(szFilename1, image1);
	//cvSaveImage(szFilename2, image2);
}

void drawLastResult(string strImage1, string strImage2,
					const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const vector<int>& matcher, const char* szFilename1, const char* szFilename2)
{
	// plot
	// for debug
	cv::Mat image1 = cv::imread(strImage1);
	int nwidth1 = image1.cols;
	int nheight1 = image1.rows;

	cv::Mat image2 = cv::imread(strImage2);
	int nwidth2 = image2.cols;
	int nheight2 = image2.rows;


	// create a pic
	//IplImage* pic1 = cvCreateImage( cvSize(nwidth1, nheight1), 8, 3 );
	//IplImage* pic2 = cvCreateImage( cvSize(nwidth2, nheight2), 8, 3 );

	int n = (int)matcher.size() / 2;
	CvRNG rng;
	for (int i = 0;i < n;++i)
	{
		//cv::Scalar color = random_color(&rng);
		cv::Scalar color = CV_RGB(255, 0, 0);
		char strId[256];
		sprintf_s(strId, "%d", i+1);
		int thickness = 2;
		drawResultLine(image1, lineList1[matcher[2*i]], color, strId, 0, thickness);
		drawResultLine(image2, lineList2[matcher[2*i+1]], color, strId, 0, thickness);
	}
	cv::imwrite(szFilename1, image1);
	cv::imwrite(szFilename2, image2);
	//cvSaveImage(szFilename1, image1);
	//cvSaveImage(szFilename2, image2);
}

void drawLastResult(const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const vector<int>& matcher, const char* szFilename1, const char* szFilename2)
{
	// plot
	// for debug
	double minX = DBL_MAX;
	double minY = DBL_MAX;
	double maxX = 0;
	double maxY = 0;

	for (int i = 0;i < (int)lineList1.size();++i)
	{
		minX = (minX > lineList1[i].pt1.x) ? lineList1[i].pt1.x : minX;
		minY = (minY > lineList1[i].pt1.y) ? lineList1[i].pt1.y : minY;
		minX = (minX > lineList1[i].pt2.x) ? lineList1[i].pt2.x : minX;
		minY = (minY > lineList1[i].pt2.y) ? lineList1[i].pt2.y : minY;
		maxX = (maxX < lineList1[i].pt1.x) ? lineList1[i].pt1.x : maxX;
		maxY = (maxY < lineList1[i].pt1.y) ? lineList1[i].pt1.y : maxY;
		maxX = (maxX < lineList1[i].pt2.x) ? lineList1[i].pt2.x : maxX;
		maxY = (maxY < lineList1[i].pt2.y) ? lineList1[i].pt2.y : maxY;
	}
	for (int i = 0;i < (int)lineList2.size();++i)
	{
		minX = (minX > lineList2[i].pt1.x) ? lineList2[i].pt1.x : minX;
		minY = (minY > lineList2[i].pt1.y) ? lineList2[i].pt1.y : minY;
		minX = (minX > lineList2[i].pt2.x) ? lineList2[i].pt2.x : minX;
		minY = (minY > lineList2[i].pt2.y) ? lineList2[i].pt2.y : minY;
		maxX = (maxX < lineList2[i].pt1.x) ? lineList2[i].pt1.x : maxX;
		maxY = (maxY < lineList2[i].pt1.y) ? lineList2[i].pt1.y : maxY;
		maxX = (maxX < lineList2[i].pt2.x) ? lineList2[i].pt2.x : maxX;
		maxY = (maxY < lineList2[i].pt2.y) ? lineList2[i].pt2.y : maxY;
	}

	int offset = 100;
	int offset2 = 2.0*offset;
	int nwidth = (int)(maxX - minX + 0.5) + offset2;
	int nheight = (int)(maxY - minY + 0.5) + offset2;
	nwidth = (nwidth < offset2) ? offset2 : nwidth;
	nheight = (nheight < offset2) ? offset2 : nheight;
	// create a pic
	cv::Mat image1(cv::Size(nwidth, nheight),CV_8UC3);
	cv::Mat image2(cv::Size(nwidth, nheight),CV_8UC3);
	cv::floodFill(image1, cvPoint(nwidth/2+offset, nheight/2+offset), CV_RGB(255, 255, 255));
	cv::floodFill(image2, cvPoint(nwidth/2+offset, nheight/2+offset), CV_RGB(255, 255, 255));

	int n = (int)matcher.size() / 2;
	CvRNG rng;
	srand( (unsigned int)time(0) );
	for (int i = 0;i < n;++i)
	{
		cv::Scalar color = random_color(&rng);
		char strId[256];
		sprintf_s(strId, "%d", i+1);
		drawResultLine(image1, lineList1[matcher[2*i]], color, strId, 100, 2);
		drawResultLine(image2, lineList2[matcher[2*i+1]], color, strId, 100, 2);
	}
	cv::imwrite(szFilename1, image1);
	cv::imwrite(szFilename2, image2);
}

void drawPolarResult(const std::vector<cvline_polar>& lineList1, const std::vector<cvline_polar>& lineList2, const char* szFilename)
{
	// plot
	// for debug
	double minTheta = DBL_MAX;
	double minRho = DBL_MAX;
	double maxTheta = -DBL_MAX;
	double maxRho = -DBL_MAX;

	for (int i = 0;i < (int)lineList1.size();++i)
	{
		minTheta = (minTheta > lineList1[i].theta) ? lineList1[i].theta : minTheta;
		minRho = (minRho > lineList1[i].rho) ? lineList1[i].rho : minRho;
		maxTheta = (maxTheta < lineList1[i].theta) ? lineList1[i].theta : maxTheta;
		maxRho = (maxRho < lineList1[i].rho) ? lineList1[i].rho : maxRho;
	}
	for (int i = 0;i < (int)lineList2.size();++i)
	{
		minTheta = (minTheta > lineList2[i].theta) ? lineList2[i].theta : minTheta;
		minRho = (minRho > lineList2[i].rho) ? lineList2[i].rho : minRho;
		maxTheta = (maxTheta < lineList2[i].theta) ? lineList2[i].theta : maxTheta;
		maxRho = (maxRho < lineList2[i].rho) ? lineList2[i].rho : maxRho;
	}

	int offset = 100;
	int offset2 = 2.0*offset;
	// imagesize = (720 + 100 * 2 ) * (720 + 100 * 2 )
	int nwidth = 920;
	int nheight = 920;
	double scale_Theta = 720.0 / (maxTheta - minTheta + DBL_EPSILON);
	double scale_Rho = 720.0 / (maxRho - minRho + DBL_EPSILON);

	// create a pic
	IplImage* pic = cvCreateImage( cvSize(nwidth, nheight), 8, 3 );

	int nCount = 0;
	for (int i = 0;i < (int)lineList1.size();++i)
	{
		cvCircle(pic, cvPoint(cvRound(offset + (lineList1[i].rho - minRho)*scale_Rho), cvRound(offset + 360 + lineList1[i].theta*360.0/PI)), 5, CV_RGB(255, 0, 0));
	}
	for (int i = 0;i < (int)lineList2.size();++i)
	{
		cvCircle(pic, cvPoint(cvRound(offset + (lineList2[i].rho - minRho)*scale_Rho), cvRound(offset + 360 + lineList2[i].theta*360.0/PI)), 5, CV_RGB(0, 255, 0));
	}
	cvSaveImage(szFilename, pic);
	cvReleaseImage(&pic);
}

cv::Size getNewSize(int iRows, int iCols, cv::Mat transMat)
{
	double pos[8] = {0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
	pos[3] = double(iRows);
	pos[4] = double(iCols);
	pos[6] = double(iCols);
	pos[7] = double(iRows);


	//double *m = new double [transMat.rows * transMat.cols];
	vector<double> m(transMat.rows * transMat.cols);
	for(int i = 0;i < transMat.rows;i++)
	{
		for(int j = 0;j < transMat.cols;j++)
		{
			m[i * transMat.cols + j] = transMat.at<double>(i, j);
		}
	}

	double newPos[8];
	for(int i = 0;i < 4;i++)
	{
		double element_x = pos[i * 2 + 0];
		double element_y = pos[i * 2 + 1];
		//float newX = m[3 * 1 + 1] * element_x + m[3 * 0 + 1] * element_y - m[3 * 0 + 2] * m[3 * 1 + 1] + m[3 * 0 + 1] * m[3 * 1 + 2];
		//newX = newX / (m[3 * 0 + 0] * m[3 * 1 + 1] - m[3 * 0 + 1] * m[3 * 1 + 0]);

		//float newY = m[3 * 1 + 0] * element_x + m[3 * 0 + 0] * element_y - m[3 * 0 + 2] * m[3 * 1 + 0] + m[3 * 0 + 0] * m[3 * 1 + 2];
		//newY = newY / (m[3 * 0 + 1] * m[3 * 1 + 0] - m[3 * 0 + 0] * m[3 * 1 + 1]);

		double newX = m[3 * 0 + 0] * element_x + m[3 * 0 + 1] * element_y + m[3 * 0 + 2];
		double newY = m[3 * 1 + 0] * element_x + m[3 * 1 + 1] * element_y + m[3 * 1 + 2];

		newPos[i * 2 + 0] = newX;
		newPos[i * 2 + 1] = newY;
	}
	double element_x = newPos[0];
	double element_y = newPos[1];
	double minX = element_x;
	double maxX = element_x;
	double minY = element_y;
	double maxY = element_y;
	for(int i = 1;i < 4;i++)
	{
		double element_x = newPos[i * 2 + 0];
		double element_y = newPos[i * 2 + 1];
		if(element_x < minX)
		{
			minX = element_x;
		}
		if(element_x > maxX)
		{
			maxX = element_x;
		}
		if(element_y < minY)
		{
			minY = element_y;
		}
		if(element_y > maxY)
		{
			maxY = element_y;
		}
	}

	int newCols = int(fabs(maxX) + 0.5);
	int newRows = int(fabs(maxY) + 0.5);
	//int newRows = int(max(fabs(maxX), fabs(minX)) + 0.5);
	//int newCols = int(max(fabs(maxY), fabs(minY)) + 0.5);
	//delete []m;
	cv::Size newSize(newCols, newRows);
	return newSize;
}

cv::Size getPerspectiveNewSize(int iRows, int iCols, cv::Mat transMat)
{
	double pos[8] = {0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
	pos[3] = double(iRows);
	pos[4] = double(iCols);
	pos[6] = double(iCols);
	pos[7] = double(iRows);


	//double *m = new double [transMat.rows * transMat.cols];
	vector<double> m(transMat.rows * transMat.cols);
	for(int i = 0;i < transMat.rows;i++)
	{
		for(int j = 0;j < transMat.cols;j++)
		{
			m[i * transMat.cols + j] = transMat.at<double>(i, j);
		}
	}

	double newPos[8];
	for(int i = 0;i < 4;i++)
	{
		double element_x = pos[i * 2 + 0];
		double element_y = pos[i * 2 + 1];
		
		double newX = (m[3 * 0 + 0] * element_x + m[3 * 0 + 1] * element_y + m[3 * 0 + 2])/
			(m[3 * 2 + 0] * element_x + m[3 * 2 + 1] * element_y + m[3 * 2 + 2]);
		double newY = (m[3 * 1 + 0] * element_x + m[3 * 1 + 1] * element_y + m[3 * 1 + 2])/
			(m[3 * 2 + 0] * element_x + m[3 * 2 + 1] * element_y + m[3 * 2 + 2]);

		newPos[i * 2 + 0] = newX;
		newPos[i * 2 + 1] = newY;
	}
	double element_x = newPos[0];
	double element_y = newPos[1];
	double minX = element_x;
	double maxX = element_x;
	double minY = element_y;
	double maxY = element_y;
	for(int i = 1;i < 4;i++)
	{
		double element_x = newPos[i * 2 + 0];
		double element_y = newPos[i * 2 + 1];
		if(element_x < minX)
		{
			minX = element_x;
		}
		if(element_x > maxX)
		{
			maxX = element_x;
		}
		if(element_y < minY)
		{
			minY = element_y;
		}
		if(element_y > maxY)
		{
			maxY = element_y;
		}
	}

	int newCols = int(fabs(maxX) + 0.5);
	int newRows = int(fabs(maxY) + 0.5);
	//int newRows = int(max(fabs(maxX), fabs(minX)) + 0.5);
	//int newCols = int(max(fabs(maxY), fabs(minY)) + 0.5);
	//delete []m;
	cv::Size newSize(newCols, newRows);
	return newSize;
}
Eigen::VectorXd affine_ransac(vector<cvline_polar> matched_ref_lines_list,
	vector<cvline_polar> matched_src_lines_list, vector<int> &inliers)
{
	// RANSAC detect outliers
	auto_ptr< estimators::Solver<mat,vec> > ptrSolver(
		new estimators::lineAffineSolver<mat,vec>);

	//-- Create input data
	int nTotal = (int)matched_ref_lines_list.size();
	mat dataPoints(nTotal, 6);
	for(int i = 0;i < nTotal;++i)
	{
		dataPoints(i, 0) = matched_src_lines_list[i].rho;
		dataPoints(i, 1) = matched_src_lines_list[i].theta;
		dataPoints(i, 2) = matched_ref_lines_list[i].pt1.x;
		dataPoints(i, 3) = matched_ref_lines_list[i].pt1.y;
		dataPoints(i, 4) = matched_ref_lines_list[i].pt2.x;
		dataPoints(i, 5) = matched_ref_lines_list[i].pt2.y;
	}
	
	vector<vec> models;
	ransac::Ransac_Handler ransac_fun_Handler;
	bool result = ransac::Ransac_RobustEstimator
		(
		dataPoints, // the input data
		estimators::affineSolver<mat,vec>::extractor, // How select sampled point from indices
		dataPoints.n_rows,  // the number of putatives data
		*(ptrSolver.get()),  // compute the underlying model given a sample set
		estimators::affineSolver<mat,vec>::defaultEvaluator,  // the function to evaluate a given model
		//Ransac Object that contain function:
		// CandidatesSelector, Sampler and TerminationFunction
		ransac_fun_Handler, // the basic ransac object
		1000,  // the maximum rounds for RANSAC routine
		inliers, // inliers to the final solution
		models, // models array that fit input data
		0.95, // the confidence want to achieve at the end
		15.0
		);

	int nParameters = (int)models[0].n_rows;
	Eigen::VectorXd parameters(nParameters);
	for (int i = 0;i < nParameters;++i)
	{
		cout<<models[0](i)<<endl;
		parameters(i) = models[0](i);
	}
	cout<<endl;
	return parameters;
}