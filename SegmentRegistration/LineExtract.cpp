#include "LineExtract.h"
#include "lsd.h"
#include "func.h"
#include "LSWMS.h"

#include "quick_selection.h"
//
//IplImage* get_lines(IplImage* img,vector<cvline_polar>& vec_lines)
//{
//    //to grey
//    //IplImage* grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
//    //cvCvtColor(img, grey, CV_BGR2GRAY);
//
//    image_double image;
//    ntuple_list out;
//    unsigned int x,y,i,j;
//    image = new_image_double(img->width,img->height);
//    for(x=0;x</*grey*/img->width;x++)
//    for(y=0;y</*grey*/img->height;y++)
//    {
//      CvScalar s= cvGet2D(/*grey*/img,y,x);
//      double pix= s.val[0];
//      image->data[ x + y * image->xsize ]= pix;
//    }
//
//    /* call LSD */
//    out = lsd(image);
//    //out= lsd_scale(image,1);
//
//    /* print output */
//    //printf("%u line segments found:\n",out->size);
//    //vector<Line> vec;
//    for(i=0;i<out->size;i++)
//    {
//      //for(j=0;j<out->dim;j++)
//      {
//        //printf("%f ",out->values[ i * out->dim + j ]);
//
//		  cvline_polar outLine;
//		  fPoint line[2];
//		  line[0].x = out->values[ i * out->dim + 0];
//		  line[0].y = out->values[ i * out->dim + 1];
//		  line[1].x = out->values[ i * out->dim + 2];
//		  line[1].y = out->values[ i * out->dim + 3];
//
//		  outLine = line2polar(line);
//		  double len = segment_length(outLine);
//		  if(len > 5.0)
//          /*vec*/vec_lines.push_back(outLine);
//      }
//      //printf("\n");
//    }
//
//    IplImage* black= cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
//    cvZero(black);
//    for(int i=0;i<vec_lines.size();++i)
//    {
//        //if(vec[i].x1==vec[i].x2||vec[i].y1==vec[i].y2)
//        cvLine(black,cvPoint(vec_lines[i].pt1.x,vec_lines[i].pt1.y),cvPoint(vec_lines[i].pt2.x,vec_lines[i].pt2.y),CV_RGB(255,255,255),1, CV_AA);
//    }
//    /*cvNamedWindow("img", 0);
//    cvShowImage("img", img);*/
//    //cvSaveImage("lines_detect.png",black/*img*/);
//    /* free memory */
//    //cvReleaseImage(&grey);
//    free_image_double(image);
//    free_ntuple_list(out);
//
//    return black;
//}

bool isBoundary(IplImage *img, cvline_polar line, cv::Scalar bgColor/* = CV_RGB(0,0,0)*/)
{
	int X = img->width;
	int Y = img->height;
	int band = img->nChannels;
	int step = img->widthStep;
	uchar * data = (uchar *)img->imageData;

	fPoint pt1 = line.pt1;
	fPoint pt2 = line.pt2;

	int k = 3;
	// 外扩k个像素
	cvline_polar l1 = segment_move_vertical(line, (double)k);
	cvline_polar l2 = segment_move_vertical(line, -(double)k);

	double len = (int)(sqrt((l1.pt1.x-l1.pt2.x)*(l1.pt1.x-l1.pt2.x)+(l1.pt1.y-l1.pt2.y)*(l1.pt1.y-l1.pt2.y)) + 0.5);
	int i;
	for (i = 0;i <= len;++i)
	{
		double lammda = i / (double)len;
		int x = (int)(lammda*(l1.pt1.x)+(1.0-lammda)*l1.pt2.x+0.5);
		int y = (int)(lammda*(l1.pt1.y)+(1.0-lammda)*l1.pt2.y+0.5);
		if (x <0 || x>= X || y < 0 || y >=Y)
		{
			continue;
		}

		bool same_color = true;
		for (int k = 0;k < band;++k)
		{
			if (bgColor.val[k] != data[ y*step + x*band + k])
			{
				same_color = false;
				break;
			}
		}
		if (!same_color)
		{
			break;
		}
	}
	if (i > len)
	{
		return true;
	}

	len = (int)(sqrt((l2.pt1.x-l2.pt2.x)*(l2.pt1.x-l2.pt2.x)+(l2.pt1.y-l2.pt2.y)*(l2.pt1.y-l2.pt2.y)) + 0.5);
	for (i = 0;i <= len;++i)
	{
		double lammda = i / (double)len;
		int x = (int)(lammda*(l2.pt1.x)+(1.0-lammda)*l2.pt2.x+0.5);
		int y = (int)(lammda*(l2.pt1.y)+(1.0-lammda)*l2.pt2.y+0.5);
		if (x <0 || x>= X || y < 0 || y >=Y)
		{
			continue;
		}

		bool same_color = true;
		for (int k = 0;k < band;++k)
		{
			if (bgColor[k] != img->imageData[ y*step + x*band + k])
			{
				same_color = false;
				break;
			}
		}
		if (!same_color)
		{
			break;
		}
	}
	if (i > len)
	{
		return true;
	}

	return false;
}


bool isBoundary(cv::Mat img, cvline_polar line, cv::Scalar bgColor/* = CV_RGB(0,0,0)*/)
{
	int X = img.cols;
	int Y = img.rows;
	int band = img.channels();
	int step = img.step;
	uchar * data = (uchar *)img.data;

	fPoint pt1 = line.pt1;
	fPoint pt2 = line.pt2;

	int k = 5;
	// 外扩k个像素
	cvline_polar l1 = segment_move_vertical(line, (double)k);
	cvline_polar l2 = segment_move_vertical(line, -(double)k);

	double len = (int)(sqrt((l1.pt1.x-l1.pt2.x)*(l1.pt1.x-l1.pt2.x)+(l1.pt1.y-l1.pt2.y)*(l1.pt1.y-l1.pt2.y)) + 0.5);
	int i;
	for (i = 0;i <= len;++i)
	{
		double lammda = i / (double)len;
		int x = (int)(lammda*(l1.pt1.x)+(1.0-lammda)*l1.pt2.x+0.5);
		int y = (int)(lammda*(l1.pt1.y)+(1.0-lammda)*l1.pt2.y+0.5);
		if (x <0 || x>= X || y < 0 || y >=Y)
		{
			continue;
		}

		bool same_color = true;
		for (int k = 0;k < band;++k)
		{
			if (bgColor.val[k] != data[ y*step + x*band + k])
			{
				same_color = false;
				break;
			}
		}
		if (!same_color)
		{
			break;
		}
	}
	if (i > len)
	{
		return true;
	}

	len = (int)(sqrt((l2.pt1.x-l2.pt2.x)*(l2.pt1.x-l2.pt2.x)+(l2.pt1.y-l2.pt2.y)*(l2.pt1.y-l2.pt2.y)) + 0.5);
	for (i = 0;i <= len;++i)
	{
		double lammda = i / (double)len;
		int x = (int)(lammda*(l2.pt1.x)+(1.0-lammda)*l2.pt2.x+0.5);
		int y = (int)(lammda*(l2.pt1.y)+(1.0-lammda)*l2.pt2.y+0.5);
		if (x <0 || x>= X || y < 0 || y >=Y)
		{
			continue;
		}

		bool same_color = true;
		for (int k = 0;k < band;++k)
		{
			if (bgColor[k] != img.data[ y*step + x*band + k])
			{
				same_color = false;
				break;
			}
		}
		if (!same_color)
		{
			break;
		}
	}
	if (i > len)
	{
		return true;
	}

	return false;
}

IplImage* get_lines(IplImage* img, vector<cvline_polar>& vec_lines, double scale/* = 0.8*/, int seletion_num/* = 200*/)
{
    //to grey
    //IplImage* grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    //cvCvtColor(img, grey, CV_BGR2GRAY);

	double * image;
	double * out;

	int X = img->width;
	int Y = img->height;


	/* create a simple image: left half black, right half gray */
	image = (double *) malloc( X * Y * sizeof(double) );
	if( image == NULL )
	{
		fprintf(stderr,"error: not enough memory\n");
		exit(EXIT_FAILURE);
	}

    for(int i=0;i</*grey*/img->width;i++)
	{
		for(int j=0;j</*grey*/img->height;j++)
		{
		  CvScalar s= cvGet2D(/*grey*/img,j,i);
		  double pix= s.val[0];
		  image[ i + j * X ]= pix;
		}
	}

    /* call LSD */
	int n;
	//out = lsd(&n, image, X, Y);
	out = lsd_scale(&n, image, X, Y, scale);

	// x1,y1,x2,y2,width,p,-log10(NFA);
    //out= lsd_scale(image,1);

    /* print output */
    //printf("%u line segments found:\n",out->size);
    //vector<Line> vec;

	double minLen, minNFA;
	//double *nfas = new double[n];
	vector<double> nfas;
	minLen = 5.0;
	vector<int> notBoundary;
	for(int i=0;i<n;i++)
	{
		cvline_polar outLine;
		fPoint line[2];
		line[0].x = out[ i * 7 + 0];
		line[0].y = out[ i * 7 + 1];
		line[1].x = out[ i * 7 + 2];
		line[1].y = out[ i * 7 + 3];

		outLine = line2polar(line);
		if (isBoundary(img, outLine, cvScalar((0,0,0))))
		{
			//continue;
		}
		notBoundary.push_back(i);
		double len = segment_length(outLine);
		if(len >= minLen) nfas.push_back(-out[i*7+6]);
	}

	minNFA = -quick_select(&nfas[0], 0, (int)nfas.size()-1, seletion_num);
	//for(int i=0;i<(int)nfas.size();i++)
	for(int m=0;m<(int)notBoundary.size();m++)
    {
		int i = notBoundary[m];
      //for(j=0;j<out->dim;j++)
      {
        //printf("%f ",out->values[ i * out->dim + j ]);

		  cvline_polar outLine;
		  fPoint line[2];
		  line[0].x = out[ i * 7 + 0];
		  line[0].y = out[ i * 7 + 1];
		  line[1].x = out[ i * 7 + 2];
		  line[1].y = out[ i * 7 + 3];

		  //double w = out[i*7+4];
		  //double p = out[i*7+5];
		  //double NFA = out[i*7+6];

		  outLine = line2polar(line);

		  if(out[i*7+6] >= minNFA)
		  {
			  vec_lines.push_back(outLine);
		  }
      }
      //printf("\n");
    }

    IplImage* black= cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    cvZero(black);
    for(int i=0;i<(int)vec_lines.size();++i)
    {
        //if(vec[i].x1==vec[i].x2||vec[i].y1==vec[i].y2)
        cvLine(black,cvPoint(vec_lines[i].pt1.x,vec_lines[i].pt1.y),cvPoint(vec_lines[i].pt2.x,vec_lines[i].pt2.y),CV_RGB(255,255,255),2, CV_AA);
	}

	/* free memory */
	free( (void *) image );
	free( (void *) out );

    return black;
}

IplImage* get_lines(const char* strImage, vector<cvline_polar>& vec_lines, double scale/* = 0.8*/, int seletion_num/* = 200*/)
{
    //to grey
    //IplImage* grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    //cvCvtColor(img, grey, CV_BGR2GRAY);
	
	IplImage* img = cvLoadImage(strImage, 0 );
	if (!img)
	{
		cout<<"open image failed for \""<<strImage<<"\"!"<<endl;
		return NULL;
	}

	double * image;
	double * out;

	int X = img->width;
	int Y = img->height;


	/* create a simple image: left half black, right half gray */
	image = (double *) malloc( X * Y * sizeof(double) );
	if( image == NULL )
	{
		fprintf(stderr,"error: not enough memory\n");
		exit(EXIT_FAILURE);
	}

    for(int i=0;i</*grey*/img->width;i++)
	{
		for(int j=0;j</*grey*/img->height;j++)
		{
		  CvScalar s= cvGet2D(/*grey*/img,j,i);
		  double pix= s.val[0];
		  image[ i + j * X ]= pix;
		}
	}

    /* call LSD */
	int n;
	//out = lsd(&n, image, X, Y);
	out = lsd_scale(&n, image, X, Y, scale);

	cout<<"all detected lines:"<<n<<endl;

	// x1,y1,x2,y2,width,p,-log10(NFA);
    //out= lsd_scale(image,1);

    /* print output */
    //printf("%u line segments found:\n",out->size);
    //vector<Line> vec;

	double minLen, minNFA;
	//double *nfas = new double[n];
	vector<double> nfas;
	minLen = 5.0;
	vector<int> notBoundary;
	for(int i=0;i<n;i++)
	{
		cvline_polar outLine;
		fPoint line[2];
		line[0].x = out[ i * 7 + 0];
		line[0].y = out[ i * 7 + 1];
		line[1].x = out[ i * 7 + 2];
		line[1].y = out[ i * 7 + 3];

		outLine = line2polar(line);
		if (isBoundary(img, outLine, cvScalar((0,0,0))))
		{
			continue;
		}
		notBoundary.push_back(i);
		double len = segment_length(outLine);
		if(len >= minLen) nfas.push_back(-out[i*7+6]);
	}
	if (nfas.size() < 3)
	{
		return NULL;
	}
	minNFA = -quick_select(&nfas[0], 0, (int)nfas.size()-1, seletion_num);
	//for(int i=0;i<(int)nfas.size();i++)
	for(int m=0;m<(int)notBoundary.size();m++)
    {
		int i = notBoundary[m];
      //for(j=0;j<out->dim;j++)
      {
        //printf("%f ",out->values[ i * out->dim + j ]);

		  cvline_polar outLine;
		  fPoint line[2];
		  line[0].x = out[ i * 7 + 0];
		  line[0].y = out[ i * 7 + 1];
		  line[1].x = out[ i * 7 + 2];
		  line[1].y = out[ i * 7 + 3];

		  //double w = out[i*7+4];
		  //double p = out[i*7+5];
		  //double NFA = out[i*7+6];

		  outLine = line2polar(line);

		  if(out[i*7+6] >= minNFA)
		  {
			  vec_lines.push_back(outLine);
		  }
      }
      //printf("\n");
    }

    //IplImage* black= cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    //cvZero(black);
	IplImage* line_img = cvLoadImage(strImage, 1 );
    for(int i=0;i<(int)vec_lines.size();++i)
    {
        //if(vec[i].x1==vec[i].x2||vec[i].y1==vec[i].y2)
        cvLine(line_img,cvPoint(vec_lines[i].pt1.x,vec_lines[i].pt1.y),cvPoint(vec_lines[i].pt2.x,vec_lines[i].pt2.y),CV_RGB(255,0,0),2, CV_AA);
	}

	/* free memory */
	free( (void *) image );
	free( (void *) out );
	cvReleaseImage(&img);

    return line_img;
}

bool GetImGray(IplImage* im, double x, double y, double &gray)
{
	if (x < 1 || x > (im->width-2) || y < 1 || y > (im->height-2))
		return false;
	int x1 = (int)x;
	int y1 = (int)y;
	int x2 = x1 + 1;
	int y2 = y1 + 1;
	int step = im->widthStep;
	int _x1 = x1 < 0 ? 0 : x1;
	int _y1 = y1 < 0 ? 0 : y1;
	int _x2 = x2 > (im->width-1) ? (im->width-1) : x2;
	int _y2 = y2 > (im->height-1) ? (im->height-1) : y2;
	gray = (x2 - x) * (y2 - y) * ((BYTE*)im->imageData)[_y1*step+_x1]/1.0 + 
		(x - x1) * (y2 - y) * ((BYTE*)im->imageData)[_y1*step+_x2]/1.0 + 
		(x2 - x) * (y - y1) * ((BYTE*)im->imageData)[_y2*step+_x1]/1.0 + 
		(x - x1) * (y - y1) * ((BYTE*)im->imageData)[_y2*step+_x2]/1.0;
	return true;
}


bool lp_LineDetect(const char* strImage, std::vector<Line>& lineList, double scale/* = 0.8*/, int seletion_num/* = 200*/)
{
	//to grey
	//IplImage* grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
	//cvCvtColor(img, grey, CV_BGR2GRAY);

	IplImage* img = cvLoadImage(strImage, 0 );
	if (!img)
	{
		cout<<"open image failed for \""<<strImage<<"\"!"<<endl;
		return NULL;
	}

	double * image;
	double * out;

	int X = img->width;
	int Y = img->height;


	/* create a simple image: left half black, right half gray */
	image = (double *) malloc( X * Y * sizeof(double) );
	if( image == NULL )
	{
		fprintf(stderr,"error: not enough memory\n");
		exit(EXIT_FAILURE);
	}

	for(int i=0;i</*grey*/img->width;i++)
	{
		for(int j=0;j</*grey*/img->height;j++)
		{
			CvScalar s= cvGet2D(/*grey*/img,j,i);
			double pix= s.val[0];
			image[ i + j * X ]= pix;
		}
	}

	/* call LSD */
	int n;
	//out = lsd(&n, image, X, Y);
	out = lsd_scale(&n, image, X, Y, scale);

	cout<<"all detected lines:"<<n<<endl;


	double minLen, minNFA;
	//double *nfas = new double[n];
	vector<double> nfas;
	minLen = 5.0;
	vector<int> notBoundary;
	for(int i=0;i<n;i++)
	{
		cvline_polar outLine;
		fPoint line[2];
		line[0].x = out[ i * 7 + 0];
		line[0].y = out[ i * 7 + 1];
		line[1].x = out[ i * 7 + 2];
		line[1].y = out[ i * 7 + 3];

		outLine = line2polar(line);
		if (isBoundary(img, outLine, cvScalar((0,0,0))))
		{
			continue;
		}
		notBoundary.push_back(i);
		double len = segment_length(outLine);
		if(len >= minLen) nfas.push_back(-out[i*7+6]);
	}

	minNFA = -quick_select(&nfas[0], 0, (int)nfas.size()-1, seletion_num);

	IplImage* smoothed_im = cvCreateImage(cvSize(img->width,img->height),IPL_DEPTH_8U,1);
	cvSmooth(img,smoothed_im);
	cvReleaseImage(&img);
	//for(int i=0;i<(int)nfas.size();i++)
	for(int m=0;m<(int)notBoundary.size();m++)
	{
		int i = notBoundary[m];
		double x1 = out[ i * 7 + 0];
		double y1 = out[ i * 7 + 1];
		double x2 = out[ i * 7 + 2];
		double y2 = out[ i * 7 + 3];
		double len = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
		len = sqrt(len);
		if(len <= 20) continue;
		double a = y1 - y2;
		double b = -(x1 - x2);
		double c = x1 * y2 - x2 * y1;
		vector<CvPoint2D32f> Pts;
		int nPts = 0;
		CvPoint2D32f p;
		p.x = x1;
		p.y = y1;
		Pts.push_back(p);
		nPts++;
		if ( abs((a+0.0000000000001) / (b+0.0000000000001)) < 1 )
		{
			double iterval = (x2 - x1) / len;
			while (1)
			{
				p.x += iterval;
				if( (iterval > 0 && p.x > x2) || (iterval < 0 && p.x < x2) ) break;
				p.y = -(a * p.x + c) / b;
				Pts.push_back(p);
				nPts++;
			}
		}
		else
		{
			double iterval = (y2 - y1) / len;
			while (1)
			{
				p.y += iterval;
				if( (iterval > 0 && p.y > y2) || (iterval < 0 && p.y < y2) ) break;
				p.x = -(b * p.y + c) / a;
				Pts.push_back(p);
				nPts++;
			}
		}
		double line_grad_x, line_grad_y;
		int nValid = 0;
		line_grad_x = line_grad_y = 0;
		for (int j = 0;j < nPts;j++)
		{
			double gray2;
			double gray1;
			if( !GetImGray(smoothed_im, Pts[j].x + 1, Pts[j].y, gray2) ) continue;
			if( !GetImGray(smoothed_im, Pts[j].x - 1, Pts[j].y, gray1) ) continue;
			line_grad_x += (gray2 - gray1) / 2;
			if( !GetImGray(smoothed_im, Pts[j].x, Pts[j].y + 1, gray2) ) continue;
			if( !GetImGray(smoothed_im, Pts[j].x, Pts[j].y - 1, gray1) ) continue;
			line_grad_y += (gray2 - gray1) / 2;
			nValid++;
		}
		if( nValid == 0) continue;
		line_grad_x /= nValid;
		line_grad_y /= nValid;
		double ExProd = (Pts[nPts-1].x - Pts[0].x) * line_grad_y - line_grad_x * (Pts[nPts-1].y - Pts[0].y);

		Line lp_line;
		if (ExProd > 0)
		{
			lp_line.StartPt = Pts[0];
			lp_line.EndPt = Pts[nPts-1];
		}
		else
		{
			lp_line.StartPt = Pts[nPts-1];
			lp_line.EndPt = Pts[0];
		}
		lp_line.length = len;
		lp_line.Center.x = (x1 + x2) / 2;
		lp_line.Center.y = (y1 + y2) / 2;
		lp_line.para_a = a;
		lp_line.para_b = b;
		lp_line.para_c = c;

		if(out[i*7+6] >= minNFA)
		{
			lineList.push_back(lp_line);
		}
	}


	/* free memory */
	free( (void *) image );
	free( (void *) out );
	cvReleaseImage(&smoothed_im);
	return true;
}

IplImage* get_lines_lswms(const char* strImage, vector<cvline_polar>& vec_lines, int numMaxLSegs/* = 200*/)
{
	// Images
	cv::Mat inputImg, imgGRAY;	
	inputImg = cv::imread(strImage);
	if(inputImg.empty())
		return NULL;

	int width = 0, height = 0;
	width = inputImg.cols;
	height = inputImg.rows;

	cv::Size procSize = cv::Size(width, height);
	printf("Input image: %s, Size (%d x %d)\n", strImage, width, height);

	bool playMode = false;
	// ---------------------------
	// Create and init LSWMS
	int R = 3;
	LSWMS lswms(procSize, R, numMaxLSegs, playMode);
	if(numMaxLSegs==0)
		printf("LSWMS object created: R=%d\n\n", R);
	else
		printf("LSWMS object created: R=%d, numMaxLSegs=%d\n\n", R, numMaxLSegs);

	// Color Conversion
	if(inputImg.channels() == 3)
	{
		cv::cvtColor(inputImg, imgGRAY, CV_BGR2GRAY);	
	}
	else
	{
		inputImg.copyTo(imgGRAY);
	}
	std::vector<LSEG> lSegs;
	std::vector<double> errors;
	// Process LSWMS

	lswms.run(inputImg, lSegs, errors);


	int n = (int)lSegs.size();
	vector<int> notBoundary;
	for(int i=0;i<n;i++)
	{
		cvline_polar outLine;
		fPoint line[2];
		line[0].x = lSegs[i][0].x;
		line[0].y = lSegs[i][0].y;
		line[1].x = lSegs[i][1].x;
		line[1].y = lSegs[i][1].y;

		outLine = line2polar(line);
		if (isBoundary(inputImg, outLine, cvScalar((0,0,0))))
		{
			continue;
		}
		notBoundary.push_back(i);
		double len = segment_length(outLine);
	}

	for(int m=0;m<(int)notBoundary.size();m++)
	{
		int i = notBoundary[m];
		{

			cvline_polar outLine;
			fPoint line[2];
			line[0].x = lSegs[i][0].x;
			line[0].y = lSegs[i][0].y;
			line[1].x = lSegs[i][1].x;
			line[1].y = lSegs[i][1].y;


			outLine = line2polar(line);
			vec_lines.push_back(outLine);
		}
	}

	//IplImage* black= cvCreateImage(procSize, IPL_DEPTH_8U, 1);
	IplImage* line_img = cvLoadImage(strImage, 1 );
	for(int i=0;i<(int)vec_lines.size();++i)
	{
		cvLine(line_img,cvPoint(vec_lines[i].pt1.x,vec_lines[i].pt1.y),cvPoint(vec_lines[i].pt2.x,vec_lines[i].pt2.y),CV_RGB(255,0,0),2, CV_AA);
	}

	return line_img;
}