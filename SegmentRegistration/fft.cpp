#include "fft.h"

void fft2(IplImage *src, IplImage *dst)
{   //ʵ�����鲿
	IplImage *image_Re = 0, *image_Im = 0, *Fourier = 0;
	//   int i, j;
	image_Re = cvCreateImage(cvGetSize(src), IPL_DEPTH_64F, 1);  //ʵ��
	//Imaginary part
	image_Im = cvCreateImage(cvGetSize(src), IPL_DEPTH_64F, 1);  //�鲿
	//2 channels (image_Re, image_Im)
	Fourier = cvCreateImage(cvGetSize(src), IPL_DEPTH_64F, 2);
	// Real part conversion from u8 to 64f (double)
	cvConvertScale(src, image_Re);
	// Imaginary part (zeros)
	cvZero(image_Im);
	// Join real and imaginary parts and stock them in Fourier image
	cvMerge(image_Re, image_Im, 0, 0, Fourier);

	// Application of the forward Fourier transform
	cvDFT(Fourier, dst, CV_DXT_FORWARD);
	cvReleaseImage(&image_Re);
	cvReleaseImage(&image_Im);
	cvReleaseImage(&Fourier);
}

void fft2shift(IplImage *src, IplImage *dst)
{
	IplImage *image_Re = 0, *image_Im = 0;
	int nRow, nCol, i, j, cy, cx;
	double scale, shift, tmp13, tmp24;
	image_Re = cvCreateImage(cvGetSize(src), IPL_DEPTH_64F, 1);
	//Imaginary part
	image_Im = cvCreateImage(cvGetSize(src), IPL_DEPTH_64F, 1);
	cvSplit( src, image_Re, image_Im, 0, 0 );
	//����ԭ���������˹����ͼ����p123
	// Compute the magnitude of the spectrum Mag = sqrt(Re^2 + Im^2)
	//���㸵��Ҷ��
	cvPow( image_Re, image_Re, 2.0);
	cvPow( image_Im, image_Im, 2.0);
	cvAdd( image_Re, image_Im, image_Re);
	cvPow( image_Re, image_Re, 0.5 );
	//�����任����ǿ�Ҷȼ�ϸ��(���ֱ任ʹ��խ���ͻҶ�����ͼ��ֵӳ��
	//һ������ֵ������ɼ�������˹����ͼ����p62)
	// Compute log(1 + Mag);
	cvAddS( image_Re, cvScalar(1.0), image_Re ); // 1 + Mag
	cvLog( image_Re, image_Re ); // log(1 + Mag)

	//Rearrange the quadrants of Fourier image so that the origin is at the image center
	nRow = src->height;
	nCol = src->width;
	cy = nRow/2; // image center
	cx = nCol/2;
	//CV_IMAGE_ELEMΪOpenCV����ĺ꣬������ȡͼ�������ֵ����һ���־��ǽ������ı任
	for( j = 0; j < cy; j++ ){
		for( i = 0; i < cx; i++ ){
			//���Ļ���������ݳ��Ŀ���жԽǽ���
			tmp13 = CV_IMAGE_ELEM( image_Re, double, j, i);
			CV_IMAGE_ELEM( image_Re, double, j, i) = CV_IMAGE_ELEM(
				image_Re, double, j+cy, i+cx);
			CV_IMAGE_ELEM( image_Re, double, j+cy, i+cx) = tmp13;

			tmp24 = CV_IMAGE_ELEM( image_Re, double, j, i+cx);
			CV_IMAGE_ELEM( image_Re, double, j, i+cx) =
				CV_IMAGE_ELEM( image_Re, double, j+cy, i);
			CV_IMAGE_ELEM( image_Re, double, j+cy, i) = tmp24;
		}
	}
	//��һ�����������Ԫ��ֵ��һΪ[0,255]
	//[(f(x,y)-minVal)/(maxVal-minVal)]*255
	double minVal = 0, maxVal = 0;
	// Localize minimum and maximum values
	cvMinMaxLoc( image_Re, &minVal, &maxVal );
	// Normalize image (0 - 255) to be observed as an u8 image
	scale = 255/(maxVal - minVal);
	shift = -minVal * scale;
	cvConvertScale(image_Re, dst, scale, shift);
	cvReleaseImage(&image_Re);
	cvReleaseImage(&image_Im);
}