#include "PointMatch.h"
#include "lineConstant.h"
#include <levmar.h>
#include <misc.h>
#include "ecm.h"
//#include "RANSAC_nD_Plane.h"
////#include <cv.h>
////#include <cxcore.h>
////#include <highgui.h>

PointMatch::PointMatch()
{
}

fPoint PointMatch::transModel(double x, double y)
{
	// 0	 1	 2	 3	 4	 6
	//a0	b0	a1	b1	a2	b2
	double xx = modelParameters[0] + modelParameters[1] * x + modelParameters[2] * y;
	double yy = modelParameters[3] + modelParameters[4] * x + modelParameters[5] * y;
	return fPoint(xx, yy);
}

fPoint PointMatch::transModel(fPoint pt)
{
	// 0	 1	 2	 3	 4	 6
	//a0	b0	a1	b1	a2	b2
	double xx = modelParameters[0] + modelParameters[1] * pt.x + modelParameters[2] * pt.y;
	double yy = modelParameters[3] + modelParameters[4] * pt.x + modelParameters[5] * pt.y;
	return fPoint(xx, yy);
}

fPoint PointMatch::pointPairDistance(PointPair pointPair)
{
	fPoint newDpt = transModel(pointPair.point);
	return fPoint(pointPair.point_prime.x-newDpt.x, pointPair.point_prime.y-newDpt.y);
}

double PointMatch::huber(double v, double c)
{
	//if(fabs(v) < c)
	//{
	//	return v*v;
	//}
	//else
	//{
	//	return c*(2*fabs(v) - c);
	//}
	double sign = v / (fabs(v) + DBL_EPSILON);
	if(fabs(v) < c)
	{
		return v;
	}
	else
	{
		//return sign*c*(2*fabs(v) - c);
		return sign*c;
	}
}

void PointMatch::funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	DataStruct *pDtStruct = (DataStruct*)adata;
	PointMatch *pModel = pDtStruct->pThis;
	//vector<cvline_polar> src_lines = pModel->m_src_lines;
	//vector<cvline_polar> ref_lines = pModel->m_ref_lines;
	const vector<PointPair>& pointPairList = pDtStruct->pointPairList;

	int nobs = static_cast<int>(pointPairList.size());
	int np = nparameter;

	int pos = 0;
	int i;
	for(i=0;i<np;++i)
	{
		pModel->modelParameters[i] = param[i]; //do not update right now
	}
	//pModel->setParameter(1, param[0]);
	//pModel->setParameter(2, param[1]);
	//pModel->setParameter(4, param[2]);
	//pModel->setParameter(5, param[3]);

	// º∆À„»®æÿ’Û
	//double max_dist = 0.0;
	double sigma = 0.0;
	for(i = 0;i < nobs;++i)
	{
		fPoint res = pModel->pointPairDistance(pointPairList[i]);
		sigma += res.x*res.x + res.y*res.y;
	}
	sigma = sqrt(sigma / (nobs*2 - 1));
	double c = 2 * sigma;
	//dlevmar_chol();

	for(i = 0;i < nobs;++i)
	{
		//fPoint line[2];
		//line[0] = pModel->transModel(ref_lines[i].x1, ref_lines[i].y1);
		//line[1] = pModel->transModel(ref_lines[i].x2, ref_lines[i].y2);

		//cvline_polar newLine = line2polar(line);
		fPoint res = pModel->pointPairDistance(pointPairList[i]);
		//double dis = sqrt(res.x*res.x + res.y*res.y);
		//double wght = 1.0 - dis / max_dist;
		//if(wght < 0.5) wght = 0.0;
		hx[pos++] = -huber(res.x, c);//-res.x;
		hx[pos++] = -huber(res.y, c);//-res.y;
		//hx[pos++] = -res.x;
		//hx[pos++] = -res.y;

		//fPoint dis = minDistanceFromLineset(newLine, src_lines);
		//fPoint dis = minPolarDistanceFromLineset(newLine, src_lines);
		//hx[c++] = -dis.x;
		//hx[c++] = -dis.y;
		//hx[c++] = -minPolarDistanceFromLineset(newLine, src_lines);
	}
}

fPoint PointMatch::getForwardDeriv( int p , fPoint pt, double pstep_scale)
{
	static double den = 0.5 / pstep_scale;
	fPoint pt1, pt2;
	double middle = modelParameters[p];
	modelParameters[p] = middle + pstep_scale;
	pt1 = transModel(pt.x, pt.y);
	modelParameters[p] = middle - pstep_scale;
	pt2 = transModel(pt.x, pt.y);
	double x = pt1.x - pt2.x;
	double y = pt1.y - pt2.y;
	x = x*den;
	y = y*den;
	modelParameters[p] = middle;
	return fPoint(x,y);
}

void PointMatch::jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata)
{
	DataStruct *pDtStruct = (DataStruct*)adata;
	PointMatch *pModel = pDtStruct->pThis;
	const vector<PointPair>& pointPairList = pDtStruct->pointPairList;

	int nobs = static_cast<int>(pointPairList.size());
	int np = nparameter;

	int pos = 0;
	int i;
	for(i=0;i<np;++i)
	{
		pModel->modelParameters[i] = param[i]; //do not update right now
	}

	vector<fPoint> imDerp(np);
	double pstep_scale = 1e-4;
	int c = 0;
	for(i = 0;i < nobs;++i)
	{
		fPoint res = pModel->pointPairDistance(pointPairList[i]);
		//double dis = sqrt(res.x*res.x + res.y*res.y);
		//double wght = 1.0 - dis / max_dist;
		//if(wght < 0.5) wght = 0.0;
		//hx[pos++] = -huber(res.x, c);//-res.x;
		//hx[pos++] = -huber(res.y, c);//-res.y;
		for(int p=0;p<np;++p)
		{
			imDerp[p] = pModel->getForwardDeriv( p , pointPairList[i].point, pstep_scale);
		}
		for(int p=0;p<np;++p)
		{
			j[c++] = imDerp[p].x;
		}
		for(int p=0;p<np;++p)
		{
			j[c++] = imDerp[p].y;
		}
	}

}

void PointMatch::registration(const vector<fPoint>& src_points, const vector<fPoint>& ref_points)
{
	modelParameters.clear();
	modelParameters.push_back(0.0f);
	modelParameters.push_back(1.0f);
	modelParameters.push_back(0.0f);
	modelParameters.push_back(0.0f);
	modelParameters.push_back(0.0f);
	modelParameters.push_back(1.0f);

	double convergence_epsilon = 1.0e-6;
	double delta = 0.0;
	vector<double> oldModelParameters = modelParameters;
	int max_iteration = 20;
	int iIter = 0;
	do 
	{
		vector<PointPair> pointPairList;
		for(int i = 0;i < static_cast<int>(ref_points.size());++i)
		{
			PointPair pointPair = findNearestPointPair(ref_points[i], src_points);
			pointPairList.push_back(pointPair);
		}
		findBestModel(pointPairList);
		delta = 0.0;
		for(int i = 0;i < static_cast<int>(modelParameters.size());++i)
		{
			delta += (modelParameters[i] - oldModelParameters[i]) * (modelParameters[i] - oldModelParameters[i]);
		}
		delta = sqrt(delta) / modelParameters.size();
		oldModelParameters = modelParameters;
		iIter++;
	} while (delta > convergence_epsilon && iIter < max_iteration);
}

void PointMatch::findBestModel(std::vector<PointPair> pointPairList)
{
	// RANSAC detect outliers
	auto_ptr< estimators::Solver<mat,vec> > ptrSolver(
		new estimators::affineSolver<mat,vec>);

	//-- Create input data
	mat dataPoints(pointPairList.size(), 4);// = "0 0; 1 1; 2 2; 3 3";
	for(int i = 0;i < (int)pointPairList.size()-1;++i)
	{
		dataPoints(i, 0) = pointPairList[i].point.x;
		dataPoints(i, 1) = pointPairList[i].point.y;
		dataPoints(i, 2) = pointPairList[i].point_prime.x;
		dataPoints(i, 3) = pointPairList[i].point_prime.y;
	}

	vector<int> inliers;
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
		0.95 // the confidence want to achieve at the end
		);

	modelParameters.clear();
	for (int i = 0;i < (int)models[0].n_rows;++i)
	{
		modelParameters.push_back(models[0](i));
	}
	drawResult(pointPairList, inliers, "pic\\finished1.png");

	// LM Method
	//modelParameters.clear();
	//modelParameters.push_back(0.0f);
	//modelParameters.push_back(1.0f);
	//modelParameters.push_back(0.0f);
	//modelParameters.push_back(0.0f);
	//modelParameters.push_back(0.0f);
	//modelParameters.push_back(1.0f);

	//int nparam = modelParameters.size();
	//int nobs = pointPairList.size() * 2;
	////get current adjustment (between -1 and 1 normally) and convert to ColumnVector
	//std::vector<double> cparm(nparam);

	//for(int i=0;i<nparam;++i)
	//{
	//	cparm[i] = modelParameters[i];
	//}

	//double *p = &cparm[0];
	//double *x = new double[nobs];
	//for(int i=0; i<nobs; i++) x[i]=0.0;

	//double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	//opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	//opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	////build Least Squares initial normal equation
	//// don't waste memory, add samples one at a time

	////double errors;
	////dlevmar_chkjac(funcErrorEquation, jacErrorEquation, p, nparam, nobs, &optStruct, &errors);

	///* Rosenbrock function */
	////int ret = dlevmar_der(funcErrorEquation, jacErrorEquation, p, x, nparam, nobs, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian

	//DataStruct dtStruct;
	//dtStruct.pThis = this;
	//dtStruct.pointPairList = pointPairList;
	////int ret = dlevmar_dif(funcErrorEquation, p, x, nparam, nobs, 1000, opts, info, NULL, NULL, &dtStruct);  // no Jacobian
	//int ret = dlevmar_der(funcErrorEquation, jacErrorEquation, p, x, nparam, nobs, 1000, opts, info, NULL, NULL, &dtStruct); // with analytic Jacobian
}

void PointMatch::registration(std::vector<PointPair> pointPairList)
{
	modelParameters.clear();
	modelParameters.push_back(0.0f);
	modelParameters.push_back(1.0f);
	modelParameters.push_back(0.0f);
	modelParameters.push_back(0.0f);
	modelParameters.push_back(0.0f);
	modelParameters.push_back(1.0f);

	int nparam = modelParameters.size();
	int nobs = pointPairList.size() * 2;
	//get current adjustment (between -1 and 1 normally) and convert to ColumnVector
	std::vector<double> cparm(nparam);

	for(int i=0;i<nparam;++i)
	{
		cparm[i] = modelParameters[i];
	}

	double *p = &cparm[0];
	double *x = new double[nobs];
	for(int i=0; i<nobs; i++) x[i]=0.0;

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	//build Least Squares initial normal equation
	// don't waste memory, add samples one at a time

	//double errors;
	//dlevmar_chkjac(funcErrorEquation, jacErrorEquation, p, nparam, nobs, &optStruct, &errors);

	/* Rosenbrock function */
	//int ret = dlevmar_der(funcErrorEquation, jacErrorEquation, p, x, nparam, nobs, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian

	DataStruct dtStruct;
	dtStruct.pThis = this;
	dtStruct.pointPairList = pointPairList;
	//int ret = dlevmar_dif(funcErrorEquation, p, x, nparam, nobs, 1000, opts, info, NULL, NULL, &dtStruct);  // no Jacobian
	int ret = dlevmar_der(funcErrorEquation, jacErrorEquation, p, x, nparam, nobs, 1000, opts, info, NULL, NULL, &dtStruct); // with analytic Jacobian
}

PointPair PointMatch::findNearestPointPair(const fPoint &basePoint, const vector<fPoint> &pointList)
{
	int bestIndex = 0;
	double nearestDistance = DBL_MAX;
	for(int i = 0;i < static_cast<int>(pointList.size());++i)
	{
		PointPair pointPair;
		pointPair.point = basePoint;
		pointPair.point_prime = pointList[i];

		fPoint dis = pointPairDistance(pointPair);
		double d = sqrt(dis.x*dis.x + dis.y*dis.y);
		if(d < nearestDistance)
		{
			nearestDistance = d;
			bestIndex = i;
		}
	}

	PointPair pointPair;
	pointPair.point = basePoint;
	pointPair.point_prime = pointList[bestIndex];

	return pointPair;
}

void PointMatch::drawResult(const std::vector<PointPair>& pointPairList, const vector<int>& inliers, const char* szFilename)
{
	// plot
	// for debug
	double minX = DBL_MAX;
	double minY = DBL_MAX;
	double maxX = 0;
	double maxY = 0;
	int num = (int)inliers.size();
	for (int i = 0;i < num;++i)
	{
		minX = (minX > pointPairList[inliers[i]].point.x) ? pointPairList[inliers[i]].point.x : minX;
		minY = (minY > pointPairList[inliers[i]].point.y) ? pointPairList[inliers[i]].point.y : minY;
		minX = (minX > pointPairList[inliers[i]].point_prime.x) ? pointPairList[inliers[i]].point_prime.x : minX;
		minY = (minY > pointPairList[inliers[i]].point_prime.y) ? pointPairList[inliers[i]].point_prime.y : minY;
		maxX = (maxX < pointPairList[inliers[i]].point.x) ? pointPairList[inliers[i]].point.x : maxX;
		maxY = (maxY < pointPairList[inliers[i]].point.y) ? pointPairList[inliers[i]].point.y : maxY;
		maxX = (maxX < pointPairList[inliers[i]].point_prime.x) ? pointPairList[inliers[i]].point_prime.x : maxX;
		maxY = (maxY < pointPairList[inliers[i]].point_prime.y) ? pointPairList[inliers[i]].point_prime.y : maxY;
	}

	int offset = 100;
	int offset2 = 2.0*offset;
	int nwidth = (int)(maxX - minX + 0.5) + offset2;
	int nheight = (int)(maxY - minY + 0.5) + offset2;
	nwidth = (nwidth < offset2) ? offset2 : nwidth;
	nheight = (nheight < offset2) ? offset2 : nheight;
	// create a pic
	IplImage* pic = cvCreateImage( cvSize(nwidth, nheight), 8, 3 );

	int nCount = 0;
	for(int i = 1; i < num;i++)
	{
		fPoint transedPt = transModel(pointPairList[inliers[i]].point);
		cvCircle(pic, cvPoint(cvRound(pointPairList[inliers[i]].point_prime.x) + offset, cvRound(pointPairList[inliers[i]].point_prime.y) + offset), 5, CV_RGB(255, 0, 0));
		cvCircle(pic, cvPoint(cvRound(transedPt.x) + offset, cvRound(transedPt.y) + offset), 5, CV_RGB(0, 255, 0));
		cvLine(pic, cvPoint(cvRound(pointPairList[inliers[i]].point_prime.x) + offset, cvRound(pointPairList[inliers[i]].point_prime.y) + offset),
			cvPoint(cvRound(transedPt.x) + offset, cvRound(transedPt.y) + offset), CV_RGB(255,255,255));
		nCount++;
	}
	cvSaveImage(szFilename, pic);
}