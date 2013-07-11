#include "ecm.h"

bool ecmpr::solve()
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;
	initialization();

	arma::vec oldParameters = m_parameters;
	double convergence_epsilon = 2.0e-4;
	double delta;
	int max_iteration = 200;
	int iTime = 0;
	do 
	{
		double minX = DBL_MAX;
		double minY = DBL_MAX;
		double maxX = 0;
		double maxY = 0;
		//if(iTime > 0)
		//{
		//	classification();
		//	parameters_refine();
		//}

		for (int i = 0;i < nX;++i)
		{
			arma::vec transedPt = forward(m_modelData.row(i).t());
			m_transData.row(i) = transedPt.t();
			minX = (minX > transedPt[0]) ? transedPt[0] : minX;
			minY = (minY > transedPt[1]) ? transedPt[1] : minY;
			maxX = (maxX < transedPt[0]) ? transedPt[0] : maxX;
			maxY = (maxY < transedPt[1]) ? transedPt[1] : maxY;
		}
		for (int j = 0;j < nY;++j)
		{
			minX = (minX > m_observationData(j, 0)) ?  m_observationData(j, 0) : minX;
			minY = (minY >  m_observationData(j, 1)) ? m_observationData(j, 1) : minY;
			maxX = (maxX <  m_observationData(j, 0)) ?  m_observationData(j, 0) : maxX;
			maxY = (maxY < m_observationData(j, 1)) ? m_observationData(j, 1) : maxY;
		}
		// plot
		// for debug
		if(iTime%1 == 0)
		{
			char szFilename[56];
			sprintf_s(szFilename, "pic\\%03d.png", iTime+1);
			int offset = 100;
			int nwidth = (int)(maxX - minX + 0.5) + 2 * offset;
			int nheight = (int)(maxY - minY + 0.5) + 2 * offset;
			// create a pic
			IplImage* pic = cvCreateImage( cvSize(nwidth, nheight), 8, 3 );
			for(int i = 1; i < nX;i++)
			{
				//cvPoint(cvRound(m_transData(i, 0)) + offset, cvRound(m_transData(i, 1)) + offset, CV_RGB(255, 0, 0) );
				cvCircle(pic, cvPoint(cvRound(m_transData(i, 0)) + offset, cvRound(m_transData(i, 1)) + offset), 3, CV_RGB(255, 0, 0));
			}
			for(int j = 1; j < nY;j++)
			{
				//cvPoint(cvRound(m_observationData(j, 0)) + offset, cvRound(m_observationData(j, 1)) + offset, CV_RGB(0, 255, 0) );
				cvCircle(pic, cvPoint(cvRound(m_observationData(j, 0)) + offset, cvRound(m_observationData(j, 1)) + offset), 3, CV_RGB(0, 255, 0));
			}
			cvSaveImage(szFilename, pic);


			if (fexists(m_sourceImageFile))
			{
				char szImageFilename[56];
				sprintf_s(szImageFilename, "pic\\image_%03d.png", iTime+1);
				vector<double> param(get_PARAMETER_NUM());
				for (int i = 0;i < get_PARAMETER_NUM();++i) param[i] = oldParameters[i];
				affine_warp_Image< vector<double> >(m_sourceImageFile, szImageFilename, param);
			}
		}

		double tmp = m_delta_2;
		updateDelta_2(iTime == 0);
		cout<<"iteration:\t"<<iTime+1<<endl;
		cout<<"rms:\t"<<m_delta_2<<endl;
		cout<<oldParameters<<endl;
		if (m_delta_2 < 1.0)
		{
			break;
			cout<<"rms epsilon reached!"<<endl;
		}
		if (fabs(m_delta_2 - tmp) < 0.5)
		{
			break;
			cout<<"delta rms epsilon reached!"<<endl;
		}
		// E-step
		E_Step();

		// CM-step A
		// estimate the trans parameters
		parameters_optimization();
		delta = norm2_Vec(oldParameters - m_parameters);
		oldParameters = m_parameters;
		// CM-step B
		// estimate the covariance matrices
		//updateCovariances();
		//m_delta_2 *= 0.9;
		//cout<<m_parameters<<endl;
	} while (++iTime < max_iteration && fabs(delta) > convergence_epsilon);
	if(iTime == max_iteration)
	{
		cout<<"Warning: the max number of iteration reached, and the results may be not accurate."<<endl;
	}

	//// classification
	classification();

	// plot
	// for debug
	double minX = DBL_MAX;
	double minY = DBL_MAX;
	double maxX = 0;
	double maxY = 0;
	for (int i = 0;i < nX;++i)
	{
		minX = (minX > m_modelData(i, 0)) ? m_modelData(i, 0) : minX;
		minY = (minY > m_modelData(i, 1)) ? m_modelData(i, 1) : minY;
		maxX = (maxX < m_modelData(i, 0)) ? m_modelData(i, 0) : maxX;
		maxY = (maxY < m_modelData(i, 1)) ? m_modelData(i, 1) : maxY;
	}
	for (int j = 0;j < nY;++j)
	{
		minX = (minX > m_observationData(j, 0)) ?  m_observationData(j, 0) : minX;
		minY = (minY >  m_observationData(j, 1)) ? m_observationData(j, 1) : minY;
		maxX = (maxX <  m_observationData(j, 0)) ?  m_observationData(j, 0) : maxX;
		maxY = (maxY < m_observationData(j, 1)) ? m_observationData(j, 1) : maxY;
	}
	char szFilename[56];
	sprintf_s(szFilename, "pic\\finished.png");
	int offset = 100;
	int nwidth = (int)(maxX - minX + 0.5) + 2 * offset;
	int nheight = (int)(maxY - minY + 0.5) + 2 * offset;
	// create a pic
	IplImage* pic = cvCreateImage( cvSize(nwidth, nheight), 8, 3 );

	for(int i = 1; i < nX;i++)
	{
		//cvPoint(cvRound(m_transData(i, 0)) + offset, cvRound(m_transData(i, 1)) + offset, CV_RGB(255, 0, 0) );
		cvCircle(pic, cvPoint(cvRound(m_transData(i, 0)) + offset, cvRound(m_transData(i, 1)) + offset), 3, CV_RGB(255, 0, 0));
	}
	for(int j = 1; j < nY;j++)
	{
		//cvPoint(cvRound(m_observationData(j, 0)) + offset, cvRound(m_observationData(j, 1)) + offset, CV_RGB(0, 255, 0) );
		cvCircle(pic, cvPoint(cvRound(m_observationData(j, 0)) + offset, cvRound(m_observationData(j, 1)) + offset), 3, CV_RGB(0, 255, 0));
	}
	int nCount = 0;
	for(int j = 1; j < nY;j++)
	{
		if(m_z[j] == nX) continue;
		//cvPoint(cvRound(m_observationData(j, 0)) + offset, cvRound(m_observationData(j, 1)) + offset, CV_RGB(0, 255, 0) );
		//cvCircle(pic, cvPoint(cvRound(m_transData(m_z[j], 0)) + offset, cvRound(m_transData(m_z[j], 1)) + offset), 2, CV_RGB(255, 0, 0));
		//cvCircle(pic, cvPoint(cvRound(m_observationData(j, 0)) + offset, cvRound(m_observationData(j, 1)) + offset), 2, CV_RGB(0, 255, 0));
		cvLine(pic, cvPoint(cvRound(m_transData(m_z[j], 0)) + offset, cvRound(m_transData(m_z[j], 1)) + offset), cvPoint(cvRound(m_observationData(j, 0)) + offset, cvRound(m_observationData(j, 1)) + offset), CV_RGB(255,255,255));
		nCount++;
	}
	cvSaveImage(szFilename, pic);

	cout<<"EM based point set registration finished after "<<iTime<<" times of iterations."<<endl;
	cout<<"********************************************************************************"<<endl<<endl;
	return true;
}

void ecmpr::updateDelta_2(bool bfirstTime)
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	double sum = 0.0;
	for(int i = 0;i < nX;++i)
	{
		m_delta_2List[i] = 0.0;
		for (int j = 0;j < nY;++j)
		{
			double dist = arma::norm(m_transData.row(i) - m_observationData.row(j), 2);
			double tmp = dist * dist;
			sum += tmp * m_alpha(j,i);
			m_delta_2List[i] += tmp * m_alpha(j,i);
		}
		m_delta_2List[i] /= m_lambda[i] * ndim;
	}
	//if (bfirstTime) m_delta_2 = sum / (nX * nY * ndim);
	//else m_delta_2 = sum / (nY * ndim);
	m_delta_2 = sum / (m_lammda_sum * ndim);
}


bool ecmpr::classification()
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	// classification
	for (int j = 0;j < nY;++j)
	{
		int maxIndex = 0;
		double maxValue = 0.0;
		for (int i = 0;i < nX + 1;++i)
		{
			if(m_alpha(j, i) > maxValue)
			{
				maxValue = m_alpha(j, i);
				maxIndex = i;
			}
		}
		m_z[j] = maxIndex;
	}
	return true;
}

bool ecmpr::parameters_refine()
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	// RANSAC detect outliers
	auto_ptr< estimators::Solver<mat,vec> > ptrSolver(
		new estimators::affineSolver<mat,vec>);

	//-- Create input data
	int nCount = 0;
	for (int j = 0;j < nY;++j)
	{
		if(m_z[j] < nX) nCount++;
	}
	mat dataPoints(nCount, 4);// = "0 0; 1 1; 2 2; 3 3";
	nCount = 0;
	for(int j = 0;j < nY;++j)
	{
		if(m_z[j] == nX) continue;
		dataPoints(nCount, 0) = m_modelData(m_z[j], 0);
		dataPoints(nCount, 1) = m_modelData(m_z[j], 1);
		dataPoints(nCount, 2) = m_observationData(j, 0);
		dataPoints(nCount, 3) = m_observationData(j, 1);
		nCount++;
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
	m_parameters = models[0];
	return true;
}

bool ecmpr::initialization()
{
	initialParameters();
	
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;
	// initialize the covariance matrices
	for (int i = 0;i < nX;++i)
	{
		arma::mat covariance(ndim, ndim);
		covariance.eye();
		m_covariances.push_back(covariance);
	}

	//// initialize delta
	//m_delta_2 = 0.0;
	//for (int i = 0;i < nX;++i)
	//{
	//	arma::vec transedPt = forward(m_modelData.row(i).t());
	//	for (int j = 0;j < nY;++j)
	//	{
	//		double dist = arma::norm(transedPt.t() - m_observationData.row(j), 2);
	//		m_delta_2 += dist * dist;
	//	}
	//}
	//m_delta_2 /= ndim * nX * nY;

	//set the size of posterior matrix
	m_alpha.set_size(nY, nX + 1);
	double w = 1.0 / (nX + 1);
	m_alpha.fill(w);

	// set the size of virtual observation matrix
	m_w.set_size(nX, ndim);

	// set the weight vector of virtual observations
	m_lambda.set_size(nX);


	m_w.fill(0.0);
	m_lammda_sum = 0.0;
	for (int i = 0;i < nX;++i)
	{
		m_lambda[i] = 0;
		for (int j = 0;j < nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			m_w.row(i) += m_alpha(j, i) * m_observationData.row(j);
		}
		//m_w.row(i) /= (m_lambda[i] + FLOAT_EPSILON);
		m_lammda_sum += m_lambda[i];
	}

	//
	m_z.resize(nY);
	//m_z.set_size(nY);

	//
	m_transData.set_size(nX, ndim);

	//
	m_delta_2List.resize(nX, 0.0);

	return true;
}

bool ecmpr::E_Step()
{
	// E-step evaluate the posteriors
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	//double outlier = 2 / (r * r + FLOAT_EPSILON);
	//double outlier = pow(2.0*PI*m_delta_2, ndim / 2) * nX / nY;
	//double outlier_dist = 20.0;
	//double r = 12.0;
	//double r = 3.0 * sqrt(m_delta_2);

	double outlier_dist = 3.0 * sqrt(m_delta_2);
	//double r = 3.8;
	//double r = 16.0;
	//double r = max(outlier_dist, 5.0);
	double r = outlier_dist;
	double beta = -log(0.03) / log(outlier_dist + 1.0);
	//double outlier = 2.0 * sqrt(m_delta_2) / (r * r + DBL_EPSILON);
	//double outlier = 2.0 / (r * r + DBL_EPSILON);
	//double outlier = 2.0 * m_delta_2 / (r * r + DBL_EPSILON);
	double outlier = 2.0 / 9.0;
	for(int j = 0;j < nY;++j)
	{
		double sum = 0.0;
		vector<double> tmp(nX);
		double minDist = DBL_MAX;
		for (int i = 0;i < nX;++i)
		{
			//r = 3.0 * sqrt(m_delta_2List[i]);
			//outlier = 2.0 * sqrt(m_delta_2List[i]) / (r * r + DBL_EPSILON);
			//arma::vec calcPos = forward(m_modelData.row(i).t());
			//double dist = MahalanobisDistance(m_observationData.row(j).t(), calcPos, m_covariances[i]);
			double dist = arma::norm(m_observationData.row(j).t() - m_transData.row(i).t(), 2);
			tmp[i] = exp(-dist*dist / (2.0 * m_delta_2+ DBL_EPSILON));///sqrt(det(m_covariances[i]));
			
			if(dist < minDist) minDist = dist;
			sum += tmp[i];
		}
		//if(minDist < outlier_dist)
		//{
		//	// inlier
		//	//outlier = 0.0;
		//	//outlier = 2.0 / (r * r + DBL_EPSILON) * exp(-(outlier_dist-minDist)*(outlier_dist-minDist)/(minDist*minDist+DBL_EPSILON));
		//	//outlier = 2.0 / (r * r + DBL_EPSILON);// * exp(-(outlier_dist-minDist)*(outlier_dist-minDist)/(minDist*minDist+DBL_EPSILON));
		//	//outlier = exp(-(outlier_dist-minDist)*(outlier_dist-minDist)/(minDist*minDist+DBL_EPSILON));
		//	//outlier = exp(-(outlier_dist-minDist)*(outlier_dist-minDist)/(minDist*minDist+DBL_EPSILON));
		//	//outlier = exp(-(outlier_dist-minDist)/(minDist+DBL_EPSILON) / (2.0*m_delta_2));
		//	//outlier = 2.0 / (r * r + DBL_EPSILON) * (1.0 - exp(-minDist / (double)(outlier_dist - minDist + DBL_EPSILON)  / (2.0*m_delta_2)));
		//	outlier = 2.0 / (r * r + DBL_EPSILON) * (1.0 - exp(-minDist*minDist / (double)(2.0*(outlier_dist - minDist)*(outlier_dist - minDist) + DBL_EPSILON)));
		//	//outlier = minDist / outlier_dist;
		//}
		//else
		//{
		//	// outlier
		//	//outlier = pow(minDist, -2.0) * sqrt(2.0*PI*m_delta_2);
		//	//outlier = 1.0 * sqrt(2.0*PI*m_delta_2);
		//	//outlier = (1.0 - pow(outlier_dist, -beta));
		//	outlier = 2.0 / (r * r + DBL_EPSILON);
		//	//outlier = 1.0;
		//}
		//outlier = (1.0 - pow((minDist + 1.0), -beta));
		//outlier = beta * pow((minDist + 1.0), -beta - 1) * sqrt(2.0*PI*m_delta_2);
		//outlier = (1.0 - maxTmp);
		sum += outlier;
		for(int i = 0;i < nX;++i)
		{
			m_alpha(j, i) = tmp[i] / (sum + DBL_EPSILON);
		}
		m_alpha(j, nX) = outlier / (sum + DBL_EPSILON);
	}

	//// calculate the weights each observations
	//vector<double> weight(nY);
	//double minOutline = arma::min(m_alpha.col(nX));
	//double mu = arma::mean(m_alpha.col(nX));
	////double rms = 0.0;
	////for (int j = 0;j < nY;++j)
	////{
	////	rms += (m_alpha(j, nX) - mu) * (m_alpha(j, nX) - mu);
	////	//weight[j] = exp(1.0 - m_alpha(j, nX) / (minOutline + FLOAT_EPSILON));
	////	//weight[j] = minOutline / (m_alpha(j, nX) + FLOAT_EPSILON);
	////}
	////rms = sqrt((rms / (nY - 1)));
	//double c = 1.0 * mu;
	//for (int j = 0;j < nY;++j)
	//{
	//	if(fabs(m_alpha(j, nX)) < c)
	//	{
	//		weight[j] = 1.0;//m_alpha(j, nX);
	//	}
	//	else
	//	{
	//		weight[j] = c/m_alpha(j, nX);//c*(2*fabs(m_alpha(j, nX)) - c);
	//	}
	//}

	// calculate the weights of virtual observations
	// calculate the virtual observation matrix
	m_w.fill(0.0);
	m_lammda_sum = 0.0;
	for (int i = 0;i < nX;++i)
	{
		m_lambda[i] = 0;
		for (int j = 0;j < nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			m_w.row(i) += m_alpha(j, i) * m_observationData.row(j);
		}
		m_lammda_sum += m_lambda[i];
		//m_w.row(i) /= (m_lambda[i] + FLOAT_EPSILON);
	}
	return true;
}

void ecmpr::funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	ecmpr *pThis = (ecmpr*)adata;
	int nX = (int)pThis->m_modelData.n_rows;
	int nY = (int)pThis->m_observationData.n_rows;
	int ndim = (int)pThis->m_modelData.n_cols;

	int pos = 0;
	int i;
	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
	}

	for(i = 0;i < nequation;++i)
	{
		arma::vec outPt = pThis->forward(pThis->m_modelData.row(i).t());
		double dist = MahalanobisDistance(pThis->m_w.row(i).t(), outPt, pThis->m_covariances[i]);
		double res = sqrt(pThis->m_lambda[i] * dist);
		hx[pos++] = res*100.0;
	}
	//for(i = 0;i < nX;++i)
	//{
	//	arma::vec outPt = pThis->forward(pThis->m_modelData.row(i).t());
	//	for (int j = 0;j < nY;++j)
	//	{
	//		double dist = MahalanobisDistance(pThis->m_observationData.row(j).t(), outPt, pThis->m_covariances[i]);
	//		double res = sqrt(pThis->m_alpha(j, i) * dist);
	//		hx[pos++] = res;
	//	}
	//}
}
void ecmpr::jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata)
{
	ecmpr *pThis = (ecmpr*)adata;

	int nX = (int)pThis->m_modelData.n_rows;
	int nY = (int)pThis->m_observationData.n_rows;
	int ndim = (int)pThis->m_modelData.n_cols;

	int pos = 0;
	int i;
	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
	}

	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;
	int c = 0;
	for(i = 0;i < nequation;++i)
	{
		for(int p=0;p<nparameter;++p)
		{
			double middle = pThis->m_parameters[p];
			pThis->m_parameters[p] = middle + pstep_scale;
			arma::vec outPt1 = pThis->forward(pThis->m_modelData.row(i).t());
			double dist1 = MahalanobisDistance(pThis->m_w.row(i).t(), outPt1, pThis->m_covariances[i]);
			double res1 = sqrt(pThis->m_lambda[i] * dist1);

			pThis->m_parameters[p] = middle - pstep_scale;
			arma::vec outPt2 = pThis->forward(pThis->m_modelData.row(i).t());
			double dist2 = MahalanobisDistance(pThis->m_w.row(i).t(), outPt2, pThis->m_covariances[i]);
			double res2 = sqrt(pThis->m_lambda[i] * dist2);

			j[c++] = (res1 - res2)*den*100.0;
			pThis->m_parameters[p] = middle;
		}
	}

}

bool ecmpr::parameters_optimization()
{

	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	int nparam = get_PARAMETER_NUM();
	std::vector<double> cparm(nparam);

	for(int i=0;i<nparam;++i)
	{
		cparm[i] = m_parameters[i];
	}

	arma::vec outParameters(nparam);
	double *p = &cparm[0];
	double *x = new double[nX];
	for(int i=0; i<nX; i++) x[i]=0.0;

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	//int ret = dlevmar_dif(funcErrorEquation, p, x, nparam, nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
	int ret = dlevmar_der(funcErrorEquation, jacErrorEquation, p, x, nparam, nX, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian
	
	return true;
}

bool ecmpr::updateCovariances()
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	for(int i = 0;i < nX;++i)
	{
		m_covariances[i].fill(0.0);
		for(int j = 0;j < nY;++j)
		{
			arma::vec calcPos = forward(m_modelData.row(i).t());
			arma::rowvec delta = m_observationData.row(j) - calcPos.t();
			m_covariances[i] += m_alpha(j, i) * delta.t() * delta;
		}
		m_covariances[i] /= (m_lambda[i] + DBL_EPSILON);
	}
	
	return true;
}