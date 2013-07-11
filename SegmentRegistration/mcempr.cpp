#include "mcempr.h"
#include <time.h>

bool mcempr::solve()
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;
	initialization();
	
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
		if(iTime%5 == 0)
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
		}

		// E-step
		E_Step();

		// CM-step A
		// estimate the trans parameters
		arma::vec oldParameters = m_parameters;
		parameters_optimization();
		delta = norm2_Vec(oldParameters - m_parameters);

		// CM-step B
		// estimate the covariance matrices
		//updateCovariances();
		//m_delta_2 *= 0.9;
		cout<<m_parameters<<endl;
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

double mcempr::acceptance_ratio(const correspondence_struct& c, const correspondence_struct& c_prime)
{
	double result = 1.0;
	for (int j = 0;j < (int)c.forward_map.size();++j)
	{
		result *= (m_alpha(j, c_prime.forward_map[j]))/(m_alpha(j, c.forward_map[j]) + DBL_EPSILON);
	}
	return result;
}

correspondence_struct mcempr::getValidCorrespondence()
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	// valid correspondence vector
	// indices start from 0
	// nX indicates invisible or outlier
	correspondence_struct correspond;
	correspond.forward_map.resize(nY, nX);
	correspond.inverse_map.resize(nX, -1);

	// randomly choose the start index
	// initialize random seed:
	srand ( time(NULL) );
	// generate start index:
	int start_pos = rand() % nY;

	vector<int> candidate_list(nX+1);
	for (int i = 0;i < nX+1;++i) candidate_list[i] = i;

	for (int j = 0;j < nY;++j)
	{
		int real_pose = (j + start_pos) % nY;
		int maxIndex = 0;
		double maxValue = 0.0;
		if(candidate_list.size() < 2) break;
		for (int i = 0;i < (int)candidate_list.size();++i)
		{
			if(m_alpha(real_pose, i) > maxValue)
			{
				maxValue = m_alpha(real_pose, candidate_list[i]);
				maxIndex = i;
			}
		}
		correspond.forward_map[real_pose] = candidate_list[maxIndex];

		// if it is not outlier
		if(candidate_list[maxIndex] != nX)
		{
			correspond.inverse_map[candidate_list[maxIndex]] = real_pose;
			candidate_list.erase(candidate_list.begin()+maxIndex);
		}
	}
	return correspond;
}

bool mcempr::MetropolisHastings()
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	int nSample = 1000;
	int nDiscard = 50;
	correspondence_struct current_correspondence = getValidCorrespondence();

	arma::Mat<int> accumulateMatrix(nY, nX + 1);
	accumulateMatrix.fill(0);
	int iCount = 0;
	do 
	{
		//double f1 = getCorrepondenceProbility(current_correspondence);
		correspondence_struct new_correspondence = MarkovNext(current_correspondence);
		//double f2 = getCorrepondenceProbility(new_correspondence);
		double acceptance = acceptance_ratio(current_correspondence, new_correspondence);
		if(acceptance < 0.99999)
		{
			// reject
			continue;
		}
		else
		{
			// accept
			current_correspondence = new_correspondence;
			if (iCount >= nDiscard)
			{
				for (int k = 0;k < (int)current_correspondence.forward_map.size();++k)
				{
					accumulateMatrix(k, current_correspondence.forward_map[k]) += 1;
				}
			}
			iCount++;
		}
	} while (iCount < nSample);

	// update the corresponding weights
	int num = nSample - nDiscard;
	for (int j = 0;j < nY;++j)
	{
		for(int i = 0;i < nX;++i)
		{
			m_alpha(j, i) = accumulateMatrix(j, i) / (double)num;
		}
	}

	return true;
}

double mcempr::getCorrepondenceProbility(const correspondence_struct& correspond)
{
	double result = 1.0;
	for (int j = 0;j < (int)correspond.forward_map.size();++j)
	{
		result *= m_alpha(j, correspond.forward_map[j]);
	}
	return result;
}

correspondence_struct mcempr::MarkovNext(const correspondence_struct& correspond)
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	correspondence_struct newCorrespond = correspond;

	// randomly choose the start index
	// initialize random seed:
	srand ( time(NULL) );
	// generate start index:
	int oldSource = rand() % nY;

	vector<int> openList(nX + 1);
	for (int i = 0;i < nX+1;++i) openList[i] = i;
	do 
	{
		int maxIndex = 0;
		double maxValue = 0.0;
		for (int i = 0;i < (int)openList.size();++i)
		{
			if(m_alpha(oldSource, i) > maxValue)
			{
				maxValue = m_alpha(oldSource, openList[i]);
				maxIndex = i;
			}
		}

		int oldTarget = newCorrespond.forward_map[oldSource];
		int newTarget = openList[maxIndex];
		if(openList[maxIndex] == nX || newCorrespond.inverse_map[newTarget] == -1)
		{
			// outlier
			if (oldTarget != nX) newCorrespond.inverse_map[oldTarget] = -1;
			newCorrespond.forward_map[oldSource] = newTarget;
			break;
		}
		else
		{
			int newSource = newCorrespond.inverse_map[newTarget];
			if (oldTarget != nX) newCorrespond.inverse_map[oldTarget] = -1;
			newCorrespond.inverse_map[newTarget] = oldSource;
			newCorrespond.forward_map[oldSource] = newTarget;
			oldSource = newSource;

			openList.erase(openList.begin() + maxIndex);
		}
	} while ((int)openList.size() > 0);

	return newCorrespond;
}


bool mcempr::classification()
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

bool mcempr::parameters_refine()
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

bool mcempr::initialization()
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

	// initialize delta
	m_delta_2 = 0.0;
	for (int i = 0;i < nX;++i)
	{
		arma::vec transedPt = forward(m_modelData.row(i).t());
		for (int j = 0;j < nY;++j)
		{
			double dist = arma::norm(transedPt.t() - m_observationData.row(j), 2);
			m_delta_2 += dist * dist;
		}
	}
	m_delta_2 /= ndim * nX * nY;

	//set the size of posterior matrix
	m_alpha.set_size(nY, nX + 1);

	// set the size of virtual observation matrix
	m_w.set_size(nX, ndim);

	// set the weight vector of virtual observations
	m_lambda.set_size(nX);

	//
	m_z.resize(nY);
	//m_z.set_size(nY);

	//
	m_transData.set_size(nX, ndim);

	return true;
}

bool mcempr::E_Step()
{
	// E-step evaluate the posteriors
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	double r = 50.0;
	//double outlier = 2 / (r * r + FLOAT_EPSILON);
	//double outlier = pow(2.0*PI*m_delta_2, ndim / 2) * nX / nY;
	double outlier = 2.0 * sqrt(m_delta_2) / (r * r + DBL_EPSILON);
	for(int j = 0;j < nY;++j)
	{
		double sum = 0.0;
		vector<double> tmp(nX);
		for (int i = 0;i < nX;++i)
		{
			//arma::vec calcPos = forward(m_modelData.row(i).t());
			//double dist = MahalanobisDistance(m_observationData.row(j).t(), calcPos, m_covariances[i]);
			double dist = arma::norm(m_observationData.row(j).t() - m_transData.row(i).t(), 2);
			tmp[i] = exp(-dist*dist / (2.0 * m_delta_2));///sqrt(det(m_covariances[i]));
			sum += tmp[i];
		}
		sum += outlier;
		for(int i = 0;i < nX;++i)
		{
			m_alpha(j, i) = tmp[i] / (sum + DBL_EPSILON);
		}
		m_alpha(j, nX) = outlier / (sum + DBL_EPSILON);
	}

	MetropolisHastings();

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
	for (int i = 0;i < nX;++i)
	{
		m_lambda[i] = 0;
		for (int j = 0;j < nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			m_w.row(i) += m_alpha(j, i) * m_observationData.row(j);
		}
		//m_w.row(i) /= (m_lambda[i] + FLOAT_EPSILON);
	}
	return true;
}

void mcempr::funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	mcempr *pThis = (mcempr*)adata;
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
void mcempr::jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata)
{
	mcempr *pThis = (mcempr*)adata;

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

bool mcempr::parameters_optimization()
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

bool mcempr::updateCovariances()
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