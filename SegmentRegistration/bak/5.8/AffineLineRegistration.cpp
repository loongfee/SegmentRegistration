#include "AffineLineRegistration.h"
#include "LineMatch.h"


bool AffineLineRegistration::solve()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;
	initialization();

	double convergence_epsilon = 2.0e-4;
	double delta;
	int max_iteration = 200;
	int iTime = 0;

	Eigen::VectorXd oldParameters = getParameters();
	do 
	{
		// compute transformed straight lines
		updateTransLines();

		//plot
		// for debug
		if(iTime%1 == 0)
		{
			char szFilename[56];
			//sprintf_s(szFilename, "pic\\%03d.jpg", iTime+1);
			sprintf_s(szFilename, "pic\\line-%03d.png", iTime+1);
			drawResult(m_reference_lines_polar, m_trans_lines_polar, szFilename);
			//drawPolarResult(m_observe_lines_polar, m_trans_lines_polar, szFilename);

			if (fexists(m_sourceImageFile))
			{
				char szImageFilename[56];
				//sprintf_s(szImageFilename, "pic\\image_%03d.jpg", iTime+1);
				sprintf_s(szImageFilename, "pic\\image-%03d.png", iTime+1);
				vector<double> param(get_PARAMETER_NUM());
				for (int i = 0;i < get_PARAMETER_NUM();++i) param[i] = oldParameters[i];
				affine_warp_Image< vector<double> >(m_sourceImageFile, szImageFilename, param);
			}
		}

		double tmp = m_delta_2;
		updateDeltaSquare();
		cout<<"iteration:\t"<<iTime+1<<endl;
		cout<<"rms:\t"<<m_delta_2<<endl;
		cout<<oldParameters<<endl;
		if (iTime > 0 && m_delta_2 > tmp)
		{
			cout<<"Warning: the model is getting worse!"<<endl;
			break;
		}
		if (m_delta_2 < 1.0)
		{
			cout<<"rms epsilon reached!"<<endl;
			break;
		}
		if (fabs(m_delta_2 - tmp) < 1e-8)
		{
			cout<<"delta rms epsilon reached!"<<endl;
			break;
		}

		// E-step
		E_Step();

		// CM-step A
		// estimate the trans parameters
		parameters_optimization();
		Eigen::VectorXd tmpParameters = getParameters();
		delta = (oldParameters - tmpParameters).norm();
		oldParameters = tmpParameters;

	} while (++iTime < max_iteration && fabs(delta) > convergence_epsilon);
	if(iTime == max_iteration)
	{
		cout<<"Warning: the max number of iteration reached, and the results may be not accurate."<<endl;
	}

	//for (int jj = 0;jj < nY;++jj)
	//{
	//	int maxIndex = 0;
	//	double maxValue = 0.0;
	//	for (int ii = 0;ii < nX+1;++ii)
	//	{
	//		if(m_alpha(jj, ii) > maxValue)
	//		{
	//			maxValue = m_alpha(jj, ii);
	//			maxIndex = ii;
	//		}
	//	}
	//	if (46 == maxIndex)
	//	{
	//		cout<<"debug.."<<endl;
	//	}
	//	cvline_polar l_obs = m_observe_lines_polar[jj];
	//	cvline_polar l_trans = m_trans_lines_polar[46];
	//	fPoint center_obs((l_obs.pt1.x+l_obs.pt2.x)*0.5, (l_obs.pt1.y+l_obs.pt2.y)*0.5);
	//	fPoint center_trans((l_trans.pt1.x+l_trans.pt2.x)*0.5, (l_trans.pt1.y+l_trans.pt2.y)*0.5);
	//	double cen_dis = (center_obs.x-center_trans.x)*(center_obs.x-center_trans.x) + (center_obs.y-center_trans.y)*(center_obs.y-center_trans.y);
	//	if (5.0 > cen_dis)
	//	{
	//		cout<<"debug.."<<endl;
	//	}
	//}
	//for (int ii = 0;ii < nX+1;++ii)
	//{
	//	cout<<ii<<"  "<<m_alpha(85, ii)<<endl;
	//}

	//correspondence_struct current_correspondence = assignment_one();
	//m_z = current_correspondence.forward_map;
	classification();
	int nRemoved = 0;
	double delta_ = sqrt(m_delta_2);
	for(int j = 1; j < nY;j++)
	{
		if(m_z[j] == nX) continue;
		//// remove outliers
		//fPoint dist = segment2polarline(m_trans_lines_polar[m_z[j]], m_reference_lines_polar[j]);
		//double dist2 = sqrt(dist.x*dist.x + dist.y*dist.y);
		//cvline_polar l_obs = m_reference_lines_polar[j];
		//cvline_polar l_trans = m_trans_lines_polar[m_z[j]];
		//fPoint center_obs((l_obs.pt1.x+l_obs.pt2.x)*0.5, (l_obs.pt1.y+l_obs.pt2.y)*0.5);
		//fPoint center_trans((l_trans.pt1.x+l_trans.pt2.x)*0.5, (l_trans.pt1.y+l_trans.pt2.y)*0.5);
		//double cen_dis = (center_obs.x-center_trans.x)*(center_obs.x-center_trans.x) + (center_obs.y-center_trans.y)*(center_obs.y-center_trans.y);
		//if (cen_dis > 1000.0)
		//{
		//	nRemoved++;
		//	continue;
		//}
		//if (dist2 > 2*delta_)
		//{	
		//	nRemoved++;
		//	continue;
		//}
		m_inliers.push_back(j);
		m_inliers.push_back(m_z[j]);
	}

	//// RANSAC detect outliers
	//auto_ptr< estimators::Solver<mat,vec> > ptrSolver(
	//	new estimators::lineAffineSolver<mat,vec>);

	////-- Create input data
	//int nTotal = (int)m_inliers.size() / 2;
	//mat dataPoints(nTotal, 6);
	//for(int i = 0;i < nTotal;++i)
	//{
	//	dataPoints(i, 0) = m_observe_lines_polar[m_inliers[2*i]].rho;
	//	dataPoints(i, 1) = m_observe_lines_polar[m_inliers[2*i]].theta;
	//	dataPoints(i, 2) = m_model_lines_polar[m_inliers[2*i+1]].pt1.x;
	//	dataPoints(i, 3) = m_model_lines_polar[m_inliers[2*i+1]].pt1.y;
	//	dataPoints(i, 4) = m_model_lines_polar[m_inliers[2*i+1]].pt2.x;
	//	dataPoints(i, 5) = m_model_lines_polar[m_inliers[2*i+1]].pt2.y;
	//}
	//
	//vector<int> inliers;
	//vector<vec> models;

	//ransac::Ransac_Handler ransac_fun_Handler;
	//bool result = ransac::Ransac_RobustEstimator
	//	(
	//	dataPoints, // the input data
	//	estimators::affineSolver<mat,vec>::extractor, // How select sampled point from indices
	//	dataPoints.n_rows,  // the number of putatives data
	//	*(ptrSolver.get()),  // compute the underlying model given a sample set
	//	estimators::affineSolver<mat,vec>::defaultEvaluator,  // the function to evaluate a given model
	//	//Ransac Object that contain function:
	//	// CandidatesSelector, Sampler and TerminationFunction
	//	ransac_fun_Handler, // the basic ransac object
	//	1000,  // the maximum rounds for RANSAC routine
	//	inliers, // inliers to the final solution
	//	models, // models array that fit input data
	//	0.95, // the confidence want to achieve at the end
	//	15.0
	//	);

	//for (int i = 0;i < (int)models[0].n_rows;++i)
	//{
	//	cout<<models[0](i)<<endl;
	//	setAdjustableParameter(i, models[0](i));
	//}
	//cout<<endl;

	
	int nTotal = (int)m_inliers.size() / 2;
	*m_fs<<"iteration time:"<<iTime<<endl;
	*m_fs<<"final covariance:"<<m_delta_2<<endl;
	*m_fs<<"after assignment:"<<nTotal+nRemoved<<endl;
	*m_fs<<"after outlier elimination:"<<nTotal<<endl;
	vector<int> inliers(nTotal);
	//mat dataPoints(nTotal, 6);
	for(int i = 0;i < nTotal;++i)
	{
		inliers[i] = i;
	}

	//int n = (int)m_inliers.size() / 2;
	int n = (int)inliers.size();
	////output precision report
	FILE* pf;
	//pf = fopen("pic\\report.txt", "w+");
	//fprintf(pf, "%4s%15s%15s%15s%15s%15s%15s%15s\n", "ID", "cen_dis", "e1", "e2", "x1", "y1", "x2", "y2");
	////fReport.open("pic\\report.txt", std::ios_base::out);
	////fReport<<"ID\t e1\t e1\t x1\t y1\t x2\t y2\n";
	//for (int p = 0;p < n;++p)
	//{
	//	int j = inliers[p];
	//	fPoint dist = segment2polarline(m_trans_lines_polar[m_inliers[2*j+1]], m_reference_lines_polar[m_inliers[2*j]]);
	//	//fReport<<j+1
	//	//	<<"\t"<<dist.x
	//	//	<<"\t"<<dist.y
	//	//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt1.x
	//	//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt1.y
	//	//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt2.x
	//	//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt2.y
	//	//	<<endl;
	//	cvline_polar l_obs = m_reference_lines_polar[m_inliers[2*j]];
	//	cvline_polar l_trans = m_trans_lines_polar[m_inliers[2*j+1]];
	//	fPoint center_obs((l_obs.pt1.x+l_obs.pt2.x)*0.5, (l_obs.pt1.y+l_obs.pt2.y)*0.5);
	//	fPoint center_trans((l_trans.pt1.x+l_trans.pt2.x)*0.5, (l_trans.pt1.y+l_trans.pt2.y)*0.5);
	//	double cen_dis = (center_obs.x-center_trans.x)*(center_obs.x-center_trans.x) + (center_obs.y-center_trans.y)*(center_obs.y-center_trans.y);
	//	fprintf(pf, "%4d%15f%15f%15f%15f%15f%15f%15f\n", j+1, cen_dis, dist.x, dist.y, m_reference_lines_polar[m_inliers[2*j]].pt1.x,
	//		m_reference_lines_polar[m_inliers[2*j]].pt1.y, m_reference_lines_polar[m_inliers[2*j]].pt2.x, m_reference_lines_polar[m_inliers[2*j]].pt2.y);
	//}
	//fclose(pf);

	// control line segments
	pf = fopen("pic\\lines.txt", "w+");
	for (int p = 0;p < n;++p)
	{
		int j = inliers[p];
		fprintf(pf, "%4d%15f%15f%15f%15f%15f\n", j+1, m_source_lines_polar[m_inliers[2*j+1]].pt1.x, m_source_lines_polar[m_inliers[2*j+1]].pt1.y, 
			m_reference_lines_polar[m_inliers[2*j]].pt1.x, m_reference_lines_polar[m_inliers[2*j]].pt1.y);
		fprintf(pf, "%4d%15f%15f%15f%15f%15f\n", j+1, m_source_lines_polar[m_inliers[2*j+1]].pt2.x, m_source_lines_polar[m_inliers[2*j+1]].pt2.y,
			m_reference_lines_polar[m_inliers[2*j]].pt2.x, m_reference_lines_polar[m_inliers[2*j]].pt2.y);
	}
	fclose(pf);
	//fReport.close();


	// draw mathes
	//vector<cvline_polar> model_ls;
	//vector<cvline_polar> observe_ls;
	m_matched_src_ls.clear();
	m_matched_ref_ls.clear();
	//int n = (int)m_inliers.size() / 2;
	for (int p = 0;p < n;++p)
	{
		int j = inliers[p];
		m_matched_src_ls.push_back(m_source_lines_polar[m_inliers[2*j+1]]);
		m_matched_ref_ls.push_back(m_reference_lines_polar[m_inliers[2*j]]);
	}

	//cv::Mat img_observe = cv::imread(m_sourceImageFile);
	//cv::Mat img_model = cv::imread(m_referenceImageFile);
	//cv::Mat img_matches;
	//drawLineMatches(img_model, m_matched_ref_ls, img_observe, m_matched_src_ls,
	//	img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
	//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, 2);

	//cv::imwrite("pic\\matches.png", img_matches);

	vector<double> param(get_PARAMETER_NUM());
	for (int i = 0;i < get_PARAMETER_NUM();++i) param[i] = oldParameters[i];
	//drawLastResult(m_observeImageFile, m_modelImageFile, param, m_observe_lines_polar, m_trans_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	drawLastResult(m_sourceImageFile, m_referenceImageFile, m_reference_lines_polar, m_source_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	//drawLastResult(m_observe_lines_polar, m_trans_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	cout<<"EM based point set registration finished after "<<iTime<<" times of iterations."<<endl;
	cout<<"********************************************************************************"<<endl<<endl;

	return true;
}


bool AffineLineRegistration::classification()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

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

correspondence_struct AffineLineRegistration::assignment_multi()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	// valid correspondence vector
	// indices start from 0
	// nX indicates invisible or outlier
	correspondence_struct correspond;
	correspond.forward_map.resize(nY, nX);
	correspond.inverse_map.resize(nX, -1);

	vector<int> openList(nY);
	for (int i = 0;i < nY;++i) openList[i] = i;

	// randomly choose the start index
	// initialize random seed:
	srand ( time(NULL) );
	// generate start index:
	//int j = rand() % openList.size();
	int j = 0;
	vector<int> closeList;
	while((int)openList.size() > 0)
	{
		int real_pose = j % (int)openList.size();
		//if (27 == openList[real_pose])
		//{
		//	cout<<"debug.."<<endl;
		//}

		int maxIndex = 0;
		double maxValue = 0.0;
		for (int i = 0;i < nX+1;++i)
		{
			if(m_alpha(openList[real_pose], i) > maxValue)
			{
				maxValue = m_alpha(openList[real_pose], i);
				maxIndex = i;
			}
		}

		if(maxIndex == nX)
		{
			closeList.push_back(openList[real_pose]);
			openList.erase(openList.begin()+real_pose);
			continue;
		}

		correspond.forward_map[openList[real_pose]] = maxIndex;
		correspond.inverse_map[maxIndex] = openList[real_pose];
		closeList.push_back(openList[real_pose]);
		openList.erase(openList.begin()+real_pose);
	}

	return correspond;
}

correspondence_struct AffineLineRegistration::assignment_one()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	// valid correspondence vector
	// indices start from 0
	// nX indicates invisible or outlier
	correspondence_struct correspond;
	correspond.forward_map.resize(nY, nX);
	correspond.inverse_map.resize(nX, -1);

	vector<int> openList(nY);
	for (int i = 0;i < nY;++i) openList[i] = i;

	// randomly choose the start index
	// initialize random seed:
	srand ( time(NULL) );
	// generate start index:
	//int j = rand() % openList.size();
	int j = 0;
	vector<int> closeList;
	while((int)openList.size() > 0)
	{
		int real_pose = j % (int)openList.size();
		//if (27 == openList[real_pose])
		//{
		//	cout<<"debug.."<<endl;
		//}


		int maxIndex = nX;
		double maxValue = 0.0;
		for (int i = 0;i < nX+1;++i)
		{
			if(m_alpha(openList[real_pose], i) > maxValue)
			{
				maxValue = m_alpha(openList[real_pose], i);
				maxIndex = i;
			}
		}
		//if (46 == maxIndex)
		//{
		//	cout<<"debug.."<<endl;
		//}

		if(maxIndex == nX)
		{
			closeList.push_back(openList[real_pose]);
			openList.erase(openList.begin()+real_pose);
			continue;
		}
		// 
		if(correspond.inverse_map[maxIndex] != -1)
		{
			// assigned
			int old_observation = correspond.inverse_map[maxIndex];
			cvline_polar old_line = m_reference_lines_polar[old_observation];
			cvline_polar new_line = m_reference_lines_polar[openList[real_pose]];
			cvline_polar trans_line = m_trans_lines_polar[maxIndex];
			fPoint center_old((old_line.pt1.x+old_line.pt2.x)*0.5, (old_line.pt1.y+old_line.pt2.y)*0.5);
			fPoint center_new((new_line.pt1.x+new_line.pt2.x)*0.5, (new_line.pt1.y+new_line.pt2.y)*0.5);
			fPoint center_trans((trans_line.pt1.x+trans_line.pt2.x)*0.5, (trans_line.pt1.y+trans_line.pt2.y)*0.5);
			double dis_old = (center_old.x-center_trans.x)*(center_old.x-center_trans.x) + (center_old.y-center_trans.y)*(center_old.y-center_trans.y);
			double dis_new = (center_new.x-center_trans.x)*(center_new.x-center_trans.x) + (center_new.y-center_trans.y)*(center_new.y-center_trans.y);
			if(dis_new < dis_old)
			{
				correspond.forward_map[openList[real_pose]] = maxIndex;
				correspond.inverse_map[maxIndex] = openList[real_pose];
				closeList.push_back(openList[real_pose]);
				openList.erase(openList.begin()+real_pose);

				correspond.forward_map[old_observation] = nX;
				openList.push_back(old_observation);
				//j++;
			}
			else
			{
				int maxIndex = nX;
				double maxValue = 0.0;
				for (int i = 0;i < nX+1;++i)
				{
					if(correspond.inverse_map[i] == -1 && m_alpha(openList[real_pose], i) > maxValue)
					{
						maxValue = m_alpha(openList[real_pose], i);
						maxIndex = i;
					}
				}
				if(maxIndex == nX)
				{
					closeList.push_back(openList[real_pose]);
					openList.erase(openList.begin()+real_pose);
					continue;
				}

				correspond.forward_map[openList[real_pose]] = maxIndex;
				correspond.inverse_map[maxIndex] = openList[real_pose];
				closeList.push_back(openList[real_pose]);
				openList.erase(openList.begin()+real_pose);
			}
		}
		else
		{
			// not assigned
			correspond.forward_map[openList[real_pose]] = maxIndex;
			correspond.inverse_map[maxIndex] = openList[real_pose];
			closeList.push_back(openList[real_pose]);
			openList.erase(openList.begin()+real_pose);
		}
	}

	return correspond;
}

void AffineLineRegistration::updateDeltaSquare()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	double sum = 0.0;
	for(int i = 0;i < nX;++i)
	{
		for (int j = 0;j < nY;++j)
		{
			double dist = (m_refData.row(j).transpose() - m_transData.row(i).transpose()).norm();
			sum += dist*dist*m_alpha(j,i);
		}
	}

	m_delta_2 = sum / (m_lammda_sum * ndim);
	if(m_delta_2 < 1.0) m_delta_2 = 1.0;
}


void AffineLineRegistration::initAdjustableParameters()
{
	m_parameters.resize(NUM_ADJUSTABLE_PARAMS);
	int numParams = (int)m_parameters.size();

	double a0,a1,a2,b0,b1,b2;
	a0 = m_parameters[a0_param] = 0.0;	//a0
	a1 = m_parameters[a1_param] = 1.0;	//a1
	a2 = m_parameters[a2_param] = 0.0;	//a2
	b0 = m_parameters[b0_param] = 0.0;	//b0
	b1 = m_parameters[b1_param] = 0.0;	//b1
	b2 = m_parameters[b2_param] = 1.0;	//b2
}

bool AffineLineRegistration::initialization()
{
	initAdjustableParameters();
	m_trans_lines_polar.resize(m_nX);
	m_refData.resize(m_nY	, 2);
	for (int i = 0;i < m_nY;++i)
	{
		double rho = m_reference_lines_polar[i].rho;
		double theta = m_reference_lines_polar[i].theta;
		double x = -sin(theta) * (rho + DBL_EPSILON);
		double y = cos(theta) * (rho + DBL_EPSILON);
		m_refData(i,0) = x;
		m_refData(i,1) = y;
	}

	// initialize delta
	//// debug
	//m_delta_2 = 0.0;
	//for(int i = 0;i < m_nX;++i)
	//{
	//	for (int j = 0;j < m_nY;++j)
	//	{
	//		fPoint dist = distance2(m_model_lines_polar[i], m_observe_lines_polar[j]);
	//		m_delta_2 += dist.x*dist.x + dist.y*dist.y;
	//	}
	//}
	//m_delta_2 /= (m_nX * m_nY * m_ndim);

	////m_delta_2 = 0.0;
	////for (int i = 0;i < m_nX;++i)
	////{
	////	cvline_polar transed_segment = forward(m_model_lines_polar[i]);
	////	for (int j = 0;j < m_nY;++j)
	////	{
	////		double delta_rho = m_observe_lines_polar[j].rho - transed_segment.rho;
	////		double delta_theta = m_observe_lines_polar[j].theta - transed_segment.theta;
	////		double dist_2 = delta_rho * delta_rho + scale_angle* delta_theta * delta_theta;
	////		//double dist = segment_distance(m_observe_lines_polar[j], transed_segment);
	////		m_delta_2 += dist_2;
	////	}
	////}
	//m_delta_2 /= m_ndim * m_nX * m_nY;

	//set the size of posterior matrix
	m_alpha.resize(m_nY, m_nX + 1);
	double w = 1.0 / (m_nX + 1);
	m_alpha.fill(w);

	// set the size of virtual observation matrix
	m_w.resize(m_nX, 2);

	// set the weight vector of virtual observations
	m_lambda.resize(m_nX);

	m_lammda_sum = 0.0;
	for (int i = 0;i < m_nX;++i)
	{
		m_lambda[i] = 0;
		double theta = 0.0;
		double rho = 0.0;
		for (int j = 0;j < m_nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			m_w.row(i) += m_alpha(j, i) * m_refData.row(j);
		}
		m_lammda_sum += m_lambda[i];
	}

	//
	m_z.resize(m_nY);
	//m_z.set_size(nY);

	//
	m_transData.resize(m_nX, 2);

	return true;
}

Eigen::Vector2d AffineLineRegistration::forward(const cvline_polar& l)const
{
	Eigen::VectorXd parameters = getParameters();
	double rho = l.rho;
	double sin_theta = sin(l.theta);
	double cos_theta = cos(l.theta);
	double a0 = m_parameters[a0_param];
	double a1 = m_parameters[a1_param];
	double a2 = m_parameters[a2_param];
	double b0 = m_parameters[b0_param];
	double b1 = m_parameters[b1_param];
	double b2 = m_parameters[b2_param];
	Eigen::Vector2d newPt;
	double m = (a1*a1+b1*b1)*cos_theta*cos_theta + 2*(a1*a2+b1*b2)*sin_theta*cos_theta + (a2*a2+b2*b2)*sin_theta*sin_theta;
	newPt(0) = (a0*a1+b0*b1)*cos_theta*cos_theta + (a0*a2+b0*b2)*sin_theta*cos_theta + (a1*a2+b1*b2)*rho*cos_theta + (a2*a2+b2*b2)*rho*sin_theta;
	newPt(1) = (a0*a2+b0*b2)*sin_theta*sin_theta + (a0*a1+b0*b1)*sin_theta*cos_theta - (a1*a2+b1*b2)*rho*sin_theta - (a1*a1+b1*b1)*rho*cos_theta;
	newPt(0) /= -(m+DBL_EPSILON);
	newPt(1) /= -(m+DBL_EPSILON);
	return newPt;
};

Eigen::VectorXd AffineLineRegistration::getParameters() const
{
	int num = get_PARAMETER_NUM();
	Eigen::VectorXd parameters(num);
	for (int i = 0;i < num;++i)
	{
		parameters[i] = m_parameters[i];
	}
	return parameters;
}

bool AffineLineRegistration::E_Step()
{
	// E-step evaluate the posteriors

	double outlier = 2.0 / 9.0;
	for(int j = 0;j < m_nY;++j)
	{
		double sum = 0.0;
		vector<double> tmp(m_nX);
		double n_ = 1.0 / (double)m_nX;
		for (int i = 0;i < m_nX;++i)
		{
			double dist = (m_refData.row(j).transpose() - m_transData.row(i).transpose()).norm();
			tmp[i] = exp(-dist*dist / (2.0 * m_delta_2 + DBL_EPSILON));
			sum += tmp[i];
		}
		
		//cout<<sum<<"\t"<<outlier<<endl;
		sum += outlier;
		for(int i = 0;i < m_nX;++i)
		{
			m_alpha(j, i) = tmp[i] / (sum + DBL_EPSILON);
		}
		m_alpha(j, m_nX) = outlier / (sum + DBL_EPSILON);
	}

	m_w.fill(0.0);
	m_lammda_sum = 0.0;
	for (int i = 0;i < m_nX;++i)
	{
		m_lambda[i] = 0;
		for (int j = 0;j < m_nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			m_w.row(i) += m_alpha(j,i) * m_refData.row(j);
		}
		m_lammda_sum += m_lambda[i];
	}
	return true;
}


void AffineLineRegistration::updateTransLines()
{
	for (int i = 0;i < m_nX;++i)
	{
		m_transData.row(i) = forward(m_source_lines_polar[i]);

		fPoint line[2];
		line[0].x = m_parameters[0] + m_parameters[1] * m_source_lines_polar[i].pt1.x + m_parameters[2] * m_source_lines_polar[i].pt1.y;
		line[0].y = m_parameters[3] + m_parameters[4] * m_source_lines_polar[i].pt1.x + m_parameters[5] * m_source_lines_polar[i].pt1.y;
		line[1].x = m_parameters[0] + m_parameters[1] * m_source_lines_polar[i].pt2.x + m_parameters[2] * m_source_lines_polar[i].pt2.y;
		line[1].y = m_parameters[3] + m_parameters[4] * m_source_lines_polar[i].pt2.x + m_parameters[5] * m_source_lines_polar[i].pt2.y;
		m_trans_lines_polar[i] = line2polar(line);
	}
}

bool AffineLineRegistration::parameters_optimization()
{
	//return linear_optimization();
	// levmar
	//return levmar_optimization();

	// spase LM
	//return spaseLM_optimization();

	// alglib improved LM
	//return alglib_optimization();

	// Eigen
	return eigen_levmar_optimization();

}

bool AffineLineRegistration::linear_optimization()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	VectorXd b(6);
	MatrixXd A(6, 6);
	b.fill(0.0);
	A.fill(0.0);

	int c = 0;
	for(int i = 0;i < (int)nX;++i)
	{
		double sin_theta = sin(m_source_lines_polar[i].theta);
		double cos_theta = cos(m_source_lines_polar[i].theta);
		double rho = m_source_lines_polar[i].rho;
		
		for (int j = 0;j < nY;++j)
		{
			double ref_sin_theta = sin(m_reference_lines_polar[j].theta);
			double ref_cos_theta = cos(m_reference_lines_polar[j].theta);
			double ref_rho = m_reference_lines_polar[j].rho;

			double len = segment_length(m_reference_lines_polar[j]);
			double weight = m_alpha(j,i);

			vector<double> coeff(6);
			coeff[0] = 0.0;
			coeff[1] = -ref_sin_theta * cos_theta;
			coeff[2] = -ref_sin_theta * sin_theta;
			coeff[3] = 0.0;
			coeff[4] = ref_cos_theta * cos_theta;
			coeff[5] = ref_cos_theta * sin_theta;

			double y = 0.0;
			for(int p1=0;p1<6;++p1)
			{        
				b[p1] += coeff[p1] * y * weight;
				for(int p2=0;p2<6;++p2)
				{
					A(p1,p2) += coeff[p1] * coeff[p2] * weight;
				}
			}

			coeff[0] = -ref_sin_theta;
			coeff[1] = ref_sin_theta *sin_theta*rho ;
			coeff[2] = -ref_sin_theta *cos_theta*rho;
			coeff[3] = ref_cos_theta;
			coeff[4] = -ref_cos_theta *sin_theta*rho ;
			coeff[5] = ref_cos_theta *cos_theta*rho;
			 y = ref_rho;
			for(int p1=0;p1<6;++p1)
			{        
				b[p1] += coeff[p1] * rho * weight;
				for(int p2=0;p2<6;++p2)
				{
					A(p1,p2) += coeff[p1] * coeff[p2] * weight;
				}
			}
		}
	}

	Eigen::MatrixXd damper = MatrixXd::Identity(6, 6);
	//A = A + damper * 1.0E-6;
	double tmp = A.determinant();
	if ( fabs(tmp) > DBL_EPSILON )
	{
		Eigen::VectorXd parameter = A.inverse() * b;
		for(int i=0;i<get_PARAMETER_NUM();++i) m_parameters[i] = parameter[i];
		return true;
	}
	else
	{
		cout << "ERROR : cannot find a solution" << endl;
		return false;
	}
	return true;
}


bool AffineLineRegistration::eigen_levmar_optimization()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nparam = get_PARAMETER_NUM();
	
	int info;
	//double fnorm, covfac;
	Eigen::VectorXd x(nparam);
	for(int i=0;i<nparam;++i) x[i] = m_parameters[i];
	
	// do the computation
	lmder_functor_affine functor(nparam, 2*nX);
	functor.pData = this;
	Eigen::LevenbergMarquardt<lmder_functor_affine> lm(functor);
	info = lm.minimize(x);
	
	for(int i=0;i<nparam;++i) m_parameters[i] = x[i];

	return true;
}



int lmder_functor_affine::operator()(const VectorXd &x, VectorXd &fvec) const
{
	AffineLineRegistration *pThis = (AffineLineRegistration*)pData;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int nparameter = pThis->get_PARAMETER_NUM();

	for(int i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = x[i];
	}

	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;

	int c = 0;
	for (int i = 0;i < nX;++i)
	{
		Eigen::Vector2d outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		fvec[c++] = -pThis->m_lambda[i]*(pThis->m_w(i, 0) - outPt(0));
		fvec[c++] = -pThis->m_lambda[i]*(pThis->m_w(i, 1) - outPt(1));
	}
	return 0;
}

//int lmder_functor_affine::df(const VectorXd &x, MatrixXd &fjac) const
//{
//	AffineLineRegistration *pThis = (AffineLineRegistration*)pData;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	for(int i=0;i<nparameter;++i)
//	{
//		pThis->m_parameters[i] = x[i];
//	}
//
//	int c = 0;
//
//	for (int i = 0;i < nX;++i)
//	{
//		Eigen::Vector2d pt = pThis->m_w.row(i);
//		Eigen::Vector2d outPt = pThis->forward(pThis->m_source_lines_polar[i]);
//
//		double a0 = pThis->m_parameters[pThis->a0_param];
//		double a1 = pThis->m_parameters[pThis->a1_param];
//		double a2 = pThis->m_parameters[pThis->a2_param];
//		double b0 = pThis->m_parameters[pThis->b0_param];
//		double b1 = pThis->m_parameters[pThis->b1_param];
//		double b2 = pThis->m_parameters[pThis->b2_param];
//
//		double rho = pThis->m_source_lines_polar[i].rho;
//		double sin_theta = sin(pThis->m_source_lines_polar[i].theta);
//		double cos_theta = cos(pThis->m_source_lines_polar[i].theta);
//
//		double m = (a1*a1+b1*b1)*cos_theta*cos_theta + 2*(a1*a2+b1*b2)*sin_theta*cos_theta + (a2*a2+b2*b2)*sin_theta*sin_theta;
//		int p = 0;
//		// a0
//		fjac(c,p) = pThis->m_lambda[i]*(-a1*cos_theta*cos_theta-a2*sin_theta*cos_theta) / (m+DBL_EPSILON);
//		fjac(c+1,p++) = pThis->m_lambda[i]*(-a2*sin_theta*sin_theta-a1*sin_theta*cos_theta) / (m+DBL_EPSILON);
//
//		// a1
//		fjac(c,p) = pThis->m_lambda[i]*pt(1) / m;
//		fjac(c+1,p++) = pThis->m_lambda[i]*0.0;
//
//		// a2
//		fjac(c,p) = pThis->m_lambda[i]*0.0;
//		fjac(c+1,p++) = pThis->m_lambda[i]*pt(0) / m;
//
//		// b0
//		fjac(c,p) = pThis->m_lambda[i]*0.0;
//		fjac(c+1,p++) = pThis->m_lambda[i]*pt(1) / m;
//
//		// b1
//		fjac(c,p) = -pThis->m_lambda[i]*(c1*pt(0)+c2*pt(1))*pt(0) / (m*m);
//		fjac(c+1,p++) = -pThis->m_lambda[i]*(c3*pt(0)+c4*pt(1))*pt(0) / (m*m);
//
//		// b2
//		fjac(c,p) = -pThis->m_lambda[i]*(c1*pt(0)+c2*pt(1))*pt(1) / (m*m);
//		fjac(c+1,p++) = -pThis->m_lambda[i]*(c3*pt(0)+c4*pt(1))*pt(1) / (m*m);
//
//		c += 2;
//	}
//	return 0;
//}

int lmder_functor_affine::df(const VectorXd &x, MatrixXd &fjac) const
{
	AffineLineRegistration *pThis = (AffineLineRegistration*)pData;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int nparameter = pThis->get_PARAMETER_NUM();

	for(int i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = x[i];
	}


	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;
	int c = 0;
	for (int i = 0;i < nX;++i)
	{
		Eigen::Vector2d pt = pThis->m_w.row(i);
		for(int p=0;p<nparameter;++p)
		{
			double middle = pThis->m_parameters[p];
			pThis->m_parameters[p] = middle + pstep_scale;
			Eigen::Vector2d outPt1 = pThis->forward(pThis->m_source_lines_polar[i]);

			pThis->m_parameters[p] = middle - pstep_scale;
			Eigen::Vector2d outPt2 = pThis->forward(pThis->m_source_lines_polar[i]);

			pThis->m_parameters[p] = middle;

			double derivative_x = (outPt1(0) - outPt2(0)) * den;
			double derivative_y = (outPt1(1) - outPt2(1)) * den;

			fjac(c,p) = pThis->m_lambda[i] * derivative_x;
			fjac(c+1,p) = pThis->m_lambda[i] * derivative_y;
		}
		c += 2;
	}
	return 0;
}