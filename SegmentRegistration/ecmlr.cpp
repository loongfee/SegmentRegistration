#include "ecmlr.h"
#include "LineMatch.h"
#include "SpectralMatching.h"

ecmlr::ecmlr( const vector<cvline_polar>& source_lines_polar, const vector<cvline_polar>& reference_lines_polar)
{
	// merge the colinear straight lines
	//double merge_threshold = 1.0;
	//for (int i = 0;i < (int)model_lines_polar.size();++i)
	//{
	//	bool is_exist = false;
	//	for (int j = 0;j < (int)m_model_lines_polar.size();++j)
	//	{
	//		fPoint dis = segment2polarline(model_lines_polar[i], m_model_lines_polar[j]);
	//		if((dis.x*dis.x + dis.y*dis.y) < merge_threshold)
	//		{
	//			// merge
	//			m_model_lines_polar[j] = mergeline(m_model_lines_polar[j], model_lines_polar[i]);
	//			is_exist = true;
	//			break;
	//		}
	//	}
	//	if(!is_exist)
	//	{
	//		m_model_lines_polar.push_back(model_lines_polar[i]);
	//	}
	//}

	//for (int i = 0;i < (int)observe_lines_polar.size();++i)
	//{
	//	bool is_exist = false;
	//	for (int j = 0;j < (int)m_observe_lines_polar.size();++j)
	//	{
	//		fPoint dis = segment2polarline(observe_lines_polar[i], m_observe_lines_polar[j]);
	//		if((dis.x*dis.x + dis.y*dis.y) < merge_threshold)
	//		{
	//			// merge
	//			m_observe_lines_polar[j] = mergeline(m_observe_lines_polar[j], observe_lines_polar[i]);
	//			is_exist = true;
	//			break;
	//		}
	//	}
	//	if(!is_exist)
	//	{
	//		m_observe_lines_polar.push_back(observe_lines_polar[i]);
	//	}
	//}

	m_source_lines_polar = source_lines_polar;
	//for (int i = 0;i < (int)m_model_lines_polar.size();++i)
	//{
	//	normalize_segment(m_model_lines_polar[i]);
	//}
	m_reference_lines_polar = reference_lines_polar;
	//for (int i = 0;i < (int)m_observe_lines_polar.size();++i)
	//{
	//	normalize_segment(m_observe_lines_polar[i]);
	//}

	m_nX = (int)m_source_lines_polar.size();
	m_nY = (int)m_reference_lines_polar.size();
	m_ndim = 2;
	//initAdjustableParameters();
}

bool ecmlr::solve()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;
	initialization();

	double convergence_epsilon = 2.0e-4;
	double delta;
	int max_iteration = 200;
	int iTime = 0;

	string foldName = m_sourceImageFile.substr(0, m_sourceImageFile.find_last_of('\\'));

	Eigen::VectorXd oldParameters = getParameters();
	do 
	{
		// compute transformed straight lines
		updateTransLines();

		//plot
		// for debug
		//if(iTime%1 == 0)
		//{
		//	char szFilename[56];
		//	//sprintf_s(szFilename, "pic\\%03d.jpg", iTime+1);
		//	sprintf_s(szFilename, "pic\\line-%03d.png", iTime+1);
		//	drawResult(m_observe_lines_polar, m_trans_lines_polar, szFilename);
		//	//drawPolarResult(m_observe_lines_polar, m_trans_lines_polar, szFilename);

		//	if (fexists(m_modelImageFile))
		//	{
		//		char szImageFilename[56];
		//		//sprintf_s(szImageFilename, "pic\\image_%03d.jpg", iTime+1);
		//		sprintf_s(szImageFilename, "pic\\image-%03d.png", iTime+1);
		//		vector<double> param(get_PARAMETER_NUM());
		//		for (int i = 0;i < get_PARAMETER_NUM();++i) param[i] = oldParameters[i];
		//		affine_warp_Image< vector<double> >(m_modelImageFile, szImageFilename, param);
		//	}
		//}

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
		if (fabs(m_delta_2 - tmp) < 0.01)
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

		// CM-step B
		// estimate the covariance matrices
		//updateCovariances();
		//m_delta_2 *= 0.9;
	} while (++iTime < max_iteration && fabs(delta) > convergence_epsilon);
	if(iTime == max_iteration)
	{
		cout<<"Warning: the max number of iteration reached, and the results may be not accurate."<<endl;
	}
	//MetropolisHastings();
	//// classification
	//classification();

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
	//double delta_ = sqrt(m_delta_2);
	//for(int j = 1; j < nY;j++)
	//{
	//	if(m_z[j] == nX) continue;
	//	// remove outliers
	//	fPoint dist = segment2polarline(m_trans_lines_polar[m_z[j]], m_reference_lines_polar[j]);
	//	double dist2 = sqrt(dist.x*dist.x + dist.y*dist.y);
	//	cvline_polar l_obs = m_reference_lines_polar[j];
	//	cvline_polar l_trans = m_trans_lines_polar[m_z[j]];
	//	fPoint center_obs((l_obs.pt1.x+l_obs.pt2.x)*0.5, (l_obs.pt1.y+l_obs.pt2.y)*0.5);
	//	fPoint center_trans((l_trans.pt1.x+l_trans.pt2.x)*0.5, (l_trans.pt1.y+l_trans.pt2.y)*0.5);
	//	double cen_dis = (center_obs.x-center_trans.x)*(center_obs.x-center_trans.x) + (center_obs.y-center_trans.y)*(center_obs.y-center_trans.y);
	//	if (cen_dis > 1000.0)
	//	{
	//		nRemoved++;
	//		continue;
	//	}
	//	if (dist2 > 2*delta_)
	//	{	
	//		nRemoved++;
	//		continue;
	//	}
	//	m_inliers.push_back(j);
	//	m_inliers.push_back(m_z[j]);
	//}

	
	//SpectralMatching Matcher;
	//m_inliers.clear();
	//Matcher.Matching(m_reference_lines_polar,m_trans_lines_polar,m_inliers);


	//// RANSAC detect outliers
	//auto_ptr< estimators::Solver<mat,vec> > ptrSolver(
	//	new estimators::lineAffineSolver<mat,vec>);

	////-- Create input data
	//int nTotal = (int)m_inliers.size() / 2;
	//mat dataPoints(nTotal, 6);
	//for(int i = 0;i < nTotal;++i)
	//{
	//	dataPoints(i, 0) = m_reference_lines_polar[m_inliers[2*i]].rho;
	//	dataPoints(i, 1) = m_reference_lines_polar[m_inliers[2*i]].theta;
	//	dataPoints(i, 2) = m_source_lines_polar[m_inliers[2*i+1]].pt1.x;
	//	dataPoints(i, 3) = m_source_lines_polar[m_inliers[2*i+1]].pt1.y;
	//	dataPoints(i, 4) = m_source_lines_polar[m_inliers[2*i+1]].pt2.x;
	//	dataPoints(i, 5) = m_source_lines_polar[m_inliers[2*i+1]].pt2.y;
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
	//	10.0
	//	);
	//
	//*m_fs<<"after outlier elimination:"<<(int)inliers.size()<<endl;
	////for (int i = 0;i < (int)models[0].n_rows;++i)
	////{
	////	cout<<models[0](i)<<endl;
	////	setAdjustableParameter(i, models[0](i));
	////}
	////cout<<endl;

	
	int nTotal = (int)m_inliers.size() / 2;
	*m_fs<<"iteration time:"<<iTime<<endl;
	*m_fs<<"final covariance:"<<m_delta_2<<endl;
	*m_fs<<"after assignment:"<<nTotal+nRemoved<<endl;
	*m_fs<<"after outlier elimination:"<<nTotal<<endl;
	*m_fs<<"parameters:"<<endl<<getParameters()<<endl;
	vector<int> inliers(nTotal);
	//mat dataPoints(nTotal, 6);
	for(int i = 0;i < nTotal;++i)
	{
		inliers[i] = i;
	}

	//output precision report
	FILE* pf;
	pf = fopen((foldName+"\\report.txt").c_str(), "w+");
	fprintf(pf, "%4s%15s%15s%15s%15s%15s%15s%15s\n", "ID", "cen_dis", "e1", "e2", "x1", "y1", "x2", "y2");
	//fReport.open("pic\\report.txt", std::ios_base::out);
	//fReport<<"ID\t e1\t e1\t x1\t y1\t x2\t y2\n";
	//int n = (int)m_inliers.size() / 2;
	int n = (int)inliers.size();
	for (int p = 0;p < n;++p)
	{
		int j = inliers[p];
		fPoint dist = segment2polarline(m_trans_lines_polar[m_inliers[2*j+1]], m_reference_lines_polar[m_inliers[2*j]]);
		//fReport<<j+1
		//	<<"\t"<<dist.x
		//	<<"\t"<<dist.y
		//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt1.x
		//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt1.y
		//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt2.x
		//	<<"\t"<<m_observe_lines_polar[m_inliers[j]].pt2.y
		//	<<endl;
		cvline_polar l_obs = m_reference_lines_polar[m_inliers[2*j]];
		cvline_polar l_trans = m_trans_lines_polar[m_inliers[2*j+1]];
		fPoint center_obs((l_obs.pt1.x+l_obs.pt2.x)*0.5, (l_obs.pt1.y+l_obs.pt2.y)*0.5);
		fPoint center_trans((l_trans.pt1.x+l_trans.pt2.x)*0.5, (l_trans.pt1.y+l_trans.pt2.y)*0.5);
		double cen_dis = (center_obs.x-center_trans.x)*(center_obs.x-center_trans.x) + (center_obs.y-center_trans.y)*(center_obs.y-center_trans.y);
		fprintf(pf, "%4d%15f%15f%15f%15f%15f%15f%15f\n", j+1, cen_dis, dist.x, dist.y, m_reference_lines_polar[m_inliers[2*j]].pt1.x,
			m_reference_lines_polar[m_inliers[2*j]].pt1.y, m_reference_lines_polar[m_inliers[2*j]].pt2.x, m_reference_lines_polar[m_inliers[2*j]].pt2.y);
	}
	fclose(pf);

	// control line segments
	pf = fopen((foldName+"\\lines.txt").c_str(), "w+");
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

	//cv::Mat img_observe = cv::imread(m_referenceImageFile);
	//cv::Mat img_model = cv::imread(m_sourceImageFile);
	//cv::Mat img_matches;
	//drawLineMatches(img_model, m_matched_ref_ls, img_observe, m_matched_src_ls,
	//	img_matches, cv::Scalar::all(-1), cv::Scalar::all(-1),
	//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, 2);

	//cv::imwrite((foldName+"\\matches.png").c_str(), img_matches);

	vector<double> param(get_PARAMETER_NUM());
	for (int i = 0;i < get_PARAMETER_NUM();++i) param[i] = oldParameters[i];
	//drawLastResult(m_observeImageFile, m_modelImageFile, param, m_observe_lines_polar, m_trans_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	//drawLastResult(m_sourceImageFile, m_referenceImageFile, m_reference_lines_polar, m_source_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	//drawLastResult(m_observe_lines_polar, m_trans_lines_polar, m_inliers, "pic\\result1.jpg", "pic\\result2.jpg");
	cout<<"EM based point set registration finished after "<<iTime<<" times of iterations."<<endl;
	cout<<"********************************************************************************"<<endl<<endl;
	return true;
}


bool ecmlr::classification()
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

double ecmlr::acceptance_ratio(const correspondence_struct& c, const correspondence_struct& c_prime)
{
	double result = 1.0;
	for (int j = 0;j < (int)c.forward_map.size();++j)
	{
		double test = m_alpha(34, 91);
		double a_prime = m_alpha(j, c_prime.forward_map[j]);
		double a = m_alpha(j, c.forward_map[j]);
		result *= (m_alpha(j, c_prime.forward_map[j]))/(m_alpha(j, c.forward_map[j]) + DBL_EPSILON);
	}
	return result;
}

correspondence_struct ecmlr::getValidCorrespondence()
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
			if(m_alpha(real_pose, candidate_list[i]) > maxValue)
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
	//for (int i = 0;i < nX+1;++i)
	//{
	//	cout<<m_alpha(2, i)<<endl;
	//}
	//cout<<endl;

	//for (int i = 0;i < (int)correspond.forward_map.size();++i)
	//{
	//	cout<<m_alpha(i, correspond.forward_map[i])<<endl;
	//}
	//cout<<endl;
	return correspond;
}

correspondence_struct ecmlr::assignment_multi()
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

correspondence_struct ecmlr::assignment_one()
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


//correspondence_struct ecmlr::assignment_one2one()
//{
//	int nX = m_nX;
//	int nY = m_nY;
//	int ndim = m_ndim;
//
//	// valid correspondence vector
//	// indices start from 0
//	// nX indicates invisible or outlier
//	correspondence_struct correspond;
//	correspond.forward_map.resize(nY, nX);
//	correspond.inverse_map.resize(nX, -1);
//
//	vector<int> openList(nY);
//	for (int i = 0;i < nY;++i) openList[i] = i;
//
//	// randomly choose the start index
//	// initialize random seed:
//	srand ( time(NULL) );
//	// generate start index:
//	//int j = rand() % openList.size();
//	int j = 0;
//	vector<int> closeList;
//	while((int)openList.size() > 0)
//	{
//		int real_pose = j % (int)openList.size();
//		//if (27 == openList[real_pose])
//		//{
//		//	cout<<"debug.."<<endl;
//		//}
//
//		int maxIndex = nX;
//		double maxValue = 0.0;
//		for (int i = 0;i < nX+1;++i)
//		{
//			if(m_alpha(openList[real_pose], i) > maxValue)
//			{
//				maxValue = m_alpha(openList[real_pose], i);
//				maxIndex = i;
//			}
//		}
//		//if (46 == maxIndex)
//		//{
//		//	cout<<"debug.."<<endl;
//		//}
//
//		if(maxIndex == nX)
//		{
//			closeList.push_back(openList[real_pose]);
//			openList.erase(openList.begin()+real_pose);
//			continue;
//		}
//		// 
//		if(correspond.inverse_map[maxIndex] != -1)
//		{
//			// assigned
//			int old_observation = correspond.inverse_map[maxIndex];
//			cvline_polar old_line = m_observe_lines_polar[old_observation];
//			cvline_polar new_line = m_observe_lines_polar[openList[real_pose]];
//			cvline_polar trans_line = m_trans_lines_polar[maxIndex];
//			fPoint center_old((old_line.pt1.x+old_line.pt2.x)*0.5, (old_line.pt1.y+old_line.pt2.y)*0.5);
//			fPoint center_new((new_line.pt1.x+new_line.pt2.x)*0.5, (new_line.pt1.y+new_line.pt2.y)*0.5);
//			fPoint center_trans((trans_line.pt1.x+trans_line.pt2.x)*0.5, (trans_line.pt1.y+trans_line.pt2.y)*0.5);
//			double dis_old = (center_old.x-center_trans.x)*(center_old.x-center_trans.x) + (center_old.y-center_trans.y)*(center_old.y-center_trans.y);
//			double dis_new = (center_new.x-center_trans.x)*(center_new.x-center_trans.x) + (center_new.y-center_trans.y)*(center_new.y-center_trans.y);
//			if(dis_new < dis_old)
//			{
//				correspond.forward_map[openList[real_pose]] = maxIndex;
//				correspond.inverse_map[maxIndex] = openList[real_pose];
//				closeList.push_back(openList[real_pose]);
//				openList.erase(openList.begin()+real_pose);
//
//				correspond.forward_map[old_observation] = nX;
//				openList.push_back(old_observation);
//				//j++;
//			}
//			else
//			{
//				int maxIndex = nX;
//				double maxValue = 0.0;
//				for (int i = 0;i < nX+1;++i)
//				{
//					if(correspond.inverse_map[i] == -1 && m_alpha(openList[real_pose], i) > maxValue)
//					{
//						maxValue = m_alpha(openList[real_pose], i);
//						maxIndex = i;
//					}
//				}
//				if(maxIndex == nX)
//				{
//					closeList.push_back(openList[real_pose]);
//					openList.erase(openList.begin()+real_pose);
//					continue;
//				}
//
//				correspond.forward_map[openList[real_pose]] = maxIndex;
//				correspond.inverse_map[maxIndex] = openList[real_pose];
//				closeList.push_back(openList[real_pose]);
//				openList.erase(openList.begin()+real_pose);
//			}
//		}
//		else
//		{
//			// not assigned
//			correspond.forward_map[openList[real_pose]] = maxIndex;
//			correspond.inverse_map[maxIndex] = openList[real_pose];
//			closeList.push_back(openList[real_pose]);
//			openList.erase(openList.begin()+real_pose);
//		}
//	}
//
//	return correspond;
//}

bool ecmlr::classification1()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nSample = 1000;
	int nDiscard = 50;
	correspondence_struct current_correspondence = getValidCorrespondence();

	Eigen::MatrixXi accumulateMatrix(nY, nX + 1);
	accumulateMatrix.fill(0);
	int iCount = 0;
	do 
	{
		//double f1 = getCorrepondenceProbility(current_correspondence);
		correspondence_struct new_correspondence = MarkovNext(current_correspondence);
		//double f2 = getCorrepondenceProbility(new_correspondence);
		double acceptance = acceptance_ratio(current_correspondence, new_correspondence);
		if(acceptance < 0.999)
		{
			// reject
			continue;
		}
		else
		{
			// accept
			current_correspondence = new_correspondence;
			iCount++;
		}
	} while (iCount < nSample);

	m_z = current_correspondence.forward_map;

	return true;
}

bool ecmlr::MetropolisHastings()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nSample = 10000;
	int nDiscard = 50;
	correspondence_struct current_correspondence = getValidCorrespondence();

	Eigen::MatrixXi accumulateMatrix(nY, nX + 1);
	accumulateMatrix.fill(0);
	int iCount = 0;
	do 
	{
		//double f1 = getCorrepondenceProbility(current_correspondence);
		correspondence_struct new_correspondence = MarkovNext(current_correspondence);
		//double f2 = getCorrepondenceProbility(new_correspondence);
		double acceptance = acceptance_ratio(current_correspondence, new_correspondence);
		if(acceptance < 0.999)
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

	//m_alpha.fill(0.0);
	//for (int j = 0;j < nY;++j)
	//{
	//	m_alpha(j, current_correspondence.forward_map[j]) = 1.0;
	//}

	return true;
}

double ecmlr::getCorrepondenceProbility(const correspondence_struct& correspond)
{
	double result = 1.0;
	for (int j = 0;j < (int)correspond.forward_map.size();++j)
	{
		result *= m_alpha(j, correspond.forward_map[j]);
	}
	return result;
}

correspondence_struct ecmlr::MarkovNext(const correspondence_struct& correspond)
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

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
			if(m_alpha(oldSource, openList[i]) > maxValue)
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


void ecmlr::updateDeltaSquare()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	double sum = 0.0;
	for(int i = 0;i < nX;++i)
	{
		for (int j = 0;j < nY;++j)
		{
			//fPoint dist = distance2(m_trans_lines_polar[i], m_observe_lines_polar[j]);
			fPoint dist = segment2polarline(m_trans_lines_polar[i], m_reference_lines_polar[j]);
			//fPoint dist;
			//dist.x = point2polarline(m_trans_lines_polar[i].pt1, m_source_lines_polar[j]);
			//dist.y = point2polarline(m_trans_lines_polar[i].pt2, m_source_lines_polar[j]);
			sum += (dist.x*dist.x + dist.y*dist.y)*m_alpha(j,i);
		}
	}

	//m_delta_2 = sum / (nY * ndim);
	//m_delta_2 = sum / nY;
	m_delta_2 = sum / (m_lammda_sum * ndim);
	if(m_delta_2 < 1.0) m_delta_2 = 1.0;
}

bool ecmlr::initialization()
{
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
	m_w.resize(m_nX);

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
			theta += m_alpha(j, i) * m_reference_lines_polar[j].theta;
			rho += m_alpha(j, i) * m_reference_lines_polar[j].rho;
		}
		m_lammda_sum += m_lambda[i];
		m_w[i].theta = theta;
		m_w[i].rho = rho;
	}

	//
	m_z.resize(m_nY);
	//m_z.set_size(nY);

	//
	m_trans_lines_polar.resize(m_nX);
	//m_transData.set_size(nX, ndim);

	return true;
}

bool ecmlr::setMatched()
{
	initAdjustableParameters();
	if (m_nY != m_nX)
	{
		return false;
	}

	//set the size of posterior matrix
	m_alpha.resize(m_nY, m_nX + 1);
	double w = 1.0 / (m_nX + 1);
	m_alpha.fill(0.0);
	for (int i = 0;i < m_nY;++i)
	{
		m_alpha(i,i) = 1.0;
	}

	// set the size of virtual observation matrix
	m_w.resize(m_nX);

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
			theta += m_alpha(j, i) * m_reference_lines_polar[j].theta;
			rho += m_alpha(j, i) * m_reference_lines_polar[j].rho;
		}
		m_lammda_sum += m_lambda[i];
		m_w[i].theta = theta;
		m_w[i].rho = rho;
	}

	//
	m_z.resize(m_nY);
	//m_z.set_size(nY);

	//
	m_trans_lines_polar.resize(m_nX);
	//m_transData.set_size(nX, ndim);

	return true;
}

Eigen::VectorXd ecmlr::getParameters() const
{
	int num = get_PARAMETER_NUM();
	Eigen::VectorXd parameters(num);
	for (int i = 0;i < num;++i)
	{
		parameters[i] = m_parameters[i];
	}
	return parameters;
}

bool ecmlr::E_Step()
{
	// E-step evaluate the posteriors

	//double outlier = 2.0 * m_delta_2 / (r * r + FLOAT_EPSILON);
	//double outlier = pow(2.0*PI*m_delta_2, ndim / 2) * nX / nY;
	//double outlier_dist = 3.0 * sqrt(m_delta_2);
	//double r = 5.0;
	//double r = outlier_dist;
	//double r = 2.0 / 9.0;
	//double outlier = 3.0 * m_delta_2 / (r * r + DBL_EPSILON);
	double outlier = 2.0 / 9.0;
	//double outlier = 2.0 / 4.0;
	for(int j = 0;j < m_nY;++j)
	{
		double sum = 0.0;
		vector<double> tmp(m_nX);
		double n_ = 1.0 / (double)m_nX;
		for (int i = 0;i < m_nX;++i)
		{
			//fPoint dist = distance2(m_observe_lines_polar[j], m_trans_lines_polar[i]);
			fPoint dist = segment2polarline(m_trans_lines_polar[i], m_reference_lines_polar[j]);
			//fPoint dist;
			//dist.x = point2polarline(m_trans_lines_polar[i].pt1, m_observe_lines_polar[j]);
			//dist.y = point2polarline(m_trans_lines_polar[i].pt2, m_observe_lines_polar[j]);
			double dist_2 = dist.x*dist.x + dist.y*dist.y;

			//////////////////////////////////////////////////////////////////////////
			double centerx1 = (m_trans_lines_polar[i].pt1.x + m_trans_lines_polar[i].pt2.x)*0.5;
			double centery1 = (m_trans_lines_polar[i].pt1.y + m_trans_lines_polar[i].pt2.y)*0.5;
			double centerx2 = (m_reference_lines_polar[j].pt1.x + m_reference_lines_polar[j].pt2.x)*0.5;
			double centery2 = (m_reference_lines_polar[j].pt1.y + m_reference_lines_polar[j].pt2.y)*0.5;
			//fPoint shift = fPoint(centerx1-centerx2, centery1-centery2);
			double shift = sqrt((centerx1-centerx2)*(centerx1-centerx2) + (centery1-centery2)*(centery1-centery2));
			//shift = max(dist.x,dist.y)*exp(-1.0/(shift*shift+DBL_EPSILON));
			shift = 0.0*exp(-m_delta_2/(shift*shift+DBL_EPSILON));
			//////////////////////////////////////////////////////////////////////////
			tmp[i] = exp(-(dist_2 + shift*shift) / (2.0 * m_delta_2));
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

	//MetropolisHastings();
	// calculate the weights of virtual observations
	// calculate the virtual observation matrix
	//m_w.fill(0.0);
	m_lammda_sum = 0.0;
	for (int i = 0;i < m_nX;++i)
	{
		m_lambda[i] = 0;
		double theta = 0.0;
		double rho = 0.0;
		double x1 = 0.0;
		double y1 = 0.0;
		double x2 = 0.0;
		double y2 = 0.0;

		for (int j = 0;j < m_nY;++j)
		{
			m_lambda[i] += m_alpha(j, i);
			theta += m_alpha(j, i) * m_reference_lines_polar[j].theta;
			rho += m_alpha(j, i) * m_reference_lines_polar[j].rho;
			x1 += m_alpha(j, i) * m_reference_lines_polar[j].pt1.x;
			y1 += m_alpha(j, i) * m_reference_lines_polar[j].pt1.y;
			x2 += m_alpha(j, i) * m_reference_lines_polar[j].pt2.x;
			y2 += m_alpha(j, i) * m_reference_lines_polar[j].pt2.y;
		}
		m_lammda_sum += m_lambda[i];
		fPoint line[2];
		line[0] = fPoint(x1, y1);
		line[1] = fPoint(x2, y2);
		m_w[i] = line2polar(line);
		//m_w[i].theta = theta;
		//m_w[i].rho = rho;
	}
	return true;
}

void ecmlr::funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	ecmlr *pThis = (ecmlr*)adata;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int pos = 0;
	int i;

	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
	}

	for(i = 0;i < nequation;++i)
	{
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		double res = pThis->m_lambda[i] * outPt.rho * outPt.rho;
		res += scale_angle * pThis->m_lambda[i] * outPt.theta * outPt.theta;
		res -= 2.0 * outPt.rho * pThis->m_w[i].rho;
		res -= scale_angle * 2.0 * outPt.theta * pThis->m_w[i].theta;
		//double dist = segment_distance(pThis->m_w[i], outPt);
		hx[pos++] = res;
	}
}

void ecmlr::jacErrorEquation(double *param, double *jac, int nparameter, int nequation, void *adata)
{
	ecmlr *pThis = (ecmlr*)adata;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

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
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		double temp_theta = 2.0 * (pThis->m_lambda[i] * outPt.theta - pThis->m_w[i].theta);
		double temp_rho = 2.0 * (pThis->m_lambda[i] * outPt.rho - pThis->m_w[i].rho);
		for(int p=0;p<nparameter;++p)
		{
			double middle = pThis->m_parameters[p];
			pThis->m_parameters[p] = middle + pstep_scale;
			cvline_polar outLine1 = pThis->forward(pThis->m_source_lines_polar[i]);

			pThis->m_parameters[p] = middle - pstep_scale;
			cvline_polar outLine2 = pThis->forward(pThis->m_source_lines_polar[i]);

			pThis->m_parameters[p] = middle;

			double derivative_theta = (outLine1.theta - outLine2.theta) * den;
			double derivative_rho = (outLine1.rho - outLine2.rho) * den;

			jac[c++] = derivative_rho * temp_rho + scale_angle * derivative_theta * temp_theta;
		}
	}
}


//void  ecmlr::alglib_function_hess(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, real_2d_array &hess, void *ptr)
//{
//	ecmlr *pThis = (ecmlr*)ptr;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	for(int i=0;i<nparameter;++i)
//	{
//		pThis->setAdjustableParameter(i, x[i], true);
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//	int c = 0;
//
//
//	for (int i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i]);
//		for (int j = 0;j < nY;++j)
//		{
//			fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
//
//			fi[c] = pThis->m_alpha(j,i) * dist.x;
//			fi[c+1] = pThis->m_alpha(j,i) * dist.y;
//
//			for(int p=0;p<nparameter;++p)
//			{
//				double middle = pThis->getAdjustableParameter(p);
//				pThis->setAdjustableParameter(p, middle + pstep_scale, true);
//				cvline_polar outLine1 = pThis->forward(pThis->m_model_lines_polar[i]);
//				fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
//
//				pThis->setAdjustableParameter(p, middle - pstep_scale, true);
//				cvline_polar outLine2 = pThis->forward(pThis->m_model_lines_polar[i]);
//				fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
//
//				pThis->setAdjustableParameter(p, middle, true);
//
//				double derivative_x = (dist1.x - dist2.x) * den;
//				double derivative_y = (dist1.y - dist2.y) * den;
//
//				jac[c][p] = pThis->m_alpha(j,i) * derivative_x;
//				jac[c+1][p] = pThis->m_alpha(j,i) * derivative_y;
//			}
//			c += 2;
//		}
//	}
//}


//double ecmlr::distance(const cvline_polar& l1, const cvline_polar& l2, Eigen::VectorXd parameters)
//{
//	double sin_theta = sin(l1.theta);
//	double cos_theta = cos(l1.theta);
//	double rho = l1.rho;
//
//	double sin_theta_prime = sin(l2.theta);
//	double cos_theta_prime = cos(l2.theta);
//	double rho_prime = l2.rho;
//
//	Eigen::VectorXd coeff(6);
//	coeff[0] = -sin_theta_prime;
//	coeff[1] = 0.0;
//	coeff[2] = -rho * sin_theta_prime;
//	coeff[3] = cos_theta_prime;
//	coeff[4] = 0.0;
//	coeff[5] = rho * cos_theta_prime;
//	double dist1 = arma::dot(coeff.t(), m_parameters) - rho_prime;
//
//	coeff[0] = 0.0;
//	coeff[1] = cos_theta * sin_theta_prime;
//	coeff[2] = sin_theta * sin_theta_prime;
//	coeff[3] = 0.0;
//	coeff[4] = -cos_theta * cos_theta_prime;
//	coeff[5] = -sin_theta * cos_theta_prime;
//	double dist2 = arma::dot(coeff.t(), m_parameters);
//
//	return dist1*dist1 + scale_angle * dist2*dist2;
//}

void ecmlr::setParameters(const vector<double> & parameters)
{
	m_parameters = parameters;
}

fPoint ecmlr::distance2(const cvline_polar& l1, const cvline_polar& l2)
{
	double centerx1 = (l1.pt1.x + l1.pt2.x)*0.5;
	double centery1 = (l1.pt1.y + l1.pt2.y)*0.5;
	double centerx2 = (l2.pt1.x + l2.pt2.x)*0.5;
	double centery2 = (l2.pt1.y + l2.pt2.y)*0.5;
	//fPoint shift = fPoint(centerx1-centerx2, centery1-centery2);
	double shift = sqrt((centerx1-centerx2)*(centerx1-centerx2) + (centery1-centery2)*(centery1-centery2));
	//double shift_weight = 1.0E-2;
	//double shift_weight = 0.08;

	fPoint dist;
	// dist1
	dist.x = distance(l1, l2);
	// dist2
	dist.y = distance(l2, l1);

	//shift = max(dist.x,dist.y)*exp(-1.0/(shift*shift+DBL_EPSILON));
	//shift = 3.0*exp(-m_delta_2/(shift*shift+DBL_EPSILON));
	shift = 0.0*exp(-1.0/(shift+DBL_EPSILON));
	//if((dist.x*dist.x+dist.y*dist.y) < 4.0)
	//{
	dist.x += shift;
	dist.y += shift;
	//}

	//dist.x = fabs(l1.rho - l2.rho);
	//dist.y = scale_angle * fabs(l1.theta - l2.theta);
	////dist.y = dist.y > lammda_rho*dist.x ? lammda_rho*dist.x : dist.y;
	////dist.y = fabs(tan(l1.theta) - tan(l2.theta));
	return dist;
}

double ecmlr::point2polarline(const fPoint& p, const cvline_polar& l)
{
	double dist = -p.x * sin(l.theta) + p.y * cos(l.theta) - l.rho;
	return dist;
}

fPoint ecmlr::segment2polarline(const cvline_polar& segment, const cvline_polar& line)
{
	return distance2(line, segment);
	fPoint dist(point2polarline(segment.pt1, line), point2polarline(segment.pt2, line));
	return dist;
}

double ecmlr::distance(const cvline_polar& l1, const cvline_polar& l2)
{
	//return (point2polarline(l1.pt1, l2) + point2polarline(l1.pt2, l2))*0.5;
	double d1 = point2polarline(l1.pt1, l2);
	double d2 = point2polarline(l1.pt2, l2);
	return sqrt(d1*d1+d2*d2);
	//return (fabs(d1)+fabs(d2))*0.5;
	// dist1
	fPoint gpt1 = l2.pt1;
	fPoint gpt2 = l2.pt2;

	double hemline = sqrt((gpt1.x - gpt2.x)*(gpt1.x - gpt2.x) + (gpt1.y - gpt2.y)*(gpt1.y - gpt2.y));
	double k1 = (gpt1.y - gpt2.y) / (hemline + DBL_EPSILON);
	double k2 = (gpt1.x - gpt2.x) / (hemline + DBL_EPSILON);

	//for the first endpoint of a line
	double dist1 = fabs((gpt2.x * gpt1.y - gpt1.x * gpt2.y) / (hemline + DBL_EPSILON) - k1 * l1.pt1.x + k2 * l1.pt1.y);

	//for the second endpoint of a line
	double dist2 = fabs((gpt2.x * gpt1.y - gpt1.x * gpt2.y) / (hemline + DBL_EPSILON) - k1 * l1.pt2.x + k2 * l1.pt2.y);

	// dist2
	return (dist1+dist2)/2.0;
}

void ecmlr::updateTransLines()
{
	for (int i = 0;i < m_nX;++i)
	{
		cvline_polar transedPt = forward(m_source_lines_polar[i]);
		m_trans_lines_polar[i] = transedPt;
	}
}

//// initial point: p[i] =  -1.2 for i odd, 1.0 for i even, minimum at (1, 1, ..., 1)
//void ecmlr::sparseLM_func(double *p, double *hx, int m, int n, void *adata)
//{
//	ecmlr *pThis = (ecmlr*)adata;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	for(int i=0;i<nparameter;++i)
//	{
//		pThis->setAdjustableParameter(i, p[i], true);
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//
//	int c = 0;
//	for (int i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i]);
//		for (int j = 0;j < nY;++j)
//		{
//			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
//			fPoint dist = pThis->segment2polarline(outPt, pThis->m_observe_lines_polar[j]);
//
//			hx[c++] = -pThis->m_alpha(j,i) * dist.x;
//			hx[c++] = -pThis->m_alpha(j,i) * dist.y;
//		}
//	}
//}
//
//
//void ecmlr::sparseLM_jac(double *p, struct splm_crsm *jac, int m, int n, void *adata)
//{
//	ecmlr *pThis = (ecmlr*)adata;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	for(int i=0;i<nparameter;++i)
//	{
//		pThis->setAdjustableParameter(i, p[i], true);
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//	int c = 0;
//
//
//	int icount = 0;
//	for (int i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i]);
//		for (int j = 0;j < nY;++j)
//		{
//			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
//			fPoint dist = pThis->segment2polarline(outPt, pThis->m_observe_lines_polar[j]);
//
//			vector<double> tmp1;
//			vector<double> tmp2;
//			for(int p=0;p<nparameter;++p)
//			{
//				double middle = pThis->getAdjustableParameter(p);
//				pThis->setAdjustableParameter(p, middle + pstep_scale, true);
//				cvline_polar outLine1 = pThis->forward(pThis->m_model_lines_polar[i]);
//				//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
//				fPoint dist1 = pThis->segment2polarline(outLine1, pThis->m_observe_lines_polar[j]);
//
//				pThis->setAdjustableParameter(p, middle - pstep_scale, true);
//				cvline_polar outLine2 = pThis->forward(pThis->m_model_lines_polar[i]);
//				//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
//				fPoint dist2 = pThis->segment2polarline(outLine2, pThis->m_observe_lines_polar[j]);
//
//				pThis->setAdjustableParameter(p, middle, true);
//
//				double derivative_x = (dist1.x - dist2.x) * den;
//				double derivative_y = (dist1.y - dist2.y) * den;
//				tmp1.push_back(pThis->m_alpha(j,i) * derivative_x);
//				tmp2.push_back(pThis->m_alpha(j,i) * derivative_y);
//
//				//jac[c][p] = pThis->m_alpha(j,i) * derivative_x;
//				//jac[c+1][p] = pThis->m_alpha(j,i) * derivative_y;
//			}
//			jac->rowptr[(i*nY+j)*2]=c;
//			for(int p=0;p<nparameter;++p)
//			{
//				if (fabs(tmp1[p]) > DBL_EPSILON)
//				{
//					jac->val[c]=tmp1[p]; jac->colidx[c++]=icount;
//				}
//				icount++;
//			}
//			jac->rowptr[(i*nY+j)*2+1]=c;
//			for(int p=0;p<nparameter;++p)
//			{
//				if (fabs(tmp2[p]) > DBL_EPSILON)
//				{
//					jac->val[c]=tmp2[p]; jac->colidx[c++]=icount;
//				}
//				icount++;
//			}
//		}
//	}
//	jac->rowptr[nX*nY*2]=c;
//}

//void  ecmlr::alglib_function_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr)
//{
//	ecmlr *pThis = (ecmlr*)ptr;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	for(int i=0;i<nparameter;++i)
//	{
//		pThis->setAdjustableParameter(i, x[i], true);
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//
//	int c = 0;
//	for (int i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i]);
//		for (int j = 0;j < nY;++j)
//		{
//			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
//			fPoint dist = pThis->segment2polarline(outPt, pThis->m_observe_lines_polar[j]);
//
//			fi[c++] = pThis->m_alpha(j,i) * dist.x;
//			fi[c++] = pThis->m_alpha(j,i) * dist.y;
//		}
//	}
//	//for (int i = 0;i < nX;++i)
//	//{
//	//	cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i]);
//	//	fPoint dist = pThis->segment2polarline(outPt, pThis->m_w[i], pThis->m_lambda[i]);
//	//	//fPoint dist = pThis->distance2(pThis->m_w[i], outPt);
//
//	//	fi[c++] = -dist.x;
//	//	fi[c++] = -dist.y;
//	//}
//}
//
//void  ecmlr::alglib_function_jac(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr)
//{
//	ecmlr *pThis = (ecmlr*)ptr;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	for(int i=0;i<nparameter;++i)
//	{
//		pThis->setAdjustableParameter(i, x[i], true);
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//	int c = 0;
//
//
//	for (int i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i]);
//		for (int j = 0;j < nY;++j)
//		{
//			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
//			fPoint dist = pThis->segment2polarline(outPt, pThis->m_observe_lines_polar[j]);
//
//			fi[c] = pThis->m_alpha(j,i) * dist.x;
//			fi[c+1] = pThis->m_alpha(j,i) * dist.y;
//
//			for(int p=0;p<nparameter;++p)
//			{
//				double middle = pThis->getAdjustableParameter(p);
//				pThis->setAdjustableParameter(p, middle + pstep_scale, true);
//				cvline_polar outLine1 = pThis->forward(pThis->m_model_lines_polar[i]);
//				//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
//				fPoint dist1 = pThis->segment2polarline(outLine1, pThis->m_observe_lines_polar[j]);
//
//				pThis->setAdjustableParameter(p, middle - pstep_scale, true);
//				cvline_polar outLine2 = pThis->forward(pThis->m_model_lines_polar[i]);
//				//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
//				fPoint dist2 = pThis->segment2polarline(outLine2, pThis->m_observe_lines_polar[j]);
//
//				pThis->setAdjustableParameter(p, middle, true);
//
//				double derivative_x = (dist1.x - dist2.x) * den;
//				double derivative_y = (dist1.y - dist2.y) * den;
//
//				jac[c][p] = pThis->m_alpha(j,i) * derivative_x;
//				jac[c+1][p] = pThis->m_alpha(j,i) * derivative_y;
//			}
//			c += 2;
//		}
//	}
//	//for (int i = 0;i < nX;++i)
//	//{
//	//	cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i]);
//
//	//	//fPoint dist = pThis->distance2(pThis->m_w[i], outPt);
//	//	fPoint dist = pThis->segment2polarline(outPt, pThis->m_w[i]);
//
//	//	fi[c] = pThis->m_lambda[i] * dist.x;
//	//	fi[c+1] = pThis->m_lambda[i] * dist.y;
//
//	//	for(int p=0;p<nparameter;++p)
//	//	{
//	//		double middle = pThis->getAdjustableParameter(p);
//	//		pThis->setAdjustableParameter(p, middle + pstep_scale, true);
//	//		cvline_polar outLine1 = pThis->forward(pThis->m_model_lines_polar[i]);
//	//		//fPoint dist1 = pThis->distance2(pThis->m_w[i], outLine1);
//	//		fPoint dist1 = pThis->segment2polarline(outLine1, pThis->m_w[i], pThis->m_lambda[i]);
//
//	//		pThis->setAdjustableParameter(p, middle - pstep_scale, true);
//	//		cvline_polar outLine2 = pThis->forward(pThis->m_model_lines_polar[i]);
//	//		//fPoint dist2 = pThis->distance2(pThis->m_w[i], outLine2);
//	//		fPoint dist2 = pThis->segment2polarline(outLine2, pThis->m_w[i], pThis->m_lambda[i]);
//
//	//		pThis->setAdjustableParameter(p, middle, true);
//
//	//		double derivative_x = (dist1.x - dist2.x) * den;
//	//		double derivative_y = (dist1.y - dist2.y) * den;
//
//	//		jac[c][p] = derivative_x;
//	//		jac[c+1][p] = derivative_y;
//	//	}
//	//	c += 2;
//	//}
//}

//void ecmlr::function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr) 
//{
//	ecmlr *pThis = (ecmlr*)ptr;
//	int nX = (int)pThis->m_nX;
//	int nY = (int)pThis->m_nY;
//	int ndim = (int)pThis->m_ndim;
//
//	int pos = 0;
//	int i;
//
//	int nparameter = pThis->get_PARAMETER_NUM();
//
//	func = 0.0;
//	for(i=0;i<nparameter;++i)
//	{
//		pThis->m_parameters[i] = x[i];
//		grad[i] = 0.0;
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//	int c = 0;
//
//	for(i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i]);
//		double res = pThis->m_lambda[i] * outPt.rho * outPt.rho;
//		res += scale_angle * pThis->m_lambda[i] * outPt.theta * outPt.theta;
//		res -= 2.0 * outPt.rho * pThis->m_w[i].rho;
//		res -= scale_angle * 2.0 * outPt.theta * pThis->m_w[i].theta;
//		//double dist = segment_distance(pThis->m_w[i], outPt);
//		func += res;
//
//		double temp_theta = 2.0 * (pThis->m_lambda[i] * outPt.theta - pThis->m_w[i].theta);
//		double temp_rho = 2.0 * (pThis->m_lambda[i] * outPt.rho - pThis->m_w[i].rho);
//		for(int p=0;p<nparameter;++p)
//		{
//			double middle = pThis->m_parameters[p];
//			pThis->m_parameters[p] = middle + pstep_scale;
//			cvline_polar outLine1 = pThis->forward(pThis->m_model_lines_polar[i]);
//
//			pThis->m_parameters[p] = middle - pstep_scale;
//			cvline_polar outLine2 = pThis->forward(pThis->m_model_lines_polar[i]);
//
//			pThis->m_parameters[p] = middle;
//
//			double derivative_theta = (outLine1.theta - outLine2.theta) * den;
//			double derivative_rho = (outLine1.rho - outLine2.rho) * den;
//
//			grad[p] += derivative_rho * temp_rho + scale_angle * derivative_theta * temp_theta;
//		}
//	}
//	double tmp = func;
//}

//bool ecmlr::parameters_optimization()
//{
//
//	int nX = m_nX;
//	int nY = m_nY;
//	int ndim = m_ndim;
//
//	int nparam = get_PARAMETER_NUM();
//	std::vector<double> cparm(nparam);
//
//	for(int i=0;i<nparam;++i)
//	{
//		cparm[i] = m_parameters[i];
//	}
//
//	//Eigen::VectorXd outParameters(nparam);
//	//double *p = &cparm[0];
//	//double *x = new double[nX];
//	//for(int i=0; i<nX; i++) x[i]=0.0;
//
//	//double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
//	//opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
//	//opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 
//
//	////int ret = dlevmar_dif(funcErrorEquation, p, x, nparam, nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
//	//int ret = dlevmar_der(funcErrorEquation, jacErrorEquation, p, x, nparam, nX, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian
//
//
//	nlopt::opt opt(nlopt::LD_VAR1, nparam);
//
//	//std::vector<double> lb(nparam);
//	//lb[0] = -HUGE_VAL; lb[1] = 0;
//	//opt.set_lower_bounds(lb);
//
//	opt.set_min_objective(myfunc, this);
//
//	//my_constraint_data data[2] = { {2,0}, {-1,1} };
//	//opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
//	//opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);
//
//	opt.set_xtol_rel(1e-4);
//
//	std::vector<double> x(nparam);
//	for(int i=0;i<nparam;++i)
//	{
//		x[i] = m_parameters[i];
//	}
//	double minf;
//	nlopt::result result = opt.optimize(x, minf);
//
//	for(int i=0;i<nparam;++i)
//	{
//		m_parameters[i] = x[i];
//	}
//
//	double sum = 0.0;
//	for (int i = 0;i < nX;++i)
//	{
//		cvline_polar outPt = forward(m_model_lines_polar[i]);
//		double a1 = m_lambda[i]*(outPt.rho*outPt.rho + scale_angle*outPt.theta*outPt.theta);
//		double a2 = 2.0 * (outPt.rho*m_w[i].rho + scale_angle*outPt.theta*m_w[i].theta);
//		sum += a1 - a2;
//	}
//	for(int j = 0;j < nY;++j)
//	{
//		double a3 = (1.0 - m_alpha(j,nX)) * (m_observe_lines_polar[j].rho*m_observe_lines_polar[j].rho + scale_angle*m_observe_lines_polar[j].theta*m_observe_lines_polar[j].theta);
//		sum += a3;
//	}
//	//Eigen::VectorXd delta = b - A * m_parameters;
//	cout<<"rms:\t"<<sum<<endl;
//	m_delta_2 = sum / (nY * ndim);
//
//	return true;
//}


void ecmlr::levmar_function_fvec(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	ecmlr *pThis = (ecmlr*)adata;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int pos = 0;
	int i;
	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
		//pThis->m_parameters[i] = param[i];
	}

	int c = 0;
	//fstream ofs;
	//ofs.open("b.txt", ios_base::out);
	//for(i = 0;i < nX;++i)
	//{
	//	cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);
	//	double delta_theta = pThis->m_w[i].theta - pThis->m_lambda[i] * outPt.theta;
	//	double delta_rho = pThis->m_w[i].rho - pThis->m_lambda[i] * outPt.rho;
	//	hx[c++] = delta_rho;
	//	hx[c++] = scale_angle * delta_theta;
	//	ofs<<-delta_rho<<endl;
	//	ofs<<-scale_angle * delta_theta<<endl;
	//}
	//ofs.close();

	//for (int i = 0;i < nX;++i)
	//{
	//	cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);
	//	for (int j = 0;j < nY;++j)
	//	{
	//		double delta_theta = pThis->m_observe_lines_polar[j].theta - outPt.theta;
	//		double delta_rho = pThis->m_observe_lines_polar[j].rho - outPt.rho;
	//		//func += pThis->m_alpha(j,i) * (delta_rho*delta_rho + scale_angle*delta_theta*delta_theta);
	//		hx[c++] = -pThis->m_alpha(j,i) * delta_rho;
	//		hx[c++] = -pThis->m_alpha(j,i) * scale_angle * delta_theta;
	//	}
	//}

	for (int i = 0;i < nX;++i)
	{
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		for (int j = 0;j < nY;++j)
		{
			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
			fPoint dist = pThis->segment2polarline(outPt, pThis->m_reference_lines_polar[j]);

			hx[c++] = pThis->m_alpha(j,i) * dist.x;
			hx[c++] = pThis->m_alpha(j,i) * dist.y;
		}
	}
}

void ecmlr::levmar_function_jac(double *param, double *jac, int nparameter, int nequation, void *adata)
{
	ecmlr *pThis = (ecmlr*)adata;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int pos = 0;
	int i;
	for(i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = param[i];
		//pThis->m_parameters[i] = param[i];
	}

	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;
	int c = 0;

	for (int i = 0;i < nX;++i)
	{
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		for (int j = 0;j < nY;++j)
		{
			for(int p=0;p<nparameter;++p)
			{
				double middle = pThis->m_parameters[p];
				pThis->m_parameters[p] = middle + pstep_scale;
				cvline_polar outLine1 = pThis->forward(pThis->m_source_lines_polar[i]);
				//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
				fPoint dist1 = pThis->segment2polarline(outLine1, pThis->m_reference_lines_polar[j]);

				pThis->m_parameters[p] = middle - pstep_scale;
				cvline_polar outLine2 = pThis->forward(pThis->m_source_lines_polar[i]);
				//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
				fPoint dist2 = pThis->segment2polarline(outLine2, pThis->m_reference_lines_polar[j]);

				pThis->m_parameters[p] = middle;

				double derivative_x = (dist1.x - dist2.x) * den;
				double derivative_y = (dist1.y - dist2.y) * den;

				jac[c*nparameter+p] = pThis->m_alpha(j,i) * derivative_x;
				jac[(c+1)*nparameter+p] = pThis->m_alpha(j,i) * derivative_y;
			}
			c += 2;
		}
	}

	////vector< vector<double> > test;
	//fstream ofs;
	//ofs.open("A.txt", ios_base::out);
	//for(i = 0;i < nX;++i)
	//{
	//	vector<double> t1(nparameter);
	//	vector<double> t2(nparameter);
	//	for(int p=0;p<nparameter;++p)
	//	{
	//		double middle = pThis->m_parameters[p];
	//		pThis->m_parameters[p] = middle + pstep_scale;
	//		cvline_polar outLine1 = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);

	//		pThis->m_parameters[p] = middle - pstep_scale;
	//		cvline_polar outLine2 = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);

	//		pThis->m_parameters[p] = middle;

	//		double derivative_theta = (outLine1.theta - outLine2.theta) * den;
	//		double derivative_rho = (outLine1.rho - outLine2.rho) * den;

	//		jac[c*nparameter + p] = pThis->m_lambda[i] * derivative_rho;// * parameter_scale[p];
	//		jac[(c+1)*nparameter + p] = pThis->m_lambda[i] * scale_angle * derivative_theta;// * parameter_scale[p];

	//		t1[p] = pThis->m_lambda[i] * derivative_rho;// * parameter_scale[p];
	//		t2[p] = pThis->m_lambda[i] * scale_angle * derivative_theta;// * parameter_scale[p];
	//	}
	//	//test.push_back(t1);
	//	//test.push_back(t2);
	//	for (int p = 0;p < nparameter;++p)
	//	{
	//		ofs<<t1[p]<<"\t";
	//	}
	//	ofs<<endl;
	//	for (int p = 0;p < nparameter;++p)
	//	{
	//		ofs<<t2[p]<<"\t";
	//	}
	//	ofs<<endl;
	//	c += 2;
	//}
	//ofs.close();
	////cout<<"test"<<endl;
	////for (int i = 0;i < nX;++i)
	////{
	////	//cvline_polar outPt = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);
	////	for(int p=0;p<nparameter;++p)
	////	{
	////		double middle = pThis->m_parameters[p];
	////		pThis->m_parameters[p] = middle + pstep_scale;
	////		cvline_polar outLine1 = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);

	////		pThis->m_parameters[p] = middle - pstep_scale;
	////		cvline_polar outLine2 = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);

	////		pThis->m_parameters[p] = middle;

	////		double derivative_theta = (outLine1.theta - outLine2.theta) * den;
	////		double derivative_rho = (outLine1.rho - outLine2.rho) * den;

	////		jac[c*nparameter + p] = pThis->m_alpha(j,i) * derivative_rho;
	////		jac[(c+1)*nparameter + p] = pThis->m_alpha(j,i) * scale_angle * derivative_theta;
	////	}
	////	for (int j = 0;j < nY;++j)
	////	{
	////		//double delta_theta = pThis->m_observe_lines_polar[j].theta - outPt.theta;
	////		//double delta_rho = pThis->m_observe_lines_polar[j].rho - outPt.rho;

	////		for(int p=0;p<nparameter;++p)
	////		{
	////			double middle = pThis->m_parameters[p];
	////			pThis->m_parameters[p] = middle + pstep_scale;
	////			cvline_polar outLine1 = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);

	////			pThis->m_parameters[p] = middle - pstep_scale;
	////			cvline_polar outLine2 = pThis->forward(pThis->m_model_lines_polar[i], pThis->m_parameters);

	////			pThis->m_parameters[p] = middle;

	////			double derivative_theta = (outLine1.theta - outLine2.theta) * den;
	////			double derivative_rho = (outLine1.rho - outLine2.rho) * den;

	////			jac[c*nparameter + p] = pThis->m_alpha(j,i) * derivative_rho;
	////			jac[(c+1)*nparameter + p] = pThis->m_alpha(j,i) * scale_angle * derivative_theta;
	////		}
	////		c += 2;
	////	}
	////}
}

bool ecmlr::parameters_optimization()
{
	// levmar
	//return levmar_optimization();

	// spase LM
	//return spaseLM_optimization();

	// alglib improved LM
	//return alglib_optimization();

	// Eigen
	return eigen_levmar_optimization();


	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nparam = get_PARAMETER_NUM();
	std::vector<double> cparm(nparam);

	for(int i=0;i<nparam;++i) cparm[i] = m_parameters[i];

	//Eigen::VectorXd outParameters(nparam);
	//// lm
	//double *p = &cparm[0];
	////double *x = new double[2*nX];
	////vector<double> x(2*nX, 0.0);
	////for(int i=0; i<nX*nY; i++) x[i]=0.0;

	//double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	//opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	//opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	////int ret = dlevmar_dif(levmar_function_fvec, &cparm[0], &x[0], nparam, 2*nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
	//int ret = dlevmar_der(levmar_function_fvec, levmar_function_jac, p, NULL, nparam, 2*nX*nY, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian
	////delete []x;

	
	return true;
}

//bool ecmlr::alglib_optimization()
//{
//	int nX = m_nX;
//	int nY = m_nY;
//	int ndim = m_ndim;
//
//	int nparam = get_PARAMETER_NUM();
//	std::vector<double> cparm(nparam);
//
//	for(int i=0;i<nparam;++i) cparm[i] = getAdjustableParameter(i);
//
//	//Eigen::VectorXd outParameters(nparam);
//	//// lm
//	//double *p = &cparm[0];
//	////double *x = new double[2*nX];
//	////vector<double> x(2*nX, 0.0);
//	////for(int i=0; i<nX*nY; i++) x[i]=0.0;
//
//	//double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
//	//opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
//	//opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 
//
//	////int ret = dlevmar_dif(levmar_function_fvec, &cparm[0], &x[0], nparam, 2*nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
//	//int ret = dlevmar_der(levmar_function_fvec, levmar_function_jac, p, NULL, nparam, 2*nX*nY, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian
//	////delete []x;
//
//
//	// alglib
//	real_1d_array x;
//	x.setcontent(nparam, &cparm[0]);
//	double epsg = 0.0000000001;
//	double epsf = 0;
//	double epsx = 0;
//	ae_int_t maxits = 0;
//	minlmstate state;
//	minlmreport rep;
//
//	minlmcreatevj(2*nX*nY, x, state);
//	minlmsetcond(state, epsg, epsf, epsx, maxits);
//	//minlmsetacctype(state, 1);
//	alglib::minlmoptimize(state, alglib_function_fvec, alglib_function_jac, NULL, this);
//	minlmresults(state, x, rep);
//	for(int i=0;i<nparam;++i) setAdjustableParameter(i, x[i], true);//m_parameters[i] = x[i];
//	return true;
//}

//bool ecmlr::spaseLM_optimization()
//{
//	int nX = m_nX;
//	int nY = m_nY;
//	int ndim = m_ndim;
//
//	int nparam = get_PARAMETER_NUM();
//	std::vector<double> cparm(nparam);
//
//	for(int i=0;i<nparam;++i) cparm[i] = getAdjustableParameter(i);
//
//	double opts[SPLM_OPTS_SZ], info[SPLM_INFO_SZ];
//	double *p = &cparm[0];
//	int ret;
//
//	opts[0]=SPLM_INIT_MU; opts[1]=SPLM_STOP_THRESH; opts[2]=SPLM_STOP_THRESH;
//	opts[3]=SPLM_STOP_THRESH;
//	opts[4]=SPLM_DIFF_DELTA; // relevant only if finite difference approximation to Jacobian is used
//	opts[5]=SPLM_CHOLMOD; // use CHOLMOD
//	//opts[5]=SPLM_PARDISO; // use PARDISO
//
//	ret=sparselm_dercrs(sparseLM_func, sparseLM_jac, p, NULL, nparam, 0, 2*nX*nY, 2*nparam*nX*nY, -1, 1000, opts, info, this); // CRS Jacobian
//	for(int i=0;i<nparam;++i) setAdjustableParameter(i, p[i], true);//m_parameters[i] = x[i];
//
//	return true;
//}

bool ecmlr::levmar_optimization()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	int nparam = get_PARAMETER_NUM();
	std::vector<double> cparm(nparam);

	for(int i=0;i<nparam;++i) cparm[i] = m_parameters[i];

	// lm
	double *p = &cparm[0];
	//double *x = new double[2*nX];
	//vector<double> x(2*nX, 0.0);
	//for(int i=0; i<nX*nY; i++) x[i]=0.0;

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	//int ret = dlevmar_dif(levmar_function_fvec, &cparm[0], &x[0], nparam, 2*nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
	int ret = dlevmar_der(levmar_function_fvec, levmar_function_jac, p, NULL, nparam, 2*nX*nY, 1000, opts, info, NULL, NULL, this); // with analytic Jacobian
	//delete []x;

	for(int i=0;i<nparam;++i) m_parameters[i] = p[i];
	return true;
}

bool ecmlr::eigen_levmar_optimization()
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
	lmder_functor functor(nparam, 2*nX*nY);
	functor.pData = this;
	Eigen::LevenbergMarquardt<lmder_functor> lm(functor);
	info = lm.minimize(x);
	
	for(int i=0;i<nparam;++i) m_parameters[i] = x[i];
	return true;
}

//bool ecmlr::ceres_optimization()
//{
//	int nX = m_nX;
//	int nY = m_nY;
//	int ndim = m_ndim;
//
//	int nparam = get_PARAMETER_NUM();
//	std::vector<double> cparm(nparam);
//
//	for(int i=0;i<nparam;++i) cparm[i] = getAdjustableParameter(i);
//	
//	double *p = &cparm[0];
//	ceres::Problem problem;
//	// The problem object takes ownership of the newly allocated
//	// SimpleCostFunction and uses it to optimize the value of x.
//	problem.AddResidualBlock(new LineMatchCost(this), NULL, p);
//	// Run the solver!
//	Solver::Options options;
//	options.max_num_iterations = 10;
//	options.linear_solver_type = ceres::DENSE_QR;
//	options.minimizer_progress_to_stdout = true;
//	Solver::Summary summary;
//	Solve(options, &problem, &summary);
//	std::cout << summary.BriefReport() << "\n";
//	std::cout << "x : 5.0 -> " << p << "\n";
//	return 0;
//
//	for(int i=0;i<nparam;++i) setAdjustableParameter(i, p[i], true);
//	return true;
//}


int lmder_functor::operator()(const VectorXd &x, VectorXd &fvec) const
{
	ecmlr *pThis = (ecmlr*)pData;
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
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		for (int j = 0;j < nY;++j)
		{
			fPoint dist = pThis->segment2polarline(outPt, pThis->m_reference_lines_polar[j]);

			fvec[c++] = pThis->m_alpha(j,i) * dist.x;
			fvec[c++] = pThis->m_alpha(j,i) * dist.y;
		}
	}
	return 0;
}

int lmder_functor::df(const VectorXd &x, MatrixXd &fjac) const
{
	ecmlr *pThis = (ecmlr*)pData;
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
		cvline_polar outPt = pThis->forward(pThis->m_source_lines_polar[i]);
		for (int j = 0;j < nY;++j)
		{
			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
			fPoint dist = pThis->segment2polarline(outPt, pThis->m_reference_lines_polar[j]);

			for(int p=0;p<nparameter;++p)
			{
				double middle = pThis->m_parameters[p];
				pThis->m_parameters[p] = middle + pstep_scale;
				cvline_polar outLine1 = pThis->forward(pThis->m_source_lines_polar[i]);
				//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
				fPoint dist1 = pThis->segment2polarline(outLine1, pThis->m_reference_lines_polar[j]);

				pThis->m_parameters[p] = middle - pstep_scale;
				cvline_polar outLine2 = pThis->forward(pThis->m_source_lines_polar[i]);
				//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
				fPoint dist2 = pThis->segment2polarline(outLine2, pThis->m_reference_lines_polar[j]);

				pThis->m_parameters[p] = middle;

				double derivative_x = (dist1.x - dist2.x) * den;
				double derivative_y = (dist1.y - dist2.y) * den;

				fjac(c,p) = pThis->m_alpha(j,i) * derivative_x;
				fjac(c+1,p) = pThis->m_alpha(j,i) * derivative_y;
			}
			c += 2;
		}
	}
	return 0;
}
