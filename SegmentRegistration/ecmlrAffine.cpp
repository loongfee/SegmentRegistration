#include "ecmlrAffine.h"
#include "alglib/interpolation.h"

using namespace alglib;

static const char*  a0_kw				= "a0";
static const char*  a1_kw      = "a1";
static const char*  a2_kw       = "a2";
static const char*  b0_kw       = "b0";
static const char*  b1_kw         = "b1";
static const char*  b2_kw         = "b2";
static const char* PARAM_NAMES[] ={"a0",
									"a1",
									"a2",
									"b0",
									"b1",
									"b2"};
static const char* PARAM_UNITS[] ={"meters",
									"meters",
									"meters",
									"meters",
									"degrees",
									"degrees"};
cvline_polar ecmlrAffine::forward(const cvline_polar& l, Eigen::VectorXd parameters)const
{
	cvline_polar outLine;
	fPoint line[2];
	line[0].x = parameters[0] + parameters[1] * l.pt1.x + parameters[2] * l.pt1.y;
	line[0].y = parameters[3] + parameters[4] * l.pt1.x + parameters[5] * l.pt1.y;
	line[1].x = parameters[0] + parameters[1] * l.pt2.x + parameters[2] * l.pt2.y;
	line[1].y = parameters[3] + parameters[4] * l.pt2.x + parameters[5] * l.pt2.y;
	outLine = line2polar(line);
	return outLine;
}

cvline_polar ecmlrAffine::forward(const cvline_polar& l)const
{
	Eigen::VectorXd parameters = getParameters();
	return forward(l, parameters);
}


void ecmlrAffine::initAdjustableParameters()
{
	m_parameters.resize(NUM_ADJUSTABLE_PARAMS);
	int numParams = (int)m_parameters.size();

	m_parameters[a0_param] = 0.0;	//a0
	m_parameters[a1_param] = 1.0;	//a1
	m_parameters[a2_param] = 0.0;	//a2
	m_parameters[b0_param] = 0.0;	//b0
	m_parameters[b1_param] = 0.0;	//b1
	m_parameters[b2_param] = 1.0;	//b2
}

bool ecmlrAffine::parameters_optimization()
{
	//return ecmlr::parameters_optimization();
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
		double x1 = m_source_lines_polar[i].pt1.x;
		double y1 = m_source_lines_polar[i].pt1.y;
		double x2 = m_source_lines_polar[i].pt2.x;
		double y2 = m_source_lines_polar[i].pt2.y;
		
		for (int j = 0;j < nY;++j)
		{
			double sin_theta = sin(m_reference_lines_polar[j].theta);
			double cos_theta = cos(m_reference_lines_polar[j].theta);
			double rho = m_reference_lines_polar[j].rho;
			double len = segment_length(m_reference_lines_polar[j]);
			double weight = m_alpha(j,i);

			vector<double> coeff(6);
			coeff[0] = -sin_theta;
			coeff[1] = -sin_theta * x1;
			coeff[2] = -sin_theta * y1;
			coeff[3] = cos_theta;
			coeff[4] = cos_theta * x1;
			coeff[5] = cos_theta * y1;

			for(int p1=0;p1<6;++p1)
			{        
				b[p1] += coeff[p1] * rho * weight;
				for(int p2=0;p2<6;++p2)
				{
					A(p1,p2) += coeff[p1] * coeff[p2] * weight;
				}
			}

			coeff[0] = -sin_theta;
			coeff[1] = -sin_theta * x2;
			coeff[2] = -sin_theta * y2;
			coeff[3] = cos_theta;
			coeff[4] = cos_theta * x2;
			coeff[5] = cos_theta * y2;
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

//bool ecmlrAffine::parameters_optimization()
//{
//	//return ecmlr::parameters_optimization();
//	int nX = m_nX;
//	int nY = m_nY;
//	int ndim = m_ndim;
//
//	VectorXd b(6);
//	MatrixXd A(6, 6);
//	b.fill(0.0);
//	A.fill(0.0);
//	//b = candidates.col(6);
//	real_2d_array fmatrix;
//	real_1d_array y;
//	real_1d_array w;
//	fmatrix.setlength(2*nX*nY, 6);
//	y.setlength(2*nX*nY);
//	w.setlength(2*nX*nY);
//	int c = 0;
//	for(int i = 0;i < (int)nX;++i)
//	{
//		double x1 = m_model_lines_polar[i].pt1.x;
//		double y1 = m_model_lines_polar[i].pt1.y;
//		double x2 = m_model_lines_polar[i].pt2.x;
//		double y2 = m_model_lines_polar[i].pt2.y;
//		
//		for (int j = 0;j < nY;++j)
//		{
//			double sin_theta = sin(m_observe_lines_polar[j].theta);
//			double cos_theta = cos(m_observe_lines_polar[j].theta);
//			double rho = m_observe_lines_polar[j].rho;
//			double len = segment_length(m_observe_lines_polar[j]);
//			double weight = m_alpha(j,i);
//
//			vector<double> coeff(6);
//			coeff[0] = -sin_theta;
//			coeff[1] = -sin_theta * x1;
//			coeff[2] = -sin_theta * y1;
//			coeff[3] = cos_theta;
//			coeff[4] = cos_theta * x1;
//			coeff[5] = cos_theta * y1;
//			fmatrix(c, 0) = -sin_theta;
//			fmatrix(c, 1) = -sin_theta * x1;
//			fmatrix(c, 2) = -sin_theta * y1;
//			fmatrix(c, 3) = cos_theta;
//			fmatrix(c, 4) = cos_theta * x1;
//			fmatrix(c, 5) = cos_theta * y1;
//			y(c) = rho;
//			w(c++) = weight;
//
//			for(int p1=0;p1<6;++p1)
//			{        
//				b[p1] += coeff[p1] * rho * weight;
//				for(int p2=0;p2<6;++p2)
//				{
//					A(p1,p2) += coeff[p1] * coeff[p2] * weight;
//				}
//				//for(int p2=p1;p2>=0;--p2)
//				//{
//				//	A(p1,p2) += coeff[p1] * coeff[p2] * weight;
//				//}
//			}
//
//			coeff[0] = -sin_theta;
//			coeff[1] = -sin_theta * x2;
//			coeff[2] = -sin_theta * y2;
//			coeff[3] = cos_theta;
//			coeff[4] = cos_theta * x2;
//			coeff[5] = cos_theta * y2;
//			fmatrix(c, 0) = -sin_theta;
//			fmatrix(c, 1) = -sin_theta * x2;
//			fmatrix(c, 2) = -sin_theta * y2;
//			fmatrix(c, 3) = cos_theta;
//			fmatrix(c, 4) = cos_theta * x2;
//			fmatrix(c, 5) = cos_theta * y2;
//			y(c) = rho;
//			w(c++) = weight;
//			for(int p1=0;p1<6;++p1)
//			{        
//				b[p1] += coeff[p1] * rho * weight;
//				for(int p2=0;p2<6;++p2)
//				{
//					A(p1,p2) += coeff[p1] * coeff[p2] * weight;
//				}
//				//for(int p2=p1;p2>=0;--p2)
//				//{
//				//	A(p1,p2) += coeff[p1] * coeff[p2] * weight;
//				//}
//			}
//		}
//	}
//
//	//ae_int_t info;
//	//real_1d_array p;
//	//lsfitreport rep;
//	//lsfitlinearw(y, w, fmatrix, info, p, rep);
//	//cout<<A<<endl;
//	Eigen::MatrixXd damper = MatrixXd::Identity(6, 6);
//	//A = A + damper * 1.0E-6;
//	double tmp = A.determinant();
//	if ( fabs(tmp) > DBL_EPSILON )
//	//if ( ::solve(parameter, A, b) )
//	{
//		Eigen::VectorXd parameter = A.inverse() * b;
//		for(int i=0;i<get_PARAMETER_NUM();++i) setAdjustableParameter(i, parameter[i], true);
//		//for(int i=0;i<get_PARAMETER_NUM();++i) setAdjustableParameter(i, p[i], true);
//		return true;
//	}
//	else
//	{
//		cout << "ERROR : cannot find a solution" << endl;
//		return false;
//	}
//
//	//int nparam = get_PARAMETER_NUM();
//	//std::vector<double> cparm(nparam);
//
//	//for(int i=0;i<nparam;++i) cparm[i] = getAdjustableParameter(i);
//
//	//// alglib
//	//real_1d_array x;
//	//x.setcontent(nparam, &cparm[0]);
//	//double epsg = 0.0000000001;
//	//double epsf = 0;
//	//double epsx = 0;
//	//ae_int_t maxits = 0;
//	//minlmstate state;
//	//minlmreport rep;
//
//	//minlmcreatefgh(x, state);
//	//minlmsetcond(state, epsg, epsf, epsx, maxits);
//	//alglib::minlmoptimize(state, alglib_function_func, alglib_function_grad, alglib_function_hess, NULL, this);
//	//minlmresults(state, x, rep);
//
//	//for(int i=0;i<nparam;++i) setAdjustableParameter(i, x[i], true);
//
//	//int nX = m_nX;
//	//int nY = m_nY;
//	//int ndim = m_ndim;
//
//	//vec b(6);
//	//mat A(6, 6);
//	//b.fill(0.0);
//	//A.fill(0.0);
//	////b = candidates.col(6);
//	////for(int i = 0;i < nX;++i)
//	////{
//	////	double sin_theta = sin(m_model_lines_polar[i].theta * m_lambda[i]);
//	////	double cos_theta = cos(m_model_lines_polar[i].theta * m_lambda[i]);
//	////	double rho = m_model_lines_polar[i].rho * m_lambda[i];
//
//	////	double sin_theta_prime = sin(m_w[i].theta);
//	////	double cos_theta_prime = cos(m_w[i].theta);
//	////	double rho_prime = m_w[i].rho;
//
//	////	A(2*i + 0, 0) = -sin_theta_prime;
//	////	A(2*i + 0, 2) = -rho * sin_theta_prime;
//	////	A(2*i + 0, 3) = cos_theta_prime;
//	////	A(2*i + 0, 5) = rho * cos_theta_prime;
//
//	////	A(2*i + 1, 1) = cos_theta * sin_theta_prime * scale_angle;
//	////	A(2*i + 1, 2) = sin_theta * sin_theta_prime * scale_angle;
//	////	A(2*i + 1, 4) = -cos_theta * cos_theta_prime * scale_angle;
//	////	A(2*i + 1, 5) = -sin_theta * cos_theta_prime * scale_angle;
//
//	////	b(2*i) = rho_prime;
//	////	b(2*i+1) = 0.0;
//	////}
//
//	//// firstly, calculate a0, a2, b0, b2
//	//for(int i = 0;i < nX;++i)
//	//{
//	//	//double sin_theta = sin(m_model_lines_polar[i].theta);
//	//	//double cos_theta = cos(m_model_lines_polar[i].theta);
//	//	//double rho = m_model_lines_polar[i].rho;
//	//	double x1 = m_model_lines_polar[i].x1;
//	//	double y1 = m_model_lines_polar[i].y1;
//	//	double x2 = m_model_lines_polar[i].x2;
//	//	double y2 = m_model_lines_polar[i].y2;
//	//	for (int j = 0;j < nY;++j)
//	//	{
//	//		double sin_theta_prime = sin(m_observe_lines_polar[j].theta);
//	//		double cos_theta_prime = cos(m_observe_lines_polar[j].theta);
//	//		double rho_prime = m_observe_lines_polar[j].rho;
//
//	//		vector<double> coeff(6);
//	//		coeff[0] = -sin_theta_prime;
//	//		coeff[1] = -sin_theta_prime * x1;
//	//		coeff[2] = -sin_theta_prime * y1;
//	//		coeff[3] = cos_theta_prime;
//	//		coeff[4] = cos_theta_prime * x1;
//	//		coeff[5] = cos_theta_prime * y1;
//	//		for(int p1=0;p1<6;++p1)
//	//		{        
//	//			b[p1] += coeff[p1] * rho_prime * m_alpha(j,i);
//	//			for(int p2=p1;p2<6;++p2)
//	//			{
//	//				A(p1,p2) += coeff[p1] * coeff[p2] * m_alpha(j,i);
//	//			}
//	//		}
//
//	//		coeff[0] = -sin_theta_prime;
//	//		coeff[1] = -sin_theta_prime * x2;
//	//		coeff[2] = -sin_theta_prime * y2;
//	//		coeff[3] = cos_theta_prime;
//	//		coeff[4] = cos_theta_prime * x2;
//	//		coeff[5] = cos_theta_prime * y2;
//	//		for(int p1=0;p1<6;++p1)
//	//		{        
//	//			b[p1] += coeff[p1] * rho_prime * m_alpha(j,i);
//	//			for(int p2=p1;p2<6;++p2)
//	//			{
//	//				A(p1,p2) += coeff[p1] * coeff[p2] * m_alpha(j,i);
//	//			}
//	//		}
//	//	}
//	//}
//	//Eigen::MatrixXd damper(6, 6);
//	//damper.eye();
//	//A = A + damper * 1.0E-6;
//	//double tmp = det(A);
//	//if ( fabs(tmp) > DBL_EPSILON )
//	//{
//	//	Eigen::VectorXd parameters = arma::inv(A) * b;
//	//	for (int i = 0;i < NUM_ADJUSTABLE_PARAMS;++i)
//	//	{
//	//		setAdjustableParameter(i, parameters[i], true);
//	//	}
//
//	//	//double sum = 0.0;
//	//	//for(int i = 0;i < nX;++i)
//	//	//{
//	//	//	for (int j = 0;j < nY;++j)
//	//	//	{
//	//	//		fPoint dist = distance2(m_model_lines_polar[i], m_observe_lines_polar[j]);
//	//	//		sum += (dist.x*dist.x + dist.y*dist.y)*m_alpha(j,i);
//	//	//	}
//	//	//}
//	//	//// debug
//	//	//cout<<"rms:\t"<<sum<<endl;
//	//	//m_delta_2 = sum / (nY * ndim);
//	//	return true;		
//	//}
//	//else
//	//{
//	//	cout << "ERROR : cannot find a solution" << endl;
//	//	return false;
//	//}
//
//	//double tmp = det(A);
//	//if ( fabs(tmp) > DBL_EPSILON )
//	//{
//	//	Eigen::VectorXd tmp1 = arma::inv(A) * b;
//
//	//	// secondly, calculate a1, b1
//	//	b.set_size(2);
//	//	A.set_size(2, 2);
//	//	b.fill(0.0);
//	//	A.fill(0.0);
//	//	for(int i = 0;i < nX;++i)
//	//	{
//	//		double sin_theta = sin(m_model_lines_polar[i].theta * m_lambda[i]);
//	//		double cos_theta = cos(m_model_lines_polar[i].theta * m_lambda[i]);
//	//		double rho = m_model_lines_polar[i].rho * m_lambda[i];
//
//	//		for (int j = 0;j < nY;++j)
//	//		{
//	//			double sin_theta_prime = sin(m_observe_lines_polar[j].theta);
//	//			double cos_theta_prime = cos(m_observe_lines_polar[j].theta);
//	//			double rho_prime = m_observe_lines_polar[j].rho;
//
//	//			vector<double> coeff(2);
//	//			coeff[0] = cos_theta * sin_theta_prime;
//	//			coeff[1] = -cos_theta * cos_theta_prime;
//
//	//			for(int p1=0;p1<2;++p1)
//	//			{        
//	//				b[p1] += coeff[p1] * rho_prime * m_alpha(j,i);
//	//				for(int p2=p1;p2<2;++p2)
//	//				{
//	//					A(p1,p2) += coeff[p1] * coeff[p2] * m_alpha(j,i);
//	//				}
//	//			}
//	//		}
//	//	}
//	//	double tmp = det(A);
//	//	if ( fabs(tmp) > DBL_EPSILON )
//	//	{
//	//		Eigen::VectorXd tmp2 = arma::inv(A) * b;
//	//		//m_parameters.set_size(6);
//	//		setAdjustableParameter(0, tmp1[0], true);
//	//		setAdjustableParameter(1, tmp2[0], true);
//	//		setAdjustableParameter(2, tmp1[1], true);
//	//		setAdjustableParameter(3, tmp1[2], true);
//	//		setAdjustableParameter(4, tmp2[1], true);
//	//		setAdjustableParameter(5, tmp1[3], true);
//	//		//m_parameters[0] = tmp1[0];	// a0
//	//		//m_parameters[1] = tmp2[0];	// a1
//	//		//m_parameters[2] = tmp1[1];	// a2
//	//		//m_parameters[3] = tmp1[2];	// b0
//	//		//m_parameters[4] = tmp2[1];	// b1
//	//		//m_parameters[5] = tmp1[3];	// b2
//
//	//		// debug
//	//		double sum = 0.0;
//	//		for(int i = 0;i < nX;++i)
//	//		{
//	//			for (int j = 0;j < nY;++j)
//	//			{
//	//				fPoint dist = distance2(m_model_lines_polar[i], m_observe_lines_polar[j]);
//	//				sum += (dist.x*dist.x + dist.y*dist.y)*m_alpha(j,i);
//	//			}
//	//		}
//	//		cout<<"rms:\t"<<sum<<endl;
//	//		m_delta_2 = sum / (nX * ndim);
//	//		return true;
//	//	}
//	//	else
//	//	{
//	//		cout << "ERROR : cannot find a solution" << endl;
//	//		return false;
//	//	}
//
//	//}
//	//else
//	//{
//	//	cout << "ERROR : cannot find a solution" << endl;
//	//	return false;
//	//}
//	return true;
//}

void ecmlrAffine::updateDeltaSquare()
{
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	double sum = 0.0;
	for(int i = 0;i < nX;++i)
	{
		for (int j = 0;j < nY;++j)
		{
			double wgt = m_alpha(j,i);
			//fPoint dist = distance2(m_trans_lines_polar[i], m_observe_lines_polar[j]);
			fPoint dist = segment2polarline(m_trans_lines_polar[i], m_reference_lines_polar[j]);
			//fPoint dist;
			//dist.x = point2polarline(m_trans_lines_polar[i].pt1, m_observe_lines_polar[j]);
			//dist.y = point2polarline(m_trans_lines_polar[i].pt2, m_observe_lines_polar[j]);
			sum += (dist.x*dist.x + dist.y*dist.y)*wgt;
		}
	}

	m_delta_2 = sum / (m_lammda_sum * ndim);
	//if(m_delta_2 < 1.0) m_delta_2 = 1.0;
}