#include "ecmlrPerspective.h"
#include "alglib/interpolation.h"

using namespace alglib;

cvline_polar ecmlrPerspective::forward(const cvline_polar& l, Eigen::VectorXd parameters)const
{
	cvline_polar outLine;
	fPoint line[2];
	line[0].x = (m_parameters[M11_param] * l.pt1.x + m_parameters[M12_param] * l.pt1.y + m_parameters[M13_param])/
		(m_parameters[M31_param] * l.pt1.x + m_parameters[M32_param] * l.pt1.y + 1.0);
	line[0].y =  (m_parameters[M21_param] * l.pt1.x + m_parameters[M22_param] * l.pt1.y + m_parameters[M23_param])/
		(m_parameters[M31_param] * l.pt1.x + m_parameters[M32_param] * l.pt1.y + 1.0);
	line[1].x =  (m_parameters[M11_param] * l.pt2.x + m_parameters[M12_param] * l.pt2.y + m_parameters[M13_param])/
		(m_parameters[M31_param] * l.pt2.x + m_parameters[M32_param] * l.pt2.y + 1.0);
	line[1].y = (m_parameters[M21_param] * l.pt2.x + m_parameters[M22_param] * l.pt2.y + m_parameters[M23_param])/
		(m_parameters[M31_param] * l.pt2.x + m_parameters[M32_param] * l.pt2.y + 1.0);
	outLine = line2polar(line);
	return outLine;
}

cvline_polar ecmlrPerspective::forward(const cvline_polar& l)const
{
	Eigen::VectorXd parameters = getParameters();
	return forward(l, parameters);
}

void ecmlrPerspective::initAdjustableParameters()
{
	m_parameters.resize(NUM_ADJUSTABLE_PARAMS);
	int numParams = (int)m_parameters.size();

	//m_parameters[M11_param] = 1.0;	
	//m_parameters[M12_param] = 0.0;
	//m_parameters[M13_param] = 0.0;	
	//m_parameters[M21_param] = 0.0;
	//m_parameters[M22_param] = 1.0;
	//m_parameters[M23_param] = 0.0;	
	//m_parameters[M31_param] = 0.0;
	//m_parameters[M32_param] = 0.0;

	m_parameters[M11_param] = 1.857;	
	m_parameters[M12_param] = 0.2857;
	m_parameters[M13_param] = 14.2857;	
	m_parameters[M21_param] = 0.142857;
	m_parameters[M22_param] = 1.142857;
	m_parameters[M23_param] = 28.857;	
	m_parameters[M31_param] = 0.0;
	m_parameters[M32_param] = 0.0;
}


bool ecmlrPerspective::parameters_optimization()
{
	return eigen_levmar_optimization();
	return ecmlr::parameters_optimization();
	int nX = m_nX;
	int nY = m_nY;
	int ndim = m_ndim;

	
	real_2d_array fmatrix;
	real_1d_array y;
	real_1d_array w;
	fmatrix.setlength(2*nX*nY, 8);
	y.setlength(2*nX*nY);
	w.setlength(2*nX*nY);
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

			fmatrix(c, 0) = -sin_theta * x1;
			fmatrix(c, 1) = -sin_theta * y1;
			fmatrix(c, 2) = -sin_theta;
			fmatrix(c, 3) = cos_theta * x1;
			fmatrix(c, 4) = cos_theta * y1;
			fmatrix(c, 5) = cos_theta;
			fmatrix(c, 6) = -rho * x1;
			fmatrix(c, 7) = -rho * y1;
			y(c) = rho;
			w(c++) = weight;
			
			fmatrix(c, 0) = -sin_theta * x2;
			fmatrix(c, 1) = -sin_theta * y2;
			fmatrix(c, 2) = -sin_theta;
			fmatrix(c, 3) = cos_theta * x2;
			fmatrix(c, 4) = cos_theta * y2;
			fmatrix(c, 5) = cos_theta;
			fmatrix(c, 6) = -rho * x2;
			fmatrix(c, 7) = -rho * y2;
			y(c) = rho;
			w(c++) = weight;
		}
	}

	ae_int_t info;
	real_1d_array p;
	lsfitreport rep;
	lsfitlinearw(y, w, fmatrix, info, p, rep);

	
	for(int i=0;i<get_PARAMETER_NUM();++i) m_parameters[i] = p[i];
	return true;
}

//bool ecmlrPerspective::parameters_optimization()
//{
//	//return ecmlr::parameters_optimization();
//	int nX = m_nX;
//	int nY = m_nY;
//	int ndim = m_ndim;
//
//	VectorXd b(8);
//	MatrixXd A(8, 8);
//	b.fill(0.0);
//	A.fill(0.0);
//
//	int c = 0;
//	for(int i = 0;i < (int)nX;++i)
//	{
//		double x1 = m_reference_lines_polar[i].pt1.x;
//		double y1 = m_reference_lines_polar[i].pt1.y;
//		double x2 = m_reference_lines_polar[i].pt2.x;
//		double y2 = m_reference_lines_polar[i].pt2.y;
//
//		for (int j = 0;j < nY;++j)
//		{
//			double sin_theta = sin(m_source_lines_polar[j].theta);
//			double cos_theta = cos(m_source_lines_polar[j].theta);
//			double rho = m_source_lines_polar[j].rho;
//			double len = segment_length(m_source_lines_polar[j]);
//			double weight = m_alpha(j,i);
//
//			vector<double> coeff(8);
//			coeff[0] = -sin_theta * x1;
//			coeff[1] = -sin_theta * y1;
//			coeff[2] = -sin_theta;
//			coeff[3] = cos_theta * x1;
//			coeff[4] = cos_theta * y1;
//			coeff[5] = cos_theta;
//			coeff[6] = -rho * x1;
//			coeff[7] = -rho * y1;
//
//			for(int p1=0;p1<8;++p1)
//			{        
//				b[p1] += coeff[p1] * rho * weight;
//				for(int p2=0;p2<8;++p2)
//				{
//					A(p1,p2) += coeff[p1] * coeff[p2] * weight;
//				}
//			}
//
//			coeff[0] = -sin_theta * x2;
//			coeff[1] = -sin_theta * y2;
//			coeff[2] = -sin_theta;
//			coeff[3] = cos_theta * x2;
//			coeff[4] = cos_theta * y2;
//			coeff[5] = cos_theta;
//			coeff[6] = -rho * x2;
//			coeff[7] = -rho * y2;
//			for(int p1=0;p1<8;++p1)
//			{        
//				b[p1] += coeff[p1] * rho * weight;
//				for(int p2=0;p2<8;++p2)
//				{
//					A(p1,p2) += coeff[p1] * coeff[p2] * weight;
//				}
//			}
//		}
//	}
//
//	Eigen::MatrixXd damper = MatrixXd::Identity(8, 8);
//	//A = A + damper * 1.0E-6;
//	double tmp = A.determinant();
//	if ( fabs(tmp) > DBL_EPSILON )
//	{
//		Eigen::VectorXd parameter = A.inverse() * b;
//		for(int i=0;i<get_PARAMETER_NUM();++i) m_parameters[i] = parameter[i];
//		return true;
//	}
//	else
//	{
//		cout << "ERROR : cannot find a solution" << endl;
//		return false;
//	}
//	return true;
//}

void ecmlrPerspective::updateDeltaSquare()
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


bool ecmlrPerspective::eigen_levmar_optimization()
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
	lmder_functor_perspective functor(nparam, 2*nX*nY);
	functor.pData = this;
	Eigen::LevenbergMarquardt<lmder_functor_perspective> lm(functor);
	info = lm.minimize(x);
	
	for(int i=0;i<nparam;++i) m_parameters[i] = x[i];
	return true;
}

int lmder_functor_perspective::operator()(const VectorXd &x, VectorXd &fvec) const
{
	ecmlrPerspective *pThis = (ecmlrPerspective*)pData;
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
	double M11 = pThis->m_parameters[ecmlrPerspective::M11_param];
	double M12 = pThis->m_parameters[ecmlrPerspective::M12_param];
	double M13 = pThis->m_parameters[ecmlrPerspective::M13_param];
	double M21 = pThis->m_parameters[ecmlrPerspective::M21_param];
	double M22 = pThis->m_parameters[ecmlrPerspective::M22_param];
	double M23 = pThis->m_parameters[ecmlrPerspective::M23_param];
	double M31 = pThis->m_parameters[ecmlrPerspective::M31_param];
	double M32 = pThis->m_parameters[ecmlrPerspective::M32_param];
	for (int i = 0;i < nX;++i)
	{
		double x1 = pThis->m_source_lines_polar[i].pt1.x;
		double y1 = pThis->m_source_lines_polar[i].pt1.y;
		double x2 = pThis->m_source_lines_polar[i].pt2.x;
		double y2 = pThis->m_source_lines_polar[i].pt2.y;

		double a1 = M11*x1+M12*y1+M13;
		double a2 = M21*x1+M22*y1+M23;
		double a3 = M31*x1+M32*y1+1.0;
		double b1 = M11*x2+M12*y2+M13;
		double b2 = M21*x2+M22*y2+M23;
		double b3 = M31*x2+M32*y2+1.0;
		for (int j = 0;j < nY;++j)
		{
			double sin_theta = sin(pThis->m_reference_lines_polar[j].theta);
			double cos_theta = cos(pThis->m_reference_lines_polar[j].theta);
			double rho = pThis->m_reference_lines_polar[j].rho;

			fvec[c++] = (-sin_theta*a1/a3+cos_theta*a2/a3-rho) * pThis->m_alpha(j,i);
			fvec[c++] = (-sin_theta*b1/b3+cos_theta*b2/b3-rho) * pThis->m_alpha(j,i);
		}
	}
	return 0;
}

int lmder_functor_perspective::df(const VectorXd &x, MatrixXd &fjac) const
{
	ecmlrPerspective *pThis = (ecmlrPerspective*)pData;
	int nX = (int)pThis->m_nX;
	int nY = (int)pThis->m_nY;
	int ndim = (int)pThis->m_ndim;

	int nparameter = pThis->get_PARAMETER_NUM();

	for(int i=0;i<nparameter;++i)
	{
		pThis->m_parameters[i] = x[i];
	}

	int c = 0;

	double M11 = pThis->m_parameters[ecmlrPerspective::M11_param];
	double M12 = pThis->m_parameters[ecmlrPerspective::M12_param];
	double M13 = pThis->m_parameters[ecmlrPerspective::M13_param];
	double M21 = pThis->m_parameters[ecmlrPerspective::M21_param];
	double M22 = pThis->m_parameters[ecmlrPerspective::M22_param];
	double M23 = pThis->m_parameters[ecmlrPerspective::M23_param];
	double M31 = pThis->m_parameters[ecmlrPerspective::M31_param];
	double M32 = pThis->m_parameters[ecmlrPerspective::M32_param];
	for (int i = 0;i < nX;++i)
	{
		double x1 = pThis->m_source_lines_polar[i].pt1.x;
		double y1 = pThis->m_source_lines_polar[i].pt1.y;
		double x2 = pThis->m_source_lines_polar[i].pt2.x;
		double y2 = pThis->m_source_lines_polar[i].pt2.y;

		double a1 = M11*x1+M12*y1+M13;
		double a2 = M21*x1+M22*y1+M23;
		double a3 = M31*x1+M32*y1+1.0;

		double b1 = M11*x2+M12*y2+M13;
		double b2 = M21*x2+M22*y2+M23;
		double b3 = M31*x2+M32*y2+1.0;
		for (int j = 0;j < nY;++j)
		{
			double sin_theta = sin(pThis->m_reference_lines_polar[j].theta);
			double cos_theta = cos(pThis->m_reference_lines_polar[j].theta);
			double rho = pThis->m_reference_lines_polar[j].rho;

			// M11
			fjac(c,0) = -sin_theta*x1/a3*pThis->m_alpha(j,i);
			fjac(c+1,0) = -sin_theta*x2/b3*pThis->m_alpha(j,i);

			// M12
			fjac(c,1) = -sin_theta*y1/a3*pThis->m_alpha(j,i);
			fjac(c+1,1) = -sin_theta*y2/b3*pThis->m_alpha(j,i);

			// M13
			fjac(c,2) = -sin_theta/a3*pThis->m_alpha(j,i);
			fjac(c+1,2) = -sin_theta/b3*pThis->m_alpha(j,i);

			// M21
			fjac(c,3) = cos_theta*x1/a3*pThis->m_alpha(j,i);
			fjac(c+1,3) = cos_theta*x2/b3*pThis->m_alpha(j,i);

			// M22
			fjac(c,4) = cos_theta*y1/a3*pThis->m_alpha(j,i);
			fjac(c+1,4) = cos_theta*y2/b3*pThis->m_alpha(j,i);

			// M23
			fjac(c,5) = cos_theta/a3;
			fjac(c+1,5) = cos_theta/b3;

			// M31
			fjac(c,6) = (sin_theta*a1-cos_theta*a2)*x1/(a3*a3)*pThis->m_alpha(j,i);
			fjac(c+1,6) = (sin_theta*b1-cos_theta*b2)*x2/(b3*b3)*pThis->m_alpha(j,i);

			// M32
			fjac(c,7) = (sin_theta*a1-cos_theta*a2)*y1/(a3*a3)*pThis->m_alpha(j,i);
			fjac(c+1,7) = (sin_theta*b1-cos_theta*b2)*y2/(b3*b3)*pThis->m_alpha(j,i);

			c += 2;
		}
	}
	return 0;
}
