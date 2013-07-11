#include "mcemprAffine.h"

arma::vec mcemprAffine::forward(const arma::vec& original, arma::vec parameters)const
{
	arma::vec outPt(2);
	outPt(0) = parameters[0] + parameters[1] * original[0] + parameters[2] * original[1];
	outPt(1) = parameters[3] + parameters[4] * original[0] + parameters[5] * original[1];
	return outPt;
}

bool mcemprAffine::initialParameters()
{
	m_parameters.set_size(get_PARAMETER_NUM());
	m_parameters[0] = 0.0;
	m_parameters[1] = 1.0;
	m_parameters[2] = 0.0;
	m_parameters[3] = 0.0;
	m_parameters[4] = 0.0;
	m_parameters[5] = 1.0;
	return true;
}

bool mcemprAffine::parameters_optimization()
{
	int nX = (int)m_modelData.n_rows;
	int nY = (int)m_observationData.n_rows;
	int ndim = (int)m_modelData.n_cols;

	//arma::mat weightMatrix(nX * 2, nX * 2);
	//weightMatrix.fill(0.0);
	//weightMatrix.diag() = m_lambda;


	//vec b(nX * 2);
	//mat A(nX * 2, 6);
	////b = candidates.col(6);
	//for(int i = 0;i < (int)nX;++i)
	//{
	//	A(2*i, 0) = 1.0;
	//	A(2*i, 1) = m_modelData(i, 0);
	//	A(2*i, 2) = m_modelData(i, 1);
	//	A(2*i, 3) = 0.0;
	//	A(2*i, 4) = 0.0;
	//	A(2*i, 5) = 0.0;
	//	b(2*i) = m_w(i, 0);

	//	A(2*i+1, 0) = 0.0;
	//	A(2*i+1, 1) = 0.0;
	//	A(2*i+1, 2) = 0.0;
	//	A(2*i+1, 3) = 1.0;
	//	A(2*i+1, 4) = m_modelData(i, 0);
	//	A(2*i+1, 5) = m_modelData(i, 1);
	//	b(2*i+1) = m_w(i, 1);

	//	weightMatrix(2*i, 2*i) = m_lambda(i);
	//	weightMatrix(2*i+1, 2*i+1) = m_lambda(i);
	//}
	vec b(6);
	mat A(6, 6);
	b.fill(0.0);
	A.fill(0.0);
	//b = candidates.col(6);
	for(int i = 0;i < (int)nX;++i)
	{
		double x = m_modelData(i, 0);
		double y = m_modelData(i, 1);
		double xx = m_w(i, 0);
		double yy = m_w(i, 1);
		A(0,0) += m_lambda(i);
		A(0,1) += m_lambda(i) * x;
		A(0,2) += m_lambda(i) * y;
		A(1,0) += m_lambda(i) * x;
		A(1,1) += m_lambda(i) * x * x;
		A(1,2) += m_lambda(i) * y * x;
		A(2,0) += m_lambda(i) * y;
		A(2,1) += m_lambda(i) * x * y;
		A(2,2) += m_lambda(i) * y * y;

		A(3,3) += m_lambda(i);
		A(3,4) += m_lambda(i) * x;
		A(3,5) += m_lambda(i) * y;
		A(4,3) += m_lambda(i) * x;
		A(4,4) += m_lambda(i) * x * x;
		A(4,5) += m_lambda(i) * y * x;
		A(5,3) += m_lambda(i) * y;
		A(5,4) += m_lambda(i) * x * y;
		A(5,5) += m_lambda(i) * y * y;

		//b(0) += m_lambda(i) * xx;
		//b(1) += m_lambda(i) * x * xx;
		//b(2) += m_lambda(i) * y * xx;
		//b(3) += m_lambda(i) * yy;
		//b(4) += m_lambda(i) * x * yy;
		//b(5) += m_lambda(i) * y * yy;
		b(0) += xx;
		b(1) += x * xx;
		b(2) += y * xx;
		b(3) += yy;
		b(4) += x * yy;
		b(5) += y * yy;
	}

	//arma::mat AA = A.t() * weightMatrix * A;
	//arma::mat bb = A.t() * weightMatrix * b;
	double tmp = det(A);
	if ( fabs(tmp) > DBL_EPSILON )
	//if ( ::solve(m_parameters, A, b) )
	{
		m_parameters = arma::inv(A) * b;
		// debug
		double sum = 0.0;
		for (int i = 0;i < nX;++i)
		{
			arma::vec outPt = forward(m_modelData.row(i).t(), m_parameters);
			double a1 = m_lambda[i]*(outPt(0)*outPt(0) + outPt(1) * outPt(1));
			double a2 = 2.0 * (outPt(0)*m_w(i,0) + outPt(1)*m_w(i,1));
			sum += a1 - a2;
			//sum += m_lambda[i]*norm2_Vec(m_w.row(i).t() - outPt);
		}
		for(int j = 0;j < nY;++j)
		{
			double a3 = (1.0 - m_alpha(j,nX)) * (m_observationData(j,0)*m_observationData(j,0) + m_observationData(j,1)*m_observationData(j,1));
			sum += a3;
		}
		//arma::vec delta = b - A * m_parameters;
		m_delta_2 = sum / (nY * ndim);
		m_delta_2 = m_delta_2 < 1.0 ? 1.0 : m_delta_2;
		cout<<"rms:\t"<<m_delta_2<<endl;
		//m_delta_2 = 0.90;
		return true;
	}
	else
	{
		cout << "ERROR : cannot find a solution" << endl;
		return false;
	}

	return true;
}