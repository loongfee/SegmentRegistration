#include "LineMatchCost.h"

bool LineMatchCost::Evaluate(double const* const* parameters,
	double* residuals,
	double** jacobians) const
 {
	const double x = parameters[0][0];
	residuals[0] = 10 - x;
	// f(x) = 10 − x
	// Compute the Jacobian if asked for.
	if (jacobians != NULL && jacobians[0] != NULL) {
	jacobians[0][0] = -1;
	}
	return true;

	int nX = (int)m_pEcmlr->m_nX;
	int nY = (int)m_pEcmlr->m_nY;
	int ndim = (int)m_pEcmlr->m_ndim;

	int nparameter = m_pEcmlr->get_PARAMETER_NUM();

	for(int i=0;i<nparameter;++i)
	{
		m_pEcmlr->setAdjustableParameter(i, parameters[0][i], true);
	}

	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;
	int c = 0;


	for (int i = 0;i < nX;++i)
	{
		cvline_polar outPt = m_pEcmlr->forward(m_pEcmlr->m_model_lines_polar[i]);
		for (int j = 0;j < nY;++j)
		{
			//fPoint dist = pThis->distance2(pThis->m_observe_lines_polar[j], outPt);
			fPoint dist = m_pEcmlr->segment2polarline(outPt, m_pEcmlr->m_observe_lines_polar[j]);

			residuals[c] = m_pEcmlr->m_alpha(j,i) * dist.x;
			residuals[c+1] = m_pEcmlr->m_alpha(j,i) * dist.y;

			for(int p=0;p<nparameter;++p)
			{
				double middle = m_pEcmlr->getAdjustableParameter(p);
				m_pEcmlr->setAdjustableParameter(p, middle + pstep_scale, true);
				cvline_polar outLine1 = m_pEcmlr->forward(m_pEcmlr->m_model_lines_polar[i]);
				//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
				fPoint dist1 = m_pEcmlr->segment2polarline(outLine1, m_pEcmlr->m_observe_lines_polar[j]);

				m_pEcmlr->setAdjustableParameter(p, middle - pstep_scale, true);
				cvline_polar outLine2 = m_pEcmlr->forward(m_pEcmlr->m_model_lines_polar[i]);
				//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
				fPoint dist2 = m_pEcmlr->segment2polarline(outLine2, m_pEcmlr->m_observe_lines_polar[j]);

				m_pEcmlr->setAdjustableParameter(p, middle, true);

				double derivative_x = (dist1.x - dist2.x) * den;
				double derivative_y = (dist1.y - dist2.y) * den;

				jacobians[c][p] = m_pEcmlr->m_alpha(j,i) * derivative_x;
				jacobians[c+1][p] = m_pEcmlr->m_alpha(j,i) * derivative_y;
			}
			c += 2;
		}
	}
}