#pragma once

#include "ecm.h"

class ecmprAffine : public ecmpr
{
public:
	ecmprAffine(const arma::mat & modelData, const arma::mat & observationData)
		:ecmpr(modelData, observationData){};
	~ecmprAffine(){};

	enum{PARAMETER_NUMBER = 6};

	virtual bool parameters_optimization();
	virtual arma::vec forward(const arma::vec&, arma::vec parameters)const;

	virtual int get_PARAMETER_NUM() const{return PARAMETER_NUMBER;};
	virtual bool initialParameters();
	virtual void updateDelta_2(bool bfirstTime);
protected:
};