#pragma once

#include "mcempr.h"

class mcemprAffine : public mcempr
{
public:
	mcemprAffine(const arma::mat & modelData, const arma::mat & observationData)
		:mcempr(modelData, observationData){};
	~mcemprAffine(){};

	enum{PARAMETER_NUMBER = 6};

	virtual bool parameters_optimization();
	virtual arma::vec forward(const arma::vec&, arma::vec parameters)const;

	virtual int get_PARAMETER_NUM() const{return PARAMETER_NUMBER;};
	virtual bool initialParameters();
protected:
};