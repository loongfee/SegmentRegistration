//*******************************************************************
//
// License:  See top level LICENSE.txt file.
//
// Author: Garrett Potts (gpotts@imagelinks.com)
//
//*************************************************************************
// $Id: rspfAdjustableParameterInfo.h 9968 2006-11-29 14:01:53Z gpotts $
#ifndef AdjustableParameterInfo_HEADER
#define AdjustableParameterInfo_HEADER
#include <iostream>
#include <float.h>

class AdjustableParameterInfo
{
public:
   friend std::ostream& operator <<(std::ostream& out, const AdjustableParameterInfo& data);
   
   AdjustableParameterInfo()
      : theParameter(0.0),
        theSigma(0.0),
	theCenter(0.0),
      theLockFlag(false)
      {
      }
   AdjustableParameterInfo(const AdjustableParameterInfo& rhs)
      :theParameter(rhs.theParameter),
      theSigma(rhs.theSigma),
      theCenter(rhs.theCenter),
      theLockFlag(rhs.theLockFlag)
      {
      }
   double getParameter()const
      {
         return theParameter;
      }
   void setParameter(double parameter)
      {
         if(!theLockFlag)
         {
            theParameter = parameter;
         }
      }
   double getSigma()const
      {
         return theSigma;
      }
   void setSigma(double sigma)
      {
         if(!theLockFlag)
         {
            theSigma = sigma;
         }
      }

	void setCenter(double center);
	double getCenter()const;

	void setOffset(double value);
   
  /*!
   * will return theCenter + theSigma*theParameter
   */
	double computeOffset()const;

	void setLockFlag(bool flag)
	{
		theLockFlag = flag;
	}
	bool getLockFlag()const
	{
		return theLockFlag;
	}
protected:
	double        theParameter;
	double        theSigma;
	double        theCenter;
	bool          theLockFlag;
};

#endif
