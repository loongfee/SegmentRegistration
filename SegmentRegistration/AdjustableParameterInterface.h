//*******************************************************************
//
// License:  See top level LICENSE.txt file.
//
// Author: Garrett Potts (gpotts@imagelinks.com)
//
//*************************************************************************
// $Id: rspfAdjustableParameterInterface.h 9968 2006-11-29 14:01:53Z gpotts $
#ifndef AdjustableParameterInterface_HEADER
#define AdjustableParameterInterface_HEADER
#include "lineConstant.h"
#include <vector>
#include "AdjustmentInfo.h"
#include "AdjustableParameterInfo.h"

class AdjustableParameterInterface
{
public:
   AdjustableParameterInterface();
   AdjustableParameterInterface(const AdjustableParameterInterface& rhs);
   virtual ~AdjustableParameterInterface(){}
   void newAdjustment(int numberOfParameters=0);
   void setCurrentAdjustment(int adjustmentIndex, bool notify=false);
   void eraseAdjustment(bool notify);
   void eraseAdjustment(int idx, bool notify);
   virtual void initAdjustableParameters();
   void resetAdjustableParameters(bool notify=false);
   void copyAdjustment(int idx, bool notify);
   void copyAdjustment(bool notify = false);
   
   /*!
    * Will copy the adjustment but will set the new center to the
    * applied current center plus the application of the adjustment
    *
    */
   void keepAdjustment(int idx, bool createCopy);
   void keepAdjustment(bool createCopy=true);
   

   const AdjustableParameterInterface& operator = (const AdjustableParameterInterface& rhs);
   void removeAllAdjustments();
   virtual int getNumberOfAdjustableParameters()const;
   double       getAdjustableParameter(int idx)const;
   void         setAdjustableParameter(int idx, double value,
                                       bool notify=false);
   void         setAdjustableParameter(int idx,
                                       double value,
                                       double sigma,
                                       bool notify=false);
   double       getParameterSigma(int idx)const;
   void         setParameterSigma(int idx,
                                  double value,
                                  bool notify=false);

   void           setParameterCenter(int idx,
                                     double center,
                                     bool notify = false);
   double        getParameterCenter(int idx)const;
   double        computeParameterOffset(int idx)const;
   void          setParameterOffset(int idx,
                                    double value,
                                    bool notify = false);
   
   
   bool isParameterLocked(int idx)const;

   void setParameterLockFlag(int idxParam, bool flag);
   bool getParameterLockFlag(int idx)const;

   void lockAllParametersCurrentAdjustment();
   void unlockAllParametersCurrentAdjustment();

   void lockAllParameters(int idxAdjustment);
   void unlockAllParameters(int idxAdjustment);
   
   void resizeAdjustableParameterArray(int numberOfParameters);

   void setAdjustment(const AdjustmentInfo& adj, bool notify=false);
   void setAdjustment(int idx, const AdjustmentInfo& adj, bool notify=false);
   
   void addAdjustment(const AdjustmentInfo& adj, bool notify);
   void getAdjustment(AdjustmentInfo& adj);
   void getAdjustment(int idx, AdjustmentInfo& adj);
   
   int getNumberOfAdjustments()const;
   int getCurrentAdjustmentIdx()const;

   
   void setDirtyFlag(bool flag=true);
   void setAllDirtyFlag(bool flag = true);
   bool hasDirtyAdjustments()const;
   
private:
   std::vector<AdjustmentInfo> theAdjustmentList;
   int                     theCurrentAdjustment;
   
public:
   virtual void adjustableParametersChanged();
   
};

#endif
