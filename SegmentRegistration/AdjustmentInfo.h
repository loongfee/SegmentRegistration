//*******************************************************************
//
// License:  See top level LICENSE.txt file.
//
// Author: Garrett Potts (gpotts@imagelinks.com)
//
//*************************************************************************
// $Id: rspfAdjustmentInfo.h 9968 2006-11-29 14:01:53Z gpotts $
#ifndef AdjustmentInfo_HEADER
#define AdjustmentInfo_HEADER
#include <vector>
#include "AdjustableParameterInfo.h"

class AdjustmentInfo
{
public:
   friend std::ostream& operator <<(std::ostream& out, const AdjustmentInfo& data);
   
   
   AdjustmentInfo(int numberOfAdjustableParameters=0);
   AdjustmentInfo(const AdjustmentInfo& rhs);
   
   void setNumberOfAdjustableParameters(int numberOfAdjustableParameters);
   int getNumberOfAdjustableParameters()const;
   bool isDirty()const;
   void setDirtyFlag(bool flag=true);
   void setLockParameterFlag(bool flag,
                             int idx);
   void keep();
   
   std::vector<AdjustableParameterInfo>& getParameterList();
   const std::vector<AdjustableParameterInfo>& getParameterList()const;
   void setLockFlag(bool flag,int idx);
   
private:
   std::vector<AdjustableParameterInfo> theParameterList;
   mutable bool                              theDirtyFlag;
};

#endif
