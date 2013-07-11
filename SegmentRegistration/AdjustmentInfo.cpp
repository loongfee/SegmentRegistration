//*******************************************************************
// Copyright (C) 2000 ImageLinks Inc.
//
// LICENSE: See top level LICENSE.txt file.
//
// Author: Garrett Potts
//
//*************************************************************************
// $Id: rspfAdjustmentInfo.cpp 15833 2009-10-29 01:41:53Z eshirschorn $
#include "AdjustmentInfo.h"

static const char* PARAM_PREFIX             = "adj_param_";
static const char* NUMBER_OF_PARAMS_KW      = "number_of_params";
static const char* DIRTY_FLAG_KW            = "dirty_flag";


std::ostream& operator <<(std::ostream& out, const AdjustmentInfo& data)
{
   int idx = 0;

   out << "\nNumber of Params: " << data.theParameterList.size()
       << "\nDirty flag:       " << data.theDirtyFlag << std::endl;

   for(idx = 0; idx < (int)data.getNumberOfAdjustableParameters(); ++idx)
   {
      out << "Param " << idx << std::endl;
      out << data.theParameterList[idx] << std::endl;
   }

   return out;
}


AdjustmentInfo::AdjustmentInfo(int numberOfAdjustableParameters)
   :theParameterList(numberOfAdjustableParameters),
    theDirtyFlag(false)
{
}

AdjustmentInfo::AdjustmentInfo(const AdjustmentInfo& rhs)
   :theParameterList(rhs.theParameterList),
    theDirtyFlag(rhs.theDirtyFlag)
{
}

void AdjustmentInfo::setNumberOfAdjustableParameters(int numberOfAdjustableParameters)
{
   std::vector<AdjustableParameterInfo> temp = theParameterList;

   theParameterList.resize(numberOfAdjustableParameters);
   if(temp.size() < numberOfAdjustableParameters)
   {
      std::copy(temp.begin(),
                temp.end(),
                theParameterList.begin());
   }
   else if(temp.size() > numberOfAdjustableParameters)
   {
      if(numberOfAdjustableParameters > 0)
      {
         std::copy(temp.begin(),
                   temp.begin()+numberOfAdjustableParameters,
                   theParameterList.begin());
      }
   }
}

int AdjustmentInfo::getNumberOfAdjustableParameters()const
{
   return (int)theParameterList.size();
}

bool AdjustmentInfo::isDirty()const
{
   return theDirtyFlag;
}

void AdjustmentInfo::setDirtyFlag(bool flag)
{
   theDirtyFlag = flag;
}

void AdjustmentInfo::setLockFlag(bool flag,
                                      int idx)
{
   if(idx < theParameterList.size())
   {
      theParameterList[idx].setLockFlag(flag);
   }
}

void AdjustmentInfo::keep()
{
   int idx = 0;

   for(idx = 0; idx < theParameterList.size();++idx)
   {
      double center = theParameterList[idx].computeOffset();
      theParameterList[idx].setParameter(0.0);
      theParameterList[idx].setCenter(center);
   }
}


std::vector<AdjustableParameterInfo>& AdjustmentInfo::getParameterList()
{
   return theParameterList;
}

const std::vector<AdjustableParameterInfo>& AdjustmentInfo::getParameterList()const
{
   return theParameterList;
}