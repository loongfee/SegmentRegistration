//*******************************************************************
// Copyright (C) 2000 ImageLinks Inc. 
//
// License:  See top level LICENSE.txt file.
//
// Author: Garrett Potts
//
//*************************************************************************
// $Id: rspfAdjustableParameterInterface.cpp 20610 2012-02-27 12:19:25Z gpotts $
#include <algorithm>
#include "AdjustableParameterInterface.h"

static const char* NUMBER_OF_ADJUSTMENTS_KW = "number_of_adjustments";
static const char* NUMBER_OF_ADJUSTMENTS_OLD_KW = "number_of_adjustements";
static const char* CURRENT_ADJUSTMENT_OLD_KW    = "current_adjustement";
static const char* CURRENT_ADJUSTMENT_KW    = "current_adjustment";
static const char* ADJUSTMENT_PREFIX        = "adjustment_";

AdjustableParameterInterface::AdjustableParameterInterface()
{
   theCurrentAdjustment = 0;
}

AdjustableParameterInterface::AdjustableParameterInterface(const AdjustableParameterInterface& rhs)
   :theAdjustmentList(rhs.theAdjustmentList),
    theCurrentAdjustment(rhs.theCurrentAdjustment)
{
}

void AdjustableParameterInterface::newAdjustment(int numberOfParameters)
{
   theAdjustmentList.push_back(AdjustmentInfo());
   if(numberOfParameters > 0)
   {
      theAdjustmentList[theAdjustmentList.size()-1].setNumberOfAdjustableParameters(numberOfParameters);
   }

   theCurrentAdjustment = (int)theAdjustmentList.size() - 1;

}

void AdjustableParameterInterface::setCurrentAdjustment(int adjustmentIdx, bool notify)
{
   if(adjustmentIdx < theAdjustmentList.size())
   {
      theCurrentAdjustment = adjustmentIdx;
      if(notify)
      {
         adjustableParametersChanged();
      }
   }
}


void AdjustableParameterInterface::initAdjustableParameters()
{
}

void AdjustableParameterInterface::resetAdjustableParameters(bool notify)
{
    if(!theAdjustmentList.size())
    {
       return;
    }
    
    int saveCurrent = theCurrentAdjustment;
    copyAdjustment();
    initAdjustableParameters();
    int numberOfAdjustables = getNumberOfAdjustableParameters();
    int idx = 0;
    
    for(idx = 0; idx < numberOfAdjustables; ++idx)
    {
       theAdjustmentList[saveCurrent].getParameterList()[idx].setParameter(theAdjustmentList[theAdjustmentList.size()-1].getParameterList()[idx].getParameter());
    }

    setCurrentAdjustment(saveCurrent);

    eraseAdjustment((int)theAdjustmentList.size()-1, false);
    
    if(notify)
    {
       adjustableParametersChanged();
    }
}

void AdjustableParameterInterface::copyAdjustment(int idx, bool notify)
{
    if(!theAdjustmentList.size())
    {
       return;
    }
    if(idx < theAdjustmentList.size())
    {
       theAdjustmentList.push_back(theAdjustmentList[idx]);

       if(idx == theCurrentAdjustment)
       {
          theCurrentAdjustment = (int)theAdjustmentList.size() - 1;
       }
       if(notify)
       {
          adjustableParametersChanged();
       }
    }
    
}

void AdjustableParameterInterface::copyAdjustment(bool notify)
{
   copyAdjustment(theCurrentAdjustment, notify);
}

void AdjustableParameterInterface::keepAdjustment(int idx,
                                                       bool createCopy)
{
    if(!theAdjustmentList.size())
    {
       return;
    }
    if(idx < theAdjustmentList.size())
    {
       if(createCopy)
       {
          copyAdjustment(idx);
       }
       theAdjustmentList[theCurrentAdjustment].keep();
    }
}

void AdjustableParameterInterface::keepAdjustment(bool createCopy)
{
   keepAdjustment(theCurrentAdjustment, createCopy);
}

const AdjustableParameterInterface& AdjustableParameterInterface::operator = (const AdjustableParameterInterface& rhs)
{
   theAdjustmentList    = rhs.theAdjustmentList;
   theCurrentAdjustment = rhs.theCurrentAdjustment;

   return *this;
}

void AdjustableParameterInterface::removeAllAdjustments()
{
   theAdjustmentList.clear();
   theCurrentAdjustment = 0;
}

int AdjustableParameterInterface::getNumberOfAdjustableParameters()const
{
   if(theAdjustmentList.size())
   {
      return theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters();
   }

   return 0;
}

void AdjustableParameterInterface::eraseAdjustment(bool notify)
{
   eraseAdjustment(theCurrentAdjustment, notify);
}

void AdjustableParameterInterface::eraseAdjustment(int idx, bool notify)
{
   if(!theAdjustmentList.size())
   {
      return;
   }
   
   if(theCurrentAdjustment == idx)
   {
      theAdjustmentList.erase(theAdjustmentList.begin() + theCurrentAdjustment);
      if(theCurrentAdjustment >= theAdjustmentList.size())
      {
         if(theAdjustmentList.size() < 1)
         {
            theCurrentAdjustment = 0;
         }
         else
         {
            theCurrentAdjustment = (int)theAdjustmentList.size() - 1;
         }
         
      }
      
      if(notify)
      {
         adjustableParametersChanged();
      }
   }
   else if(idx < theAdjustmentList.size())
   {
      theAdjustmentList.erase(theAdjustmentList.begin() + idx);
      if(theAdjustmentList.size() < 1)
      {
         theCurrentAdjustment = 0;
      }
      else
      {
         if(theCurrentAdjustment > idx)
         {
            --theCurrentAdjustment;
            if(notify)
            {
               adjustableParametersChanged();
            }
         }
      }
      if(notify)
      {
         adjustableParametersChanged();
      }
   }
}

double AdjustableParameterInterface::getAdjustableParameter(int idx)const
{
   if(theAdjustmentList.size())
   {
      if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
      {
         return theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].getParameter();
      }
   }
   
   return 0.0;
}

void AdjustableParameterInterface::setAdjustableParameter(int idx, double value, double sigma, bool notify)
{
   if(!theAdjustmentList.size())
   {
      return;
   }
   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].setParameter(value);
      theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].setSigma(sigma);
      if(notify)
      {
         adjustableParametersChanged();
      }
   }
   
}

void AdjustableParameterInterface::setAdjustableParameter(int idx, double value, bool notify)
{
   if(!theAdjustmentList.size())
   {
      return;
   }
   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].setParameter(value);

      if(notify)
      {
         adjustableParametersChanged();
      }
   }
}

double AdjustableParameterInterface::getParameterSigma(int idx)const
{
   if(theAdjustmentList.size())
   {
      if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
      {
         return theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].getSigma();
      }
   }

   return 0.0;
}

void AdjustableParameterInterface::setParameterSigma(int idx, double value, bool notify)
{
   if(!theAdjustmentList.size())
   {
      return;
   }
   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].setSigma(value);
      if(notify)
      {
         adjustableParametersChanged();
      }
   }
}

void  AdjustableParameterInterface::setParameterCenter(int idx, double center, bool notify)
{
   if(!theAdjustmentList.size())
   {
      return;
   }

   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
     theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].setCenter(center);

	 if(notify)
	 {
	   adjustableParametersChanged();
	 }
   }
}

double AdjustableParameterInterface::getParameterCenter(int idx)const
{
   if(!theAdjustmentList.size())
   {
      return 0.0;
   }
   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      return theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].getCenter();
   }

   return 0.0;
}

double   AdjustableParameterInterface::computeParameterOffset(int idx)const
{
   if(!theAdjustmentList.size())
   {
      return 0.0;
   }
   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      return theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].computeOffset();
   }

   return 0.0;
}


bool AdjustableParameterInterface::isParameterLocked(int idx)const
{
   if(!theAdjustmentList.size())
   {
      return false;
   }
   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      return theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].getLockFlag();
   }

   return false;
   
}

void AdjustableParameterInterface::setParameterLockFlag(int idxParam, bool flag)
{
   if(!theAdjustmentList.size())
   {
      return;
   }
   if(idxParam < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      theAdjustmentList[theCurrentAdjustment].getParameterList()[idxParam].setLockFlag(flag);
   }
}

bool AdjustableParameterInterface::getParameterLockFlag(int idx)const
{
   if(!theAdjustmentList.size())
   {
      return false;
   }
   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      return theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].getLockFlag();
   }

   return false;
}

void AdjustableParameterInterface::lockAllParametersCurrentAdjustment()
{
   lockAllParameters(theCurrentAdjustment);
}

void AdjustableParameterInterface::unlockAllParametersCurrentAdjustment()
{
   unlockAllParameters(theCurrentAdjustment);
}

void AdjustableParameterInterface::lockAllParameters(int idxAdjustment)
{
   if(idxAdjustment < getNumberOfAdjustments())
   {
      int idx = 0;
      int n   = theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters();
      
      for(idx = 0; idx < n; ++idx)
      {
         theAdjustmentList[idxAdjustment].getParameterList()[idx].setLockFlag(true);
      }
   }
}

void AdjustableParameterInterface::unlockAllParameters(int idxAdjustment)
{
   if(idxAdjustment < getNumberOfAdjustments())
   {
      int idx = 0;
      int n   = theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters();
      
      for(idx = 0; idx < n; ++idx)
      {
         theAdjustmentList[idxAdjustment].getParameterList()[idx].setLockFlag(false);
      }
   }
}


void AdjustableParameterInterface::setParameterOffset(int idx,
                                                           double value,
                                                           bool notify)
{
//    double center   = getParameterCenter(idx);
//    double sigma    = getParameterSigma(idx);
//    double minValue = center - sigma;
//    double maxValue = center + sigma;
//    double x = 0.0;
   
//    if(sigma != 0.0)
//    {
//       x = (value - center)/sigma;
      
//       value = center + x*sigma;
      
//       if(value < minValue)
//       {
//          x = -1;
//       }
//       else if(value >maxValue)
//       {
//          x = 1.0;
//       }
//       setAdjustableParameter(idx, x, false);
//    }
   
   if(!theAdjustmentList.size())
   {
      return;
   }
   if(idx < theAdjustmentList[theCurrentAdjustment].getNumberOfAdjustableParameters())
   {
      theAdjustmentList[theCurrentAdjustment].getParameterList()[idx].setOffset(value);
      if(notify)
      {
         adjustableParametersChanged();
      }
   }
}

void AdjustableParameterInterface::resizeAdjustableParameterArray(int numberOfParameters)
{
   if(!theAdjustmentList.size())
   {
      newAdjustment(numberOfParameters);
      return;
   }

   theAdjustmentList[theCurrentAdjustment].setNumberOfAdjustableParameters(numberOfParameters);
}

void AdjustableParameterInterface::setAdjustment(const AdjustmentInfo& adj, bool notify)
{
   setAdjustment(theCurrentAdjustment, adj, notify);
}

void AdjustableParameterInterface::setAdjustment(int idx, const AdjustmentInfo& adj, bool notify)
{
   if(idx < getNumberOfAdjustments())
   {
      theAdjustmentList[(int)idx] = adj;
      if(notify)
      {
         adjustableParametersChanged();
      }
   }
}


void AdjustableParameterInterface::addAdjustment(const AdjustmentInfo& adj, bool notify)
{
   theAdjustmentList.push_back(adj);
   if(notify)
   {
      adjustableParametersChanged();
   }
}

void AdjustableParameterInterface::getAdjustment(AdjustmentInfo& adj)
{
   getAdjustment(theCurrentAdjustment,  adj);
}

void AdjustableParameterInterface::getAdjustment(int idx, AdjustmentInfo& adj)
{
   adj.setNumberOfAdjustableParameters(0);

   if(idx < getNumberOfAdjustments())
   {
      adj = theAdjustmentList[(int)idx];
   }
}

int AdjustableParameterInterface::getNumberOfAdjustments()const
{
   return (int)theAdjustmentList.size();
}

int AdjustableParameterInterface::getCurrentAdjustmentIdx()const
{
   return theCurrentAdjustment;
}

void AdjustableParameterInterface::setDirtyFlag(bool flag)
{
   if(theAdjustmentList.size() > 0)
   {
      theAdjustmentList[theCurrentAdjustment].setDirtyFlag(flag);
   }
}

void AdjustableParameterInterface::setAllDirtyFlag(bool flag)
{
   int idx = 0;
   
   for(idx = 0; idx < theAdjustmentList.size(); ++idx)
   {
      theAdjustmentList[idx].setDirtyFlag(flag);
   }
}

bool AdjustableParameterInterface::hasDirtyAdjustments()const
{
   int idx = 0;
      
   for(idx = 0; idx < theAdjustmentList.size(); ++idx)
   {
      if(theAdjustmentList[idx].isDirty())
      {
         return true;
      }
   }

   return false;
}
void AdjustableParameterInterface::adjustableParametersChanged()
{
}
