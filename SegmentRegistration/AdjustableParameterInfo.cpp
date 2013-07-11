//*******************************************************************
// Copyright (C) 2000 ImageLinks Inc. 
//
// License:  See top level LICENSE.txt file.
//
// Author: Garrett Potts
//
//*************************************************************************
// $Id: rspfAdjustableParameterInfo.cpp 11347 2007-07-23 13:01:59Z gpotts $
#include <sstream>
#include <algorithm>
#include "AdjustableParameterInfo.h"

// static const char* PARAM_NAME_KW       = "name";
// static const char* PARAM_UNITS_KW      = "units";
static const char* PARAM_KW            = "parameter";
static const char* PARAM_SIGMA_KW      = "sigma";
static const char* PARAM_CENTER_KW     = "center";
static const char* PARAM_LOCK_FLAG_KW  = "lock_flag";

std::ostream& operator <<(std::ostream& out, const AdjustableParameterInfo& data)
{
   out << "center:      " << data.theCenter <<  std::endl
       << "parameter:   " << data.theParameter << std::endl
       << "sigma:       " << data.theSigma << std::endl
       << std::endl
       << "locked:       " << (data.theLockFlag?"true":"false") << std::endl;
   
   return out;
}

void AdjustableParameterInfo::setCenter(double center)
{
   if(!theLockFlag)
   {
      theCenter = center;
   }
}

double AdjustableParameterInfo::getCenter()const
{
  return theCenter;
}

double AdjustableParameterInfo::computeOffset()const
{
  return theCenter + theSigma*theParameter;
}

void AdjustableParameterInfo::setOffset(double value)
{
   if(!theLockFlag)
   {
      double minValue = theCenter - theSigma;
      double maxValue = theCenter + theSigma;
      double x = 0.0;
      
      if(std::abs(theSigma) > DBL_EPSILON)
      {
         x = (value - theCenter)/theSigma;
         
         value = theCenter + x*theSigma;
         
         if(value < minValue)
         {
            x = -1;
         }
         else if(value > maxValue)
         {
            x = 1.0;
         }
         theParameter = x;
      }
   }
}