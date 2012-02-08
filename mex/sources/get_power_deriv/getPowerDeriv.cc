/*
** mex file for matlab getPowerDeriv.m.
**
** Copyright (C) 2012 Dynare Team
**
** This file is part of Dynare.
**
** Dynare is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** Dynare is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
**
** AUTHOR(S): stephane DOT adjemian AT univ DASH lemans DOT fr
**/


#include <math.h>
#include <dynmex.h>
#include "mex.h"

#ifdef _MSC_VER
# define nearbyint(x) (fabs((x)-floor(x)) < fabs((x)-ceil(x)) ? floor(x) : ceil(x))
#endif

double getPowerDeriv(double *x, double *p, int k)
{
  if ( fabs(*x) < 1e-12 && *p > 0 && k >= *p && fabs(*p-nearbyint(*p)) < 1e-12 )
    return 0.0;
  else
    {
      //int i = 0;
      double dxp = pow(*x, *p-k);
      for (int i=0; i<k; i++)
        dxp *= (*p)--;
      return dxp;
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*
  ** INPUTS:
  ** prhs[0]    [double]    scalar, x.
  ** prhs[1]    [double]    scalar, p.
  ** prhs[2]    [integer]   scalar, k.
  **
  ** OUTPUTS:
  ** plhs[0]    [double]    scalar, The k-th derivative of x^p.
  ** plhs[1]    [double]    scalar, info variable (equal to zero if call to getPowerDeriv is successfull).
  */
  if ( !( nrhs==3)  )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("getPowerDeriv:: Three input arguments are required (x,p and k)!");
    }
  if (nlhs>2)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("getPowerDeriv:: I only return one or two output arguments!");
    }
  if ( !(mxGetN(prhs[0])==1 && (mxGetM(prhs[0]))==1) )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("getPowerDeriv:: First argument, x, must be a double scalar!");
    }
  if ( !(mxGetN(prhs[1])==1 && (mxGetM(prhs[1]))==1) )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("getPowerDeriv:: Second argument, p, must be a double scalar!");
    }
  if ( !(mxGetN(prhs[2])==1 && (mxGetM(prhs[2]))==1) )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("getPowerDeriv:: Third argument, k, must be an integer scalar!");
    }
  double *x = mxGetPr(prhs[0]);
  double *p = mxGetPr(prhs[1]);
  int k = (int) mxGetScalar(prhs[2]);
  plhs[0] = mxCreateDoubleScalar(getPowerDeriv(x,p,k));
  if (nlhs>1)
    {
      plhs[1] = mxCreateDoubleScalar(0.0);
    }
}
