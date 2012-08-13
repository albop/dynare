/*
 * Copyright (C) 2010-2012 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#if defined(_WIN32) || defined(__CYGWIN32__)
# ifndef NOMINMAX
#  define NOMINMAX // Do not define "min" and "max" macros
# endif
# include <windows.h>
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

#include <string>
#include <Eigen/Core>

using namespace Eigen;

#include "ts_exception.h"

// <model>_Dynamic DLL pointer
typedef void (*DynamicFn)
(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state,
 int it_, double *residual, double *g1, double *g2, double *g3);

/**
 * creates pointer to Dynamic function inside <model>_dynamic.dll
 * and handles calls to it.
 **/
class DynamicModelDLL
{
private:
  DynamicFn Dynamic; // pointer to the Dynamic function in DLL
#if defined(_WIN32) || defined(__CYGWIN32__)
  HINSTANCE dynamicHinstance;  // DLL instance pointer in Windows
#else
  void *dynamicHinstance; // and in Linux or Mac
#endif

public:
  // construct and load Dynamic model DLL
  DynamicModelDLL(const std::string &dynamicDllFile) throw (TSException);
  virtual ~DynamicModelDLL();

  //! evaluate Dynamic model DLL
  template<class Derived1, class Derived2>
  void eval(const VectorXd &y, const MatrixXd &x, const PlainObjectBase<Derived1> &modParams,
            PlainObjectBase<Derived2> &ySteady, VectorXd &residual, MatrixXd *g1, MatrixXd *g2,
            MatrixXd *g3) throw (TSException)
  {
    Dynamic(y.data(), x.data(), 1, modParams.data(), ySteady.data(), 0, residual.data(),
	    g1 == NULL ? NULL : g1->data(), g2 == NULL ? NULL : g2->data(), g3 == NULL ? NULL : g3->data());
  };
};
