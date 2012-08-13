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

///////////////////////////////////////////////////////////
//  ComputeModelSolution.h
//  Implementation of the Class ModelSolution
//  Created on:      15-Jan-2010 07:37:47
///////////////////////////////////////////////////////////

#if !defined(ModelSolution_5ADFF920_9C74_46f5_9FE9_88AD4D4BBF19__INCLUDED_)
#define ModelSolution_5ADFF920_9C74_46f5_9FE9_88AD4D4BBF19__INCLUDED_

#include <Eigen/Core>

using namespace Eigen;

#include "DecisionRules.hh"
#include "dynamic_dll.hh"

/**
 * compute the steady state (2nd stage), and
 * computes first order approximation
 *
 */
class ModelSolution
{

public:
  ModelSolution(const std::string &dynamicDllFile, ptrdiff_t n_endo, ptrdiff_t n_exo, const std::vector<ptrdiff_t> &zeta_fwrd_arg,
                const std::vector<ptrdiff_t> &zeta_back_arg, const std::vector<ptrdiff_t> &zeta_mixed_arg,
                const std::vector<ptrdiff_t> &zeta_static_arg, double qz_criterium);
  virtual ~ModelSolution() {};
  template <class Derived1, class Derived2>
  void compute(PlainObjectBase<Derived1> &steadyState, const PlainObjectBase<Derived2> &deepParams,
               MatrixXd &ghx, MatrixXd &ghu) throw (DecisionRules::BlanchardKahnException,
                                                    GeneralizedSchurDecomposition::GSDException)
  {
    // compute Steady State
    ComputeSteadyState(steadyState, deepParams);

    // then get jacobian and

    ComputeModelSolution(steadyState, deepParams, ghx, ghu);

  }

private:
  const ptrdiff_t n_endo;
  const ptrdiff_t n_exo;
  const ptrdiff_t n_jcols; // Num of Jacobian columns
  std::vector<ptrdiff_t> zeta_fwrd_mixed, zeta_back_mixed;
  MatrixXd jacobian;
  VectorXd residual;
  MatrixXd Mx;
  DecisionRules decisionRules;
  DynamicModelDLL dynamicDLLp;
  VectorXd llXsteadyState;
  //Matrix jacobian;
  template <class Derived1, class Derived2>
  void ComputeModelSolution(PlainObjectBase<Derived1> &steadyState, const PlainObjectBase<Derived2> &deepParams,         
			    MatrixXd &ghx, MatrixXd &ghu) 
    throw (DecisionRules::BlanchardKahnException, GeneralizedSchurDecomposition::GSDException)
  {
    // set extended Steady State

    for (ptrdiff_t i = 0; i < (ptrdiff_t) zeta_back_mixed.size(); i++)
      llXsteadyState(i) = steadyState(zeta_back_mixed[i]);

    for (ptrdiff_t i = 0; i < n_endo; i++)
      llXsteadyState(zeta_back_mixed.size() + i) = steadyState(i);

    for (ptrdiff_t i = 0; i < (ptrdiff_t) zeta_fwrd_mixed.size(); i++)
      llXsteadyState(zeta_back_mixed.size() + n_endo + i) = steadyState(zeta_fwrd_mixed[i]);

    //get jacobian
    dynamicDLLp.eval(llXsteadyState, Mx, deepParams, steadyState, residual, &jacobian, NULL, NULL);

    //compute rules
    decisionRules.compute(jacobian, ghx, ghu);
  }
  template <class Derived1, class Derived2>
  void ComputeSteadyState(PlainObjectBase<Derived1> &steadyState, const PlainObjectBase<Derived2> &deepParams)
  {
    // does nothig for time being.
  }


};

#endif // !defined(5ADFF920_9C74_46f5_9FE9_88AD4D4BBF19__INCLUDED_)
