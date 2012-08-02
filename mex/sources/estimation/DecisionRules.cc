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

#include <cassert>

#include <algorithm>

#include "DecisionRules.hh"

DecisionRules::DecisionRules(ptrdiff_t n_arg, ptrdiff_t p_arg,
                             const std::vector<ptrdiff_t> &zeta_fwrd_arg,
                             const std::vector<ptrdiff_t> &zeta_back_arg,
                             const std::vector<ptrdiff_t> &zeta_mixed_arg,
                             const std::vector<ptrdiff_t> &zeta_static_arg,
                             double qz_criterium) :
  n(n_arg), p(p_arg), zeta_fwrd(zeta_fwrd_arg), zeta_back(zeta_back_arg),
  zeta_mixed(zeta_mixed_arg), zeta_static(zeta_static_arg),
  n_fwrd(zeta_fwrd.size()), n_back(zeta_back.size()),
  n_mixed(zeta_mixed.size()), n_static(zeta_static.size()),
  n_back_mixed(n_back+n_mixed), n_fwrd_mixed(n_fwrd+n_mixed),
  n_dynamic(n-n_static),
  S(n, n_static),
  A(n, n_back_mixed + n + n_fwrd_mixed),
  D(n_fwrd + n_back + 2*n_mixed, n_fwrd + n_back + 2*n_mixed),
  E(n_fwrd + n_back + 2*n_mixed, n_fwrd + n_back + 2*n_mixed),
  Z_prime(n_fwrd + n_back + 2*n_mixed, n_fwrd + n_back + 2*n_mixed),
  QR(n, n_static, n_back_mixed + n + n_fwrd_mixed),
  GSD(n_fwrd + n_back + 2*n_mixed, qz_criterium),
  Solver1(n_fwrd_mixed, n_fwrd_mixed),
  Solver2(n_back_mixed, n_back_mixed),
  Solver3(n_static, n_static),
  Solver4(n, n),
  Z21(n_fwrd_mixed, n_back_mixed),
  g_y_back(n_back_mixed, n_back_mixed),
  g_y_fwrd(n_fwrd_mixed, n_fwrd_mixed),
  g_y_static(n_static, n_back_mixed),
  A0s(n_static, n_static),
  A0d(n_static, n_dynamic),
  g_y_dynamic(n_dynamic, n_back_mixed),
  g_u_tmp1(n, n_back_mixed),
  g_u_tmp2(n, n)
{
  assert(n == n_back + n_fwrd + n_mixed + n_static);

  set_union(zeta_fwrd.begin(), zeta_fwrd.end(),
            zeta_mixed.begin(), zeta_mixed.end(),
            back_inserter(zeta_fwrd_mixed));
  set_union(zeta_back.begin(), zeta_back.end(),
            zeta_mixed.begin(), zeta_mixed.end(),
            back_inserter(zeta_back_mixed));
  set_union(zeta_back_mixed.begin(), zeta_back_mixed.end(),
            zeta_fwrd.begin(), zeta_fwrd.end(),
            back_inserter(zeta_dynamic));

  // Compute beta_back and pi_back
  for (ptrdiff_t i = 0; i < n_back_mixed; i++)
    if (find(zeta_mixed.begin(), zeta_mixed.end(), zeta_back_mixed[i])
        == zeta_mixed.end())
      pi_back.push_back(i);
    else
      beta_back.push_back(i);

  // Compute beta_fwrd and pi_fwrd
  for (ptrdiff_t i = 0; i < n_fwrd_mixed; i++)
    if (find(zeta_mixed.begin(), zeta_mixed.end(), zeta_fwrd_mixed[i])
        == zeta_mixed.end())
      pi_fwrd.push_back(i);
    else
      beta_fwrd.push_back(i);
}

void
DecisionRules::compute(const MatrixXd &jacobian, MatrixXd &g_y, MatrixXd &g_u) throw (BlanchardKahnException, GeneralizedSchurDecomposition::GSDException)
{
  assert(jacobian.rows() == n
         && jacobian.cols() == (n_back_mixed + n + n_fwrd_mixed + p));
  assert(g_y.rows() == n && g_y.cols() == n_back_mixed);
  assert(g_u.rows() == n && g_u.cols() == p);

  // Construct S, perform QR decomposition and get A = Q*jacobian
  for (ptrdiff_t i = 0; i < n_static; i++)
    S.col(i) = jacobian.col(n_back_mixed+zeta_static[i]);

  A = jacobian.block(0, 0, n, n_back_mixed + n + n_fwrd_mixed);
  QR.computeAndLeftMultByQ(S, "T", A);

  // Construct matrix D
  D.setZero();
  for (ptrdiff_t i = 0; i < n_mixed; i++)
    D(n - n_static + i, beta_back[i]) = 1.0;
  for (ptrdiff_t j = 0; j < n_back_mixed; j++)
    D.col(j).head(n-n_static) = A.col(n_back_mixed + zeta_back_mixed[j]).tail(n-n_static);
  
  D.block(0, n_back_mixed, n - n_static, n_fwrd_mixed) = A.block(n_static, n_back_mixed + n, n - n_static, n_fwrd_mixed);

  // Construct matrix E
  E.setZero();
  for (ptrdiff_t i = 0; i < n_mixed; i++)
    E(n - n_static + i, n_back_mixed + beta_fwrd[i]) = 1.0;
  E.block(0, 0, n - n_static, n_back_mixed) = -A.block(n_static, 0, n - n_static, n_back_mixed);
  for (ptrdiff_t j = 0; j < n_fwrd; j++)
    E.col(n_back_mixed + pi_fwrd[j]).head(n-n_static) = -A.col(n_back_mixed + zeta_fwrd_mixed[pi_fwrd[j]]).tail(n - n_static);

  // Perform the generalized Schur
  ptrdiff_t sdim;
  GSD.compute(E, D, Z_prime, sdim);

  if (n_back_mixed != sdim)
    throw BlanchardKahnException(true, n_fwrd_mixed, n_fwrd + n_back + 2*n_mixed - sdim);

  // Compute DR for forward variables w.r. to endogenous
  Solver1.compute(Z_prime.block(n_back_mixed, n_back_mixed, n_fwrd_mixed, n_fwrd_mixed).transpose());
  
  if (!Solver1.isInvertible())
    throw BlanchardKahnException(false, n_fwrd_mixed, n_fwrd + n_back + 2*n_mixed - sdim);

  g_y_fwrd = -Solver1.solve(Z_prime.block(0, n_back_mixed, n_back_mixed, n_fwrd_mixed).transpose());

  for (ptrdiff_t i = 0; i < n_fwrd_mixed; i++)
    g_y.row(zeta_fwrd_mixed[i]) = g_y_fwrd.row(i);

  // Compute DR for backward variables w.r. to endogenous
  Solver2.compute(Z_prime.block(0, 0, n_back_mixed, n_back_mixed));
  g_y_back.noalias() = E.block(0, 0, n_back_mixed, n_back_mixed) * Solver2.solve(MatrixXd::Identity(n_back_mixed, n_back_mixed)); // S_{11} * Z'_{11}^{-1}
  
  Solver2.compute(D.block(0, 0, n_back_mixed, n_back_mixed));
  g_y_back = Z_prime.block(0, 0, n_back_mixed, n_back_mixed) * Solver2.solve(g_y_back);
  
  // TODO: avoid to copy mixed variables again, rather test it...
  for (ptrdiff_t i = 0; i < n_back_mixed; i++)
    g_y.row(zeta_back_mixed[i]) = g_y_back.row(i);

  // Compute DR for static variables w.r. to endogenous
  for (ptrdiff_t i = 0; i < n_dynamic; i++)
    {
      g_y_dynamic.row(i) = g_y.row(zeta_dynamic[i]);
      A0d.col(i) = A.col(n_back_mixed + zeta_dynamic[i]).head(n_static);
    }

  for (ptrdiff_t i = 0; i < n_static; i++)
    A0s.col(i) = A.col(n_back_mixed + zeta_static[i]).head(n_static);
  
  Solver3.compute(A0s);
  
  g_y_static = -Solver3.solve(A.block(0, n_back_mixed + n, n_static, n_fwrd_mixed)*g_y_fwrd*g_y_back
                              + A0d*g_y_dynamic + A.block(0, 0, n_static, n_back_mixed));

  for (ptrdiff_t i = 0; i < n_static; i++)
    g_y.row(zeta_static[i]) = g_y_static.row(i);
  
  // Compute DR for all endogenous w.r. to shocks
  g_u_tmp1 = jacobian.block(0, n_back_mixed + n, n, n_fwrd_mixed) * g_y_fwrd;
  
  g_u_tmp2 = jacobian.block(0, n_back_mixed, n, n);
  for (ptrdiff_t i = 0; i < n_back_mixed; i++)
    g_u_tmp2.col(zeta_back_mixed[i]) += g_u_tmp1.col(i);
  Solver4.compute(g_u_tmp2);
  g_u = -Solver4.solve(jacobian.block(0, n_back_mixed + n + n_fwrd_mixed, n, p));
}

std::ostream &
operator<<(std::ostream &out, const DecisionRules::BlanchardKahnException &e)
{
  if (e.order)
    out << "The Blanchard-Kahn order condition is not satisfied: you have " << e.n_fwrd_vars << " forward variables for " << e.n_explosive_eigenvals << " explosive eigenvalues";
  else
    out << "The Blanchard Kahn rank condition is not satisfied";
  return out;
}

