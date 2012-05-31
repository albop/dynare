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

#include <iostream>

#include "QRDecomposition.hh"

int
main(int argc, char **argv)
{
  size_t m = 4, n = 3;
  MatrixXd S(m, n), Q(m, m), A(m, n), S2(m, n), B;
  QRDecomposition QRD(m, n, m);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      S(i, j) = i*n + j + 1;

  std::cout << "Matrix to be decomposed:" << std::endl << S << std::endl;

  Q.setIdentity();

  S2 = S;
  QRD.computeAndLeftMultByQ(S2, "N", Q);

  std::cout << "Q =" << std::endl << Q << std::endl;

  B = Q.transpose()*Q;
  
  std::cout << "Q'*Q =" << std::endl << B << std::endl;

  S2 = S2.triangularView<Upper>();

  std::cout << "R =" << std::endl << S2 << std::endl;

  A = Q*S2;

  std::cout << "Q*R =" << std::endl << A << std::endl;

  // Check values
  assert((B - MatrixXd::Identity(m,m)).lpNorm<Infinity>() < 1e-4);

  assert((A-S).lpNorm<Infinity>() < 1e-4);
}
