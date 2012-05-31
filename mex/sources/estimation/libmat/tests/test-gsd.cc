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

#include "GeneralizedSchurDecomposition.hh"

int
main(int argc, char **argv)
{
  Matrix3d D, E;
  D << 1, 2, 3,
    4, 5, 6,
    7, 8, 9;
  
  E << 1, -3, 4,
    -7, 9, 1,
    -3, 4, 0;

  std::cout << "D =" << std::endl << D << std::endl;
  std::cout << "E =" << std::endl << E << std::endl;

  GeneralizedSchurDecomposition GSD(3, 1.00001);

  Matrix3d S, T, Z;
  ptrdiff_t sdim;

  GSD.compute(D, E, S, T, Z, sdim);

  std::cout << "S =" << std::endl << S << std::endl;
  std::cout << "T =" << std::endl << T << std::endl;
  std::cout << "Z =" << std::endl << Z << std::endl;

  Vector3d eig_real, eig_cmplx;
  GSD.getGeneralizedEigenvalues(eig_real, eig_cmplx);

  std::cout << "Real part of generalized eigenvalues: " << std::endl << eig_real << std::endl;
  std::cout << "Complex part of generalized eigenvalues: " << std::endl << eig_cmplx << std::endl;

  return 0;
}
