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

#include <algorithm> // For std::min()

#include <Eigen/Core>

using namespace Eigen;

#include <dynlapack.h>

/* Note that Eigen provides built-in QR decomposers (classes *HouseHolderQR).
   However, we need to have our own because we need the ability to
   left-multiply arbitrary matrices by Q, and this is not offered in a
   convenient way by the Eigen classes. */
class QRDecomposition
{
private:
  const ptrdiff_t rows, cols, mind;
  //! Number of columns of the matrix to be left-multiplied by Q
  const ptrdiff_t cols2;
  lapack_int lwork, lwork2;
  double *work, *work2, *tau;
public:
  /*!
    \todo Replace heuristic choice for workspace size by a query to determine the optimal size
    \param[in] rows_arg Number of rows of the matrix to decompose
    \param[in] cols_arg Number of columns of the matrix to decompose
    \param[in] cols2_arg Number of columns of the matrix to be multiplied by Q
  */
  QRDecomposition(ptrdiff_t rows_arg, ptrdiff_t cols_arg, ptrdiff_t cols2_arg);
  virtual ~QRDecomposition();
  //! Performs the QR decomposition of a matrix, and left-multiplies another matrix by Q
  /*!
    \param[in,out] A On input, the matrix to be decomposed. On output, equals to the output of dgeqrf
    \param[in] trans Specifies whether Q should be transposed before the multiplication, either "T" or "N"
    \param[in,out] C The matrix to be left-multiplied by Q, modified in place
  */
  template<typename Derived1, typename Derived2>
  void computeAndLeftMultByQ(PlainObjectBase<Derived1> &A, const char *trans,
                             PlainObjectBase<Derived2> &C);
};

template<typename Derived1, typename Derived2>
void
QRDecomposition::computeAndLeftMultByQ(PlainObjectBase<Derived1> &A, const char *trans,
                                       PlainObjectBase<Derived2> &C)
{
  assert(A.rows() == rows && A.cols() == cols);
  assert(C.rows() == rows && C.cols() == cols2);

  lapack_int m = rows, n = cols, lda = A.outerStride();
  lapack_int info;
  dgeqrf(&m, &n, A.data(), &lda, tau, work, &lwork, &info);
  assert(info == 0);

  n = cols2;
  lapack_int k = mind, ldc = C.outerStride();
  dormqr("L", trans, &m, &n, &k, A.data(), &lda, tau, C.data(), &ldc,
         work2, &lwork2, &info);
  assert(info == 0);
}
