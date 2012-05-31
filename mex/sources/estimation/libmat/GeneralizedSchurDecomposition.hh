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

#include <dynlapack.h>

#include <Eigen/Core>

using namespace Eigen;

class GeneralizedSchurDecomposition
{
private:
  const ptrdiff_t n;
  const double criterium;
  lapack_int lwork;
  double *alphar, *alphai, *beta, *vsl, *work;
  lapack_int *bwork;
  static double criterium_static;
  static lapack_int selctg(const double *alphar, const double *alphai, const double *beta);
public:
  class GSDException
  {
  public:
    const lapack_int info, n;
    GSDException(lapack_int info_arg, lapack_int n_arg) : info(info_arg), n(n_arg) {};
  };
  //! \todo Replace heuristic choice for workspace size by a query to determine the optimal size
  GeneralizedSchurDecomposition(ptrdiff_t n_arg, double criterium_arg);
  virtual ~GeneralizedSchurDecomposition();
  //! \todo Add a lock around the modification of criterium_static for making it thread-safe
  template<typename Mat1, typename Mat2, typename Mat3>
  void compute(PlainObjectBase<Mat1> &S, PlainObjectBase<Mat2> &T,
               PlainObjectBase<Mat3> &Z, ptrdiff_t &sdim) throw (GSDException);

  /*!
    \param[out] sdim Number of non-explosive generalized eigenvalues
  */
  template<typename Mat1, typename Mat2, typename Mat3, typename Mat4, typename Mat5>
  void compute(const MatrixBase<Mat1> &D, const MatrixBase<Mat2> &E,
               PlainObjectBase<Mat3> &S, PlainObjectBase<Mat4> &T,
               PlainObjectBase<Mat5> &Z, ptrdiff_t &sdim) throw (GSDException);
 
  template<typename Vec1, typename Vec2>
  void getGeneralizedEigenvalues(PlainObjectBase<Vec1> &eig_real,
                                 PlainObjectBase<Vec2> &eig_cmplx);
};

std::ostream &operator<<(std::ostream &out, const GeneralizedSchurDecomposition::GSDException &e);

template<typename Mat1, typename Mat2, typename Mat3>
void
GeneralizedSchurDecomposition::compute(PlainObjectBase<Mat1> &S, PlainObjectBase<Mat2> &T, PlainObjectBase<Mat3> &Z, ptrdiff_t &sdim) throw (GSDException)
{
  assert(S.rows() == n && S.cols() == n
         && T.rows() == n && T.cols() == n
         && Z.rows() == n && Z.cols() == n);

  lapack_int n2 = n;
  lapack_int info, sdim2;
  lapack_int lds = S.outerStride(), ldt = T.outerStride(), ldz = Z.outerStride();

  criterium_static = criterium;
  // Here we are forced to give space for left Schur vectors, even if we don't use them, because of a bug in dgges()
  dgges("N", "V", "S", &selctg, &n2, S.data(), &lds, T.data(), &ldt,
        &sdim2, alphar, alphai, beta, vsl, &n2, Z.data(), &ldz,
        work, &lwork, bwork, &info);

  if (info != 0)
    throw GSDException(info, n2);

  sdim = sdim2;
}

template<typename Mat1, typename Mat2, typename Mat3, typename Mat4, typename Mat5>
void
GeneralizedSchurDecomposition::compute(const MatrixBase<Mat1> &D, const MatrixBase<Mat2> &E,
                                       PlainObjectBase<Mat3> &S, PlainObjectBase<Mat4> &T, PlainObjectBase<Mat5> &Z, ptrdiff_t &sdim) throw (GSDException)
{
  assert(D.rows() == n && D.cols() == n
         && E.rows() == n && E.cols() == n);
  S = D;
  T = E;
  compute(S, T, Z, sdim);
}

template<typename Vec1, typename Vec2>
void
GeneralizedSchurDecomposition::getGeneralizedEigenvalues(PlainObjectBase<Vec1> &eig_real,
                                                         PlainObjectBase<Vec2> &eig_cmplx)
{
  assert(eig_real.size() == n && eig_cmplx.size() == n);

  double *par = alphar, *pai = alphai, *pb = beta,
    *per = eig_real.data(), *pei = eig_cmplx.data();
  while (par < alphar + n)
    {
      *per = *par / *pb;
      if (*pai == 0.0 && *pb == 0.0)
        *pei = 0.0;
      else
        *pei = *pai / *pb;
      par++;
      pai++;
      pb++;
      per += eig_real.innerStride();
      pei += eig_cmplx.innerStride();
    }
}
