/*
 * Copyright (C) 2003-2008 Ondra Kamenik
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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/cppblas.h,v 1.2 2004/11/24 20:42:52 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef CPPBLAS_H
#define CPPBLAS_H

#ifdef MATLAB
#include "mex.h"
#endif

#include "../../matlab_versions_compatibility.h"

#if defined(MATLAB) && !defined(__linux__) && !defined(OCTAVE)
#define BLAS_dgemm dgemm
#define BLAS_dgemv dgemv
#define BLAS_dtrsv dtrsv
#define BLAS_dtrmv dtrmv
#define BLAS_daxpy daxpy
#define BLAS_dcopy dcopy
#define BLAS_zaxpy zaxpy
#define BLAS_dscal dscal
#define BLAS_dtrsm dtrsm
#define BLAS_ddot  ddot
#else
#define BLAS_dgemm dgemm_
#define BLAS_dgemv dgemv_
#define BLAS_dtrsv dtrsv_
#define BLAS_dtrmv dtrmv_
#define BLAS_daxpy daxpy_
#define BLAS_dcopy dcopy_
#define BLAS_zaxpy zaxpy_
#define BLAS_dscal dscal_
#define BLAS_dtrsm dtrsm_
#define BLAS_ddot  ddot_
#endif

#if defined NO_BLAS_H
#define BLCHAR const char*
#define CONST_BLINT const int*
#define CONST_BLDOU const double*
#define BLDOU double*

extern "C" {
	void BLAS_dgemm(BLCHAR transa, BLCHAR transb, CONST_BLINT m, CONST_BLINT n,
					CONST_BLINT k, CONST_BLDOU alpha, CONST_BLDOU a, CONST_BLINT lda,
					CONST_BLDOU b, CONST_BLINT ldb, CONST_BLDOU beta,
					BLDOU c, CONST_BLINT ldc);
	void BLAS_dgemv(BLCHAR trans, CONST_BLINT m, CONST_BLINT n, CONST_BLDOU alpha,
					CONST_BLDOU a, CONST_BLINT lda, CONST_BLDOU x, CONST_BLINT incx,
					CONST_BLDOU beta, BLDOU y, CONST_BLINT incy);
	void BLAS_dtrsv(BLCHAR uplo, BLCHAR trans, BLCHAR diag, CONST_BLINT n,
					CONST_BLDOU a, CONST_BLINT lda, BLDOU x, CONST_BLINT incx);
	void BLAS_dtrmv(BLCHAR uplo, BLCHAR trans, BLCHAR diag, CONST_BLINT n,
					CONST_BLDOU a, CONST_BLINT lda, BLDOU x, CONST_BLINT incx);
	void BLAS_daxpy(CONST_BLINT n, CONST_BLDOU a, CONST_BLDOU x, CONST_BLINT incx,
					BLDOU y, CONST_BLINT incy);
	void BLAS_dcopy(CONST_BLINT n, CONST_BLDOU x, CONST_BLINT incx,
					BLDOU y, CONST_BLINT incy);
	void BLAS_zaxpy(CONST_BLINT n, CONST_BLDOU a, CONST_BLDOU x, CONST_BLINT incx,
					BLDOU y, CONST_BLINT incy);
	void BLAS_dscal(CONST_BLINT n, CONST_BLDOU a, BLDOU x, CONST_BLINT incx);
	void BLAS_dtrsm(BLCHAR side, BLCHAR uplo, BLCHAR transa, BLCHAR diag, CONST_BLINT m,
					CONST_BLINT n, CONST_BLDOU alpha, CONST_BLDOU a, CONST_BLINT lda,
					BLDOU b, CONST_BLINT ldb);
	double BLAS_ddot(CONST_BLINT n, CONST_BLDOU x, CONST_BLINT incx, CONST_BLDOU y,
					 CONST_BLINT incy);
};
#else /* NO_BLAS_H isn't defined */
#include "blas.h"
#endif

#endif /* CPPBLAS_H */

// Local Variables:
// mode:C++
// End:
