/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include <typeinfo>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iomanip>

#if OPENSMOKE_USE_MKL == 1
	#include "mkl.h"
	#include "mkl_lapacke.h"
	#if defined(_WIN32) || defined(_WIN64) 
	#else
		#include "mm_malloc.h"
	#endif
#elif OPENSMOKE_USE_OPENBLAS == 1
	#include "cblas.h"
	#include "lapacke.h"
	#include "malloc.h"
#endif

namespace OpenSMOKE
{
	#define BAND_COL(A,j) (((A->cols)[j])+(A->s_mu))
	#define BAND_COL_ELEM(col_j,i,j) (col_j[(i)-(j)])
	#define BAND_ELEM(A,i,j) ((A->cols)[j][(i)-(j)+(A->s_mu)])
	#define ROW(i,j,smu) (i-j+smu)

	#if defined(_WIN32) || defined(_WIN64) 
		// Nothing to declare
	#else
		template<typename TT>
		OpenSMOKEBandMatrix<TT>* NewBandMat(const int N, const int mu, const int ml, const int smu);
		template<typename TT>
		void bandMatVec(TT **a, const TT *x, TT *y, const int n, const int mu, const int ml, const int smu);
		template<typename TT>
		void bandMatTransposeVec(TT **a, const TT *x, TT *y, const int n, const int mu, const int ml, const int smu);
	#endif
 
	template<typename T>
	int bandGBTRS(T **a, const int n, const int smu, const int ml, const int *p, T *b);

        template<typename T>
        int bandGBTRF(T **a, const int n, const int mu, const int ml, const int smu, int *p);

        template<typename T>
        void bandCopy(T **a, T **b, const int n, const int a_smu, const int b_smu, const int copymu, const int copyml);

	template<typename T>
	OpenSMOKEBandMatrix<T>::OpenSMOKEBandMatrix(const int nEquations, const int dimBlock)
	{
		isTridiagonalBlock_ = true;

		const int nUpper = dimBlock * 2 - 1;
		const int nLower = dimBlock * 2 - 1;
		const int smu = std::min(nEquations - 1, nUpper + nLower);

		#if defined(_WIN32) || defined(_WIN64) 
		*this = *NewBandMat<T>(nEquations, nUpper, nLower, smu);
		#else
		*this = *NewBandMat<T>(nEquations, nUpper, nLower, smu);
		#endif
	}

	template<typename T>
	OpenSMOKEBandMatrix<T>::OpenSMOKEBandMatrix(const int nEquations, const int nUpper, const int nLower)
	{
		isTridiagonalBlock_ = false;
		const int smu = std::min(nEquations - 1, nUpper + nLower);

		#if defined(_WIN32) || defined(_WIN64) 
		*this = *NewBandMat<T>(nEquations, nUpper, nLower, smu);
		#else
		*this = *NewBandMat<T>(nEquations, nUpper, nLower, smu);
		#endif
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::AddIdentity()
	{
		for (int i = 0; i<M; i++)
			cols[i][s_mu] += 1.;
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::AddDiagonal(const T *d)
	{
		for (int i = 0; i<M; i++)
			cols[i][s_mu] += d[i];
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::Print(std::ostream& out)
	{
		T **a;
		a = cols;
		
		for (int i = 0; i < N; i++)
		{
			const int start = std::max((int)(0), i - ml);
			const int finish = std::min(N - 1, i + mu);
			
			for (int j = 0; j < start; j++)
				out << std::setw(12) << "";
			for (int j = start; j <= finish; j++)
				out << std::setw(12) << a[j][i - j + s_mu];
			out << std::endl;
		}
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::Scale(const double c_differential, const double c_algebraic, const int *index)
	{
		int width = mu + ml + 1;
		int ngroups = std::min(width, N);

		for (int group = 1; group <= ngroups; group++)
		{
			for (int j = group - 1; j < N; j += width)
			{
				double* col_j = BAND_COL(this, j);
				
				int i1 = std::max((int)(0), j - mu);
				int i2 = std::min(j + ml, N);
				for (int i = i1; i <= i2; i++)
				{
					if (index[i] == 1)
						BAND_COL_ELEM(col_j, i, j) *= c_differential;
					else
						BAND_COL_ELEM(col_j, i, j) *= c_algebraic;
				}
			}
		}
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::SetToZero()
	{
		const T ZERO = 0.;
		/*		
		int i, j, colSize;
		T *col_j;

		colSize = mu + ml + 1;
		for (j = 0; j<M; j++)
		{
			col_j = cols[j] + s_mu - mu;
			for (i = 0; i<colSize; i++)
				col_j[i] = ZERO;
		}
		*/
		for (int j = 0; j < ldata; j++)
			data[j] = ZERO;
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::DestroyMat()
	{
		#if (OPENSMOKE_USE_MKL == 1)
			_mm_free(data);
			data = NULL;
			_mm_free(cols);
			cols = NULL;
			_mm_free(p);
			p = NULL;
		#else
			free(data);
			data = NULL;
			free(cols);
			cols = NULL;
			free(p);
			p = NULL;
		#endif
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::CopyTo(OpenSMOKEBandMatrix<T>* B)
	{
//		#if __APPLE__
//
//			bandCopy(cols, B->cols, M, s_mu, B->s_mu, mu, ml);
//
//		#else

			#if (OPENSMOKE_USE_MKL == 1)

				cblas_dcopy(ldata, data, 1, B->data, 1);

			#elif (OPENSMOKE_USE_OPENBLAS == 1)

				cblas_dcopy(ldata, data, 1, B->data, 1);

			#else

				bandCopy(cols, B->cols, M, s_mu, B->s_mu, mu, ml);

			#endif

//		#endif
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::Scale(const double c)
	{
		#if __APPLE__

			cblas_dscal(ldata, c, data, 1);

		#elif (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)

			const int one = 1;
			const int lmat = ldata;
			dscal(&lmat, &c, data, &one);

		#else

			bandScale(c, cols, M, mu, ml, s_mu);

		#endif
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::Product(const T *v_in, T *v_out)
	{
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)

			// In principle should work, but it does not
			// I should better investigate the problem
			// For now I use the slower version
			//cblas_dgbmv(CblasColMajor, CblasNoTrans, M, M, ml, mu, 1., data, ldim, v_in, 1, 0., v_out, 1);

			bandMatVec(cols, v_in, v_out, M, mu, ml, s_mu);

		#else

			bandMatVec(cols, v_in, v_out, M, mu, ml, s_mu);

		#endif
	}

	template<typename T>
	void OpenSMOKEBandMatrix<T>::TProduct(const T *v_in, T *v_out)
	{
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)

			// In principle should work, but it does not
			// I should better investigate the problem
			// For now I use the slower version
			// cblas_dgbmv(CblasColMajor, CblasTrans, M, M, ml, mu, 1., data, ldim, v_in, 1, 0., v_out, 1);

			bandMatTransposeVec(cols, v_in, v_out, M, mu, ml, s_mu);
				
		#else

			bandMatTransposeVec(cols, v_in, v_out, M, mu, ml, s_mu);

		#endif
	}

	template<typename T>
	int OpenSMOKEBandMatrix<T>::Factorize()
	{
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)

		const int ier = LAPACKE_dgbtrf(CblasColMajor, M, M, ml, mu, data, ldim, p);

		if (ier != 0)
		{
			if (ier > 0)
			{
				std::cout << "Pivot " << ier << " is equal to zero (singular matrix). The factorization has been completed, but U is exactly singular." << std::endl;
				std::cout << "The solution cannot be calculated." << std::endl;
			}
			else
			{
				std::cout << "Parameter " << -ier << " has an illegal value" << std::endl;

				std::cout << " (1) Matrix layout: " << CblasColMajor << std::endl;
				std::cout << " (2) Rows:          " << M << std::endl;
				std::cout << " (3) Columns:       " << M << std::endl;
				std::cout << " (4) Lower band:    " << ml << std::endl;
				std::cout << " (5) Upper band:    " << mu << std::endl;
				std::cout << " (6) Data[0]:       " << data[0] << std::endl;
				std::cout << " (7) Leading dim.:  " << ldim << std::endl;
				std::cout << " (8) Pivot[0]:      " << p[0] << std::endl;
			}

			// OpenSMOKE::FatalErrorMessage("Factorizing the banded linear system: LAPACKE_dgbtrf failed");
		}

		#else
		
			const int ier = bandGBTRF(cols, M, mu, ml, s_mu, p);
			if (ier > 0)
				OpenSMOKE::FatalErrorMessage("Factorizing the banded linear system: BandGBTRF failed");

		#endif

		return ier;
	}

	template<typename T>
	int OpenSMOKEBandMatrix<T>::Solve(T *b)
	{
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)

			const int one = 1;
			const int ier = LAPACKE_dgbtrs(CblasColMajor, 'N', M, ml, mu, one, data, ldim, p, b, M);

			if (ier != 0)
			{
				std::cout << "Parameter " << -ier << " has an illegal value" << std::endl;

				std::cout << " (1) Matrix layout: " << CblasColMajor << std::endl;
				std::cout << " (2) Transpose:     " << "N" << std::endl;
				std::cout << " (3) Size:          " << M << std::endl;
				std::cout << " (4) Lower band:    " << ml << std::endl;
				std::cout << " (5) Upper band:    " << mu << std::endl;
				std::cout << " (6) #RHS:          " << one << std::endl;
				std::cout << " (7) Data[0]:       " << data[0] << std::endl;
				std::cout << " (8) Leading dim.:  " << ldim << std::endl;
				std::cout << " (9) Pivot[0]:      " << p[0] << std::endl;
				std::cout << " (10) b[0]:         " << b[0] << std::endl;
				std::cout << " (11) Size:         " << M << std::endl;

				// OpenSMOKE::FatalErrorMessage("Solving the banded linear system: LAPACKE_dgbtrs failed");
			}

		#else
			
		const int ier = bandGBTRS(cols, M, s_mu, ml, p, b);

		#endif

		return ier;
	}

	template<typename T>
	int OpenSMOKEBandMatrix<T>::Solve(const int nrhs, T *b)
	{
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)

		const int ier = LAPACKE_dgbtrs(CblasColMajor, 'N', M, ml, mu, nrhs, data, ldim, p, b, M);

		if (ier != 0)
		{
			std::cout << "Parameter " << -ier << " has an illegal value" << std::endl;

			std::cout << " (1) Matrix layout: " << CblasColMajor << std::endl;
			std::cout << " (2) Transpose:     " << "N" << std::endl;
			std::cout << " (3) Size:          " << M << std::endl;
			std::cout << " (4) Lower band:    " << ml << std::endl;
			std::cout << " (5) Upper band:    " << mu << std::endl;
			std::cout << " (6) #RHS:          " << nrhs << std::endl;
			std::cout << " (7) Data[0]:       " << data[0] << std::endl;
			std::cout << " (8) Leading dim.:  " << ldim << std::endl;
			std::cout << " (9) Pivot[0]:      " << p[0] << std::endl;
			std::cout << " (10) b[0]:         " << b[0] << std::endl;
			std::cout << " (11) Size:         " << M << std::endl;

			//OpenSMOKE::FatalErrorMessage("Solving the banded linear system: LAPACKE_dgbtrs failed");
		}

		#else

		OpenSMOKE::FatalErrorMessage("Solution of multiple rhs with banded structure requires MKL or OpenBlas");
		const int ier = bandGBTRS(cols, M, s_mu, ml, p, b);

		#endif

		return ier;
	}

	// This function must be checked!
	// Please do not use it
	template<typename T>
	int OpenSMOKEBandMatrix<T>::FactorizeAndSolve(T *b)
	{
		/*
		#if (OPENSMOKE_USE_MKL == 1 || OPENSMOKE_USE_OPENBLAS == 1)

			double* r = new double[M];
			double* c = new double[M];
			const int ldafb = 2*iml + imu + 1;

			double* afb = new double[ldafb*M];
			double* berr = new double[1];
			double* ferr = new double[1];
			double* x = new double[M];

			double rcond;
			double rpivot;
			char equed;

			const int ier = LAPACKE_dgbsvx(CblasColMajor, 'E', 'N', M, ml, mu, 1, data, ldim, afb, ldafb, p, &equed, r, c, b, M, x, M, &rcond, ferr, berr, &rpivot);
				
			if (ier != 0)
			{
				std::cout << "Equilibration type: " << equed << std::endl;
				std::cout << "1/cond:             " << rcond << std::endl;
				std::cout << "1/pivot:            " << rpivot << std::endl;

				if (ier > 0 && ier <= M)
				{
					std::cout << "Pivot " << ier << " is equal to zero (singular matrix). The factorization has been completed, but U is exactly singular." << std::endl;
					std::cout << "The solution cannot be calculated." << std::endl;
				}
				else if (ier == M + 1)
				{
					std::cout << "U is nonsingular, but the reciprocal of condition number is less than machine precision, meaning that the matrix is" << std::endl;
					std::cout << "singular working precision.Nevertheless, the solution and error bounds are computed because there are a number of " << std::endl;
					std::cout << "where the computed solution can be more accurate than the value of rcond would suggest." << std::endl;
				}
				else
				{
					std::cout << "Parameter " << -ier << " has an illegal value" << std::endl;
				}

				OpenSMOKE::FatalErrorMessage("Factorizing and solving the banded linear system: LAPACKE_dgbsvx failed");
			}

			for(unsigned int i=0;i<M;i++)
				b[i] = x[i];

		#else

			OpenSMOKE::FatalErrorMessage("Factorizing and solving the banded linear system: bandGBTRS(cols, M, s_mu, ml, p, b) not implemented");

		#endif
		*/

		const int ier = -1;
		OpenSMOKE::FatalErrorMessage("Factorizing and solving the banded linear system: bandGBTRS(cols, M, s_mu, ml, p, b) not implemented");

		return ier;
	}

	template<typename T>
	OpenSMOKEBandMatrix<T>* NewBandMat(const int N, const int mu, const int ml, const int smu)
	{
		OpenSMOKEBandMatrix<T>* A;
		int j, colSize;

		if (N <= 0) return(NULL);
		colSize = smu + ml + 1;
		
		#if (OPENSMOKE_USE_MKL == 1)
			A = NULL;
			A = (OpenSMOKEBandMatrix<T>*)_mm_malloc(sizeof *A, 64);
			if (A == NULL) return (NULL);
		#else
			A = NULL;
			A = (OpenSMOKEBandMatrix<T>*)malloc(sizeof *A);
			if (A == NULL) return (NULL);
		#endif
		
		#if (OPENSMOKE_USE_MKL == 1)
			A->data = NULL;
			A->data = (T *)_mm_malloc(N * colSize * sizeof(T),64);
			if (A->data == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL); 
			}
		#else
			A->data = NULL;
			A->data = (T *)malloc(N * colSize * sizeof(T));
			if (A->data == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}
		#endif

		#if (OPENSMOKE_USE_MKL == 1)
			A->cols = NULL;
			A->cols = (T **)_mm_malloc(N * sizeof(T *),64);
			if (A->cols == NULL)
			{
				_mm_free(A->data);
				_mm_free(A);
				A = NULL;
				return(NULL);
			}
		#else
			A->cols = NULL;
			A->cols = (T **)malloc(N * sizeof(T *));
			if (A->cols == NULL)
			{
				free(A->data);
				free(A); 
				A = NULL;
				return(NULL);
			}
		#endif

		for (j = 0; j < N; j++)
			A->cols[j] = A->data + j * colSize;

		#if (OPENSMOKE_USE_MKL == 1)
			A->p = NULL;
			A->p = (int *)_mm_malloc(N * sizeof(int),64);
			if (A->p == NULL)
			{
				_mm_free(A);
				A = NULL;
				return(NULL);
			}
		#else
			A->p = NULL;
			A->p = (int *)malloc(N * sizeof(int));
			if (A->p == NULL)
			{
				free(A);
				A = NULL;
				return(NULL);
			}
		#endif

		A->M = N;
		A->N = N;
		A->mu = mu;
		A->ml = ml;
		A->s_mu = smu;
		A->ldim = colSize;
		A->ldata = N * colSize;

		return(A);
	}

	template<typename T>
	void bandCopy(T **a, T **b, const int n, const int a_smu, const int b_smu, const int copymu, const int copyml)
	{
		int i, j, copySize;
		T *a_col_j, *b_col_j;

		copySize = copymu + copyml + 1;

		for (j = 0; j < n; j++) 
		{
			a_col_j = a[j] + a_smu - copymu;
			b_col_j = b[j] + b_smu - copymu;
			for (i = 0; i < copySize; i++)
				b_col_j[i] = a_col_j[i];
		}
	}

	template<typename T>
	void bandScale(const T c, T **a, const int n, const int mu, const int ml, const int smu)
	{
		int i, j, colSize;
		T *col_j;

		colSize = mu + ml + 1;

		for (j = 0; j < n; j++) 
		{
			col_j = a[j] + smu - mu;
			for (i = 0; i < colSize; i++)
				col_j[i] *= c;
		}
	}

	template<typename T>
	void bandAddIdentity(T **a, const int n, const int smu)
	{
        const T ONE = 1.;

		for (int j = 0; j < n; j++)
			a[j][smu] += ONE;
	}

	template<typename T>
	void bandMatVec(T **a, const T *x, T *y, const int n, const int mu, const int ml, const int smu)
	{
		int i, j, is, ie;
		T *col_j;

		for (i = 0; i<n; i++)
			y[i] = 0.0;

		for (j = 0; j<n; j++) 
		{
			col_j = a[j] + smu - mu;
			is = (0 > j - mu) ? 0 : j - mu;
			ie = (n - 1 < j + ml) ? n - 1 : j + ml;
			for (i = is; i <= ie; i++)
				y[i] += col_j[i - j + mu] * x[j];
		}
	}

	template<typename T>
	void bandMatTransposeVec(T **a, const T *x, T *y, const int n, const int mu, const int ml, const int smu)
	{
		int i, j, is, ie;
		T *col_j;

		for (i = 0; i<n; i++)
			y[i] = 0.0;

		for (j = 0; j<n; j++)
		{
			col_j = a[j] + smu - mu;
			is = (0 > j - mu) ? 0 : j - mu;
			ie = (n - 1 < j + ml) ? n - 1 : j + ml;
			for (i = is; i <= ie; i++)
				y[j] += col_j[i - j + mu] * x[i];
		}
	}

	// Factorization
	template<typename T>
	int bandGBTRF(T **a, const int n, const int mu, const int ml, const int smu, int *p)
	{
		const T ZERO = 0.;
		const T ONE = 1.;
		int c, r, num_rows;
		int i, j, k, l, storage_l, storage_k, last_col_k, last_row_k;
		T *a_c, *col_k, *diag_k, *sub_diag_k, *col_j, *kptr, *jptr;
		T max, temp, mult, a_kj;
		bool swap;

		/* zero out the first smu - mu rows of the rectangular array a */
		num_rows = smu - mu;
		if (num_rows > 0)
		{
			for (c = 0; c < n; c++)
			{
				a_c = a[c];
				for (r = 0; r < num_rows; r++)
				{
					a_c[r] = ZERO;
				}
			}
		}

		/* k = elimination step number */
		for (k = 0; k < n - 1; k++, p++)
		{
			col_k = a[k];
			diag_k = col_k + smu;
			sub_diag_k = diag_k + 1;
			last_row_k = std::min(n - 1, k + ml);

			/* find l = pivot row number */
			l = k;
			max = std::fabs(*diag_k);
			for (i = k + 1, kptr = sub_diag_k; i <= last_row_k; i++, kptr++)
			{
				if (std::fabs(*kptr) > max)
				{
					l = i;
					max = std::fabs(*kptr);
				}
			}
			storage_l = ROW(l, k, smu);
			*p = l;

			/* check for zero pivot element */
			if (col_k[storage_l] == ZERO) return(k + 1);

			/* swap a(l,k) and a(k,k) if necessary */
			if ((swap = (l != k))) {
				temp = col_k[storage_l];
				col_k[storage_l] = *diag_k;
				*diag_k = temp;
			}

			/* Scale the elements below the diagonal in         */
			/* column k by -1.0 / a(k,k). After the above swap, */
			/* a(k,k) holds the pivot element. This scaling     */
			/* stores the pivot row multipliers -a(i,k)/a(k,k)  */
			/* in a(i,k), i=k+1, ..., SUNMIN(n-1,k+ml).            */
			mult = -ONE / (*diag_k);
			for (i = k + 1, kptr = sub_diag_k; i <= last_row_k; i++, kptr++)
				(*kptr) *= mult;

			/* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., SUNMIN(n-1,k+ml) */
			/* row k is the pivot row after swapping with row l.                */
			/* The computation is done one column at a time,                    */
			/* column j=k+1, ..., SUNMIN(k+smu,n-1).                               */
			last_col_k = std::min(k + smu, n - 1);
			for (j = k + 1; j <= last_col_k; j++)
			{

				col_j = a[j];
				storage_l = ROW(l, j, smu);
				storage_k = ROW(k, j, smu);
				a_kj = col_j[storage_l];

				/* Swap the elements a(k,j) and a(k,l) if l!=k. */
				if (swap)
				{
					col_j[storage_l] = col_j[storage_k];
					col_j[storage_k] = a_kj;
				}

				/* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j) */
				/* a_kj = a(k,j), *kptr = - a(i,k)/a(k,k), *jptr = a(i,j) */

				if (a_kj != ZERO)
				{
					for (i = k + 1, kptr = sub_diag_k, jptr = col_j + ROW(k + 1, j, smu);
						i <= last_row_k;
						i++, kptr++, jptr++)
						(*jptr) += a_kj * (*kptr);
				}
			}
		}

		return -1;
	}

	// Solution 
	template<typename T>
	int bandGBTRS(T **a, const int n, const int smu, const int ml, const int *p, T *b)
	{
		const int ZERO = 0;
		int k, l, i, first_row_k, last_row_k;
		T mult, *diag_k;

		/* Solve Ly = Pb, store solution y in b */
		for (k = 0; k < n - 1; k++) 
		{
			l = p[k];
			mult = b[l];
			if (l != k) {
				b[l] = b[k];
				b[k] = mult;
			}
			diag_k = a[k] + smu;
			last_row_k = std::min(n - 1, k + ml);
			for (i = k + 1; i <= last_row_k; i++)
				b[i] += mult * diag_k[i - k];
		}

		/* Solve Ux = y, store solution x in b */
		for (k = n - 1; k >= 0; k--) 
		{
			diag_k = a[k] + smu;
			first_row_k = std::max(ZERO, k - smu);
			b[k] /= (*diag_k);
			mult = -b[k];
			for (i = first_row_k; i <= k - 1; i++)
				b[i] += mult*diag_k[i - k];
		}
		
		return 0;
	}
}
