/***************************************************************************
 *   Copyright (C) 2012 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef OpenSMOKE_OpenSMOKEBandMatrix_Hpp
#define OpenSMOKE_OpenSMOKEBandMatrix_Hpp

#include "OpenSMOKEStdInclude.h"
#include "OpenSMOKEBaseClass.h"
#include "OpenSMOKEUtilities.h"
#include "OpenSMOKEFunctions.h"

namespace OpenSMOKE
{
	//!  A class for banded or tridiagonal-block (as a special case) matrices
	/*!
		 A class for banded or tridiagonal-block (as a special case) matrices
	*/

	template<typename T>
	class OpenSMOKEBandMatrix
	{
	public:

		/**
		*@brief Constructor for a tridiagonal-block matrix
		*@param dimBlock block dimension
		*/
		OpenSMOKEBandMatrix(const int nEquations, const int dimBlock);

		/**
		*@brief Constructor for a band matrix
		*@param nUpper upper band size
		*@param nLower lower band size
		*/
		OpenSMOKEBandMatrix(const int nEquations, const int nUpper, const int nLower);

		/**
		*@brief Returns true if the matrix is tridiagonal block
		*/
		bool isTriagonalBlock() const { return isTridiagonalBlock_; }

		/**
		*@brief Returns the upper bandwidth
		*/
		inline int nUpper() const { return mu; }

		/**
		*@brief Returns the lower bandwidth
		*/
		inline int nLower() const { return ml; }

		/**
		*@brief Sets all the coefficients equal to zero
		*/
		void SetToZero();

		/**
		*@brief Add the identity matrix
		*/
		void AddIdentity();

		/**
		*@brief Adds the specified diagonal matrix
		*@param d the diagonal vector to be added
		*/
		void AddDiagonal(const T *d);

		/**
		*@brief Destroys the matrix
		*/
		void DestroyMat();

		/**
		*@brief Copies the current matrix in the B banded matrix (sizesmust be consistent)
		*@param B the matrix where to copy
		*/
		void CopyTo(OpenSMOKEBandMatrix<T>* B);

		/**
		*@brief Scales all the elements of the matrix by the same scalar
		*@param c scalar
		*/
		void Scale(const double c);

		/**
		*@brief Scales all the elements of the matrix by the two different scalars, according to the type of equation
		*@param c_differential scalar to be used by the differential equations
		*@param c_algebraic scalar to be used by the algebraic equations
		*@param index type of reaction (0=algebraic, 1=differential)
		*/
		void Scale(const double c_differential, const double c_algebraic, const int *index);	// TOIMPROVE (slow)

		/**
		*@brief Multiplies the current matrix time a vector: y = A *x 
		*@param x vector to be multiplied
		*@param y vector where to put the result
		*/
		void Product(const T *x, T *y);

		/**
		*@brief Multiplies the transpose of current matrix time a vector: y = A' * x
		*@param x vector to be multiplied
		*@param y vector where to put the result
		*/
		void TProduct(const T *x, T *y);

		/**
		*@brief LU factorization
		*/
		int Factorize();

		/**
		*@brief Solves the linear system (only after LU factorization)
		*@param b rhs of linear system
		*/
		int Solve(T *b);

		/**
		*@brief Solves the linear system (only after LU factorization)
		*@param nrhs number of right hand sides
		*@param b rhs of linear system
		*/
		int Solve(const int nrhs, T *b);

		/**
		*@brief Factorizes and solves the linear system - Please, do not use it: it is still under testing
		*@param b rhs of linear system
		*/
		int FactorizeAndSolve(T *b);

		/**
		*@brief Prints the matrix on the screen (for diagnostic purposes)
		*@param out output stream
		*/
		void Print(std::ostream& out);

	public:

		int s_mu;	//!< storage upper bandwidth, mu <= s_mu <= N - 1.
			        //!< The dgbtrf routine writes the LU factors into the storage for A.
					//!< The upper triangular factor U, however, may have an upper bandwidth as big as MIN(N - 1, mu + ml) because of
			        //!< partial pivoting.The s_mu field holds the upper bandwidth allocated for A.

		T **cols;   //!< array of pointers.cols[j] points to the first element of the j - th column of the matrix in the array data.
		T *data;	//!< pointer to a contiguous block of double variables
		
	private:

		int M;			//!< number of rows
		int N;			//!< number of columns
		int mu;		//!< upper bandwidth, 0 <= mu <= min(M, N)
		int ml;		//!< lower bandwidth, 0 <= ml <= min(M, N)
		int ldata;		//!< length of the data array = ldim*(s_mu + ml + 1)
		int ldim;		//!< leading dimension(ldim >= s_mu)
		int *p;			//!< permutation vector

		bool isTridiagonalBlock_;
		
	public:

		template<typename TT>
		friend OpenSMOKEBandMatrix<TT>* NewBandMat(const int N, const int mu, const int ml, const int smu);

		template<typename TT>
		friend void bandCopy(TT **a, TT **b, const int n, const int a_smu, const int b_smu, const int copymu, const int copyml);

		template<typename TT>
		friend void bandScale(const TT c, TT **a, const int n, const int mu, const int ml, const int smu);
		
		template<typename TT>
		friend void bandAddIdentity(TT **a, const int n, const int smu);

		template<typename TT>
		friend void bandMatVec(TT **a, const TT *x, TT *y, const int n, const int mu, const int ml, const int smu);

		template<typename TT>
		friend void bandMatTransposeVec(TT **a, const TT *x, TT *y, const int n, const int mu, const int ml, const int smu);

		template<typename TT>
		friend int bandGBTRF(TT **a, const int n, const int mu, const int ml, const int smu, int *p);

		template<typename TT>
		friend void bandGBTRS(TT **a, const int n, const int smu, const int ml, const int *p, TT *b);
	};

	typedef OpenSMOKEBandMatrix<float>		OpenSMOKEBandMatrixFloat;
	typedef OpenSMOKEBandMatrix<double>		OpenSMOKEBandMatrixDouble;
}

#include "OpenSMOKEBandMatrix.hpp"

#endif	// OpenSMOKE_OpenSMOKEBandMatrix_Hpp

