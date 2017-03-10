/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
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

#ifndef OpenSMOKE_StdInclude
#define OpenSMOKE_StdInclude

	#include <fstream>
	#include <iostream>
	#include <iomanip>
	#include <sstream>
	#include <cmath>
	#include <stdarg.h>
	#include <assert.h>
	#include <vector>
	#include <map>
    #include <stdlib.h>

	#define __OPENSMOKE_VERSION__ "0.7.0"

	#define OPENSMOKE_LONG_DOUBLE 8
	#define OPENSMOKE_DOUBLE 8

	#define OPENSMOKE_INT 4
	#define OPENSMOKE_LONG_INT 4

	#define OPENSMOKE_FATAL_ERROR_EXIT -1

	#define OPENSMOKE_SUCCESSFULL_EXIT  0

	const int ONE_INT					=	1;
	const double ONE_DOUBLE				=	1.;
	const int ZERO_INT					=	0;
	const double ZERO_DOUBLE			=	0.;
	const double OPENSMOKE_BIG			=	1.e+30;
	const float  OPENSMOKE_BIG_FLOAT	=	3.402823466e+38F;
	const float  OPENSMOKE_TINY_FLOAT	=	1.175494351e-38F;
	const double OPENSMOKE_BIG_DOUBLE	=	1.7976931348623154e+308;
	const double OPENSMOKE_TINY_DOUBLE	=	2.2250738585072014e-308;

	#if OPENSMOKE_LONG_DOUBLE == OPENSMOKE_DOUBLE
		const long double OPENSMOKE_BIG_LONG_DOUBLE		=	1.7976931348623154e+308;
		const long double OPENSMOKE_TINY_LONG_DOUBLE	=	2.2250738585072014e-308;
	#else
		const long double OPENSMOKE_BIG_LONG_DOUBLE		=	1.1E+4932L;
		const long double OPENSMOKE_TINY_LONG_DOUBLE	=	1.1E-4932L;
	#endif



	/**
	* Format of input/output files
	*/
	enum OpenSMOKE_File_Format 
	{
		OPENSMOKE_BINARY_FILE, 
		OPENSMOKE_FORMATTED_FILE 
	};

	/**
	* This function is used to check for fatal errors. If the argument is false a fatal error message
	  is reported and additional information (for developing purposes) are shown (file name and line)
	*/
	#define CheckForFatalError(tag) \
	{ \
		if (tag == false) \
		{ \
			std::cout << __FILE__ << " (" << __LINE__ << ")" << std::endl; \
			std::cout << "Press enter to continue..." << std::endl; \
			getchar(); \
			exit(OPENSMOKE_FATAL_ERROR_EXIT); \
		} \
	} \

	namespace OpenSMOKE
	{
		/**
		* Types of equations for DAE systems
		*/
		enum EquationType
		{
			EQUATION_TYPE_ALGEBRAIC,
			EQUATION_TYPE_DIFFERENTIAL
		};

		/**
		* Dense solvers
		*/
		enum DenseSolverType
		{
			SOLVER_DENSE_NONE, 
			SOLVER_DENSE_EIGEN
		};

		/**
		* Dense decompositions
		*/
		enum DenseDecompositionType 
		{ 
			DENSE_DECOMPOSITION_FULL_PIVOTING_LU, 
			DENSE_DECOMPOSITION_PARTIAL_PIVOTING_LU 
		};


		/**
		* Sparse solvers
		*/
		enum SparseSolverType 
		{ 
			SOLVER_SPARSE_NONE,
			SOLVER_SPARSE_EIGEN_SPARSE_LU,
			SOLVER_SPARSE_EIGEN_BICGSTAB,
			SOLVER_SPARSE_EIGEN_GMRES,
			SOLVER_SPARSE_EIGEN_DGMRES,
			SOLVER_SPARSE_EIGEN_PARDISO,
			SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL,
			SOLVER_SPARSE_EIGEN_UMFPACK,
			SOLVER_SPARSE_EIGEN_LIS
		};

		/**
		* Sparse preconditioners
		*/
		enum SparsePreconditionerType 
		{ 
			PRECONDITIONER_SPARSE_DIAGONAL,
			PRECONDITIONER_SPARSE_ILUT
		};
	}

#endif
