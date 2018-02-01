/*----------------------------------------------------------------------*\
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
|	License                                                               |
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

#ifndef  OpenSMOKE_CVODE_Sundials_H
#define  OpenSMOKE_CVODE_Sundials_H

#include "math/external-ode-solvers/OpenSMOKE_OdeSystemSolver.h"

#include <cvode/cvode.h>						/* prototypes for CVODE fcts., consts. */
#include <sunmatrix/sunmatrix_dense.h>			/* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_lapackdense.h>	/* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>                 /* access to CVDls interface            */
#include <nvector/nvector_serial.h>				/* serial N_Vector types, fcts., macros */

#include <sunlinsol/sunlinsol_dense.h>		    /* prototype for CVDense */
#include <sunlinsol/sunlinsol_band.h>			/* prototype for CVBand */
#include <sunlinsol/sunlinsol_lapackdense.h>    /* prototype for CVDense */
#include <sunlinsol/sunlinsol_lapackband.h>     /* prototype for CVBand */
#include <sundials/sundials_dense.h>			/* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>			/* definition of type realtype */

namespace OpenSMOKE
{
	#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
	#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */
	 
	static int check_flag(void *flagvalue, char *funcname, int opt);

	template <typename T>
	class OpenSMOKE_CVODE_Sundials : public OpenSMOKE::OpenSMOKE_OdeSystemSolver<T>
	{
		public:

			OpenSMOKE_CVODE_Sundials(T* odeSystem);

			void SetDimensions(const int n);
			void SetLapackSolver(const bool flag);
	
			void Solve(const double xf);	

			void Status() const;

			std::string Tag() const;
			int GetNumberOfSteps() const;
			int GetNumberOfFunctionEvaluations() const;
			int GetNumberOfJacobianEvaluations() const;
			int GetNumberOfLUFactorizations() const;
			int GetNumberOfNonLinearIterations() const;
			int GetLastOrderUsed() const;
			int GetNumberOfConvergenceFailures() const;
			int GetNumberOfErrorTestFailures() const;
			double GetLastStepUsed() const;

			/**
			* Default destructor
			*/
			~OpenSMOKE_CVODE_Sundials(void);

	private:

		N_Vector y0Sundials_;
		N_Vector ySundials_;
		void *cvode_mem_;

		bool firstCall_; 
		bool iUseLapack_;

		SUNMatrix A;
		SUNLinearSolver LS;

	private:

		void MemoryAllocation(const int n);
		void AnalyzeUserOptions();
	};

	/*
	 * Check function return value...
	 *   opt == 0 means SUNDIALS function allocates memory so check if
	 *            returned NULL pointer
	 *   opt == 1 means SUNDIALS function returns a flag so check if
	 *            flag >= 0
	 *   opt == 2 means function allocates memory so check if returned
	 *            NULL pointer 
	 */

	static int check_flag(void *flagvalue, const std::string funcname, int opt)
	{
		int *errflag;

		/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
		if (opt == 0 && flagvalue == NULL) 
		{
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname.c_str());
			return(1); 
		}

		/* Check if flag < 0 */
		else if (opt == 1) 
		{
			errflag = (int *) flagvalue;
			if (*errflag < 0) 
			{
				fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname.c_str(), *errflag);
				return(1);
			}
		}

		/* Check if function returned NULL pointer - no memory allocated */
		else if (opt == 2 && flagvalue == NULL) 
		{
			fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname.c_str());
			return(1); 
		}

		return(0);
	}
}

#include "OpenSMOKE_CVODE_Sundials.hpp"

#endif	// OpenSMOKE_CVODE_Sundials_H
