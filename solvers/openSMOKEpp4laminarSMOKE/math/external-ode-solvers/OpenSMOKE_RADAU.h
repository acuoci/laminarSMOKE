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

#ifndef OpenSMOKE_RADAU_H
#define OpenSMOKE_RADAU_H

#include "math/external-ode-solvers/OpenSMOKE_OdeSystemSolver.h"

namespace OpenSMOKE
{
	enum OpenSMOKE_RADAU_SolverType { solver_radau, solver_radau5 };

	template <typename T>
	class OpenSMOKE_RADAU : public OpenSMOKE::OpenSMOKE_OdeSystemSolver<T>
	{
		public:

			OpenSMOKE_RADAU(T* odeSystem);

			void SetDimensions(const int n);
	
			void Solve(const double xf);	

			void Status() const;

			void SetRadauSolver()
			{
				iSolver_ = solver_radau;
			}

			void SetRadau5Solver()
			{
				iSolver_ = solver_radau5;
			}

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
			~OpenSMOKE_RADAU(void);

	private:

		void (*ptMassMatrix_)(int *n,double *am, int *lmas,int *rpar, int *ipar);

		OpenSMOKE_RADAU_SolverType iSolver_;

		// internal variables
		int lwork_;
		int liwork_;

		// internal variables
		double *rwork_;
		int    *iwork_;

		// Jacobian options
		//int mLower_;		// mLower_=N means full matrix	0<mLower_<N means number of non zero diagonals below the main diagonal
		//int mUpper_;		//                              0<=mUpper_<N means number of non zero diagonals above the main diagonal

		// Mass matrix
		int imas_;		// 0=mass matrix equal to the I matrix	1=mass matrix must be supplied by the user
		int mlmas_;		// mlmas=N means full matrix	0<=mlmas_<N means number of non zero diagonals below the main diagonal
		int mumas_;		//                              0<=mlmas_<N means number of non zero diagonals above the main diagonal
	
		int		*ipar_;
		double	*rpar_;
		int		 idid_;	//  1  COMPUTATION SUCCESSFUL
						//  2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
						// -1  INPUT IS NOT CONSISTENT
						// -2  LARGER NMAX IS NEEDED
						// -3  STEP SIZE BECOMES TOO SMALL
						// -4  MATRIX IS REPEATEDLY SINGULAR.

	private:

		void MemoryAllocation(const int n);
		void AnalyzeUserOptions();
	};

	#if defined(_WIN32) || defined(_WIN64) 
		extern "C" {	void RADAU(int *N,
						void FCN(int*,double*,double*,double*,double*,int*),
						double *X, double *Y, double *XEND, double *H,
						double *RTOL, double *ATOL, int *ITOL,
						void JAC(int*, double*, double*, double*, int*, double*, double*),
						int *IJAC, int *MLJAC, int *MUJAC,
						void MAS(int *n,double *am, int *lmas,int *rpar, int *ipar),
						int *IMAS, int *MLMAS, int *MUMAS,
						void SOLOUT(int*,double*,double*,double*,double*,int*,int*,double*,int*,int*),
						int *IOUT,
						double *WORK, int *LWORK,int *IWORK, int *LIWORK,
						double *RPAR, int *IPAR, int *IDID); 
					}

		extern "C"  {	void RADAU5(int *N,
									void FCN(int*,double*,double*,double*,double*,int*),
									double *X, double *Y, double *XEND, double *H,
									double *RTOL, double *ATOL, int *ITOL,
									void JAC(int*, double*, double*, double*, int*, double*, double*),
									int *IJAC, int *MLJAC, int *MUJAC,
									void MAS(int *n,double *am, int *lmas,int *rpar, int *ipar),
									int *IMAS, int *MLMAS, int *MUMAS,
									void SOLOUT(int*,double*,double*,double*,double*,int*,int*,double*,int*,int*),
									int *IOUT,
									double *WORK, int *LWORK,int *IWORK, int *LIWORK,
									double *RPAR, int *IPAR, int *IDID); 
				}

		extern "C"  { double CONTRA(int *I, double *S, double *CONT, int *LRC); }
		extern "C"  { double CONTR5(int *I, double *S, double *CONT, int *LRC); }
	#else
		extern "C" {	void radau_(int *N,
						void FCN(int*,double*,double*,double*,double*,int*),
						double *X, double *Y, double *XEND, double *H,
						double *RTOL, double *ATOL, int *ITOL,
						void JAC(int*, double*, double*, double*, int*, double*, double*),
						int *IJAC, int *MLJAC, int *MUJAC,
						void MAS(int *n,double *am, int *lmas,int *rpar, int *ipar),
						int *IMAS, int *MLMAS, int *MUMAS,
						void SOLOUT(int*,double*,double*,double*,double*,int*,int*,double*,int*,int*),
						int *IOUT,
						double *WORK, int *LWORK,int *IWORK, int *LIWORK,
						double *RPAR, int *IPAR, int *IDID); 
					}

		extern "C"  {	void radau5_(int *N,
									void FCN(int*,double*,double*,double*,double*,int*),
									double *X, double *Y, double *XEND, double *H,
									double *RTOL, double *ATOL, int *ITOL,
									void JAC(int*, double*, double*, double*, int*, double*, double*),
									int *IJAC, int *MLJAC, int *MUJAC,
									void MAS(int *n,double *am, int *lmas,int *rpar, int *ipar),
									int *IMAS, int *MLMAS, int *MUMAS,
									void SOLOUT(int*,double*,double*,double*,double*,int*,int*,double*,int*,int*),
									int *IOUT,
									double *WORK, int *LWORK,int *IWORK, int *LIWORK,
									double *RPAR, int *IPAR, int *IDID); 
				}

		extern "C"  { double contra_(int *I, double *S, double *CONT, int *LRC); }
		extern "C"  { double contr5_(int *I, double *S, double *CONT, int *LRC); }
	#endif

	double ccontra(int i, double s, double *cont, int *lrc);
	double ccontr5(int i, double s, double *cont, int *lrc);

}

#include "OpenSMOKE_RADAU.hpp"

#endif	// OpenSMOKE_RADAU_H

