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

#ifndef OpenSMOKE_DASPK_H
#define OpenSMOKE_DASPK_H

#include "math/external-ode-solvers/OpenSMOKE_OdeSystemSolver.h"

namespace OpenSMOKE
{
	template <typename T>
	class OpenSMOKE_DASPK : public OpenSMOKE::OpenSMOKE_OdeSystemSolver<T>
	{
		public:

			OpenSMOKE_DASPK(T* odeSystem);

			void SetDimensions(const int n);
	
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
			~OpenSMOKE_DASPK(void);

	private:

		int info[20];
		double* dy_;
		
		// internal variables
		int lrw_;
		int liw_;

		// internal variables
		double *rwork_;
		int    *iwork_;

		// internal variables
		double *rpar_;
		int    *ipar_;

	private:

		void MemoryAllocation(const int n);
		void AnalyzeUserOptions();
	};

	#if defined(_WIN32) || defined(_WIN64) 
		extern "C" {	void DDASPK (void RES(double*,double*,double*, double*, double*, int*, double*, int*), int *NEQ, double *T, double *Y, double *YPRIME, double *TOUT,
		                int *INFO, double *RTOL, double *ATOL, int *IDID, double *RWORK, int *LRW,
						int *IWORK, int *LIW, double *RPAR, int *IPAR, 
						void JAC(double*, double*, double*, double*, double*, double*, int*), 
						void PSOL(int*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, int*, double*, int*),
						void PRINT(double *T, double *Y)
								);
}
	#else
		extern "C" {	void ddaspk_ (void RES(double*,double*,double*, double*, double*, int*, double*, int*), int *NEQ, double *T, double *Y, double *YPRIME, double *TOUT,
		                int *INFO, double *RTOL, double *ATOL, int *IDID, double *RWORK, int *LRW,
						int *IWORK, int *LIW, double *RPAR, int *IPAR, 
						void JAC(double*, double*, double*, double*, double*, double*, int*), 
						void PSOL(int*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, int*, double*, int*),
						void PRINT(double *T, double *Y)
									);
				}
	#endif

}

#include "OpenSMOKE_DASPK.hpp"

#endif	// OpenSMOKE_DASPK_H
