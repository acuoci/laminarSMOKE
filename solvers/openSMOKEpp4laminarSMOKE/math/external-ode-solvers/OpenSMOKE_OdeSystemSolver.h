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

#ifndef  OpenSMOKE_OdeSystemSolver_H
#define  OpenSMOKE_OdeSystemSolver_H

#include "math/external-ode-solvers/OpenSMOKE_OdeSystemObject.h"

namespace OpenSMOKE
{
	template <typename T>
	class OpenSMOKE_OdeSystemSolver
	{
	
		public:

			virtual void Solve(const double xf) = 0;	

			virtual void Status() const = 0;

			void SetInitialValues(const double x0, const double* y0)
			{
				memcpy(y0_, y0, n_*sizeof(double));
				x0_ = x0;
			}
	
			void SetAnalyticalJacobian(const bool flag)
			{
				if (flag == true)			iJacobian_ = 1;
				else if (flag == false)		iJacobian_ = 0;
			}

			void SetAbsoluteTolerance(const double tolerance)
			{
				absTolerance_[0] = tolerance;
			}
	
			void SetRelativeTolerance(const double tolerance)
			{
				relTolerance_[0] = tolerance;
			}

			void Solution(double* solution) const
			{
				memcpy(solution, y_, n_*sizeof(double));
			}

			double LastPoint() const 
			{ 
				return x_;
			}

			void SetMaximumNumberOfSteps(const int value)
			{
				iSetMaximumNumberOfSteps_ = true;
				maximumNumberOfSteps_ = value;
			}

			void SetMaximumOrder(const int value)
			{
				iSetMaximumOrder_ = true;
				maximumOrder_ = value;
			}

			void SetMaximumNumberOfPrints(const int value)
			{
				iSetMaximumNumberOfPrints_ = true;
				maximumNumberOfPrints_ = value;
			}

			void SetFirstStep(const double value)
			{
				iSetFirstStep_ = true;
				firstStep_ = value;
			}

			void SetMaximumStep(const double value)
			{
				iSetMaximumStep_ = true;
				maximumStep_ = value;
			}

			void SetMinimumStep(const double value)
			{
				iSetMinimumStep_ = true;
				minimumStep_ = value;
			}

			void SetVerboseOutput(const bool verbose_output)
			{
				verbose_output_ = verbose_output;
			}

			double GetCPUTime() const
			{
				return tEnd_ - tStart_;
			}

			void SetZeroConstraints(int *positions)
			{
			}

			void SetTridiagonalBlock(const int blockDimension)
			{
				mUpper_ = 2 * blockDimension - 1;
				mLower_ = 2 * blockDimension - 1;
			}

			void SetUpperHalfBandWidth(const int value)
			{
				mUpper_ = value;
			}

			void SetLowerHalfBandWidth(const int value)
			{
				mLower_ = value;
			}

			virtual std::string Tag() const						= 0;
			virtual int GetNumberOfSteps() const				= 0;
			virtual int GetNumberOfFunctionEvaluations() const	= 0;
			virtual int GetNumberOfJacobianEvaluations() const	= 0;
			virtual int GetNumberOfLUFactorizations() const		= 0;
			virtual int GetNumberOfNonLinearIterations() const  = 0;
			virtual int GetLastOrderUsed() const                = 0;
			virtual int GetNumberOfConvergenceFailures() const  = 0;
			virtual int GetNumberOfErrorTestFailures() const    = 0;
			virtual double GetLastStepUsed() const              = 0;

			//virtual ~OpenSMOKE_OdeSystemSolver_BaseClass(void)  = 0;

	protected:

		// dimension of the system
		int  n_;
	
		// compute the jacobian analytically (0=numerically, 1=analitically)
		int iJacobian_;
		
		// Tolerances
		int iTolerance_;		// 1=scalars	2=vectors
		double *relTolerance_;	// relative tolerance
		double *absTolerance_;	// absolute tolerance
		
		// initial value	
		double *y0_;	
		double *y_;

		// steps
		double x_;
		double x0_;
		double xend_;

		// bands
		int mUpper_;
		int mLower_;

		// output
		int iOutput_;

		// CPU Time
		double tStart_;
		double tEnd_;


		bool iSetMaximumNumberOfSteps_;
		int maximumNumberOfSteps_ ;

		bool iSetMaximumOrder_;
		int maximumOrder_;

		bool iSetMaximumNumberOfPrints_;
		int maximumNumberOfPrints_;

		bool iSetFirstStep_;
		double firstStep_;

		bool iSetMaximumStep_;
		double maximumStep_;

		bool iSetMinimumStep_;
		double minimumStep_;

		bool verbose_output_;

		T* odeSystem_;


	protected:

		virtual void MemoryAllocation(const int n) = 0;

		virtual void AnalyzeUserOptions() = 0;
		
		void SetDefaultValues()
		{
			mLower_		=	0;	
			mUpper_		=	0;	

			iOutput_	=	1;	// output function is enabled
			iJacobian_	=	0;	// no analytical jacobian	

			iTolerance_	=	1;	// scalar
			relTolerance_   = new double[1];
			absTolerance_   = new double[1];
			relTolerance_[0]	=	1.0e-7;
			absTolerance_[0]	=	1.e-3*relTolerance_ [0];

			iSetMaximumOrder_          = false;
			iSetMaximumNumberOfSteps_  = false;
			iSetMaximumNumberOfPrints_ = false;

			iSetFirstStep_   = false;
			iSetMaximumStep_ = false;
			iSetMinimumStep_ = false;

			maximumNumberOfPrints_ = 0;
			maximumNumberOfSteps_ = 0;
			maximumOrder_ = 0;
			maximumStep_ = 0;
			minimumStep_ = 0;
			firstStep_ = 1.e-6;

			verbose_output_ = true;

			tStart_ = 0.;
			tEnd_	= 0.;
		}

		double GetClockTime() const
		{
			return static_cast<double>(std::clock())/CLOCKS_PER_SEC;
		}
	};

}

#endif	// OpenSMOKE_OdeSystemSolver_H
