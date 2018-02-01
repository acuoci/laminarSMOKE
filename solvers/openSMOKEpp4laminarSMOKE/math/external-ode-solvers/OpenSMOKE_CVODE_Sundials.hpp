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

#include <typeinfo>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iomanip>

namespace OpenSMOKE
{
	template<typename T>
	OpenSMOKE_CVODE_Sundials<T>::OpenSMOKE_CVODE_Sundials(T* odeSystem)
	{
		firstCall_ = true;
		this->iUseLapack_ = true;

		this->SetDefaultValues();

		y0Sundials_ = NULL;
		cvode_mem_ = NULL;
		A = NULL;
		LS = NULL;

		this->odeSystem_ = odeSystem;
	}

	template <typename T>
	void OpenSMOKE_CVODE_Sundials<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_CVODE_Sundials<T>::SetLapackSolver(const bool flag)
	{
		this->iUseLapack_ = flag;
	}

	template<typename T>
	void OpenSMOKE_CVODE_Sundials<T>::MemoryAllocation(const int n)
	{	
		this->n_		=	n;					// Number of equations

		this->y0_ 	= new double[this->n_];
		this->y_ 	= new double[this->n_];
 
		y0Sundials_ = N_VNew_Serial(this->n_);
		if (check_flag((void *)y0Sundials_, std::string("N_VNew_Serial"), 0))	exit(-1);

		ySundials_ = N_VNew_Serial(this->n_);
		if (check_flag((void *)ySundials_, std::string("N_VNew_Serial"), 0))	exit(-1);
	}

	template<typename T>
	void OpenSMOKE_CVODE_Sundials<T>::AnalyzeUserOptions()
	{
		/*
		TODO
		if (this->iSetMaximumNumberOfPrints_) iwork_[6] = this->maximumNumberOfPrints_;
		if (this->iSetMaximumOrder_)			iwork_[8] = this->maximumOrder_;
		*/

		if (this->iSetMaximumNumberOfSteps_)
		{
			int flag = CVodeSetMaxNumSteps(cvode_mem_, this->maximumNumberOfSteps_);
			if (check_flag(&flag, std::string("CVodeSetMaxNumSteps"), 1)) exit(-1);
		}

		if (this->iSetMaximumStep_)
		{
			int flag = CVodeSetMaxStep(cvode_mem_, this->maximumStep_);
			if (check_flag(&flag, std::string("CVodeSetMaxStep"), 1)) exit(-1);
		}

		if (this->iSetMinimumStep_)
		{
			int flag = CVodeSetMinStep(cvode_mem_, this->minimumStep_);
			if (check_flag(&flag, std::string("CVodeSetMinStep"), 1)) exit(-1);
		}

		if (this->iSetFirstStep_)
		{
			int flag = CVodeSetInitStep(cvode_mem_, this->firstStep_);
			if (check_flag(&flag, std::string("CVodeSetInitStep"), 1)) exit(-1);
		}
	}

	template<typename T>
	void OpenSMOKE_CVODE_Sundials<T>::Solve(const double xend)
	{

		int flag;

		this->x_ = this->x0_;
		this->xend_ = xend;

		for(int i=0;i<this->n_;i++)
			NV_Ith_S(y0Sundials_,i) = this->y0_[i];

		if (firstCall_ == true)
		{
			firstCall_ = false;

			/* Call CVodeCreate to create the solver memory and specify the 
			* Backward Differentiation Formula and the use of a Newton iteration */
			cvode_mem_ = CVodeCreate(CV_BDF, CV_NEWTON);
			if (check_flag((void *)cvode_mem_, std::string("CVodeCreate"), 0)) exit(-1);

			/* Call CVodeInit to initialize the integrator memory and specify the
			* user's right hand side function in y'=f(t,y), the inital time t0, and
			* the initial dependent variable vector y0Sundials_. */
			flag = CVodeInit(cvode_mem_, this->odeSystem_->GetSystemFunctionsStatic, this->odeSystem_->GetWriteFunctionStatic, this->x0_, y0Sundials_);
			if (check_flag(&flag, std::string("CVodeInit"), 1)) exit(-1);

			/* Call CVodeSVtolerances to specify the scalar relative tolerance
			* and vector absolute tolerances */
			flag = CVodeSStolerances(cvode_mem_, this->relTolerance_[0], this->absTolerance_[0]);
			if (check_flag(&flag, std::string("CVodeSVtolerances"), 1)) exit(-1);

			/* Call Solver */
			if (this->iUseLapack_ == false)
			{
				if (this->mUpper_ == 0 && this->mLower_ == 0)
				{
					// std::cout << "CVODE Solver: Dense Jacobian (without Lapack)..." << std::endl;

					/* Create dense SUNMatrix for use in linear solves */
					A = SUNDenseMatrix(this->n_, this->n_);
					if (check_flag((void *)A, std::string("SUNDenseMatrix"), 0)) exit(-1);
					
					/* Create SUNDenseLinearSolver solver object for use by CVode */
					LS = SUNDenseLinearSolver(ySundials_, A);
					if (check_flag((void *)LS, std::string("SUNDenseLinearSolver"), 0)) exit(-1);
				}
				else
				{
					// std::cout << "CVODE Solver: Band Jacobian (without Lapack)..." << std::endl;

					/* Create banded SUNMatrix for use in linear solves -- since this will be factored,
					set the storage bandwidth to be the sum of upper and lower bandwidths */
					A = SUNBandMatrix(this->n_, this->mUpper_, this->mLower_, (this->mUpper_+this->mLower_) );
					if (check_flag((void *)A, std::string("SUNBandMatrix"), 0)) exit(-1);

					/* Create banded SUNLinearSolver object for use by CVode */
					LS = SUNBandLinearSolver(ySundials_, A);
					if (check_flag((void *)LS, std::string("SUNBandLinearSolver"), 0)) exit(-1);
				}
			}
			else
			{
				if (this->mUpper_ == 0 && this->mLower_ == 0)
				{
					// std::cout << "CVODE Solver: Dense Jacobian (with Lapack)..." << std::endl;

					/* Create dense SUNMatrix for use in linear solves */
					A = SUNDenseMatrix(this->n_, this->n_);
					if (check_flag((void *)A, std::string("SUNDenseMatrix"), 0)) exit(-1);

					/* Create SUNLapackDense solver object for use by CVode */
					LS = SUNLapackDense(ySundials_, A);
					if (check_flag((void *)LS, std::string("SUNLapackDense"), 0)) exit(-1);
				}
				else
				{
					// std::cout << "CVODE Solver: Band Jacobian (with Lapack)..." << std::endl;

					/* Create banded SUNMatrix for use in linear solves -- since this will be factored,
					set the storage bandwidth to be the sum of upper and lower bandwidths */
					A = SUNBandMatrix(this->n_, this->mUpper_, this->mLower_, (this->mUpper_ + this->mLower_));
					if (check_flag((void *)A, std::string("SUNBandMatrix"), 0)) exit(-1);

					/* Create banded SUNLapackBand solver object for use by CVode */
					LS = SUNLapackBand(ySundials_, A);
					if (check_flag((void *)LS, std::string("SUNLapackBand"), 0)) exit(-1);
				}
			}

			/* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
			flag = CVDlsSetLinearSolver(cvode_mem_, LS, A);
			if (check_flag(&flag, std::string("CVDlsSetLinearSolver"), 1)) exit(-1);
		}
		else
		{
			flag = CVodeReInit(cvode_mem_, this->x0_, y0Sundials_);
			if (check_flag(&flag, std::string("CVodeReInit"), 1)) exit(-1);
		}

		AnalyzeUserOptions();

		/* Solving */
		this->tStart_ =  this->GetClockTime();
		flag = CVode(cvode_mem_, this->xend_, ySundials_, &this->x_, CV_NORMAL);
		this->tEnd_ =  this->GetClockTime();

		this->x0_ = this->x_;
		for(int i=0;i<this->n_;i++)
			NV_Ith_S(y0Sundials_,i) = NV_Ith_S(ySundials_,i);
		for(int i=0;i<this->n_;i++)
			this->y_[i] = NV_Ith_S(ySundials_,i);
	}

	template<typename T>
	void OpenSMOKE_CVODE_Sundials<T>::Status() const
	{
		int flag;
		long int nst, nfe, nsetups, netf, nni, ncfn, nje, nfeLS, nge;
		int qcurrent, qlast;
		double hcurrent, hlast;

		flag = CVodeGetNumSteps(cvode_mem_, &nst);
		check_flag(&flag, std::string("CVodeGetNumSteps"), 1);
		flag = CVDlsGetNumJacEvals(cvode_mem_, &nje);
		check_flag(&flag, std::string("CVDlsGetNumJacEvals"), 1);
		flag = CVodeGetNumRhsEvals(cvode_mem_, &nfe);
		check_flag(&flag, std::string("CVodeGetNumRhsEvals"), 1);

		flag = CVodeGetNumLinSolvSetups(cvode_mem_, &nsetups);
		check_flag(&flag, std::string("CVodeGetNumLinSolvSetups"), 1);
		flag = CVodeGetNumErrTestFails(cvode_mem_, &netf);
		check_flag(&flag, std::string("CVodeGetNumErrTestFails"), 1);
		flag = CVodeGetNumNonlinSolvIters(cvode_mem_, &nni);
		check_flag(&flag, std::string("CVodeGetNumNonlinSolvIters"), 1);
		flag = CVodeGetNumNonlinSolvConvFails(cvode_mem_, &ncfn);
		check_flag(&flag, std::string("CVodeGetNumNonlinSolvConvFails"), 1);
		flag = CVodeGetNumGEvals(cvode_mem_, &nge);
		check_flag(&flag, std::string("CVodeGetNumGEvals"), 1);

		flag = CVDlsGetNumRhsEvals(cvode_mem_, &nfeLS);
		check_flag(&flag, std::string("CVDlsGetNumRhsEvals"), 1);

		flag = CVodeGetLastOrder(cvode_mem_, &qlast);
		check_flag(&flag, std::string("CVodeGetLastOrder"), 1);
		flag = CVodeGetCurrentOrder(cvode_mem_, &qcurrent);
		check_flag(&flag, std::string("CVodeGetCurrentOrder"), 1);
		flag = CVodeGetLastStep(cvode_mem_, &hlast);
		check_flag(&flag, std::string("CVodeGetLastStep"), 1);
		flag = CVodeGetCurrentStep(cvode_mem_, &hcurrent);
		check_flag(&flag, std::string("CVodeGetCurrentStep"), 1);


		std::cout << "CVODE Sundials Status" << std::endl;
		std::cout << " * Absolute tolerance:              " << this->absTolerance_[0]   << std::endl;	// Absolute tolerance
		std::cout << " * Relative tolerance:              " << this->relTolerance_[0]   << std::endl;	// Relative tolerance
		std::cout << " * Number of steps:                 " << nst << std::endl;	// Number of steps taken for the problem so far 
		std::cout << " * Number of function evaluations:  " << nfe << std::endl;	// Number of f evaluations for the problem so far.
		std::cout << " * Number of Jacobians:             " << nje << std::endl;	// Number of Jacobian evaluations (and of matrix LU decompositions) for the problem so far.
		std::cout << " * Last step:                       " << hlast << std::endl;	
		std::cout << " * Next  step:                      " << hcurrent << std::endl;	
		std::cout << " * Last order:                      " << qlast << std::endl;	
		std::cout << " * Next order:                      " << qcurrent << std::endl;
	}

	template<typename T>
	std::string OpenSMOKE_CVODE_Sundials<T>::Tag() const
	{
		return "CVODE";
	}

	template<typename T>
	int OpenSMOKE_CVODE_Sundials<T>::GetNumberOfSteps() const
	{
		long int nst;
		int flag = CVodeGetNumSteps(cvode_mem_, &nst);
		check_flag(&flag, std::string("CVodeGetNumSteps"), 1);

		return int(nst);
	}

	template<typename T>
	int OpenSMOKE_CVODE_Sundials<T>::GetNumberOfFunctionEvaluations() const
	{
		long int nfe;
		int flag = CVodeGetNumRhsEvals(cvode_mem_, &nfe);
		check_flag(&flag, std::string("CVodeGetNumRhsEvals"), 1);

		return int(nfe);
	}

	template<typename T>
	int OpenSMOKE_CVODE_Sundials<T>::GetNumberOfJacobianEvaluations() const
	{
		long int nje;
		int flag = CVDlsGetNumJacEvals(cvode_mem_, &nje);
		check_flag(&flag, std::string("CVDlsGetNumJacEvals"), 1);

		return int(nje);
	}

	template<typename T>
	int OpenSMOKE_CVODE_Sundials<T>::GetNumberOfLUFactorizations() const
	{
		return int(0);
	}

	template <typename T>
	int OpenSMOKE_CVODE_Sundials<T>::GetNumberOfNonLinearIterations() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_CVODE_Sundials<T>::GetLastOrderUsed() const
	{
		return 0;
	}

	template <typename T>
	double OpenSMOKE_CVODE_Sundials<T>::GetLastStepUsed() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_CVODE_Sundials<T>::GetNumberOfConvergenceFailures() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_CVODE_Sundials<T>::GetNumberOfErrorTestFailures() const
	{
		return 0;
	}

	template<typename T>
	OpenSMOKE_CVODE_Sundials<T>::~OpenSMOKE_CVODE_Sundials(void)
	{
		/* Free vectors */
		N_VDestroy_Serial(y0Sundials_);
		N_VDestroy_Serial(ySundials_);

		/* Free integrator memory */
		CVodeFree(&cvode_mem_);
		SUNLinSolFree(LS);
		SUNMatDestroy(A);

		delete[] this->y0_;
		delete[] this->y_;
	}
}
