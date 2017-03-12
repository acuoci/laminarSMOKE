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

#include "math/OpenSMOKEFunctions.h"
#include <typeinfo>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iomanip>

namespace OpenSMOKE
{
	template <typename T>
	OpenSMOKE_DLSODA<T>::OpenSMOKE_DLSODA(T* odeSystem)
	{	
		this->SetDefaultValues();

		itask_ = 4;

		istate_ = 1;		// This is the first call for a problem

		iOptions_ = 0;		// Flag indicating whether optional inputs are used (dafault = no)

		jt_ = 2;			// means an internally generated (difference quotient) full Jacobian (using NEQ extra calls to F per df/dy value).

		lrw_ = 0;
		liw_ = 0;

		this->odeSystem_ = odeSystem;
	}

	template <typename T>
	void OpenSMOKE_DLSODA<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_DLSODA<T>::MemoryAllocation(const int n)
	{
		istate_ = 1;		// This is the first call for a problem

		this->n_		=	n;					// Number of equations

		int lrn, lrs;
		//if (maxOrder == default)
		{
			lrn = 20+16*this->n_;
			if(jt_==1 || jt_==2)
				lrs = 22+9*this->n_+this->n_*this->n_;
			else if(jt_==4 || jt_==5)
				lrs = 22+10*this->n_+(2*this->mLower_+this->mUpper_)*this->n_;
		}
	/*	else
		{
			int nyh = this->n_;
			int maxordn = 12;
			int maxords = 5;
			int lmat = 0;
			if(jt_==1 || jt_==2)
				lmat = this->n_*this->n_+2;
			else if(jt_==4 || jt_==5)
				lmat = (2*this->mLower_+this->mUpper_+1)*this->n_+2;

			lrn = 20 + nyh*(maxordn+1) + 3*this->n_;
			lrs = 20 + nyh*(maxords+1) + 3*this->n_ + lmat; 
		}
		*/
		lrw_ = std::max(lrn, lrs);
		liw_ = 20+this->n_;

		
		//delete[]	this->y0_;
		//delete[]	this->y_;
		//delete[]	rwork_;
		//delete[]	iwork_;
		this->y0_ 	= new double[this->n_];
		this->y_ 	= new double[this->n_];
		rwork_ 	= new double[lrw_];
		iwork_	= new int[liw_];
		 
		// Default values
		for(int i=0; i<lrw_; i++)
			rwork_[i]	=	0.0;
		for(int i=0; i<liw_; i++)
			iwork_[i]	=	0;


		// Additional options
		iOptions_ = 1;

		iwork_[4] = 0;		// flag to generate extra printing at method switches. IXPR = 0 means no extra printing (the default).  IXPR = 1 means print data on each switch.
		iwork_[5] = 0;		// Maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
		iwork_[6] = 0;		// Maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size).  This must be positive to
							// result in a nondefault value.  The default value is 10.
		iwork_[7] = 0;		// the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
		iwork_[8] = 0;		// the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.

		rwork_[4] = 0.;		// Step size to be attempted on the first step. The default value is determined by the solver.
		rwork_[5] = 0.;		// Maximum absolute step size allowed. The default value is infinite.
		rwork_[6] = 0.;		// Minimum absolute step size allowed.  The default value is 0.  (This lower bound is not enforced on the final step before reaching TCRIT when ITASK = 4 or 5.)
	}

	template <typename T>
	void OpenSMOKE_DLSODA<T>::AnalyzeUserOptions()
	{
		if (this->iSetMaximumNumberOfSteps_ == true)
		{
			if (this->maximumNumberOfSteps_ >= 0)	iwork_[4] = this->maximumNumberOfSteps_;
			else									iwork_[4] = 0;
		}

		if (this->iSetMaximumNumberOfPrints_ == true)
		{
			if (this->maximumNumberOfPrints_ >= 0)	iwork_[5] = this->maximumNumberOfPrints_;
			else									iwork_[5] = 0;
		}

		if (this->iSetMaximumOrder_ == true)
		{
			if (this->maximumOrder_ >= 0) iwork_[8] = this->maximumOrder_;
			else                          iwork_[8] = 0;
		}


		if (this->iSetFirstStep_)
		{
			if (this->firstStep_ >= 0) rwork_[4] = this->firstStep_;
			else                       rwork_[4] = 0.;
		}

		if (this->iSetMaximumStep_)
		{
			if (this->maximumStep_ >= 0) rwork_[5] = this->maximumStep_;
			else                         rwork_[5] = 0.;
		}

		if (this->iSetMinimumStep_)
		{
			if (this->minimumStep_ >= 0) rwork_[6] = this->minimumStep_;
			else                         rwork_[6] = 0.;
		}
	}

	template <typename T>
	void OpenSMOKE_DLSODA<T>::Solve(const double xend)
	{
		istate_ = 1;		// This is the first call for a problem
		rwork_[0] = xend;

		AnalyzeUserOptions();

		this->x_ = this->x0_;
		this->xend_ = xend;
		memcpy(this->y_, this->y0_, this->n_*sizeof(double));
	
		this->tStart_ =  this->GetClockTime();

		#if defined(_WIN32) || defined(_WIN64) 
		DLSODA(	this->odeSystem_->GetSystemFunctionsStatic, &this->n_, this->y_, &this->x_, &this->xend_, 
				&this->iTolerance_, this->relTolerance_, this->absTolerance_,
				&itask_, &istate_, &iOptions_, rwork_, &lrw_, iwork_, &liw_, 
				this->odeSystem_->GetAnalyticalJacobianStatic, &jt_,
				this->odeSystem_->GetWriteFunctionStatic);
		#else
		dlsoda_(	this->odeSystem_->GetSystemFunctionsStatic, &this->n_, this->y_, &this->x_, &this->xend_, 
				&this->iTolerance_, this->relTolerance_, this->absTolerance_,
				&itask_, &istate_, &iOptions_, rwork_, &lrw_, iwork_, &liw_, 
				this->odeSystem_->GetAnalyticalJacobianStatic, &jt_,
				this->odeSystem_->GetWriteFunctionStatic);
		#endif
		this->tEnd_ =  this->GetClockTime();

		this->x0_ = this->x_;
		memcpy(this->y0_, this->y_, this->n_*sizeof(double));

	}

	template <typename T>
	void OpenSMOKE_DLSODA<T>::Status() const
	{
		std::cout << "DLSODA Status" << std::endl;
		std::cout << " * Absolute tolerance:              " << this->absTolerance_[0]   << std::endl;	// Absolute tolerance
		std::cout << " * Relative tolerance:              " << this->relTolerance_[0]   << std::endl;	// Relative tolerance
		std::cout << " * Number of steps:                 " << iwork_[10] << std::endl;	// Number of steps taken for the problem so far 
		std::cout << " * Number of function evaluations:  " << iwork_[11] << std::endl;	// Number of f evaluations for the problem so far.
		std::cout << " * Number of Jacobian evaluations:  " << iwork_[12] << std::endl;	// Number of Jacobian evaluations (and of matrix LU decompositions) for the problem so far.
		std::cout << " * Last order used:                 " << iwork_[13] << std::endl;	// Method order last used (successfully).
		std::cout << " * Next order to be used:           " << iwork_[14] << std::endl;	// The order to be attempted on the next step.
		std::cout << " * Last time step:                  " << rwork_[10] << std::endl;	// the step size in t last used (successfully)
		std::cout << " * Next time step:                  " << rwork_[11] << std::endl;	// the step size to be attempted on the next step
		std::cout << " * Current time value:              " << rwork_[12] << std::endl;	// the current value of the independent variable
		std::cout << " * Last method (1=NS, 2=S):         " << iwork_[18] << std::endl;	// the method indicator for the last successful step: 1 means Adams (nonstiff), 2 means BDF (stiff)
		std::cout << " * Current method (1=NS, 2=S):      " << iwork_[19] << std::endl;	// the method indicator for the current time: 1 means Adams (nonstiff), 2 means BDF (stiff)
	}

	template <typename T>
	std::string OpenSMOKE_DLSODA<T>::Tag() const
	{
		return "DLSODA";
	}

	template <typename T>
	int OpenSMOKE_DLSODA<T>::GetNumberOfSteps() const
	{
		return iwork_[10];
	}

	template <typename T>
	int OpenSMOKE_DLSODA<T>::GetNumberOfFunctionEvaluations() const
	{
		return iwork_[11];
	}

	template <typename T>
	int OpenSMOKE_DLSODA<T>::GetNumberOfJacobianEvaluations() const
	{
		return iwork_[12];
	}

	template <typename T>
	int OpenSMOKE_DLSODA<T>::GetNumberOfLUFactorizations() const
	{
		return iwork_[12];
	}

	template <typename T>
	int OpenSMOKE_DLSODA<T>::GetNumberOfNonLinearIterations() const
	{
		return -1;
	}

	template <typename T>
	int OpenSMOKE_DLSODA<T>::GetLastOrderUsed() const
	{
		return iwork_[13];
	}

	template <typename T>
	double OpenSMOKE_DLSODA<T>::GetLastStepUsed() const
	{
		return rwork_[10];
	}

	template <typename T>
	int OpenSMOKE_DLSODA<T>::GetNumberOfConvergenceFailures() const
	{
		return -1;
	}

	template <typename T>
	int OpenSMOKE_DLSODA<T>::GetNumberOfErrorTestFailures() const
	{
		return -1;
	}

	template <typename T>
	OpenSMOKE_DLSODA<T>::~OpenSMOKE_DLSODA(void)
	{
		delete[] this->y0_;
		delete[] this->y_;
		delete[] rwork_;
		delete[] iwork_;
	}

}
