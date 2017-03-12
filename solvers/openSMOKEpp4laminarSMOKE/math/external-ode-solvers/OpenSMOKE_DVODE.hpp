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

#include "math/OpenSMOKEFunctions.h"
#include <typeinfo>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iomanip>

namespace OpenSMOKE
{
	template <typename T>
	OpenSMOKE_DVODE<T>::OpenSMOKE_DVODE(T* odeSystem)
	{	
		this->SetDefaultValues();

		itask_ = 4;			// no overshooting of xf (rwork_[0] must be provided)			
		istate_ = 1;		// This is the first call for a problem

		iOptions_ = 0;			// Flag indicating whether optional inputs are used (dafault = no)

		mf_ = 22;

		lrw_ = 0;
		liw_ = 0;

		this->odeSystem_ = odeSystem;
	}

	template <typename T>
	void OpenSMOKE_DVODE<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_DVODE<T>::MemoryAllocation(const int n)
	{
		istate_ = 1;					// This is the first call for a problem	

		this->n_ = n;					// Number of equations

		if (mf_ == 10)					// 10: Nonstiff (Adams) method, no Jacobian used.
		{
			lrw_ = 20 + 16*this->n_;
			liw_ = 30;
		}
		else if (mf_==21 || mf_==22)	// 21: Stiff (BDF) method, user-supplied full Jacobian. 22: Stiff method, internally generated full Jacobian
		{
			lrw_ = 22+9*this->n_+2*this->n_*this->n_;
			liw_ = 30+this->n_; 
		}
		else if (mf_==24 || mf_==25)	// 24: Stiff method, user-supplied banded Jacobian. 25: Stiff method, internally generated banded Jacobian
		{
			lrw_ = 22+11*this->n_+(3*this->mLower_+2*this->mUpper_)*this->n_;
			liw_ = 30+this->n_;
		}

		//delete[]	this->y0_;
		//delete[]	this->y_;
		//delete[]	rwork_;
		//delete[]	iwork_;
		this->y0_ 	= new double[this->n_];
		this->y_ 	= new double[this->n_];
		rwork_ 		= new double[lrw_];
		iwork_		= new int[liw_];
		 
		// Default values
		for(int i=0; i<lrw_; i++)
			rwork_[i]	=	0.0;
		for(int i=0; i<liw_; i++)
			iwork_[i]	=	0;


		// Additional options
		iOptions_ = 1;

		iwork_[4] = 0;		// Maximum order to be allowed.  The default value is 12 if METH = 1, and 5 if METH = 2. (See the MF description above for METH.)
							// If MAXORD exceeds the default value, it will be reduced to the default value.  If MAXORD is changed during the problem, it may cause the current order to be reduced.
		iwork_[5] = 0;		// Maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
		iwork_[6] = 0;		// Maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size).  This must be positive to
							// result in a nondefault value.  The default value is 10.
		
		rwork_[4] = 0.;		// Step size to be attempted on the first step. The default value is determined by the solver.
		rwork_[5] = 0.;		// Maximum absolute step size allowed. The default value is infinite.
		rwork_[6] = 0.;		// Minimum absolute step size allowed.  The default value is 0.  (This lower bound is not enforced on the final step before reaching TCRIT when ITASK = 4 or 5.)
	}

	template <typename T>
	void OpenSMOKE_DVODE<T>::AnalyzeUserOptions()
	{
		if (this->iSetMaximumOrder_ == true)
		{
			if (this->maximumOrder_ >= 0)	iwork_[4] = this->maximumOrder_;
			else							iwork_[4] = 0;
		}

		if (this->iSetMaximumNumberOfSteps_ == true)
		{
			if (this->maximumNumberOfSteps_ >=0) iwork_[5] = this->maximumNumberOfSteps_;
			else                                 iwork_[5] = 0;
		}

		if (this->iSetMaximumNumberOfPrints_ == true)
		{
			if (this->maximumNumberOfPrints_ >= 0) iwork_[6] = this->maximumNumberOfPrints_;
			else                                   iwork_[6] = 0;
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
	void OpenSMOKE_DVODE<T>::Solve(const double xend)
	{
		istate_ = 1;					
		rwork_[0] = xend;	// maximum value for overshooting

		AnalyzeUserOptions();

		this->x_ = this->x0_;
		this->xend_ = xend;


		memcpy(this->y_, this->y0_, this->n_*sizeof(double));
	
		this->tStart_ =  this->GetClockTime();

		#if defined(_WIN32) || defined(_WIN64) 
			
			DVODE( this->odeSystem_->GetSystemFunctionsStatic, &this->n_, this->y_, &this->x_, &this->xend_, 
					&this->iTolerance_, this->relTolerance_, this->absTolerance_,
					&itask_, &istate_, &iOptions_, rwork_, &lrw_, iwork_, &liw_, 
					this->odeSystem_->GetAnalyticalJacobianStatic, &mf_, rpar_, ipar_, 
					this->odeSystem_->GetWriteFunctionStatic);
		#else
			dvode_( this->odeSystem_->GetSystemFunctionsStatic, &this->n_, this->y_, &this->x_, &this->xend_, 
					&this->iTolerance_, this->relTolerance_, this->absTolerance_,
					&itask_, &istate_, &iOptions_, rwork_, &lrw_, iwork_, &liw_, 
					this->odeSystem_->GetAnalyticalJacobianStatic, &mf_, rpar_, ipar_,
					this->odeSystem_->GetWriteFunctionStatic);
		#endif

		this->tEnd_ =  this->GetClockTime();

		this->x0_ = this->x_;
		memcpy(this->y0_, this->y_, this->n_*sizeof(double));

	}

	template <typename T>
	void OpenSMOKE_DVODE<T>::Status() const
	{
		std::cout << "DVODE Status:                      " << std::endl;				// Status
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
	}

	template <typename T>
	std::string OpenSMOKE_DVODE<T>::Tag() const
	{
		return "DVODE";
	}

	template <typename T>
	int OpenSMOKE_DVODE<T>::GetNumberOfSteps() const
	{
		return iwork_[10];
	}

	template <typename T>
	int OpenSMOKE_DVODE<T>::GetNumberOfFunctionEvaluations() const
	{
		return iwork_[11];
	}

	template <typename T>
	int OpenSMOKE_DVODE<T>::GetNumberOfJacobianEvaluations() const
	{
		return iwork_[12];
	}

	template <typename T>
	int OpenSMOKE_DVODE<T>:: GetNumberOfLUFactorizations() const
	{
		return iwork_[18];
	}

	template <typename T>
	int OpenSMOKE_DVODE<T>:: GetNumberOfNonLinearIterations() const
	{
		return iwork_[19];
	}

	template <typename T>
	int OpenSMOKE_DVODE<T>:: GetLastOrderUsed() const
	{
		return iwork_[13];
	}

	template <typename T>
	double OpenSMOKE_DVODE<T>:: GetLastStepUsed() const
	{
		return rwork_[10];
	}

	template <typename T>
	int OpenSMOKE_DVODE<T>:: GetNumberOfConvergenceFailures() const
	{
		return iwork_[20];
	}

	template <typename T>
	int OpenSMOKE_DVODE<T>:: GetNumberOfErrorTestFailures() const
	{
		return iwork_[21];
	}

	template <typename T>
	OpenSMOKE_DVODE<T>::~OpenSMOKE_DVODE(void)
	{
		delete[] this->y0_;
		delete[] this->y_;
		delete[] rwork_;
		delete[] iwork_;
	}

}
