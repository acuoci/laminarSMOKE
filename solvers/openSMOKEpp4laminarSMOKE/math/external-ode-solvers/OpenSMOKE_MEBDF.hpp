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
	OpenSMOKE_MEBDF<T>::OpenSMOKE_MEBDF(T* odeSystem)
	{	
		this->SetDefaultValues();

		mbdn_  = new int[4];
		masbdn_ = new int[4];

		istate_ = 1;		// This is the first call for a problem

		mf_ = 22;

		lrw_ = 0;
		liw_ = 0;

		ierr_ = 0;
		lout_ = 1;

		this->odeSystem_ = odeSystem;
	}

	template <typename T>
	void OpenSMOKE_MEBDF<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_MEBDF<T>::MemoryAllocation(const int n)
	{
		istate_ = 1;					// This is the first call for a problem	

		this->n_ = n;					// Number of equations

		masbdn_[0] = 0;					// Identity mass matrix
		masbdn_[1] = 0;					// Not used
		masbdn_[2] = 0;					// Not used
		masbdn_[3] = 0;					// Not used

		if (mf_==21 || mf_==22)	// 21: Stiff (BDF) method, user-supplied full Jacobian. 22: Stiff method, internally generated full Jacobian
		{
			lrw_ = (33+3*this->n_)*this->n_+3;
			liw_ = 14+this->n_; 
		}
		else if (mf_==23 || mf_==24)	// 23: Stiff method, user-supplied banded Jacobian. 24: Stiff method, internally generated banded Jacobian
		{
			// It is not clear from the documentations
			lrw_ = 0;
			liw_ = 0;
			
			mbdn_[0] = this->mLower_;
			mbdn_[1] = this->mUpper_;
			mbdn_[2] = this->mUpper_+this->mLower_+1;
			mbdn_[3] = 2*this->mLower_+this->mUpper_+1;
		}

	//	delete[]	this->y0_;
	//	delete[]	this->y_;
	//	delete[]	rwork_;
	//	delete[]	iwork_;
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
		iwork_[0]  = 0;
		iwork_[1]  = 0;
		iwork_[2]  = 0;

		iwork_[13] = 100000;		// Maximum number of steps
		this->maximumOrder_ = 7;	// Maximum order
		this->iTolerance_ = 2;		// Index of tolerance policy
		this->firstStep_=1.e-12;	// First step
	}

	template <typename T>
	void OpenSMOKE_MEBDF<T>::AnalyzeUserOptions()
	{
		if (this->iSetMaximumNumberOfSteps_ == true)
		{
			if (this->maximumNumberOfSteps_ >= 0)	iwork_[13] = this->maximumNumberOfSteps_;
			else									iwork_[13] = 100000;
		}

		if (this->iSetMaximumOrder_ == true)
		{
			if (this->maximumOrder_ <= 0)		this->maximumOrder_ = 7;
			else if (this->maximumOrder_ > 7)	this->maximumOrder_ = 7;
		}
	}

	template <typename T>
	void OpenSMOKE_MEBDF<T>::Solve(const double xend)
	{
		istate_ = 1;					// This is the first call for a problem	

		AnalyzeUserOptions();

		this->x_ = this->x0_;
		this->xend_ = xend;
		double xout_ = xend;
		double h0_ = this->firstStep_;
		double t0_ = this->x_;

		memcpy(this->y_, this->y0_, this->n_*sizeof(double));


		this->tStart_ =  this->GetClockTime();

		#if defined(_WIN32) || defined(_WIN64) 
			
		MEBDF( &this->n_, &t0_, &h0_, this->y_, &xout_, &this->xend_,
			&mf_, &istate_ , &lout_, &lrw_, rwork_, &liw_, iwork_, mbdn_, masbdn_, &this->maximumOrder_, 
			&this->iTolerance_, this->relTolerance_, this->absTolerance_, rpar_, ipar_,
				this->odeSystem_->GetSystemFunctionsStatic,
				this->odeSystem_->GetAnalyticalJacobianStatic,
				this->odeSystem_->GetMassMatrixMEBDFStatic, &ierr_, this->odeSystem_->GetWriteFunctionStatic);
		
		istate_ = 0;
		
		MEBDF( &this->n_, &t0_, &h0_, this->y_, &xout_, &this->xend_,
			&mf_, &istate_ , &lout_, &lrw_, rwork_, &liw_, iwork_, mbdn_, masbdn_, &this->maximumOrder_, 
			&this->iTolerance_, this->relTolerance_, this->absTolerance_, rpar_, ipar_,
				this->odeSystem_->GetSystemFunctionsStatic,
				this->odeSystem_->GetAnalyticalJacobianStatic,
				this->odeSystem_->GetMassMatrixMEBDFStatic, &ierr_, this->odeSystem_->GetWriteFunctionStatic);

		#else
		
		mebdf_( &this->n_, &t0_, &h0_, this->y_, &xout_, &this->xend_,
			&mf_, &istate_ , &lout_, &lrw_, rwork_, &liw_, iwork_, mbdn_, masbdn_, &this->maximumOrder_, 
			&this->iTolerance_, this->relTolerance_, this->absTolerance_, rpar_, ipar_,
				this->odeSystem_->GetSystemFunctionsStatic,
				this->odeSystem_->GetAnalyticalJacobianStatic,
				this->odeSystem_->GetMassMatrixMEBDFStatic,
			&ierr_, this->odeSystem_->GetWriteFunctionStatic);
		
		istate_ = 0;
		
		mebdf_( &this->n_, &t0_, &h0_, this->y_, &xout_, &this->xend_,
			&mf_, &istate_ , &lout_, &lrw_, rwork_, &liw_, iwork_, mbdn_, masbdn_, &this->maximumOrder_, 
			&this->iTolerance_, this->relTolerance_, this->absTolerance_, rpar_, ipar_,
				this->odeSystem_->GetSystemFunctionsStatic,
				this->odeSystem_->GetAnalyticalJacobianStatic,
				this->odeSystem_->GetMassMatrixMEBDFStatic,
			&ierr_, this->odeSystem_->GetWriteFunctionStatic);

		#endif

		this->tEnd_ =  this->GetClockTime();

		this->x0_ = this->x_;
		memcpy(this->y0_, this->y_, this->n_*sizeof(double));

	//	if (istate_!=0)
	//	{	
	//		std::cout << "istate: " << istate_ << std::endl;
	//		Status();
	//		getchar(); 
	//	}

	}

	template <typename T>
	void OpenSMOKE_MEBDF<T>::Status() const
	{
		std::cout << "MEBDF Status:                      " << std::endl;				// Status
		std::cout << " * Absolute tolerance:              " << this->absTolerance_[0]   << std::endl;	// Absolute tolerance
		std::cout << " * Relative tolerance:              " << this->relTolerance_[0]   << std::endl;	// Relative tolerance

		std::cout << " * Number of steps:                 " << iwork_[4] << std::endl;	// Number of steps taken for the problem so far 
		std::cout << " * Number of failed steps:          " << iwork_[5] << std::endl;	// Number of steps taken for the problem so far 

		std::cout << " * Number of function evaluations:  " << iwork_[6] << std::endl;	// Number of f evaluations for the problem so far.
		std::cout << " * Number of Jacobian evaluations:  " << iwork_[7] << std::endl;	// Number of Jacobian evaluations (and of matrix LU decompositions) for the problem so far.
		std::cout << " * Last order used:                 " << iwork_[3] << std::endl;	// Method order last used (successfully).
		std::cout << " * Next order to be used:           " << iwork_[14] << std::endl;	// The order to be attempted on the next step.
		
		std::cout << " * Last time step:                  " << rwork_[1] << std::endl;	// the step size in t last used (successfully)
	}

	template <typename T>
	std::string OpenSMOKE_MEBDF<T>::Tag() const
	{
		return "MEBDF";
	}

	template <typename T>
	int OpenSMOKE_MEBDF<T>::GetNumberOfSteps() const
	{
		return iwork_[4];
	}

	template <typename T>
	int OpenSMOKE_MEBDF<T>::GetNumberOfFunctionEvaluations() const
	{
		return iwork_[6];
	}

	template <typename T>
	int OpenSMOKE_MEBDF<T>::GetNumberOfJacobianEvaluations() const
	{
		return iwork_[7];
	}

	template <typename T>
	int OpenSMOKE_MEBDF<T>::GetNumberOfLUFactorizations() const
	{
		return iwork_[8];
	}

	template <typename T>
	int OpenSMOKE_MEBDF<T>::GetNumberOfNonLinearIterations() const
	{
		return -1;
	}

	template <typename T>
	int OpenSMOKE_MEBDF<T>::GetLastOrderUsed() const
	{
		return iwork_[3];
	}

	template <typename T>
	double OpenSMOKE_MEBDF<T>::GetLastStepUsed() const
	{
		return rwork_[1];
	}

	template <typename T>
	int OpenSMOKE_MEBDF<T>::GetNumberOfConvergenceFailures() const
	{
		return -1;
	}

	template <typename T>
	int OpenSMOKE_MEBDF<T>::GetNumberOfErrorTestFailures() const
	{
		return -1;
	}

	template <typename T>
	OpenSMOKE_MEBDF<T>::~OpenSMOKE_MEBDF(void)
	{
                delete[] this->y0_;
		delete[] this->y_;
		delete[] rwork_;
		delete[] iwork_;
		delete[] mbdn_;
		delete[] masbdn_;
	}

}
