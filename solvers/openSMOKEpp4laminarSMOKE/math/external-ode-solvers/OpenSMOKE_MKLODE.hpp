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
#include "intel_ode.h"

namespace OpenSMOKE
{
	template <typename T>
	OpenSMOKE_MKLODE<T>::OpenSMOKE_MKLODE(T* odeSystem)
	{	
		this->SetDefaultValues();

		lrw_ = 0;
		liw_ = 0;

		this->odeSystem_ = odeSystem;
	}

	template <typename T>
	void OpenSMOKE_MKLODE<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_MKLODE<T>::MemoryAllocation(const int n)
	{
		this->n_ = n;					// Number of equations

		lrw_ = max (13*this->n_,(7+2*this->n_)*this->n_);
		liw_ = this->n_;
	
	//	delete[]	this->y0_;
	//	delete[]	this->y_;
	//	delete[]	this->dy_;
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

		// Default
		for (int i=0;i<128;i++)
			ipar[i] = 0;
	}

	template <typename T>
	void OpenSMOKE_MKLODE<T>::AnalyzeUserOptions()
	{
	}

	template <typename T>
	void OpenSMOKE_MKLODE<T>::Solve(const double xend)
	{
		AnalyzeUserOptions();

		this->x_ = this->x0_;
		this->xend_ = xend;

		memcpy(this->y_, this->y0_, this->n_*sizeof(double));
		
		this->tStart_ =  this->GetClockTime();

		this->firstStep_ = 1.e-12;
		this->minimumStep_ = 1.e-12;

		int ierr;
		double scaled_absTolerance = this->absTolerance_[0]/this->relTolerance_[0];
		dodesol_mk52lfn(ipar, &this->n_, &this->x_, &this->xend_, this->y_, this->odeSystem_->GetSystemFunctionsStatic, 
						&this->firstStep_, &this->minimumStep_, this->relTolerance_, &scaled_absTolerance, rwork_, iwork_, &ierr);

		this->tEnd_ =  this->GetClockTime();

		this->x0_ = this->x_;
		memcpy(this->y0_, this->y_, this->n_*sizeof(double));
	}

	template <typename T>
	void OpenSMOKE_MKLODE<T>::Status() const
	{
		std::cout << "MKLODE Status:                      " << std::endl;				// Status
		std::cout << " * Absolute tolerance:              " << this->absTolerance_[0]   << std::endl;	// Absolute tolerance
		std::cout << " * Relative tolerance:              " << this->relTolerance_[0]   << std::endl;	// Relative tolerance
	//	std::cout << " * Number of steps:                 " << GetNumberOfSteps() << std::endl;	// Number of steps taken for the problem so far 
	//	std::cout << " * Number of function evaluations:  " << GetNumberOfFunctionEvaluations() << std::endl;	// Number of f evaluations for the problem so far.
	//	std::cout << " * Number of Jacobian evaluations:  " << GetNumberOfJacobianEvaluations() << std::endl;	// Number of Jacobian evaluations (and of matrix LU decompositions) for the problem so far.
	//	std::cout << " * Last order used:                 " << GetLastOrder() << std::endl;	// Method order last used (successfully).
	//	std::cout << " * Next order to be used:           " << GetNextOrder() << std::endl;	// The order to be attempted on the next step.
	//	std::cout << " * Last step:                       " << GetLastStep() << std::endl;	// the step size in t last used (successfully)
	//	std::cout << " * Next step:                       " << GetNextStep() << std::endl;	// the step size to be attempted on the next step
	//	std::cout << " * Current indipendent variable:    " << GetCurrentIndipendentVariable() << std::endl;	// the current value of the independent variable
	}

	template <typename T>
	std::string OpenSMOKE_MKLODE<T>::Tag() const
	{
		return "MKLODE";
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetNumberOfSteps() const
	{
		return iwork_[10];
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetNumberOfFunctionEvaluations() const
	{
		return iwork_[11];
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetNumberOfJacobianEvaluations() const
	{
		return iwork_[12];
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>:: GetNumberOfLUFactorizations() const
	{
		return iwork_[18];
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetLastOrder() const
	{
		return iwork_[7];
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetNextOrder() const
	{
		return iwork_[6];
	}

	template <typename T>
	double OpenSMOKE_MKLODE<T>::GetLastStep() const
	{
		return rwork_[6];
	}

	template <typename T>
	double OpenSMOKE_MKLODE<T>::GetNextStep() const
	{
		return rwork_[2];
	}

	template <typename T>
	double OpenSMOKE_MKLODE<T>::GetCurrentIndipendentVariable() const
	{
		return rwork_[3];
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::ErrorTestFailures() const
	{
		return iwork_[13];
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetNumberOfNonLinearIterations() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetLastOrderUsed() const
	{
		return 0;
	}

	template <typename T>
	double OpenSMOKE_MKLODE<T>::GetLastStepUsed() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetNumberOfConvergenceFailures() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_MKLODE<T>::GetNumberOfErrorTestFailures() const
	{
		return 0;
	}

	template <typename T>
	OpenSMOKE_MKLODE<T>::~OpenSMOKE_MKLODE(void)
	{
		delete[] this->y0_;
		delete[] this->y_;
		delete[] rwork_;
		delete[] iwork_;
	}

}
