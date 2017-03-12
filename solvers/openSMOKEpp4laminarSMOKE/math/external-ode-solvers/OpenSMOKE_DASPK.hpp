/*-----------------------------------------------------------------------*\
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
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
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
	OpenSMOKE_DASPK<T>::OpenSMOKE_DASPK(T* odeSystem)
	{	
		this->SetDefaultValues();

		lrw_ = 0;
		liw_ = 0;

		this->odeSystem_ = odeSystem;
	}

	template <typename T>
	void OpenSMOKE_DASPK<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_DASPK<T>::MemoryAllocation(const int n)
	{
		this->n_ = n;					// Number of equations

		lrw_ = 50 + 9*this->n_ + this->n_*this->n_;
		liw_ = 40 + this->n_;
		
		this->y0_ 	= new double[this->n_];
		this->y_ 	= new double[this->n_];
		this->dy_ 	= new double[this->n_];
		rwork_ 		= new double[lrw_];
		iwork_		= new int[liw_];
		ipar_       = new int[1];
		
		 
		// Default values
		for(int i=0; i<lrw_; i++)
			rwork_[i]	=	0.0;
		for(int i=0; i<liw_; i++)
			iwork_[i]	=	0;

		// Default
		for (int i=0;i<20;i++)
			info[i] = 0;
	}

	template <typename T>
	void OpenSMOKE_DASPK<T>::AnalyzeUserOptions()
	{
	}

	template <typename T>
	void OpenSMOKE_DASPK<T>::Solve(const double xend)
	{
		AnalyzeUserOptions();

		this->x_ = this->x0_;
		this->xend_ = xend;

		memcpy(this->y_, this->y0_, this->n_*sizeof(double));
		this->odeSystem_->GetSystemDerivativesCallBack(&this->x0_, this->y0_, this->dy_);
		
		this->tStart_ =  this->GetClockTime();

		int IDID;
		
		#if defined(_WIN32) || defined(_WIN64) 

		DDASPK (this->odeSystem_->GetSystemFunctionsStatic, &this->n_, &this->x_, this->y_, this->dy_, &this->xend_,
		        info, this->relTolerance_, this->absTolerance_, &IDID, rwork_, &lrw_, iwork_, &liw_,
				rpar_, ipar_,  this->odeSystem_->GetAnalyticalJacobianStatic, 
				this->odeSystem_->GetKrylovSolverStatic,
				this->odeSystem_->GetWriteFunctionStatic);

		#else

		ddaspk_ (this->odeSystem_->GetSystemFunctionsStatic, &this->n_, &this->x_, this->y_, this->dy_, &this->xend_,
		        info, this->relTolerance_, this->absTolerance_, &IDID, rwork_, &lrw_, iwork_, &liw_,
				rpar_, ipar_,  this->odeSystem_->GetAnalyticalJacobianStatic, 
				this->odeSystem_->GetKrylovSolverStatic,
				this->odeSystem_->GetWriteFunctionStatic);
		#endif

		this->tEnd_ =  this->GetClockTime();

		this->x0_ = this->x_;
		memcpy(this->y0_, this->y_, this->n_*sizeof(double));
	}

	template <typename T>
	void OpenSMOKE_DASPK<T>::Status() const
	{
		std::cout << "DASPK Status:                      " << std::endl;				// Status
		std::cout << " * Absolute tolerance:              " << this->absTolerance_[0]   << std::endl;	// Absolute tolerance
		std::cout << " * Relative tolerance:              " << this->relTolerance_[0]   << std::endl;	// Relative tolerance
		std::cout << " * Number of steps:                 " << GetNumberOfSteps() << std::endl;	// Number of steps taken for the problem so far 
		std::cout << " * Number of function evaluations:  " << GetNumberOfFunctionEvaluations() << std::endl;	// Number of f evaluations for the problem so far.
		std::cout << " * Number of Jacobian evaluations:  " << GetNumberOfJacobianEvaluations() << std::endl;	// Number of Jacobian evaluations (and of matrix LU decompositions) for the problem so far.
		std::cout << " * Last order used:                 " << GetLastOrderUsed() << std::endl;	// Method order last used (successfully).
		std::cout << " * Last step:                       " << GetLastStepUsed() << std::endl;	// the step size in t last used (successfully)
	}

	template <typename T>
	std::string OpenSMOKE_DASPK<T>::Tag() const
	{
		return "DASPK";
	}

	template <typename T>
	int OpenSMOKE_DASPK<T>::GetNumberOfSteps() const
	{
		return iwork_[10];
	}

	template <typename T>
	int OpenSMOKE_DASPK<T>::GetNumberOfFunctionEvaluations() const
	{
		return iwork_[11];
	}

	template <typename T>
	int OpenSMOKE_DASPK<T>::GetNumberOfJacobianEvaluations() const
	{
		return iwork_[12];
	}

	template <typename T>
	int OpenSMOKE_DASPK<T>:: GetNumberOfLUFactorizations() const
	{
		return ipar_[0];
	}

	template <typename T>
	int OpenSMOKE_DASPK<T>::GetNumberOfNonLinearIterations() const
	{
		return iwork_[18];
	}

	template <typename T>
	int OpenSMOKE_DASPK<T>::GetLastOrderUsed() const
	{
		return iwork_[7];
	}

	template <typename T>
	double OpenSMOKE_DASPK<T>::GetLastStepUsed() const
	{
		return rwork_[6];
	}

	template <typename T>
	int OpenSMOKE_DASPK<T>::GetNumberOfConvergenceFailures() const
	{
		return iwork_[15];
	}

	template <typename T>
	int OpenSMOKE_DASPK<T>::GetNumberOfErrorTestFailures() const
	{
		return iwork_[13];
	}

	template <typename T>
	OpenSMOKE_DASPK<T>::~OpenSMOKE_DASPK(void)
	{
		delete[] this->y0_;
		delete[] this->y_;
		delete[] this->dy_;
		delete[] rwork_;
		delete[] iwork_;
	}
}
