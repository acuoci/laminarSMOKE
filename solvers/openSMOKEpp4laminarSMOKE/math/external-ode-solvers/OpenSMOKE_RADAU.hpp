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
	double ccontra(int i, double s, double *cont, int *lrc)
	{
		int		I=i;
		double	S=s;

		#if defined(_WIN32) || defined(_WIN64) 
			return CONTRA(&I, &S, cont, lrc);
		#else
			return contra_(&I, &S, cont, lrc);
		#endif
	}

	template <typename T>
	double ccontr5(int i, double s, double *cont, int *lrc)
	{
		int		I=i;
		double	S=s;

		#if defined(_WIN32) || defined(_WIN64) 
			return CONTR5(&I, &S, cont, lrc);
		#else
			return contr5_(&I, &S, cont, lrc);
		#endif
	}

	template <typename T>
	OpenSMOKE_RADAU<T>::OpenSMOKE_RADAU(T* odeSystem)
	{	
		this->SetDefaultValues();

		lwork_ = 0;
		liwork_ = 0;

		
		imas_	=	0;	// identity mass matrix
		ipar_	=	0;
	
		this->firstStep_	= 1.e-5;	
		idid_				= -1;
		iSolver_			= solver_radau5;

		this->odeSystem_ = odeSystem;
	}

	template <typename T>
	void OpenSMOKE_RADAU<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_RADAU<T>::MemoryAllocation(const int n)
	{
		this->n_	=	n;
		this->mLower_	=	this->n_;		// means full matrix

		if (iSolver_ == solver_radau5)
		{
			int LJAC;
			int LMAS;
			int LE;
	
			if (this->mLower_==this->n_)	LJAC = this->n_;
			else				LJAC = this->mLower_+this->mUpper_+1;
	
			if (imas_==0)						LMAS=0;
			else if (imas_==1 && mlmas_==this->n_)	LMAS = this->n_;
			else if (mlmas_<this->n_)					LMAS = mlmas_+mumas_+1;

			if (this->mLower_==this->n_)	LE = this->n_;
			else				LE = 2*this->mLower_+this->mUpper_+1;
	
			lwork_ = this->n_*(LJAC+LMAS+3*LE+12)+20;
			liwork_ = 4*this->n_*this->n_+12*this->n_+20;
		}
		else if (iSolver_ == solver_radau)
		{
			int NSMAX = 7;
			lwork_ 	= (NSMAX+1)*this->n_*this->n_+(3*NSMAX+3)*this->n_+20;
			liwork_	= (2+(NSMAX-1)/2)*this->n_+20;
		}

	//	delete[] this->y0_;
	//	delete[] this->y_;
	//	delete[] rwork_;
	//	delete[] iwork_;

		this->y0_ 	= new double[this->n_];
		this->y_ 		= new double[this->n_];
		rwork_ 	= new double[lwork_];
		iwork_	= new int[liwork_];

		// Default values
		for(int i=0; i<20; i++)
		{
			iwork_[i]	=	0;
			rwork_[i]	=	0.0;
		}

		// Scalar Tolerances
		this->iTolerance_ = 0;
	}

	template <typename T>
	void OpenSMOKE_RADAU<T>::AnalyzeUserOptions()
	{
		if (this->verbose_output_ == false)	this->iOutput_ = 0;
		else                                this->iOutput_ = 1;

		if (this->iSetMaximumNumberOfSteps_ == true)
		{
			if (this->maximumNumberOfSteps_ >= 0)	iwork_[1] = this->maximumNumberOfSteps_;
			else									iwork_[1] = 0;
		}

		if (this->iSetMaximumStep_ == true)
		{
			if (this->maximumStep_ >= 0) rwork_[6] = this->maximumStep_;
			else                         rwork_[6] = 0.;
		}

		if (this->iSetFirstStep_ == true)
		{
			if (this->firstStep_ <= 0) this->firstStep_ = 1.e-5;
		}
	}

	template <typename T>
	void OpenSMOKE_RADAU<T>::Solve(const double xend)
	{
		AnalyzeUserOptions();

		this->x_ = this->x0_;
		this->xend_ = xend;
		memcpy(this->y_, this->y0_, this->n_*sizeof(double));
	
		this->tStart_ =  this->GetClockTime();

		#if defined(_WIN32) || defined(_WIN64) 
			if (iSolver_ == solver_radau5)
			RADAU5(	&this->n_, this->odeSystem_->GetSystemFunctionsStatic, &this->x_, this->y_, &this->xend_, 
					&this->firstStep_, &this->relTolerance_[0], &this->absTolerance_[0],&this->iTolerance_,
					this->odeSystem_->GetAnalyticalJacobianStatic, &this->iJacobian_, &this->mLower_, &this->mUpper_,
					ptMassMatrix_, &imas_, &mlmas_, &mumas_,
					this->odeSystem_->GetWriteFunctionStatic, &this->iOutput_,
					rwork_, &lwork_, iwork_, &liwork_, rpar_, ipar_, &idid_);

			else if (iSolver_ == solver_radau)
			RADAU(	&this->n_, this->odeSystem_->GetSystemFunctionsStatic, &this->x_, this->y_, &this->xend_, 
					&this->firstStep_, &this->relTolerance_[0], &this->absTolerance_[0],&this->iTolerance_,
					this->odeSystem_->GetAnalyticalJacobianStatic, &this->iJacobian_, &this->mLower_, &this->mUpper_,
					ptMassMatrix_, &imas_, &mlmas_, &mumas_,
					this->odeSystem_->GetWriteFunctionStatic, &this->iOutput_,
					rwork_, &lwork_, iwork_, &liwork_, rpar_, ipar_, &idid_);
		#else
			if (iSolver_ == solver_radau5)
			radau5_(	&this->n_, this->odeSystem_->GetSystemFunctionsStatic, &this->x_, this->y_, &this->xend_, 
					&this->firstStep_, &this->relTolerance_[0], &this->absTolerance_[0],&this->iTolerance_,
					this->odeSystem_->GetAnalyticalJacobianStatic, &this->iJacobian_, &this->mLower_, &this->mUpper_,
					ptMassMatrix_, &imas_, &mlmas_, &mumas_,
					this->odeSystem_->GetWriteFunctionStatic, &this->iOutput_,
					rwork_, &lwork_, iwork_, &liwork_, rpar_, ipar_, &idid_);

			else if (iSolver_ == solver_radau)
			radau_(	&this->n_, this->odeSystem_->GetSystemFunctionsStatic, &this->x_, this->y_, &this->xend_, 
					&this->firstStep_, &this->relTolerance_[0], &this->absTolerance_[0],&this->iTolerance_,
					this->odeSystem_->GetAnalyticalJacobianStatic, &this->iJacobian_, &this->mLower_, &this->mUpper_,
					ptMassMatrix_, &imas_, &mlmas_, &mumas_,
					this->odeSystem_->GetWriteFunctionStatic, &this->iOutput_,
					rwork_, &lwork_, iwork_, &liwork_, rpar_, ipar_, &idid_);
		#endif

		this->tEnd_ =  this->GetClockTime();

		this->x0_ = this->x_;
		memcpy(this->y0_, this->y_, this->n_*sizeof(double));

	}

	template <typename T>
	void OpenSMOKE_RADAU<T>::Status() const
	{
		std::cout << "RADAU Status:                       " << std::endl;				// Status
		std::cout << " * Absolute tolerance:              " << this->absTolerance_[0]   << std::endl;	// Absolute tolerance
		std::cout << " * Relative tolerance:              " << this->relTolerance_[0]   << std::endl;	// Relative tolerance
		std::cout << " * Number of steps:                 " << iwork_[15] << std::endl;	// Number of steps taken for the problem so far 
		std::cout << " * Number of function evaluations:  " << iwork_[13] << std::endl;	// Number of f evaluations for the problem so far.
		std::cout << " * Number of Jacobian evaluations:  " << iwork_[14] << std::endl;	// Number of Jacobian evaluations (and of matrix LU decompositions) for the problem so far.
		
		std::cout << " * Number of accepted steps:        " << iwork_[16] << std::endl;	
		std::cout << " * Number of rejected steps:        " << iwork_[17] << std::endl;	
		std::cout << " * Number of LU Factorizations:     " << iwork_[18] << std::endl;	
		std::cout << " * Number of F/B substitutions:     " << iwork_[19] << std::endl;	
	}

	template <typename T>
	std::string OpenSMOKE_RADAU<T>::Tag() const
	{
		if (iSolver_ == solver_radau5)	return "RADAU5";
		else							return "RADAU";
	}

	template <typename T>
	int OpenSMOKE_RADAU<T>::GetNumberOfSteps() const
	{
		return iwork_[15];
	}

	template <typename T>
	int OpenSMOKE_RADAU<T>::GetNumberOfFunctionEvaluations() const
	{
		return iwork_[13];
	}

	template <typename T>
	int OpenSMOKE_RADAU<T>::GetNumberOfJacobianEvaluations() const
	{
		return iwork_[14];
	}

	template <typename T>
	int OpenSMOKE_RADAU<T>:: GetNumberOfLUFactorizations() const
	{
		return iwork_[18];
	}

	template <typename T>
	int OpenSMOKE_RADAU<T>::GetNumberOfNonLinearIterations() const
	{
		return -1;
	}

	template <typename T>
	int OpenSMOKE_RADAU<T>::GetLastOrderUsed() const
	{
		return -1;
	}

	template <typename T>
	double OpenSMOKE_RADAU<T>::GetLastStepUsed() const
	{
		return this->firstStep_;
	}

	template <typename T>
	int OpenSMOKE_RADAU<T>::GetNumberOfConvergenceFailures() const
	{
		return -1;
	}

	template <typename T>
	int OpenSMOKE_RADAU<T>::GetNumberOfErrorTestFailures() const
	{
		return -1;
	}

	template <typename T>
	OpenSMOKE_RADAU<T>::~OpenSMOKE_RADAU(void)
	{
		delete[] this->y0_;
		delete[] this->y_;
		delete[] rwork_;
		delete[] iwork_;
	}
}


// -------------------------------------------------------------------------------------------------------------------	//
//														IWORK Vector													//
// -------------------------------------------------------------------------------------------------------------------	//
// Position 1
// !=0  the Jacobian is transformed to Hessemberg form (advantageous for large systems)
//      does not work for banded (mljac_<this->n_) systems or for implicit systems (imas_=1)

// Position 2
// Maximum number of steps (default 100,000)

// Position 3
// Maximum number of Newton Iterations for the solution of the implicit system in each step (default is 7)

// Position 4
// != 0 Zero starting values are used for the Newton's method (to be apply when nstep>naccept+nrejected)

// Position 5,6,7 
// Useful only for DAE of index > 1

// Position 8
// Switch for step size strategy
// 1 = model predictive controller (gustafsson) (default) (safer results)
// 2 = conventional (faster for simple problems)

// Position 9,10
// Only for a particular class of problems

// -------------------------------------------------------------------------------------------------------------------	//
//														WORK Vector														//
// -------------------------------------------------------------------------------------------------------------------	//
// Position 1
// raunding unit (default 1e-16)

// Position 2
// Safety factor for size step prediction (default 0.90)

// Position 3
// Decides if the Jacobian must be recalculated
// <0 means that the JAcobian is evaluated at each step
// 0.1 for costly jacobian evaluations
// 0.001 for small systems (default)

// Position 4
// Stopping criterion for Newtons method 
// Default min(0.03, sqrt(rtol))
// Smaller values are safer, but slower

// Position 5,6
// If work(5) < hnew/hold < work(6) the step is not changed, which means that LU decompositions can be saved
// Default (1., 1.2)
// Large, full systems: (1.0, 1.20)
// Fast and safe values: (0.99, 2.0)

// Position 7
// Maximum step size (default xend-x)

// Position 8,9
// Boundary to choose the new step size
// work(8) < hnew/hold < work(9)
// Default: (0.20, 8.0)
