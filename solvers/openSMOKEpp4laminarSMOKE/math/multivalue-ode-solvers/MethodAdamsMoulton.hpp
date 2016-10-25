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
|	License                                                           |
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

namespace OdeSMOKE
{
	template <typename ODESystemKernel>
	const unsigned int MethodAdamsMoulton<ODESystemKernel>::MIN_ORDER = 1;
	template <typename ODESystemKernel>
	const unsigned int MethodAdamsMoulton<ODESystemKernel>::MAX_ORDER = 12;
	template <typename ODESystemKernel>
	const unsigned int MethodAdamsMoulton<ODESystemKernel>::DEFAULT_MAX_CONVERGENCE_ITER = 2;
	template <typename ODESystemKernel>
	const double MethodAdamsMoulton<ODESystemKernel>::DELTA_ALFA1 = 0.075;
	template <typename ODESystemKernel>
	const double MethodAdamsMoulton<ODESystemKernel>::ALFA2 = 0.875;
	template <typename ODESystemKernel>
	const double MethodAdamsMoulton<ODESystemKernel>::DELTA_ALFA3 = 0.07;
	template <typename ODESystemKernel>
	const unsigned int MethodAdamsMoulton<ODESystemKernel>::MAX_CONVERGENCE_FAILURE = 50;

	template <typename ODESystemKernel>
	MethodAdamsMoulton<ODESystemKernel>::MethodAdamsMoulton()
	{
	}

	template <typename ODESystemKernel>
	MethodAdamsMoulton<ODESystemKernel>::~MethodAdamsMoulton()
	{
		delete[] r_;
		delete[] Ep_;
		delete[] z_;
		delete[] v_;
		delete[] alfa2_;
	}

	template <typename ODESystemKernel>
	void MethodAdamsMoulton<ODESystemKernel>::ResetMethod()
	{
		iterOrder_ = 0;
		maxOrderUsed_ = 1;
		orderInNextStep_ = 1;
		maxConvergenceIterations_ = DEFAULT_MAX_CONVERGENCE_ITER;
		p_ = 1;

		odeHStatus_ = H_STATUS_CONST;
		odeOrderStatus_ = ORDER_STATUS_CONST;
		status_ = ODE_STATUS_TO_BE_INITIALIZED;
		convergenceStatus_ = CONVERGENCE_STATUS_OK;

		// Default values
		h_ = 0.;
		min_step_size_ = 0.;
		hNextStep_ = 0.;
		hScale_ = 0.;

		// Initialize cumulative counters
		numberOfSteps_ = 0;
		numberOfFunctionCalls_ = 0;
		numberOfDecreasedSteps_ = 0;
		numberOfIncreasedSteps_ = 0;
		numberOfConvergenceFailuresForOrderMax_ = 0;

		// Initialize local counters
		iterConvergence_ = 0;
		iterConvergenceRate_ = 0;
		iterConvergenceFailure_ = 0;
	}

	template <typename ODESystemKernel>
	void MethodAdamsMoulton<ODESystemKernel>::MemoryAllocationMethod()
	{
		// Ode System Memory Allocation
		this->MemoryAllocation();

		// Safety coefficient for choosing the optimal new order
		deltaAlfa1_ = DELTA_ALFA1;
		deltaAlfa3_ = DELTA_ALFA3;
		alfa2_ = new double[MAX_ORDER + 1];
		for (unsigned int i = 0; i <= MAX_ORDER; i++)
			alfa2_[i] = ALFA2;

		// Auxiliary function
		Eigen::VectorXd gamma(MAX_ORDER + 2);
		gamma(0) = 1.;
		for (unsigned int i = 1; i <= MAX_ORDER; i++)
		{
			double sum = 0.;
			for (unsigned j = 1; j <= i; j++)
				sum -= gamma(j - 1) / double(i - j + 2);
			gamma(i) = sum;
		}

		// Allocating memory for vectors z and v
		z_ = new Eigen::VectorXd[MAX_ORDER + 3];
		v_ = new Eigen::VectorXd[MAX_ORDER + 3];
		for (unsigned int i = 0; i <= MAX_ORDER + 2; i++)
		{
			z_[i].resize(this->ne_);
			v_[i].resize(this->ne_);
		}

		// Allocating memory for vector r
		r_ = new Eigen::VectorXd[MAX_ORDER + 1];
		for (unsigned int i = 0; i <= MAX_ORDER; i++)
		{
			r_[i].resize(MAX_ORDER);
			r_[i].setZero();
		}

		// Filling vector r with constant values specific of the Adams Moulton algorithms
		// See Buzzi-Ferraris, pag 682 (equations 29.191-29.197)
		// Order 1: 1 1
		// Order 2: 1/2 1 1/2
		// Order 3: 5/12 1 3/4 1/6
		// Order 4: 3/8 1 11/12 1/3 1/24
		// Order 5: 251/720 1 25/24 35/72 5/40 1/120
		// ...
		r_[0](0) = 1.;
		r_[1](0) = 1.;
		r_[1](1) = 1.;
		r_[2](1) = 1.;
		for (unsigned int j = 3; j <= MAX_ORDER; j++)
		{
			r_[1](j - 1) = Factorial(j - 1);
			for (unsigned int i = 3; i <= j; i++)
				r_[i - 1](j - 1) = double(j - 1)*r_[i - 1](j - 1 - 1) + r_[i - 2](j - 1 - 1);
			r_[j](j - 1) = 1.;
		}
		for (unsigned int j = 2; j <= MAX_ORDER; j++)
		{
			r_[0](j - 1) = r_[0](j - 1 - 1) + gamma(j - 1);
			r_[1](j - 1) = 1.;
			for (unsigned int i = 3; i <= j + 1; i++)
				r_[i - 1](j - 1) = r_[i - 1](j - 1) / (double(i - 1)*Factorial(j - 1));
		}

		// Initialize vector of error coefficients Ep (Buzzi-Ferraris)
		Ep_ = new double[MAX_ORDER + 1];
		Ep_[0] = 1.;
		Ep_[1] = 1.;
		Ep_[2] = .5;
		Ep_[3] = 1.;
		for (unsigned int j = 3; j < MAX_ORDER; j++)
			Ep_[j + 1] = -double(Factorial(j + 2))*gamma(j + 1);

		// Memory allocation
		b_.resize(this->ne_);
		y_.resize(this->ne_);
		va_.resize(this->ne_);
		deltab_.resize(this->ne_);
		f_.resize(this->ne_);
		vb_.resize(this->ne_);

	}

	template <typename ODESystemKernel>
	unsigned int MethodAdamsMoulton<ODESystemKernel>::FindCorrection(const double t, const Eigen::VectorXd& errorWeights)
	{
		// Call the equations
		this->Equations(v_[0], t, va_);
		numberOfFunctionCalls_++;

		// 
		va_ *= h_;
		deltab_ = va_ - v_[1];
		b_ = deltab_;

		double  oldCorrectionControl;
		for (unsigned int k = 1; k <= maxConvergenceIterations_; k++)
		{
			// Estimation of the error
			const double correctionControl = r_[0](p_ - 1)*OpenSMOKE::ErrorControl(deltab_, errorWeights);

			// If the error is sufficiently small, the convergence is ok
			if (correctionControl < 1.)
			{
				convergenceStatus_ = CONVERGENCE_STATUS_OK;
				return k;
			}

			// If the new error is 2 times larger than the previous one, it is better to stop the Newton's method
			// and try to restart with a smaller step
			if (k > 1 && correctionControl > 2.*oldCorrectionControl)
			{
				convergenceStatus_ = CONVERGENCE_STATUS_FAILURE;
				return k;
			}

			// Store the error
			oldCorrectionControl = correctionControl;

			// Update the solution by adding the vector b: v0 = v0+r0*b
			va_ = r_[0][p_ - 1] * b_;
			va_ += v_[0];

			// Call the equations
			this->Equations(va_, t, f_);
			numberOfFunctionCalls_++;


			va_ = h_*f_;
			vb_ = b_;
			b_ = va_ - v_[1];
			deltab_ = b_ - vb_;
		}

		convergenceStatus_ = CONVERGENCE_STATUS_FAILURE;
		return maxConvergenceIterations_;
	}

	template <typename ODESystemKernel>
	void MethodAdamsMoulton<ODESystemKernel>::ConvergenceRate()
	{
	}

	template <typename ODESystemKernel>
	void MethodAdamsMoulton<ODESystemKernel>::WhatToDoInCaseOfConvergenceFailure(const double tInMeshPoint, double& tStabilize)
	{
		// In case of convergence failures, the step must be reduced
		hScale_ = 0.25;
		if (numberOfSteps_ < 10)
			hScale_ = 0.1;
		hScale_ = std::max(hScale_, min_step_size_ / std::fabs(h_));
		odeHStatus_ = H_STATUS_DECREASED;

		// The new step is chosen
		hNextStep_ *= hScale_;
		numberOfDecreasedSteps_++;
	}

	template <typename ODESystemKernel>
	void MethodAdamsMoulton<ODESystemKernel>::CheckConstraints(Eigen::VectorXd& y)
	{
		// Minimum constraints
		if (min_constraints_ == true)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
			if (y(i) < min_values_(i))
				y(i) = min_values_(i);
		}

		// Maximum constraints
		if (max_constraints_ == true)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
			if (y(i) > max_values_(i))
				y(i) = max_values_(i);
		}
	}

	template <typename ODESystemKernel>
	unsigned int MethodAdamsMoulton<ODESystemKernel>::CheckIllegalConstraints(const Eigen::VectorXd& y)
	{
		unsigned int number_of_illegal_constraints_ = 0;

		if (min_constraints_ == true)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
			if (y(i) < min_values_(i))
				number_of_illegal_constraints_++;
		}

		if (max_constraints_ == true)
		for (unsigned int i = 0; i < this->ne_; i++)
		{
			if (y(i) > max_values_(i))
				number_of_illegal_constraints_++;
		}

		return number_of_illegal_constraints_;
	}

	template <typename ODESystemKernel>
	void MethodAdamsMoulton<ODESystemKernel>::OdeMethodSummary(std::ostream& out)
	{
		out << std::endl;
		out << "Adams-Moulton algorithm" << std::endl;
		out << "---------------------------------------------------------------------------------------" << std::endl;
		out << "* Number of steps:                               " << numberOfSteps_ << std::endl;
		out << "* Max order used:                                " << maxOrderUsed_ << std::endl;
		out << "* Number of calls of equations:                  " << numberOfFunctionCalls_ << std::endl;
		out << "* Number of decreased step sizes:                " << numberOfDecreasedSteps_ << std::endl;
		out << "* Number of increased step sizes:                " << numberOfIncreasedSteps_ << std::endl;
		out << "* Number of convergence failures for order max:  " << numberOfConvergenceFailuresForOrderMax_ << std::endl;
		out << "---------------------------------------------------------------------------------------" << std::endl;
		out << std::endl;
	}

}
