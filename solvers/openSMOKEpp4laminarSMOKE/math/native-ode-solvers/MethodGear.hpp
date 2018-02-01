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
	const unsigned int MethodGear<ODESystemKernel>::MIN_ORDER = 1;
	template <typename ODESystemKernel>
	const unsigned int MethodGear<ODESystemKernel>::MAX_ORDER = 5;
	template <typename ODESystemKernel>
	const unsigned int MethodGear<ODESystemKernel>::DEFAULT_MAX_CONVERGENCE_ITER = 4;
	template <typename ODESystemKernel>
	const unsigned int MethodGear<ODESystemKernel>::MAX_ITERATIONS_FACTORIZATION = 20;
	template <typename ODESystemKernel>
	const unsigned int MethodGear<ODESystemKernel>::MAX_ITERATIONS_JACOBIAN = 50;
	template <typename ODESystemKernel>
	const double MethodGear<ODESystemKernel>::DELTA_ALFA1 = .05;
	template <typename ODESystemKernel>
	const double MethodGear<ODESystemKernel>::ALFA2 = .8;
	template <typename ODESystemKernel>
	const double MethodGear<ODESystemKernel>::DELTA_ALFA3 = .25;
	template <typename ODESystemKernel>
	const unsigned int MethodGear<ODESystemKernel>::MAX_CONVERGENCE_FAILURE = 50;

	template <typename ODESystemKernel>
	MethodGear<ODESystemKernel>::MethodGear()
	{
	}

	template <typename ODESystemKernel>
	MethodGear<ODESystemKernel>::~MethodGear()
	{
		delete[] r_;
		delete[] Ep_;
		delete[] z_;
		delete[] v_;
		delete[] alfa2_;
	}

	template <typename ODESystemKernel>
	void MethodGear<ODESystemKernel>::ResetMethod()
	{
		// Kernel default values
		this->ResetKernel();

		p_ = 1;
		iterOrder_ = 0;
		maxOrderUsed_ = 1;
		orderInNextStep_ = 1;
		maxConvergenceIterations_ = DEFAULT_MAX_CONVERGENCE_ITER;

		maxIterationsJacobian_ = 20 + 5 * (this->ne_ - 3);
		if (maxIterationsJacobian_ > MAX_ITERATIONS_JACOBIAN && this->ne_ < 80)
			maxIterationsJacobian_ = MAX_ITERATIONS_JACOBIAN;

		odeHStatus_ = H_STATUS_CONST;
		odeOrderStatus_ = ORDER_STATUS_CONST;
		status_ = ODE_STATUS_TO_BE_INITIALIZED;
		factorizationStatus_ = MATRIX_HAS_TO_BE_FACTORIZED;
		jacobianStatus_ = JACOBIAN_STATUS_HAS_TO_BE_CHANGED;
		jacobianType_ = JACOBIAN_TYPE_NUMERICAL;
		convergenceStatus_ = CONVERGENCE_STATUS_OK;

		// Default values
		h_ = 0.;
		min_step_size_ = 0.;
		hNextStep_ = 0.;
		hScale_ = 0.;

		// Initialize cumulative counters
		numberOfSteps_ = 0;
		numberOfFunctionCalls_ = 0;
		numberOfLinearSystemSolutions_ = 0;
		numberOfJacobians_ = 0;
		numberOfMatrixFactorizations_ = 0;
		numberOfDecreasedSteps_ = 0;
		numberOfIncreasedSteps_ = 0;
		numberOfConvergenceFailuresForOrderMax_ = 0;

		// Initialize local counters
		iterConvergence_ = 0;
		iterConvergenceRate_ = 0;
		iterConvergenceFailure_ = 0;
	}

	template <typename ODESystemKernel>
	void MethodGear<ODESystemKernel>::MemoryAllocationMethod()
	{
		// Kernel memory allocation
		this->MemoryAllocationKernel();

		// Safety coefficient for choosing the optimal new order
		deltaAlfa1_ = DELTA_ALFA1;
		deltaAlfa3_ = DELTA_ALFA3;
		alfa2_ = new double[MAX_ORDER + 1];
		for (unsigned int i = 0; i <= MAX_ORDER; i++)
			alfa2_[i] = ALFA2;

		// Initialize vector of error coefficients Ep (Buzzi-Ferraris, 29.204)
		Ep_ = new double[MAX_ORDER + 1];
		for (unsigned int j = 0; j <= MAX_ORDER; j++)
			Ep_[j] = Factorial(j);

		// Allocating memory for vector r
		r_ = new Eigen::VectorXd[MAX_ORDER + 1];
		for (unsigned int i = 0; i <= MAX_ORDER; i++)
		{
			r_[i].resize(MAX_ORDER);
			r_[i].setZero();
		}

		// Filling vector r with constant values specific of the Gear algorithms
		// See Buzzi-Ferraris, pag 683 (equations 29.199-29.203)
		// Order 1: 1 1
		// Order 2: 2/3 1 1/3
		// Order 3: 6/11 1 6/11 1/11
		// Order 4: 24/50 1 35/50 10/50 1/50
		// Order 5: 120/274 1 225/274 85/274 15/274 1/274
		{
			r_[0](0) = 1.;
			r_[1](0) = 1.;
			for (unsigned int j = 2; j <= MAX_ORDER; j++)
			{
				r_[0](j - 1) = Factorial(j);
				for (unsigned int i = 2; i <= j; i++)
					r_[i - 1](j - 1) = double(j)*r_[i - 1](j - 1 - 1) + r_[i - 2](j - 1 - 1);
				r_[j](j - 1) = 1.;
			}
			double sum = 1.;
			for (unsigned int j = 2; j <= MAX_ORDER; j++)
			{
				sum += 1. / double(j);
				for (unsigned int i = 1; i <= j + 1; i++)
					r_[i - 1](j - 1) /= (Factorial(j)*sum);
			}
		}

		// Allocating memory for vectors z and v
		z_ = new Eigen::VectorXd[MAX_ORDER + 3];
		v_ = new Eigen::VectorXd[MAX_ORDER + 3];
		for (unsigned int i = 0; i <= MAX_ORDER + 2; i++)
		{
			z_[i].resize(this->ne_);
			v_[i].resize(this->ne_);
		}

		// Vectors
		b_.resize(this->ne_);
		y_.resize(this->ne_);
		va_.resize(this->ne_);
		deltab_.resize(this->ne_);
		f_.resize(this->ne_);
		vb_.resize(this->ne_);
	}

	template <typename ODESystemKernel>
	unsigned int MethodGear<ODESystemKernel>::FindCorrection(const double t, const Eigen::VectorXd& errorWeights)
	{
		// Current solution to be updated
		vb_ = v_[0];

		// Check Constraints
		CheckConstraints(vb_);

		// Call the equations
		this->Equations(vb_, t, f_);
		numberOfFunctionCalls_++;

		// Force recalculation of Jacobian matrix (every maxIterationsJacobian_)
		if (numberOfSteps_ >= stepOfLastJacobian_ + maxIterationsJacobian_)
			jacobianStatus_ = JACOBIAN_STATUS_HAS_TO_BE_CHANGED;

		if (jacobianStatus_ == JACOBIAN_STATUS_HAS_TO_BE_CHANGED)
		{
			// The Jacobian is evaluated at the current step
			stepOfLastJacobian_ = numberOfSteps_;

			// Evaluation of the Jacobian matrix
			switch (jacobianType_)
			{
				case JACOBIAN_TYPE_CONST:
					jacobianStatus_ = JACOBIAN_STATUS_OK;
					break;

				case JACOBIAN_TYPE_USERDEFINED:
					jacobianStatus_ = JACOBIAN_STATUS_MODIFIED;
					this->UserDefinedJacobian(vb_, t);
					numberOfJacobians_++;
					factorizationStatus_ = MATRIX_HAS_TO_BE_FACTORIZED;
					break;

				case JACOBIAN_TYPE_NUMERICAL:
					jacobianStatus_ = JACOBIAN_STATUS_MODIFIED;
					this->NumericalJacobian(vb_, t, f_, h_, errorWeights, max_constraints_, max_values_);
					numberOfJacobians_++;
					factorizationStatus_ = MATRIX_HAS_TO_BE_FACTORIZED;
					break;
			}
		}

		// The hr0 value is used to build the G matrix (see Eq. 29.214)
		const double hr0 = h_*r_[0](p_ - 1);

		// If the order was changed, the matrix is factorized (of course!)
		if (factorizationStatus_ == MATRIX_FACTORIZED)
		{
			if (iterOrder_ == 0)
				factorizationStatus_ = MATRIX_HAS_TO_BE_FACTORIZED;
		}

		// The matrix is assembled and factorized
		if (factorizationStatus_ == MATRIX_HAS_TO_BE_FACTORIZED)
		{
			stepOfLastFactorization_ = numberOfSteps_;
			numberOfMatrixFactorizations_++;
			this->BuildAndFactorizeMatrixG(hr0);
			factorizationStatus_ = MATRIX_FACTORIZED;
		}

		// The firs iteration of the Newton's method is performed
		// The first guess value is assumed b=0, which means that the following linear system is solved: G*d = -v1+h*f
		// The solution d, since b=0 on first guess, is equal to the new b
		{
			deltab_ = f_*h_;
			deltab_ -= v_[1];
			this->SolveLinearSystem(deltab_);
			numberOfLinearSystemSolutions_++;
			b_ = deltab_;
		}

		// Newton's iterations
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
			va_ = b_*r_[0](p_ - 1);
			va_ += v_[0];

			// Check constraints for minimum and maximum values
			{
				// Minimum constraints
				if (min_constraints_ == true)
				{
					for (unsigned int j = 0; j < this->ne_; j++)
					if (va_(j) < min_values_(j))
					{
						va_(j) = min_values_(j);
						b_(j) = (min_values_(j) - z_[0](j)) / r_[0](p_ - 1);
					}
				}

				// Maximum constraints
				if (max_constraints_ == true)
				{
					for (unsigned int j = 0; j < this->ne_; j++)
					if (va_(j) > max_values_(j))
					{
						va_(j) = max_values_(j);
						b_(j) = (max_values_(j) - z_[0](j)) / r_[0](p_ - 1);
					}
				}
			}

			// Call the equations
			this->Equations(va_, t, f_);
			numberOfFunctionCalls_++;

			// Apply the Newton's iteration
			{
				deltab_ = f_*h_;
				deltab_ -= v_[1];
				deltab_ -= b_;
				this->SolveLinearSystem(deltab_);
				numberOfLinearSystemSolutions_++;
				b_ += deltab_;
			}
		}

		convergenceStatus_ = CONVERGENCE_STATUS_FAILURE;
		return maxConvergenceIterations_;
	}

	template <typename ODESystemKernel>
	void MethodGear<ODESystemKernel>::WhatToDoInCaseOfConvergenceFailure(const double tInMeshPoint, double& tStabilize)
	{
		// In case of convergence failures, the step must be reduced
		hScale_ = 0.25;
		if (numberOfSteps_ < 10)
			hScale_ = 0.1;
		hScale_ = std::max(hScale_, min_step_size_ / std::fabs(h_));
		odeHStatus_ = H_STATUS_DECREASED;

		// At the same time, if possible, the order must be decreased
		if (p_ > MIN_ORDER)
		{
			p_ -= 1;
			odeOrderStatus_ = ORDER_STATUS_DECREASED;
		}

		// The propensity to change the step size is decreased for the highest orders (4 and 5)
		alfa2_[4] -= 0.01;
		alfa2_[5] -= 0.02;

		// Exclude 4th order if too many convergence failures were observed
		if (alfa2_[4] < 0.)
		{
			alfa2_[4] = 0.;
			maximum_order_ = 3;
			if (p_ > 3)
			{
				p_ = 3;
				odeOrderStatus_ = ORDER_STATUS_DECREASED;
			}
		}

		// Exclude 5th order if too many convergence failures were observed
		if (alfa2_[5] < 0.)
		{
			alfa2_[5] = 0.;
			maximum_order_ = 4;
			if (p_ > 4)
			{
				p_ = 4;
				odeOrderStatus_ = ORDER_STATUS_DECREASED;
			}
		}

		if (numberOfSteps_ > 5)
		{
			numberOfConvergenceFailuresForOrderMax_++;

			// Every 5 convergence failures, the order is reset to the minimum
			if (numberOfConvergenceFailuresForOrderMax_ % 5 == 0)
			{
				if (p_ != MIN_ORDER)
				{
					p_ = MIN_ORDER;
					odeOrderStatus_ = ORDER_STATUS_DECREASED;
				}
			}

			// ???
			if (numberOfConvergenceFailuresForOrderMax_ % 10 == 0 && tStabilize != tInMeshPoint)
			{
				this->Equations(z_[0], tInMeshPoint, va_);
				numberOfFunctionCalls_++;
				z_[1] = va_*h_;
				tStabilize = tInMeshPoint;
			}
		}

		// The re-evaluation of the Jacobian matrix is forced every MAX_ITERATIONS_JACOBIAN / 10 steps
		if (numberOfSteps_ > stepOfLastJacobian_ + MAX_ITERATIONS_JACOBIAN / 10)
			jacobianStatus_ = JACOBIAN_STATUS_HAS_TO_BE_CHANGED;

		// The re-evaluation of the Jacobian matrix is forced if the number of convergence failures reached the maximum allowed
		if (numberOfSteps_ != stepOfLastJacobian_ && iterConvergenceFailure_ == MAX_CONVERGENCE_FAILURE)
			jacobianStatus_ = JACOBIAN_STATUS_HAS_TO_BE_CHANGED;

		// The new step is chosen
		hNextStep_ *= hScale_;
		numberOfDecreasedSteps_++;
	}

	template <typename ODESystemKernel>
	void MethodGear<ODESystemKernel>::ConvergenceRate()
	{
		// In case of too much large error or convergence failure
		if (iterErrorFailure_ > 0 || iterConvergenceFailure_ > 0)
			iterConvergenceRate_ = 0;

		// If the correction was find easily, with only a single iteration,
		// we do not need to apply any correction to the algorithm
		if (iterConvergence_ <= 1)
			return;

		vb_ = v_[0] - z_[0];
		this->JacobianTimesVector(vb_, &va_);
		va_ *= h_;
		va_ += z_[1];
		va_ = v_[1] - va_;
		
		#if defined(__clang__)

			double normd = 0.;
			for (unsigned int i=0;i<va_.size();i++)
				normd += std::fabs(va_(i));

                        double normf = 0.;
                        for (unsigned int i=0;i<v_[1].size();i++)
                                normf += std::fabs(v_[1](i));

		#else
			const double normd = va_.lpNorm<1>();
			const double normf = v_[1].lpNorm<1>();
		#endif

		// If more than 2 iterations were needed and the normd is much larged than the normf,
		// it is better to reduce the step and reduce the order
		if (normd > 3.*normf + .1 && iterConvergence_ > 2)
		{
			// Better to update the Jacobian matrix
			if (numberOfSteps_ > stepOfLastJacobian_ + MAX_ITERATIONS_JACOBIAN / 10)
			if (jacobianType_ != JACOBIAN_TYPE_CONST)
				jacobianStatus_ = JACOBIAN_STATUS_HAS_TO_BE_CHANGED;

			// The step is halved
			hScale_ = 0.5;
			odeHStatus_ = H_STATUS_DECREASED;
			hNextStep_ = h_*hScale_;
			numberOfDecreasedSteps_++;

			// The order is reduced
			if (p_ > 1)
			{
				orderInNextStep_ = p_ - 1;
				odeOrderStatus_ = ORDER_STATUS_DECREASED;
			}

			iterConvergenceRate_ = 0;
			iterConvergenceFailure_ = 1;

			return;
		}

		// If the normd is not so larger than the normf
		if (normd > normf + 0.1)
		{
			// The iterConvergenceRate_ is increased only when we are here
			iterConvergenceRate_++;

			// The step is decreased by a small factor
			if (iterConvergenceRate_ >= 5)
			{
				hScale_ = 0.8;
				odeHStatus_ = H_STATUS_DECREASED;
				hNextStep_ = h_*hScale_;
				numberOfDecreasedSteps_++;

				iterConvergenceRate_ = 0;
				iterConvergenceFailure_ = 1;
			}

			if (jacobianType_ == JACOBIAN_TYPE_CONST)
				return;

			// Better to re-estimate the Jacobian matrix
			if (numberOfSteps_ > stepOfLastJacobian_ + MAX_ITERATIONS_JACOBIAN / 10 && iterConvergence_ > 2)
			{
				jacobianStatus_ = JACOBIAN_STATUS_HAS_TO_BE_CHANGED;

				iterConvergenceRate_ = 0;
				iterConvergenceFailure_ = 1;
			}
		}
	}

	template <typename ODESystemKernel>
	void MethodGear<ODESystemKernel>::CheckConstraints(Eigen::VectorXd& y)
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
	unsigned int MethodGear<ODESystemKernel>::CheckIllegalConstraints(const Eigen::VectorXd& y)
	{
		unsigned int number_of_illegal_constraints_ = 0;

		if (min_constraints_ == true)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
			if (y(i) < min_values_(i))
			{
				std::cout << std::setw(15) << std::left << std::scientific << std::setprecision(3) << y(i);
				std::cout << std::setw(15) << std::left << std::scientific << std::setprecision(3) << min_values_(i);
				std::cout << std::endl;
				number_of_illegal_constraints_++;
			}
		}

		if (max_constraints_ == true)
		for (unsigned int i = 0; i < this->ne_; i++)
		{
			if (y(i) > max_values_(i))
			{
				std::cout << std::setw(15) << std::left << std::scientific << std::setprecision(3) << y(i);
				std::cout << std::setw(15) << std::left << std::scientific << std::setprecision(3) << max_values_(i);
				std::cout << std::endl;
				number_of_illegal_constraints_++;
			}
		}

		return number_of_illegal_constraints_;
	}

	template <typename ODESystemKernel>
	void MethodGear<ODESystemKernel>::OdeMethodSummary(std::ostream& out)
	{
		out << std::endl;
		out << "Gear algorithm" << std::endl;
		out << "---------------------------------------------------------------------------------------" << std::endl;
		out << "* Number of steps:                               " << numberOfSteps_ << std::endl;
		out << "* Max order used:                                " << maxOrderUsed_ << std::endl;
		out << "* Number of calls of equations:                  " << numberOfFunctionCalls_ << std::endl;
		out << "* Number of linear system solutions:             " << numberOfLinearSystemSolutions_ << std::endl;
		out << "* Number of Jacobian evaluations:                " << numberOfJacobians_ << std::endl;
		out << "* Number of matrix factorizations:               " << numberOfMatrixFactorizations_ << std::endl;
		out << "* Number of decreased step sizes:                " << numberOfDecreasedSteps_ << std::endl;
		out << "* Number of increased step sizes:                " << numberOfIncreasedSteps_ << std::endl;
		out << "* Number of convergence failures for order max:  " << numberOfConvergenceFailuresForOrderMax_ << std::endl;
		out << "---------------------------------------------------------------------------------------" << std::endl;
		out << std::endl;

		this->OdeSolverKernelSummary(out);
	}

}
