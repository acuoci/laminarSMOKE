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
	// Static consts
	template <typename Method>
	const double MultiValueSolver<Method>::DEFAULT_TOL_ABS = 1.e-10;
	template <typename Method>
	const double MultiValueSolver<Method>::DEFAULT_TOL_REL = 100.*OpenSMOKE::MachEpsFloat();
	template <typename Method>
	const double MultiValueSolver<Method>::DEFAULT_HSCALE_MAX1 = 1000.;
	template <typename Method>
	const double MultiValueSolver<Method>::DEFAULT_HSCALE_MAX2 = 10.;
	template <typename Method>
	const double MultiValueSolver<Method>::DEFAULT_HSCALE_MAX3 = 2.;
	template <typename Method>
	const unsigned int MultiValueSolver<Method>::DEFAULT_MAXIMUM_NUMBER_OF_STEPS = 500000;
	template <typename Method>
	const unsigned int MultiValueSolver<Method>::MAX_ERROR_FAILURE = 50;

	template <typename Method>
	MultiValueSolver<Method>::MultiValueSolver()
	{
		this->ne_ = 0;

		SetDefaultValues();
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetDefaultValues()
	{
		// Maximum integration order
		this->maximum_order_ = this->MAX_ORDER;

		// Maximum number of steps
		max_number_steps_ = DEFAULT_MAXIMUM_NUMBER_OF_STEPS;

		// Maximum step size
		max_step_size_ = OPENSMOKE_BIG_DOUBLE;
		user_defined_max_step_size_ = false;

		// Minimum step size
		this->min_step_size_ = 0.;
		user_defined_min_step_size_ = false;

		// First step size
		first_step_size_ = 0.;
		user_defined_first_step_size_ = false;

		// Minimum and maximum constraints
		this->max_constraints_ = false;
		this->min_constraints_ = false;

		// Direction of integration
		direction_plus_ = true;

		// Call print function
		printResults_ = true;

		// Tolerances
		abs_tolerances_scalar_ = true;
		rel_tolerances_scalar_ = true;
		abs_tolerance_ = DEFAULT_TOL_ABS;
		rel_tolerance_ = DEFAULT_TOL_REL;

		// Contraint on the maximum value of the dependent variable
		calculation_with_strict_tmax_ = false;
		tStrictMax_ = 0.;

		// Stop criteria
		stopIntegrationForSmallYPrimeNorm1_ = false;
		YPrimeNorm1_ = 0.;
	}

	template <typename Method>
	void MultiValueSolver<Method>::Reset()
	{
		hPreviousStep_ = 0.;
		hNordsieck_ = 0.;

		hMinUsed_ = OPENSMOKE_BIG_DOUBLE;
		hMaxUsed_ = 0.;
		hScaleMax_ = DEFAULT_HSCALE_MAX1;

		this->orderInNextStep_ = 1;
		fixedOrderMin_ = false;

		numberOfStepsWithoutChanges_ = 0;
		numberOfConvergenceFailure_ = 0;
		numberOfConvergenceSuccess_ = 0;
		numberOfErrorCheckFailure_ = 0;
		numberOfErrorCheckSuccess_ = 0;

		iterOfConstStepsForLargeFactorizationTime_ = 0;
		iterConstSteps_ = 0;
		iterMinimumOrder_ = 0;

		cpuTimeEquationSystem_ = 0.;
		cpuTimeToFactorize_ = 0.;
		errorStatus_ = ERROR_STATUS_OK;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetPrint(const bool flag)
	{
		printResults_ = flag;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMaximumOrder(const unsigned int maximum_order)
	{
		this->maximum_order_ = maximum_order;
		if (maximum_order > this->MAX_ORDER)
			this->maximum_order_ = this->MAX_ORDER;
		if (maximum_order < 1)
			this->maximum_order_ = this->MAX_ORDER;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetAbsoluteTolerances(const Eigen::VectorXd& abs_tolerances)
	{
		if (abs_tolerances.size() != this->ne_)
			FatalErrorMessage("The size of the absolute tolerance vector is wrong");

		abs_tolerances_.resize(this->ne_);
		abs_tolerances_ = abs_tolerances;
		abs_tolerances_scalar_ = false;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetRelativeTolerances(const Eigen::VectorXd& rel_tolerances)
	{
		if (rel_tolerances.size() != this->ne_)
			FatalErrorMessage("The size of the relative tolerance vector is wrong");

		rel_tolerances_.resize(this->ne_);
		rel_tolerances_ = rel_tolerances;
		rel_tolerances_scalar_ = false;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetAbsoluteTolerances(const double abs_tolerance)
	{
		if (abs_tolerance < 1.e-32)
			FatalErrorMessage("The user-defined absolute tolerance is too small");

		abs_tolerance_ = abs_tolerance;
		abs_tolerances_scalar_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetRelativeTolerances(const double rel_tolerance)
	{
		if (rel_tolerance < 1.e-32)
			FatalErrorMessage("The user-defined relative tolerance is too small");

		rel_tolerance_ = rel_tolerance;
		rel_tolerances_scalar_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMaximumValues(const Eigen::VectorXd& max_values)
	{
		if (max_values.size() != this->ne_)
			FatalErrorMessage("The size of the maximum value vector is wrong");

		this->max_values_.resize(this->ne_);
		this->max_values_ = max_values;
		this->max_constraints_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMaximumValues(const double max_value)
	{
		this->max_values_.resize(this->ne_);
		this->max_values_.setConstant(max_value);
		this->max_constraints_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMinimumValues(const Eigen::VectorXd& min_values)
	{
		if (min_values.size() != this->ne_)
			FatalErrorMessage("The size of the maximum value vector is wrong");

		this->min_values_.resize(this->ne_);
		this->min_values_ = min_values;
		this->min_constraints_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMinimumValues(const double min_value)
	{
		this->min_values_.resize(this->ne_);
		this->min_values_.setConstant(min_value);
		this->min_constraints_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMaximumNumberOfSteps(const unsigned int max_number_steps)
	{
		if (max_number_steps <= 0)
			FatalErrorMessage("The user-defined maximum number of steps must be at least equal to 1");

		max_number_steps_ = max_number_steps;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMaximumStepSize(const double max_step_size)
	{
		if (max_step_size <= 0)
			FatalErrorMessage("The user-defined maximum step must be larger than 0");

		this->max_step_size_ = max_step_size;
		user_defined_max_step_size_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMinimumStepSize(const double min_step_size)
	{
		if (min_step_size <= 0)
			FatalErrorMessage("The user-defined minimum step must be larger than 0");

		this->min_step_size_ = min_step_size;
		user_defined_min_step_size_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetFirstStepSize(const double first_step_size)
	{
		if (first_step_size <= 0)
			FatalErrorMessage("The user-defined first step must be larger than 0");

		this->first_step_size_ = first_step_size;
		user_defined_first_step_size_ = true;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetStopConditionMaximumYPrimeNorm1(const double YPrimeNorm1)
	{
		stopIntegrationForSmallYPrimeNorm1_ = true;
		YPrimeNorm1_ = YPrimeNorm1;
	}

	template <typename Method>
	void MultiValueSolver<Method>::UnsetStopConditions()
	{
		stopIntegrationForSmallYPrimeNorm1_ = false;
		YPrimeNorm1_ = 0.;
	}

	template <typename Method>
	void MultiValueSolver<Method>::Solution(Eigen::VectorXd& solution) const
	{
		solution = this->y_;
	}

	template <typename Method>
	void MultiValueSolver<Method>::FirstOrderDerivatives(Eigen::VectorXd& first_order_derivatives)
	{
		const double hp = (tOut_ - t_) / hNordsieck_;

		if (hp == 0.)
		{
			first_order_derivatives = this->z_[1] / hNordsieck_;
		}
		else
		{
			double hr = 1.;

			first_order_derivatives.setZero();
			for (unsigned int i = 1; i <= this->p_; i++)
			{
				hr *= hp;
				this->va_ = double(i) * this->z_[i] * hr;
				first_order_derivatives += this->va_;
			}

			first_order_derivatives /= (tOut_ - t_);
		}
	}

	template <typename Method>
	void MultiValueSolver<Method>::SecondOrderDerivatives(Eigen::VectorXd& second_order_derivatives)
	{
		const double hp = (tOut_ - t_) / hNordsieck_;

		if (hp == 0.)
		{
			second_order_derivatives = this->z_[2] * 2. / (hNordsieck_*hNordsieck_);
		}
		else
		{
			double hr = hp;

			second_order_derivatives.setZero();
			unsigned int coeff = 1;
			for (unsigned int i = 2; i <= this->p_; i++)
			{
				hr *= hp;
				this->va_ = double(coeff) * this->z_[i] * hr;
				second_order_derivatives += this->va_;
				coeff += i;
			}

			second_order_derivatives *= 2. / (tOut_ - t_) / (tOut_ - t_);
		}
	}

	template <typename Method>
	void MultiValueSolver<Method>::ErrorEstimation(Eigen::VectorXd& error_estimation)
	{
		for (unsigned int j = 0; j < this->ne_; j++)
			error_estimation(j) = this->Ep_[this->p_]*this->z_[this->p_+1](j);
	}

	template <typename Method>
	void MultiValueSolver<Method>::FatalErrorMessage(const std::string message)
	{
		std::cout << "Fatal error: " << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(-1);
	}

	template <typename Method>
	void MultiValueSolver<Method>::WarningMessage(const std::string message)
	{
		if (printResults_ == true)
			std::cout << "Warning message: " << message << std::endl;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetInitialConditions(const double t0, const Eigen::VectorXd& y0)
	{
		// Check if memory allocation is needed (first time or the number of equations is changed)
		bool memory_allocation_is_needed = false;
		if (this->ne_ != y0.size())
			memory_allocation_is_needed = true;

		// Initial value of independent variable
		tOut_ = t0;
		t_ = t0;
		t0_ = t0;
		tInMeshPoint_ = t0;
		tStabilize_ = t0;

		// Total number of equations
		this->ne_ = boost::lexical_cast<unsigned int>(y0.size());
		if (this->ne_==0)
			FatalErrorMessage("Wrong number of equations. The number of equations must be > 0");

		// Allocating local memory
		if (memory_allocation_is_needed == true)
		{
			one_over_epsilon_.resize(this->ne_);
			this->MemoryAllocationMethod();
		}

		// Reset the method class
		this->ResetMethod();

		// Initial value of unknowns
		y0_ = y0;

		// Reset counters
		Reset();
	
		// Call the equation system to track CPU time
		cpuTimeEquationSystem_ = OpenSMOKE::OpenSMOKEGetCpuTime();
		this->Equations(y0_, t0_, this->f_);
		cpuTimeEquationSystem_ = OpenSMOKE::OpenSMOKEGetCpuTime() - cpuTimeEquationSystem_;
	}

	template <typename Method>
	void MultiValueSolver<Method>::SetMaxConstraintOnIndependentVariable(const double tMax)
	{
		calculation_with_strict_tmax_ = true;
		tStrictMax_ = tMax;
	}

	template <typename Method>
	void MultiValueSolver<Method>::UnsetMaxConstraintOnIndependentVariable()
	{
		calculation_with_strict_tmax_ = false;
		tStrictMax_ = 0.;
	}

	template <typename Method>
	OdeStatus MultiValueSolver<Method>::Solve(const double tOut)
	{
		if (this->status_ < 0)
		{
			FatalErrorMessage("The requested integration cannot be performed because of previous errors");
			return this->status_;
		}

		// Check the feasibility of the constraint on the independent variable
		if (calculation_with_strict_tmax_ == true)
		{
			if ( (tStrictMax_ - tOut)*(tOut - t0_) < 0.)
			{
				this->status_ = ODE_STATUS_ILLEGAL_MAX_INDEPENDENT_VARIABLE;
				ParseOdeStatus();
				return this->status_;
			}
		}

		// Set the requested output independent variable
		tOut_ = tOut;
		if (t0_ > tOut_) direction_plus_ = false;

		// Evaluation of the initial step
		if (this->status_ == ODE_STATUS_TO_BE_INITIALIZED)
		{
			// Check the constraints on the initial conditions
			if (this->CheckIllegalConstraints(y0_) != 0)
			{
				this->status_ = ODE_STATUS_ILLEGAL_CONSTRAINTS;
				ParseOdeStatus();
				return this->status_;
			}

			// First step size
			if (user_defined_first_step_size_ == false)
				first_step_size_ = CalculateInitialStepSize();
		
			// Initialize calculations
			this->z_[0] = y0_;
			this->z_[1] = this->f_*first_step_size_;
			hPreviousStep_ = first_step_size_;
			hNordsieck_ = first_step_size_;
			this->hNextStep_ = first_step_size_;

			// Update the status
			this->status_ = ODE_STATUS_CONTINUATION;
		}
		// Continuation (i.e. the integration is restarted without changing the initial conditions)
		else
		{
			// Check the constraints on the initial conditions
			if (this->CheckIllegalConstraints(this->z_[0]) != 0)
			{
				this->status_ = ODE_STATUS_ILLEGAL_CONSTRAINTS;
				ParseOdeStatus();
				return this->status_;
			}

			// Check if the requested output time is already covered by the previous calculation
			if ( (direction_plus_ == true  && tOut_ < t_)  || 
				 (direction_plus_ == false && tOut_ > t_)  )
			{
				// The requested output time must be fall inside the last step calculated
				if ((direction_plus_ == true && tOut_ < t_ - hPreviousStep_) ||
					(direction_plus_ == false && tOut_ > t_ - hPreviousStep_))
				{
					std::cout << "Continuation procedure: illegal value of requested tOut: " << tOut_ << std::endl;

					if (direction_plus_ == true)
						std::cout << "The continuation procedure cis allowed only for tOut >= " << t_ - hPreviousStep_ << std::endl;
					else
						std::cout << "The continuation procedure cis allowed only for tOut <= " << t_ - hPreviousStep_ << std::endl;
				
					this->status_ = ODE_STATUS_ILLEGAL_CONTINUATION_REQUEST;
					ParseOdeStatus();
					return this->status_;
				}

				// Otherwise no calculations are needed. Only the interpolation is necessary
				{
					// Interpolate
					Interpolation();

					// Print
					if (printResults_ == true)
					{
						this->Equations(this->y_, tOut_, this->va_);
						this->Print(tOut_, this->y_);
					}

					return this->status_;
				}
			}

			// Check the strict constraints on the independent variable
			if (calculation_with_strict_tmax_ == true)
			{
				if ((direction_plus_ == true && t_ + this->hNextStep_ > tStrictMax_) ||
					(direction_plus_ == false && t_ + this->hNextStep_ < tStrictMax_))
				{
					this->hNextStep_ = tStrictMax_ - t_;
					this->odeHStatus_ = H_STATUS_DECREASED;
				}
			}
		}

		// Multivalue steps
		for (unsigned int k = 1; k <= max_number_steps_; k++)
		{
			if (stopOdeIntegration == true)
			{
				this->status_ = ODE_STATUS_EXCEPTION_HANDLING_STOP;
				ParseOdeStatus();
				return this->status_;
			}

			// Sets the error weight vector before each step
			CalculatesOneOverEpsilon(this->z_[0]);

			// Test for excessive accuracy requirement
			double tolSafe = OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE*OpenSMOKE::ErrorControl(this->z_[0], one_over_epsilon_);
			if (tolSafe > 1.)
			{
				this->status_ = ODE_STATUS_TOO_STRICT_TOLERANCES;
				ParseOdeStatus();
				return this->status_;
			}

			// Update the current time
			tInMeshPoint_ = t_;

			// Advances the solution of the ODE's by one integration step
			DriverStep();

			// If the norm of the first order derivatives is sufficiently small, the integration is completed with success
			if (this->status_ == ODE_STATUS_STOP_INTEGRATION_FOR_SMALL_YPRIME_NORM1)
			{
				this->y_ = this->z_[0];
				tOut_ = t_;
				return this->status_;
			}

			// Problems in integration step
			if (this->status_ <= 0)
			{
				this->y_.resize(this->ne_);
				t_ = tInMeshPoint_;
				ParseOdeStatus();
				return this->status_;
			}

			// Update the number of steps
			this->numberOfSteps_++;

			// Update the adopted maximum order (only for statistic purposes)
			if (this->maxOrderUsed_ < this->p_)
				this->maxOrderUsed_ = this->p_;

			// Update the adopted maximum step (only for statistic purposes)
			if (hMaxUsed_ < this->h_)
				hMaxUsed_ = this->h_;

			// Update the adopted minimum step (only for statistic purposes)
			if (hMinUsed_ > this->h_)
				hMinUsed_ = this->h_;

			// If the requested end time is reached exactly, the solver returns the solution
			if (tOut_ == t_)
			{
				this->y_ = this->z_[0];
				return this->status_;
			}

			// If the requested end time is overcome, the solver returns the interpolated solution
			if ((direction_plus_ == true && tOut_ < t_) || (direction_plus_ == false && tOut_ > t_))
			{
				// Interpolate
				Interpolation();

				// Print
				if (printResults_ == true)
				{
					this->Equations(this->y_, tOut_, this->va_);
					this->Print(tOut_, this->y_);
				}
				return this->status_;
			}
		}

		this->status_ = ODE_STATUS_MAX_NUMBER_OF_STEPS_REACHED;
		ParseOdeStatus();
		return this->status_;
	}

	template <typename Method>
	void MultiValueSolver<Method>::CalculatesOneOverEpsilon(const Eigen::VectorXd& w)
	{
		if (abs_tolerances_scalar_ == true && rel_tolerances_scalar_ == true)
		{
			for (unsigned int i = 0; i<this->ne_; i++)
				one_over_epsilon_(i) = abs_tolerance_ + rel_tolerance_*std::fabs(w(i));
		}
		else if (abs_tolerances_scalar_ == false && rel_tolerances_scalar_ == true)
		{
			for (unsigned int i = 0; i<this->ne_; i++)
				one_over_epsilon_(i) = abs_tolerances_(i) + rel_tolerance_*w(i);
		}
		else if (abs_tolerances_scalar_ == true && rel_tolerances_scalar_ == false)
		{
			for (unsigned int i = 0; i<this->ne_; i++)
				one_over_epsilon_(i) = abs_tolerance_ + rel_tolerances_(i)*w(i);
		}
		else if (abs_tolerances_scalar_ == false && rel_tolerances_scalar_ == false)
		{
			for (unsigned int i = 0; i<this->ne_; i++)
				one_over_epsilon_(i) = abs_tolerances_(i) + rel_tolerances_(i)*w(i);
		}

		for (unsigned int i = 0; i < this->ne_; i++)
			one_over_epsilon_(i) = 1. / one_over_epsilon_(i);
	}


	template <typename Method>
	double MultiValueSolver<Method>::CalculateInitialStepSize()
	{
		double hLowerBound;
		double hUpperBound;

		if (t0_ != 0.)
		{
			hLowerBound = 100.*OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE*std::fabs(t0_);
			hUpperBound = std::fabs(t0_);
		}
		else
		{
			hLowerBound = OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
			hUpperBound = 1.e-4;
		}

		if (calculation_with_strict_tmax_ == true)
		{
			const double delta = std::fabs(tStrictMax_ - t0_);
			if (hUpperBound > delta)
				hUpperBound = delta;
		}

		CalculatesOneOverEpsilon(y0_);
		double h0 = 0.5 / (OpenSMOKE::ErrorControl(this->f_, one_over_epsilon_) + OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT);

		if (h0 < hLowerBound)
			h0 = hLowerBound;

		if (h0 > hUpperBound)
			h0 = hUpperBound;

		// To take into account the possibility of negative direction
		h0 *= std::fabs(tOut_ - t0_) / (tOut_ - t0_);

		return h0;
	}

	template <typename Method>
	void MultiValueSolver<Method>::Interpolation()
	{
		const double hp = (tOut_ - t_) / hNordsieck_;
		double hr = 1.;

		this->y_ = this->z_[0];
		if (hp == 0.)
			return;

		for (unsigned int i = 1; i <= this->p_; i++)
		{
			hr *= hp;
			this->va_ = this->z_[i] * hr;
			this->y_ += this->va_;
		}

		// Check constraints on minimum and maximum values
		this->CheckConstraints(this->y_);
	}

	template <typename Method>
	void MultiValueSolver<Method>::ParseOdeStatus()
	{
		switch (this->status_)
		{
			case ODE_STATUS_MAX_NUMBER_OF_STEPS_REACHED:
				WarningMessage("The ODE system integration reached the maximum number of steps.");
				break;

			case ODE_STATUS_TOO_STRICT_TOLERANCES:
				WarningMessage("The required tolerances are too strict. Try to relax them!");
				break;

			case ODE_STATUS_ILLEGAL_CONTINUATION_REQUEST:
				WarningMessage("The continuation request cannot be performed because of illegal value.");
				break;

			case ODE_STATUS_MAX_NUMBER_ERRORTEST_FAILURES:
				WarningMessage("The maximum number of error test failures was reached.");
				break;

			case ODE_STATUS_MAX_NUMBER_CONVERGENCETEST_FAILURES:
				WarningMessage("The maximum number of convergence test failures was reached.");
				break;

			case ODE_STATUS_TOO_SMALL_STEP_SIZE:
				if (calculation_with_strict_tmax_ == false)
					WarningMessage("The current step is equal to the minimum allowed step.");
				break;

			case ODE_STATUS_ILLEGAL_MAX_INDEPENDENT_VARIABLE:
				WarningMessage("The requested maximum value of the indepepndent variable is illegal.");
				break;

			case ODE_STATUS_ILLEGAL_CONSTRAINTS:
				WarningMessage("The provided constraint on maximum and/or minimum values are not satisfied.");
				break;

			case ODE_STATUS_EXCEPTION_HANDLING_STOP:
				WarningMessage("The integration was stopped as requested by the user");
				break;
		}
	}


	template <typename Method>
	void MultiValueSolver<Method>::DriverStep()
	{
		this->iterConvergenceFailure_ = 0;
		this->iterErrorFailure_ = 0;
		this->convergenceStatus_ = CONVERGENCE_STATUS_OK;
		errorStatus_ = ERROR_STATUS_OK;

		// Reset the counter if the order or the step were changed
		if (this->odeOrderStatus_ != ORDER_STATUS_CONST || this->odeHStatus_ != H_STATUS_CONST)
			this->iterOrder_ = 0;

		this->p_ = this->orderInNextStep_;
		this->odeOrderStatus_ = ORDER_STATUS_CONST;

		// The scale factor is chosen according to the total number of steps
		// During the first steps, the time step can be increased a lot, because the first estimation of h0 could be excessively small
		// During the successive integration is better to reduce the possibility to increase too much the step
		hScaleMax_ = DEFAULT_HSCALE_MAX3;
		if (this->numberOfSteps_ < 10)
			hScaleMax_ = DEFAULT_HSCALE_MAX1;
		else if (this->numberOfSteps_ > 10 && this->numberOfSteps_ < 20)
			hScaleMax_ = DEFAULT_HSCALE_MAX2;

		while (1)
		{
			// If the step was modified, we have to update the z vector
			// See Buzzi-Ferraris, page 675
			if (this->odeHStatus_ != H_STATUS_CONST)
			{
				this->hScale_ = this->hNextStep_ / hNordsieck_;
				double hr = 1.;
				for (unsigned int i = 1; i <= this->p_; i++)
				{
					hr *= this->hScale_;
					this->z_[i] *= hr;
				}

				// Check Constraints
				this->CheckConstraints(this->z_[0]);

				// Since the step was changed, we reset the counter used to monitor the number of 
				// steps performed with the same order, to stabilize the algorithm, before changing the order
				this->iterOrder_ = 0;

				this->odeHStatus_ = H_STATUS_CONST;
			}

			// We update the current values of the step
			this->h_ = this->hNextStep_;
			hNordsieck_ = this->hNextStep_;

			// If the current step is too small, the integration is aborted
			if (std::fabs(this->h_) < OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT*this->min_step_size_)
			{
				this->status_ = ODE_STATUS_TOO_SMALL_STEP_SIZE;
				return;
			}

			// We update the current value of independent variable
			t_ = tInMeshPoint_ + this->h_;

			// If the current step is too small, the integration is aborted
			if (t_ == tInMeshPoint_)
			{
				this->status_ = ODE_STATUS_TOO_SMALL_STEP_SIZE;
				//if (calculation_with_strict_tmax_ == false)
				//	ParseOdeStatus();
				return;
			}

			// The user can ask to stop the simulation when the norm of first order derivative is sufficiently small
			// This is useful when the steady state solution has to be found
			if (stopIntegrationForSmallYPrimeNorm1_ == true)
			{
				// First order derivarives
				this->va_ = this->z_[1] / this->h_;

				// Norm of first order derivatives
				// See: http://wiki.ros.org/eigen/Troubleshooting
				const double norm1 = this->va_.template lpNorm<1>();

				// If the norm is sufficiently small, the integration is successfully completed
				if (norm1 < YPrimeNorm1_)
				{
					this->status_ = ODE_STATUS_STOP_INTEGRATION_FOR_SMALL_YPRIME_NORM1;
					return;
				}
			}

			// Core function of the multivalue algorithm
			// The integration is advanced over the current step: from t_ to t_+h_
			MultiValueStep();

			// If the method to find the correction failed
			if (this->convergenceStatus_ == CONVERGENCE_STATUS_FAILURE)
			{
				// Update the counters
				this->iterConvergenceFailure_++;
				numberOfConvergenceFailure_++;	// cumulative number of convergence failures

				// If too many convergence failures were observed, the integration is aborted
				if (this->iterConvergenceFailure_ >= this->MAX_CONVERGENCE_FAILURE)
				{
					this->status_ = ODE_STATUS_MAX_NUMBER_CONVERGENCETEST_FAILURES;
					return;
				}

				// If the current integration step is too small, the integration is aborted
				if (std::fabs(this->h_ - this->min_step_size_) <= OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE*std::fabs(this->min_step_size_))
				{
					this->status_ = ODE_STATUS_TOO_SMALL_STEP_SIZE;
					return;
				}

				// Find a proper solution to face the convergence failure
				// It depends on the specific multivalue method
				this->WhatToDoInCaseOfConvergenceFailure(tInMeshPoint_, tStabilize_);
			}
			else
			{
				numberOfConvergenceSuccess_++;

				// Buzzi-Ferraris eq. 29.171: hScaleCurrent = E(p)*v(p+1)/eps
				// where eps = tolAbs + tolRel*y
				double hScaleCurrent = this->Ep_[this->p_] * OpenSMOKE::ErrorControl(this->v_[this->p_ + 1], one_over_epsilon_);
			
				// The step must be decreased if the error is too large
				if (hScaleCurrent > 1.)
				{
					errorStatus_ = ERROR_STATUS_FAILURE;	// The adopted step was too large
					this->iterErrorFailure_++;					// The local number of failures is updated
					numberOfErrorCheckFailure_++;			// The cumulative number of failures is updated

					// If the local number of failures is too large, the integration is aborted
					// This means that the step was decreased too much times
					if (this->iterErrorFailure_ >= MAX_ERROR_FAILURE)
					{
						this->status_ = ODE_STATUS_MAX_NUMBER_ERRORTEST_FAILURES;
						return;
					}

					// If the current time step is too small, it cannotbe further decreased
					// Therefore the integration is aborted
					if (std::fabs(this->h_ - this->min_step_size_) <= OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE*std::fabs(this->min_step_size_))
					{
						this->status_ = ODE_STATUS_TOO_SMALL_STEP_SIZE;
						return;
					}

					// The new step must be evaluated.
					// At the same time a possible reduction of the order is evaluated

					// 1. Method of order p
					// The scaling factor is evaluated according to Eq. 29.185: pow(eps/E(p)*z(p+1), 1/p+1)
					// A safety coefficient equal to 0.8 is included
					hScaleMax_ = 1.;
					hScaleCurrent = 0.8*std::pow(1. / hScaleCurrent, 1. / double(this->p_ + 1));
					this->hScale_ = hScaleCurrent;

					// 2. Method of order p-1 (if possible)
					// The scaling factor is evaluated according to Eq. 29.184: pow(eps/E(p-1)*z(p), 1/p)
					// A safety coefficient equal to 0.8 is included
					{
						double hScaleMinus = 0.;
						if (this->p_ > this->MIN_ORDER)
						{
							hScaleMinus = this->Ep_[this->p_ - 1] * OpenSMOKE::ErrorControl(this->v_[this->p_], one_over_epsilon_);
							hScaleMinus = 1. / (hScaleMinus + OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT);
							hScaleMinus = 0.8*std::pow(hScaleMinus, 1. / double(this->p_));

							// If the step corresponding to the lower order is more convenient, the order of the method is decreased
							if (hScaleMinus > hScaleCurrent)
							{
								this->hScale_ = hScaleMinus;
								this->p_ -= 1;
								this->odeOrderStatus_ = ORDER_STATUS_DECREASED;
							}
						}
					}

					// The new step is between 0.001 and 0.5 times the previous (failed) step
					this->hScale_ = std::min(this->hScale_, 0.5);
					this->hScale_ = std::max(this->hScale_, this->min_step_size_ / std::fabs(this->h_));
					this->hScale_ = std::max(this->hScale_, 1.e-3);

					if (this->p_ == 1 && this->hScale_ < 0.1 && this->numberOfSteps_ > 10 && tInMeshPoint_ != tStabilize_)
					{
						this->Equations(this->z_[0], tInMeshPoint_, this->va_);
						this->numberOfFunctionCalls_++;

						this->z_[1] = this->va_*this->h_;
						tStabilize_ = tInMeshPoint_;
					}

					// In any case, the step is decreased
					this->odeHStatus_ = H_STATUS_DECREASED;
					this->hNextStep_ *= this->hScale_;

					// Additional safety factors (they should be better studied)
					if (cpuTimeEquationSystem_ <= 0.1  &&  3.*cpuTimeEquationSystem_<cpuTimeToFactorize_  &&  cpuTimeToFactorize_>0.1)
						this->hNextStep_ *= .9;

					// Additional safety factors (they should be better studied)
					if (100.*cpuTimeEquationSystem_ < cpuTimeToFactorize_ && cpuTimeToFactorize_ > 1.)
						this->hNextStep_ *= .9;

					// The cumulative number of decreasing step operations is updated
					this->numberOfDecreasedSteps_++;
				}

				// The step is OK because the error is sufficiently small
				// We have to check if it is convenient to change the order
				else
				{
					errorStatus_ = ERROR_STATUS_OK;
					numberOfErrorCheckSuccess_++;
				
					hNordsieck_ = this->h_;
					hPreviousStep_ = this->h_;
					this->hNextStep_ = this->h_;
					this->orderInNextStep_ = this->p_;
					this->odeOrderStatus_ = ORDER_STATUS_CONST;
					this->odeHStatus_ = H_STATUS_CONST;

					// First of all we analyze the convergence rate
					this->ConvergenceRate();

					// If the integration has to respect strictly the maximum t, we have to impose proper contraints on the next step
					if (calculation_with_strict_tmax_ == true)
					{
						if (direction_plus_ == true && t_ + this->h_ > tStrictMax_)
							this->hNextStep_ = tStrictMax_ - t_;

						if (direction_plus_ == false && t_ + this->h_ < tStrictMax_)
							this->hNextStep_ = tStrictMax_ - t_;

						if (this->h_ != this->hNextStep_)
						{
							if (std::fabs(this->hNextStep_) > std::fabs(this->h_))
								this->odeHStatus_ = H_STATUS_INCREASED;
							else
								this->odeHStatus_ = H_STATUS_DECREASED;
						}
					}

					// If failures were observed during the convergence procedure, it is better to maintain the same order
					// Only the counters are reset
					if (this->iterErrorFailure_ > 0 || this->iterConvergenceFailure_ > 0)
					{
						this->iterErrorFailure_ = 0;
						this->iterConvergenceFailure_ = 0;
					}
					// Otherwise we can try to change the order
					else
					{
						this->iterOrder_++;
						// The order is changed only if a sufficient number of steps is performed
						// In particular, for order p it is necessary to wait at least for p+1 steps in order to change the order
						if (this->iterOrder_ > this->p_ + 1)
							NewOrderNewH();
					}

					// We can update the solution on the current mesh point
					for (unsigned int i = 0; i <= this->p_ + 2; i++)
						this->z_[i] = this->v_[i];

					this->CheckConstraints(this->z_[0]);

					// Update the counter
					if (this->odeOrderStatus_ == ORDER_STATUS_CONST && this->odeHStatus_ == H_STATUS_CONST)
						numberOfStepsWithoutChanges_++;

					// Print
					if (printResults_ == true)
					{
						if (direction_plus_ == true && t_<=tOut_)
							this->Print(t_, this->z_[0]);
						if (direction_plus_ == false && t_>=tOut_)
							this->Print(t_, this->z_[0]);		
					}

					break;
				}
			}

		} // end while(1)
	}

	template <typename Method>
	void MultiValueSolver<Method>::NewOrderNewH()
	{
		double hScaleMinus = 0.;
		double hScaleCurrent = 0.;
		double hScalePlus = 0.;
		double alfa1 = 0.;
		double alfa3 = 0.;

		// Order p-1
		if (this->p_ > this->MIN_ORDER && this->odeOrderStatus_ != ORDER_STATUS_INCREASED)
		{
			// Estimate the error e corresponding to order p-1
			hScaleMinus = this->Ep_[this->p_ - 1] * OpenSMOKE::ErrorControl(this->v_[this->p_], one_over_epsilon_) + OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
			hScaleMinus = 1. / hScaleMinus;
		
			// Estimate the new step using a safetyfactor (alfa1)
			alfa1 = this->alfa2_[this->p_ - 1] - this->deltaAlfa1_;
			hScaleMinus = alfa1*std::pow(hScaleMinus, 1. / double(this->p_));
		}

		// Order p
		{
			// Estimate the error e corresponding to order p
			hScaleCurrent = this->Ep_[this->p_] * OpenSMOKE::ErrorControl(this->v_[this->p_ + 1], one_over_epsilon_) + OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
			hScaleCurrent = 1. / hScaleCurrent;

			// Estimate the new step using a safetyfactor (alfa2)
			hScaleCurrent = this->alfa2_[this->p_] * std::pow(hScaleCurrent, 1. / double(this->p_ + 1));
		}

		// Order p+1
		if (this->p_ < this->maximum_order_ &&	this->odeOrderStatus_ != ORDER_STATUS_DECREASED)
		{
			// Estimate the error e corresponding to order p+1
			hScalePlus = this->Ep_[this->p_ + 1] * OpenSMOKE::ErrorControl(this->v_[this->p_ + 2], one_over_epsilon_) + OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
			hScalePlus = 1. / hScalePlus;

			// Estimate the new step using a safetyfactor (alfa3)
            alfa3 = this->alfa2_[this->p_ + 1] - this->deltaAlfa3_;
			hScalePlus = alfa3*std::pow(hScalePlus, 1. / double(this->p_ + 2));
		}

		// Choose among the different orders
		// 1. Choosing the order p-1
		if (hScaleMinus > hScaleCurrent && hScaleMinus > hScalePlus)
		{
			this->hScale_ = (.75 / alfa1)*hScaleMinus;
			this->orderInNextStep_ = this->p_ - 1;
			this->odeOrderStatus_ = ORDER_STATUS_DECREASED;
		}
		// 2. Choosing the order p
		else if (hScaleCurrent > hScaleMinus && hScaleCurrent > hScalePlus)
		{
			this->hScale_ = (.75 / this->alfa2_[this->p_])*hScaleCurrent;
			this->orderInNextStep_ = this->p_;
			this->odeOrderStatus_ = ORDER_STATUS_CONST;
		}
		// 3. Choosing the order p+1
		else
		{
			this->hScale_ = (.75 / alfa3)*hScalePlus;
			this->orderInNextStep_ = this->p_ + 1;
			this->odeOrderStatus_ = ORDER_STATUS_INCREASED;
		}

		// If the integration has to be performed with the minimum order
		if (fixedOrderMin_ == true)
		{
			iterMinimumOrder_++;
			if (iterMinimumOrder_ < 2000)
			{
				this->hScale_ = (.75 / this->alfa2_[this->p_])*hScaleCurrent;
				this->orderInNextStep_ = this->p_;
				this->odeOrderStatus_ = ORDER_STATUS_CONST;
			}
			else
				fixedOrderMin_ = false;
		}

		// Evaluate if changing the order is really computationally and advantage
		{
			// In normal conditions the order and the step are changed only if the new step is at least 1.2 times larger than the previous
			// This is done in order to avoid too change too many times the step and factorize the Jacobian matrix
			double scl = 1.2;

			// This is a situation in which the factorization is not so heavy
			if (cpuTimeEquationSystem_ < 0.1 && cpuTimeEquationSystem_ < cpuTimeToFactorize_ / 3. && cpuTimeToFactorize_ > 0.1)
				scl = 2.;

			// In this case the factorization is very expensive and therefore it is better to keep the step and the order constant, in order
			// to avoid to factorize the Jacobian matrix
			if (cpuTimeEquationSystem_ < cpuTimeToFactorize_ / 100. && cpuTimeToFactorize_ > 1. && iterOfConstStepsForLargeFactorizationTime_ <= 100)
			{
				iterOfConstStepsForLargeFactorizationTime_++;
				scl = 12.;
				if (this->numberOfSteps_ > 10)
					hScaleMax_ = DEFAULT_HSCALE_MAX2;
			}

			// 
			if (iterOfConstStepsForLargeFactorizationTime_ >= 100)
			{
				scl = 2.;
				hScaleMax_ = DEFAULT_HSCALE_MAX3;
			}

			// This means that it is not so convenient to change the order and/or the step
			if (this->hScale_ < scl)
			{
				iterOfConstStepsForLargeFactorizationTime_ = 0;
				iterConstSteps_++;
				this->hScale_ = 1.;
				this->odeHStatus_ = H_STATUS_CONST;
				this->orderInNextStep_ = this->p_;
				this->odeOrderStatus_ = ORDER_STATUS_CONST;
				this->hNextStep_ = this->h_;
			}
			// This means that it is very convenient to chenge the step (and/or the order)
			else
			{
				iterConstSteps_ = 0;
				this->hScale_ = std::min(this->hScale_, hScaleMax_);
				this->odeHStatus_ = H_STATUS_INCREASED;
				this->hNextStep_ = this->h_*this->hScale_;
				if (this->hNextStep_ > max_step_size_)
					this->hNextStep_ = max_step_size_;
			}
		}

		// If more than 500 steps were done without changing the step size or
		// the the time step is very small it is better ho reduce the step and force the solver to use the minimum order
		if (iterConstSteps_ > 500 || (this->numberOfSteps_ % 1000 == 0 && this->numberOfSteps_ != 0 && (std::fabs(this->h_) < 1.e-7 || std::fabs(this->h_) < 1.e-5*t_)))
		{
			if (this->p_ > this->MIN_ORDER)
			{
				fixedOrderMin_ = true;
				iterMinimumOrder_ = 0;
				this->orderInNextStep_ = this->MIN_ORDER;
				this->odeOrderStatus_ = ORDER_STATUS_DECREASED;
				this->hScale_ = 0.5;
				this->hNextStep_ = this->h_*this->hScale_;
				this->odeHStatus_ = H_STATUS_DECREASED;
			}
		}


		// Apply the proper constraint if the user asked a strict control on the maximum value of the independent variable
		if (calculation_with_strict_tmax_ == true)
		{
			if (direction_plus_ == true && t_ + this->hNextStep_ > tStrictMax_)
				this->hNextStep_ = tStrictMax_ - t_;

			if (direction_plus_ == false && t_ + this->hNextStep_ < tStrictMax_)
				this->hNextStep_ = tStrictMax_ - t_;

			if (this->h_ != this->hNextStep_)
			if (std::fabs(this->hNextStep_) > std::fabs(this->h_))
				this->odeHStatus_ = H_STATUS_INCREASED;
			else
				this->odeHStatus_ = H_STATUS_DECREASED;
		}

		if (this->odeHStatus_ == H_STATUS_INCREASED)
			this->numberOfIncreasedSteps_++;

		if (this->odeHStatus_ == H_STATUS_DECREASED)
			this->numberOfDecreasedSteps_++;
	}

	template <typename Method>
	void MultiValueSolver<Method>::MultiValueStep()
	{
		// Prevision (eq. 29.149): v = Dz
		// Please note that the D matrix is not strictly needed, since only sums of vectors can be adopted as
		// suggested by Buzzi-Ferraris at page 669
		{
			for (unsigned int i = 0; i <= this->p_; i++)
			this->v_[i] = this->z_[i];

			for (int i = 0; i < int(this->p_); i++)
			for (int j = int(this->p_ - 1); j >= i; j--)
				this->v_[j] += this->v_[j + 1];
		}

		// Check constraints for minimum and maximum values
		{
			// Minimum constraints
			if (this->min_constraints_ == true)
			{
				for (unsigned int j = 0; j < this->ne_; j++)
					if (this->v_[0](j) < this->min_values_(j))
					{
						if (std::fabs(this->v_[0](j) - this->min_values_(j))*one_over_epsilon_(j) < double(this->ne_))
						{
							this->v_[0](j) = this->min_values_(j);
							this->v_[1](j) = this->min_values_(j) - this->z_[0](j);
							for (unsigned int i = 2; i <= this->p_; i++)
								this->v_[i](j) = 0.;
						}
					}
			}
		
			// Maximum constraints
			if (this->max_constraints_ == true)
			{
				for (unsigned int j = 0; j < this->ne_; j++)
				if (this->v_[0](j) > this->max_values_(j))
				{
					if (std::fabs(this->v_[0](j) - this->max_values_(j))*one_over_epsilon_(j) < double(this->ne_))
					{
						this->v_[0](j) = this->max_values_(j);
						this->v_[1](j) = this->max_values_(j) - this->z_[0](j);
						for (unsigned int i = 2; i <= this->p_; i++)
							this->v_[i](j) = 0.;
					}
				}
			}
		}

		// According to the particular multivalue (for example the Gear method) the b vector (i.e. the correction)
		// is estimated (for the Gear method this is done using the Newton's method)
		this->iterConvergence_ = this->FindCorrection(t_, one_over_epsilon_);

		// If the norm of first order derivatives is sufficiently small, the integration is completed
		if (this->status_ == ODE_STATUS_STOP_INTEGRATION_FOR_SMALL_YPRIME_NORM1)
			return;

		// If the correction was not found, the failure is managed by the DriverStep function
		if (this->convergenceStatus_ == CONVERGENCE_STATUS_FAILURE)
			return;

		// Correction to be applied to the z vector: z0 = v0 + r0*b
		this->va_ = this->b_*this->r_[0](this->p_ - 1);

		// Check constraints for minimum and maximum values
		{
			// Minimum constraints
			if (this->min_constraints_ == true)
			{
				for (unsigned int j = 0; j < this->ne_; j++)
				if (this->v_[0](j) + this->va_(j) < this->min_values_(j))
				{
					this->b_(j) = this->min_values_(j) - this->z_[0](j) - this->v_[1](j);
					this->v_[0](j) = this->min_values_(j) - this->va_(j);
					if (this->v_[0](j) < this->min_values_(j))
						this->v_[0](j) = this->min_values_(j);
				}
			}

			// Maximum constraints
			if (this->max_constraints_ == true)
			{
				for (unsigned int j = 0; j < this->ne_; j++)
				if (this->v_[0](j) + this->va_(j) > this->max_values_(j))
				{
					this->b_(j) = this->max_values_(j) - this->z_[0](j) - this->v_[1](j);
					this->v_[0](j) = this->max_values_(j) - this->va_(j);
					if (this->v_[0](j) > this->max_values_(j))
						this->v_[0](j) = this->max_values_(j);
				}
			}
		}

		// Correction to be applied to the zj vector: z1 = vj + rj*b (r1 is equal to 1)
		this->v_[0] += this->va_;
		this->v_[1] += this->b_;
		for (unsigned int i = 2; i <= this->p_; i++)
		{
			this->va_ = this->b_*this->r_[i](this->p_ - 1);
			this->v_[i] += this->va_;
		}

		// Nordsieck component p + 1 (useful for estimating the error)
		// See equation 29.169: z[p+1] ~= ( z[p](n) - z[p](n-1) ) / (p+1) 
		// Be careful: this->v_ is the updated vector (i.e. z at t(n) and z at t(n-1) )
		this->va_ = this->v_[this->p_] - this->z_[this->p_];
		this->v_[this->p_ + 1] = this->va_ / (double(this->p_) + 1.);

		// Nordsieck component p + 2 (useful for estimating the error)
		// See equation 29.169: z[p+2] ~= ( z[p+1](n) - z[p+1](n-1) ) / (p+2) 
		if (this->iterOrder_ >= this->p_)
		{
			this->va_ = this->v_[this->p_ + 1] - this->z_[this->p_ + 1];
			this->v_[this->p_ + 2] = this->va_ / (double(this->p_) + 2.);
		}
	}

	template <typename Method>
	void MultiValueSolver<Method>::OdeSummary(std::ostream& out)
	{
		if (printResults_ == true)
		{
			out << std::endl;
			out << "ODE system solution" << std::endl;
			out << "---------------------------------------------------------------------------------------" << std::endl;
			out << "* Size of the first step:                     " << first_step_size_ << std::endl;
			out << "* Minimum step size:                          " << hMinUsed_ << std::endl;
			out << "* Maximum step size:                          " << hMaxUsed_ << std::endl;
			out << "* Number of successive steps with same size:  " << numberOfStepsWithoutChanges_ << std::endl;
			out << "* Number of convergence failures:             " << numberOfConvergenceFailure_ << std::endl;
			out << "* Number of convergence successes:            " << numberOfConvergenceSuccess_ << std::endl;
			out << "* Number of failures in error checking:       " << numberOfErrorCheckFailure_ << std::endl;
			out << "* Number of successes in error checking:      " << numberOfErrorCheckSuccess_ << std::endl;
			out << "---------------------------------------------------------------------------------------" << std::endl;
			out << std::endl;

			this->OdeMethodSummary(out);
		}
	
	}
}
