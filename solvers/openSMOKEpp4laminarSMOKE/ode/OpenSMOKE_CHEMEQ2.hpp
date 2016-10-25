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
|   Copyright(C) 2016, 2015  Alberto Cuoci                                |
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
	double SIGN(const double A, const double B)
	{
		return (B >= 0.) ? std::fabs(A) : -std::fabs(A);
	}

	double MAX(const double A, const double B, const double C)
	{
		return std::max(A, std::max(B, C));
	}

	double MIN(const double A, const double B, const double C)
	{
		return std::min(A, std::min(B, C));
	}

	template <typename T>
	const double OpenSMOKE_CHEMEQ2<T>::rswitch_ = 5.965900;

	template <typename T>
	OpenSMOKE_CHEMEQ2<T>::OpenSMOKE_CHEMEQ2(T* odeSystem)
	{
		this->SetDefaultValues();

		number_system_calls_ = 0;
		number_steps_ = 0;
		iteration_ = 0;
		failures_ = 0;

		alpha_option_ = 1;
		high_level_accuracy_ = false;

		epsmin_ = 1e-2;
		sqreps_ = 0.50;
		epscl_  = 1.e2;
		epsmax_ = 10.;
		max_number_subiterations_ = 1;
		tfd_ = 1.000008;
		dtmin_ = 1.e-15;

		this->odeSystem_ = odeSystem;
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::AnalyzeUserOptions()
	{
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::MemoryAllocation(const int n)
	{
		this->n_ = n;		// Number of equations

		q_.resize(this->n_);
		d_.resize(this->n_);
		this->y0_ = new double[this->n_];
		this->y_ = new double[this->n_];
		y1_.resize(this->n_);
		ys_.resize(this->n_);
		qs_.resize(this->n_);
		scrarray_.resize(this->n_);
		rtau_.resize(this->n_);
		rtaus_.resize(this->n_);
		ymin_.resize(this->n_);
		

		for (int i = 0; i<this->n_; i++)
		{
			q_[i] = 0.;
			d_[i] = 0.;
			this->y0_[i] = 0.;
			y1_[i] = 0.;
			rtau_[i] = 0.;
			rtaus_[i] = 0.;
			ymin_[i] = 1.e-20;
			this->y_[i] = 0.;

			ys_[i] = 0.;
			qs_[i] = 0.;
			scrarray_[i] = 0.;
		}
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::SetEpsilonMinimum(const double epsmin)
	{
		if (epsmin > 0.)
		{
			epsmin_ = epsmin;
			sqreps_ = 5.*std::sqrt(epsmin_);
			epscl_ = 1. / epsmin_;
		}
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::SetEpsilonMaximum(const double epsmax)
	{
		if (epsmax > 0.)
			epsmax_ = epsmax;
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::SetMinimumStep(const double dtmin)
	{
		if (dtmin > 0.)
			dtmin_ = dtmin;
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::SetMaximumNumberSubIterations(const int max_number_subiterations)
	{
		if (max_number_subiterations > 0)
			max_number_subiterations_ = max_number_subiterations;
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::SetMinimumValue(const double ymin)
	{
		if (ymin > 0.)
		{
			for (int i = 0; i < this->n_; i++)
				ymin_[i] = ymin;
		}
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::SetHighLevelAccuracy()
	{
		high_level_accuracy_ = true;

		ym1_.resize(this->n_);
		ym2_.resize(this->n_);

		for (int i = 0; i < this->n_; i++)
		{
			ym1_[i] = 0.;
			ym2_[i] = 0.;
		}
	}

	template <typename T>
	double OpenSMOKE_CHEMEQ2<T>::GetAlpha(const double rtaui)
	{
		// One of two approximations for alpha is chosen, see (Mott et al., 2000):
		// 1) Pade b for all rtaui see (Mott et al., 2000)
		// 2) Pade a for rtaui<=rswitch,
		//    linear approximation for rtaui > rswitch

		if (alpha_option_ == 1)	// Option 1) Pade b
		{

			return (180. + rtaui*(60. + rtaui*(11. + rtaui))) /
				(360. + rtaui*(60. + rtaui*(12. + rtaui)));
		}
		else // Option 2) Pade a or linear
		{
			if (rtaui <= rswitch_)
			{
				return (840. + rtaui*(140. + rtaui*(20. + rtaui))) /
					(1680. + 40. * rtaui*rtaui);
			}
			else
			{
				return (1. - 1. / rtaui);
			}
		}
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::Solve(const double xend)
	{
		if (high_level_accuracy_ == true && max_number_subiterations_ < 3)
			ErrorMessage("Solve", "At least 3 sub-iterations are required when the high level of accuracy option is turned on.");

		this->xend_ = xend;
		const double dt_global = this->xend_ - this->x0_;

		// Current value of the independent variable relative to the start of the global timestep
		double tn = 0.;

		// Initialize and limit y to the minimum value
		for (int i = 0; i < this->n_; i++)
		{
			q_[i] = 0.;
			d_[i] = 0.;
			this->y_[i] = std::max(this->y0_[i], ymin_[i]);
		}

		// Evaluates the production and consumption contributions
		this->odeSystem_->Equations(this->y_, q_, d_, tn + this->x0_);
		number_system_calls_++;

		// Estimate the initial stepsize
		// Strongly increasing functions(q >> > d assumed here) use a step -
		// size estimate proportional to the step needed for the function to
		// reach equilibrium where as functions decreasing or in equilibrium
		// use a stepsize estimate directly proportional to the characteristic 
		// stepsize of the function.convergence of the integration
		// scheme is likely since the smallest estimate is chosen for the
		// initial stepsize.

		double scrtch = 1.e-25;
		for (int i = 0; i < this->n_; i++)
		{
			const double ascr = std::fabs(q_[i]);
			const double scr2 = SIGN(1. / this->y_[i], 0.1*epsmin_*ascr - d_[i]);
			const double scr1 = scr2*d_[i];
				  double temp = -std::fabs(ascr - d_[i])*scr2;

			// If the species is already at the minimum, disregard destruction when calculating step size
			if (this->y_[i] == ymin_[i]) temp = 0.;

			scrtch = std::max(scr1, std::max(temp, scrtch));
		}

		double dt = std::min(sqreps_ / scrtch, dt_global);

		// Loop to cover the whole integration interval
		while(1)
		{
			// Update number of time steps
			number_steps_++;

			// Independent variable at the start of the chemical timestep
			double ts = tn;

			// Starting values are stored
			for (int i = 0; i < this->n_; i++)
			{
				rtau_[i] = dt*d_[i] / this->y_[i];
				ys_[i] = this->y_[i];
				qs_[i] = q_[i];
				rtaus_[i] = rtau_[i];
			}

			// Find the predictor terms
			predictor:

				// Prediction
				for (int i = 0; i < this->n_; i++)
				{
					const double rtaui = rtau_[i];
					const double alpha = GetAlpha(rtaui);

					scrarray_[i] = (q_[i] - d_[i]) / (1.0 + alpha*rtaui);
				}

				iteration_ = 1;

				while (iteration_ <= max_number_subiterations_)
				{
					// Limit decreasing functions to their minimum values
					for (int i = 0; i < this->n_; i++)
					{
						if (high_level_accuracy_ == true)
						{
							ym2_[i] = ym1_[i];
							ym1_[i] = this->y_[i];
						}
						this->y_[i] = std::max(ys_[i] + dt*scrarray_[i], ymin_[i]);
					}

					// Removed from original algorithm so that previous, rather than first, 
					// corrector is compared to. Results in faster integration.
					if (iteration_ == 1)
					{
						// The first corrector step advances the time (tentatively) and
						// saves the initial predictor value as y1 for the timestep check later.
						tn = ts + dt;
						for (int i = 0; i < this->n_; i++)
							y1_[i] = this->y_[i];
					}

					// Evaluate the derivatives for the corrector step
					this->odeSystem_->Equations(this->y_, q_, d_, tn + this->x0_);
					number_system_calls_++;
					eps_ = 1.e-10;

					for (int i = 0; i < this->n_; i++)
					{
						// dt times average p from Eq. (37)
						const double rtaub = 0.5 * (rtaus_[i] + dt*d_[i] / this->y_[i]);
						const double alpha = GetAlpha(rtaub);

						const double qt = qs_[i] * (1. - alpha) + q_[i] * alpha;
						const double pb = rtaub / dt;

						scrarray_[i] = (qt - ys_[i] * pb) / (1.0 + alpha*rtaub);
					}

					iteration_++;
				
				} // end while 

				// Calculate new f, check for convergence, and limit decreasing functions
				// The order of operations in this loop is important
				for (int i = 0; i < this->n_; i++)
				{
					const double scr2 = std::max(ys_[i] + dt*scrarray_[i], 0.);
					double scr1 = std::fabs(scr2 - y1_[i]);
					this->y_[i] = std::max(scr2, ymin_[i]);

					if (high_level_accuracy_ == true)
					{
						ym2_[i] = ym1_[i];
						ym1_[i] = this->y_[i];
					}

					if (0.25*(ys_[i] + this->y_[i]) > ymin_[i])
					{
						scr1 = scr1 / this->y_[i];
						eps_ = std::max(0.5*(scr1 + std::min(std::fabs(q_[i] - d_[i]) / (q_[i] + d_[i] + 1.e-30), scr1)), eps_);
					}
				}
				
				eps_ *= epscl_;

				// Print diagnostics if step size becomes too small
				if (dt <= dtmin_ + 1.e-16*tn)
				{
					std::cout << " * Current step:                 " << dt << std::endl;
					std::cout << " * Current independent variable: " << tn << std::endl;
					std::cout << " * Minimum step:                 " << dtmin_ << std::endl;

					std::cout << " * Details" << std::endl;
					for (int i = 0; i < this->n_; i++)
					{
						const double dtc = epsmin_*this->y_[i] / (std::fabs(q_[i] - d_[i]) + 1.e-30);
						
						std::cout << "   " << i << " " << q_[i] << " " << d_[i] << " " << this->y_[i] << " " << dtc << std::endl;
					}

					ErrorMessage("Solve", "Time step is too small");
				}

				// Check for convergence
			
				// The following section is used for the stability check
				double stab = 0.;	// see Eqs. (52) and (53)
				if (high_level_accuracy_ == true)
				{
					stab = 0.01;
					if (max_number_subiterations_ >= 3)
					{
						for (int i = 0; i < this->n_; i++)
						{
							stab = std::max(stab, std::fabs(this->y_[i] - ym1_[i]) /
									(std::fabs(ym1_[i] - ym2_[i]) + 1.e-20*this->y_[i]));
						}
					}
				}

				// Valid step: return if delta_global has been reached
				if (eps_ <= epsmax_ && stab <= 1.)
				{
					if (dt_global <= tn*tfd_)
					{
						// Call the print function
						this->odeSystem_->GetWriteFunction(tn + this->x0_, this->y_);

						this->x0_ = this->x_;
						memcpy(this->y0_, this->y_, this->n_ * sizeof(double));

						return;
					}
				}
				// Invalid step: reset tn to ts
				else
				{
					tn = ts;
				}

				// Perform step size modifications
				// Estimate sqrt(eps) by Newton's iterations
				double rteps = 0.5*(eps_ + 1.0);
				rteps = 0.5*(rteps + eps_ / rteps);
				rteps = 0.5*(rteps + eps_ / rteps);

				double dto = dt;
				if (high_level_accuracy_ == true)
					dt = MIN(dt*(1./rteps+0.005), tfd_*(dt_global-tn), dto/(stab+0.001));
				else
					dt = std::min(dt*(1./rteps+0.005), tfd_*(dt_global - tn));

				// Begin a new step if previous step converged
				if ( (eps_ > epsmax_) || (stab > 1.))
				{
					failures_++;

					// After an unsuccessful step the initial timescales don't
					// change, but dt does, requiring rtaus to be scaled by the
					// ratio of the new and old timesteps.
					dto = dt / dto;
					for (int i = 0; i < this->n_; i++)
						rtaus_[i] *= dto;

					// Unsuccessful steps return to predictor section so that the initial source
					// terms do not get recalculated
					goto predictor;
				}

				// Successful step: get the source terms for the next step and continue
				this->odeSystem_->Equations(this->y_, q_, d_, tn + this->x0_);
				number_system_calls_++;

				// Call the print function
				this->odeSystem_->GetWriteFunction(tn + this->x0_, this->y_);
		}

		this->x0_ = this->x_;
		memcpy(this->y0_, this->y_, this->n_ * sizeof(double));
	}

	template <typename T>
	void OpenSMOKE_CHEMEQ2<T>::Status() const
	{
		std::cout << "CHEMEQ2 Status:                     " << std::endl;				// Status
		std::cout << " * Number of steps:                 " << number_steps_ << std::endl;			// Number of steps taken for the problem so far 
		std::cout << " * Number of function evaluations:  " << number_system_calls_ << std::endl;	// Number of f evaluations for the problem so far.
	}

	template <typename T>
	std::string OpenSMOKE_CHEMEQ2<T>::Tag() const
	{
		return "CHEMEQ2";
	}

	template <typename T>
	int OpenSMOKE_CHEMEQ2<T>::GetNumberOfSteps() const
	{
		return number_steps_;
	}

	template <typename T>
	int OpenSMOKE_CHEMEQ2<T>::GetNumberOfFunctionEvaluations() const
	{
		return number_system_calls_;
	}

	template <typename T>
	int OpenSMOKE_CHEMEQ2<T>::GetNumberOfJacobianEvaluations() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_CHEMEQ2<T>::GetNumberOfLUFactorizations() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_CHEMEQ2<T>::GetNumberOfNonLinearIterations() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_CHEMEQ2<T>::GetLastOrderUsed() const
	{
		return 0;
	}

	template <typename T>
	double OpenSMOKE_CHEMEQ2<T>::GetLastStepUsed() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_CHEMEQ2<T>::GetNumberOfConvergenceFailures() const
	{
		return 0;
	}

	template <typename T>
	int OpenSMOKE_CHEMEQ2<T>::GetNumberOfErrorTestFailures() const
	{
		return 0;
	}

	template <typename T>
	OpenSMOKE_CHEMEQ2<T>::~OpenSMOKE_CHEMEQ2(void)
	{
	}
}
