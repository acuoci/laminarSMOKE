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

#ifndef OpenSMOKE_CHEMEQ2_H
#define OpenSMOKE_CHEMEQ2_H

#include "math/external-ode-solvers/OpenSMOKE_OdeSystemSolver.h"

namespace OpenSMOKE
{
	//!  A class to integrate "stiff" equations similar to chemical equilibrium equations.
	/*!
	
	The source code below is taken from:
	Mott, D.R. and Oran, E.S., 2001
	"CHEMEQ2: A Solver for the Stiff Ordinary Differential Equations of Chemical Kinetics"
	Naval Research Laboratory, NRL/MR/6400-01-8553.  
	The report documentation page under distribution/availability statement states:
	"Approved for public release; distribution is unlimited."
	
	CHEMEQ2 solves a class of "stiff" ODEs associated with reactive flow problems that cannot be readily
	solved by the standard classical methods. In contrast to the original chemeq subroutine, this version uses the same
	quasi-steady-state update for every species regardless of the timescale for that species. An adaptive stepsize is chosen to
	give accurate results for the fastest changing quantity, and a stability check on the timestep is also available when the
	corrector is iterated.

	Original CHEMEQ development:
	 + originators: T.R. Young NRL 1982
	 + vax version: T.R. young NRL code 4040 May 1983
	 + workstation: G. Patnaik Berkeley Research Jun 1995

	CHEMEQ2 development: D.R. Mott NRL code 6404 may 1999

	Conversion to C++: Alberto Cuoci, Politecnico di Milano, May 2016
	*/

	template <typename T>
	class OpenSMOKE_CHEMEQ2 : public OpenSMOKE::OpenSMOKE_OdeSystemSolver<T>
	{
	public:

		/**
		*@brief Default constructor
		*@param odeSystem the ODE system to be solved
		*/
		OpenSMOKE_CHEMEQ2(T* odeSystem);

		/**
		*@brief Sets the total number of ordinary differential equations to be solved
		*@param n the total number of ordinary differential equations to be solved
		*/
		void SetDimensions(const int n);

		/**
		*@brief Sets the accuracy parameter for determining the next timestep
		*@param epsmin the accuracy parameter for determining the next timestep (default 0.01)
		*/
		void SetEpsilonMinimum(const double epsmin);

		/**
		*@brief Repeat timestep if correction is greater than epsmax*epsmin*y(i) for any i
		*@param epsmax the parameter governing the condition above (default 10.)
		*/
		void SetEpsilonMaximum(const double epsmax);

		/**
		*@brief Sets the minimum allowed step during the integration
		*@param dtmin the minimum allowed step (default 1e-15)
		*/
		void SetMinimumStep(const double dtmin);

		/**
		*@brief Sets the maximum number of corrector iterations to perform
		*@param max_number_subiterations number of corrector iterations to perform (default 1)
		*/
		void SetMaximumNumberSubIterations(const int max_number_subiterations);

		/**
		*@brief Sets the minimum values allowed for unknowns
		*@param ymin the minimum values allowed for unknowns (default 1e-20)
		*/
		void SetMinimumValue(const double ymin);

		/**
		*@brief Sets the additional checks to improve the accuracy (requires at leat 3 corrector iterations)
		*/
		void SetHighLevelAccuracy();

		/**
		*@brief Solves the ODE system
		*@param tf the final value of the independent variable to be reached
		*/
		void Solve(const double tf);

		/**
		*@brief Prints on the screen a short summary about the integration carried out
		*/
		void Status() const;

		/**
		*@brief Prints on the screen the name of the adopted ODE solver
		*/
		std::string Tag() const;

		/**
		*@brief Prints on the screen the number of steps
		*/
		int GetNumberOfSteps() const;

		/**
		*@brief Prints on the screen the number of function evaluations
		*/
		int GetNumberOfFunctionEvaluations() const;

		/**
		*@brief Prints on the screen the number of Jacobian evaluations
		*/
		int GetNumberOfJacobianEvaluations() const;

		/**
		*@brief Prints on the screen the number of LU factorizations
		*/
		int GetNumberOfLUFactorizations() const;

		/**
		*@brief Prints on the screen the number of non-linear iterations
		*/
		int GetNumberOfNonLinearIterations() const;

		/**
		*@brief Prints on the last order used
		*/
		int GetLastOrderUsed() const;

		/**
		*@brief Prints on the screen the number of convergence failures
		*/
		int GetNumberOfConvergenceFailures() const;

		/**
		*@brief Prints on the screen the number of error test failures
		*/
		int GetNumberOfErrorTestFailures() const;

		/**
		*@brief Prints on the size of the last step used
		*/
		double GetLastStepUsed() const;

		/**
		* Default destructor
		*/
		~OpenSMOKE_CHEMEQ2(void);

	private:

		int number_steps_;				//!< number of steps
		int number_system_calls_;		//!< number of calls to the ODE system
		int max_number_subiterations_;	//!< maximum number of times the corrector is applied (default 1)

		double eps_;					//!< maximum correction term, finally scaled by 1/epsmin
		int iteration_;					//!< current iteration
		int failures_;					//!< number of failures

		std::vector<double> q_;			//!< calculated production rates Eq. (1)
		std::vector<double> d_;			//!< calculated destruction rates
		std::vector<double> y1_;		//!< predicted value from Eq. (35)
		std::vector<double> ys_;		//!< initial values for the chemical time-step
		std::vector<double> qs_;		//!< initial production rates Eq. (38)
		std::vector<double> scrarray_;	//!< scratch (temporary) variable array
		std::vector<double> rtau_;		//!< ratio of timestep to timescale
		std::vector<double> rtaus_;		//!< ratio of timestep to initial timescale for current timestep
		std::vector<double> ymin_;		//!< minimum values allowed for unknowns (default 1e-20)
		
		
		std::vector<double> ym1_;	//!< vector used only when high level of accuracy option is turned on (previous corrector iterate)
		std::vector<double> ym2_;	//!< vector used only when high level of accuracy option is turned on (previous corrector iterate)

		double epsmin_;			//!< the maximum relative error allowed for convergence of the corrector step (default 0.01)
		double epsmax_;			//!< this number provides the basis for deciding weather convergence can
								//!< be achieved with out added stepsize achieved with out added stepsize
								//!< reduction. if eps / epsmin is greater than epsmx further reduction is applied (default 10)
		double sqreps_;			//!< parameter used to calculate initial timestep: sqreps_ = 5.*sqrt(epsmin_)
		double epscl_;			//!< intermediate variable used to avoid repeated divisions: epscl_ = 1./epsmin_

		double tfd_;			//!< round-off parameter used to determine when integration is complete
		double dtmin_;			//!< the smallest stepsize allowed (default 1e-15)
		

		unsigned int alpha_option_;	//!< approximation of alpha: 1=Pade b, 2=Pade a (default 1)
		

		// The accuracy-based timestep calculation can be augmented with a stability-based check when at least three corrector
		// iterations are performed. For most problems, the stability check is not needed, and eliminating the calculations
		// and logic associated with the check enhances performance.
		bool high_level_accuracy_;	//!< enables accuracy through stability based check (default false)

	private:

		/**
		*@brief Allocates the memory
		*@param n the total number of ordinary differential equations to be solved
		*/
		void MemoryAllocation(const int n);

		/**
		*@brief Checks if the provided user-defined paramets are consistent
		*/
		void AnalyzeUserOptions();

		/**
		*@brief Returns the alpha parameter
		*@param rtau the independent variable for calculating alpha
		*/
		double GetAlpha(const double rtau);

	private:

		// Value of dt/tau used to switch between Eqs. (39) and (41) when Pade(a) is used
		static const double rswitch_;	//!< constant adopted for calculation of alpha parameter (only if option 2 is enabled)
	};
}

#include "OpenSMOKE_CHEMEQ2.hpp"

#endif	// OpenSMOKE_CHEMEQ2_H
