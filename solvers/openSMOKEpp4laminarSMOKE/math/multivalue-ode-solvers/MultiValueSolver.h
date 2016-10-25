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

#ifndef MultiValueSolver_H
#define MultiValueSolver_H

#include <Eigen/Dense>
#include "math/OpenSMOKEUtilities.h"
#include "OdeSolverUtilities.h"

namespace OdeSMOKE
{
	enum OdeErrorStatus
	{
		ERROR_STATUS_FAILURE,
		ERROR_STATUS_OK
	};

	extern bool stopOdeIntegration = false;

	//!  A class to solve ODE systems using a multivalue method
	/*!
	The purpose of this class is the solution of a ODE system using a proper multivalue method.
	The class is based on the Method policy, which specifies the specific method used as multivalue algorithm. As an
	example, for stiff problems, the Method must be one of the Gear (or BDF) algorithms.
	*/

	template <typename Method>
	class MultiValueSolver : public Method
	{
	public:

		/**
		*@brief Default constructor
		*/
		MultiValueSolver();

		/**
		*@brief Set the initial conditions
		*@param t0 the initial value of the independent variable
		*@param y0 the vector containing the initial values of the dependent variables
		*/
		void SetInitialConditions(const double t0, const Eigen::VectorXd& y0);

		/**
		*@brief Solve the ODE system
		*@param tOut the final value of the independent variable
		*/
		OdeStatus Solve(const double tOut);

		/**
		*@brief Summary
		*@param out output stream
		*/
		void OdeSummary(std::ostream& out);

		/**
		*@brief Enables the possibility to print results and data on the screen
		*@param flag if true the possibility to print is turned on
		*/
		void SetPrint(const bool flag);

		/**
		*@brief Set the absolute tolerances (default 1e-12)
		*@param abs_tolerances the vector of requested absolute tolerances
		*/
		void SetAbsoluteTolerances(const Eigen::VectorXd& abs_tolerances);

		/**
		*@brief Set the absolute tolerances (default 1e-12)
		*@param abs_tolerances the requested absolute tolerances
		*/
		void SetAbsoluteTolerances(const double abs_tolerance);

		/**
		*@brief Set the relative tolerances
		*@param rel_tolerances the vector of requested relative tolerances
		*/
		void SetRelativeTolerances(const Eigen::VectorXd& rel_tolerances);

		/**
		*@brief Set the relative tolerances
		*@param rel_tolerances the requested relative tolerances
		*/
		void SetRelativeTolerances(const double rel_tolerance);

		/**
		*@brief Set the maximum allowed values for the dependent variables (default: none)
		*@param max_values the requested maximum values
		*/
		void SetMaximumValues(const Eigen::VectorXd& max_values);

		/**
		*@brief Set the maximum allowed values for the dependent variables (default: none)
		*@param max_values the requested maximum value
		*/
		void SetMaximumValues(const double max_value);

		/**
		*@brief Set the minimum allowed values for the dependent variables (default: none)
		*@param min_values the requested minimum values
		*/
		void SetMinimumValues(const Eigen::VectorXd& min_values);

		/**
		*@brief Set the minimum allowed values for the dependent variables (default: none)
		*@param min_value the requested minimum value
		*/
		void SetMinimumValues(const double min_value);

		/**
		*@brief Set the maximum number of steps
		*@param max_number_steps the maximum number of steps
		*/
		void SetMaximumNumberOfSteps(const unsigned int max_number_steps);

		/**
		*@brief Set the maximum size of the steps
		*@param  max_step_size the maximum size of the steps
		*/
		void SetMaximumStepSize(const double max_step_size);

		/**
		*@brief Set the minimum size of the steps
		*@param  max_step_size the minimum size of the steps
		*/
		void SetMinimumStepSize(const double min_step_size);

		/**
		*@brief Set the size of the first step
		*@param  initial_step_size the msize of the first step
		*/
		void SetFirstStepSize(const double initial_step_size);

		/**
		*@brief Set the maximum order which can be used during the integration
		*@param maximum_order the maximum order which can be used during the integration
		*/
		void SetMaximumOrder(const unsigned int maximum_order);

		/**
		*@brief Set a strict constraint on the maximum value of the dependent variable
		*@param tMax the maximum value of the dependent variable
		*/
		void SetMaxConstraintOnIndependentVariable(const double tMax);

		/**
		*@brief Remove the constraint on the independent variable
		*/
		void UnsetMaxConstraintOnIndependentVariable();

		/**
		*@brief The integration is stopped when the norm1 of dy/dt is smaller than the user defined amount YPrimNorm1
		*@param YPrimNorm1 the user defined threshold value
		*/
		void SetStopConditionMaximumYPrimeNorm1(const double YPrimNorm1);

		/**
		*@brief Reset the stop conditions (default: the integration is performed up to the request final value of independent variable)
		*/
		void UnsetStopConditions();

		/**
		*@brief Returns the size of the first step
		*/
		double firstStepSize() const { return first_step_size_; }

		/**
		*@brief Returns the size of the minimum step used during the integration
		*/
		double minimumStepUsed() const { return hMinUsed_; }

		/**
		*@brief Returns the size of the maximum step used during the integration
		*/
		double maximumStepUsed() const { return hMaxUsed_; }

		/**
		*@brief Returns the number of steps performed without changing neither the size, neither the order
		*/
		unsigned int numberOfStepsWithoutChanges() const { return numberOfStepsWithoutChanges_; }

		/**
		*@brief Returns the total number of failures in the Newton's convergence procedure
		*/
		unsigned int numberOfConvergenceFailure() const { return numberOfConvergenceFailure_; }

		/**
		*@brief Returns the total number of successes in the Newton's convergence procedure
		*/
		unsigned int numberOfConvergenceSuccess() const { return numberOfConvergenceSuccess_; }

		/**
		*@brief Returns the total number of failures in the local error test
		*/
		unsigned int numberOfErrorCheckFailure() const { return numberOfErrorCheckFailure_; }

		/**
		*@brief Returns the total number of successes in the local error test
		*/
		unsigned int numberOfErrorCheckSuccess() const { return numberOfErrorCheckSuccess_; }

		/**
		*@brief Returns the current solution at the end of calculations
		*@param solution the requested solution
		*/
		void Solution(Eigen::VectorXd& solution) const;

		/**
		*@brief Returns the first order derivateives at the end of calculations
		*@param first_order_derivatives the calculated first order derivatives (dy/dt)
		*/
		void FirstOrderDerivatives(Eigen::VectorXd& first_order_derivatives);

		/**
		*@brief Returns the second order derivatives at the end of calculations
		*@param second_order_derivatives the calculated second order derivatives (d2y/dt2)
		*/
		void SecondOrderDerivatives(Eigen::VectorXd& second_order_derivatives);

		/**
		*@brief Returns the estimation of the local error at the end of calculations
		*@param error_estimation the estimated error
		*/
		void ErrorEstimation(Eigen::VectorXd& error_estimation);

	private:

		/**
		*@brief Set default values
		*/
		void SetDefaultValues();

		/**
		*@brief reset the counters
		*/
		void Reset();

		/**
		*@brief Calculates the initial step size (estimation, which is rapidly adapted after some steps)
		*/
		double CalculateInitialStepSize();

		/**
		*@brief Calculates the initial step size (estimation, which is rapidly adapted after some steps)
		*/
		void CalculatesOneOverEpsilon(const Eigen::VectorXd& w);

		/**
		*@brief Interpolates to return the solution exactly for the requested final value of the dependent variable
		*/
		void Interpolation();

		/**
		*@brief Main function: solves the requested ODE system
		*       This function performs all the calculations needed to advance from the initial to the final value of the independent variable
		*/
		void DriverStep();

		/**
		*@brief Function to parse the status of the Ode solver
		*/
		void ParseOdeStatus();

		/**
		*@brief Function to estimate optimal order and step size
		*/
		void NewOrderNewH();

		/**
		*@brief Function to advance on a single step
		*/
		void MultiValueStep();

		/**
		*@brief Fatal error message
		*/
		void FatalErrorMessage(const std::string message);

		/**
		*@brief Fatal error message
		*/
		void WarningMessage(const std::string message);

	private:

		// Step size
		double first_step_size_;	//!< size of the initial step
		double hPreviousStep_;		//!< size of the previous step
		double hNordsieck_;			//!< size of Nordisieck step
		double max_step_size_;		//!< maximum size of the step
		double hMinUsed_;			//!< minimum step size reached during the integration
		double hMaxUsed_;			//!< maximum step size reached during the integration
		double hScaleMax_;			//!< maximum scale factor for increasing the step size

		bool user_defined_min_step_size_;		//!< true if the user defined the minimum step size
		bool user_defined_first_step_size_;		//!< true if the user defined the first step size
		bool user_defined_max_step_size_;		//!< true if the user defined the maximum step size

		// Independent variable
		double t_;					//!< value of independent variable
		double t0_;					//!< initial value of the independent variable
		double tInMeshPoint_;		//!< value of the independent variable at the beginning of the step
		double tStabilize_;			//!< ???

		// Cumulative counters
		unsigned int numberOfStepsWithoutChanges_;		//!< total number of steps performed without changing the order and the step size with respect to the previous step
		unsigned int numberOfConvergenceFailure_;		//!< total number of failures in finding the correction by the Method policy
		unsigned int numberOfConvergenceSuccess_;		//!< total number of successes in finding the correction by the Method policy
		unsigned int numberOfErrorCheckFailure_;		//!< total number of times in which the error after performing a step was too large
		unsigned int numberOfErrorCheckSuccess_;		//!< total number of times in which the error after performing a was ok

		// Cpu times
		double cpuTimeEquationSystem_;		//!< cpu time to solve a single system of equations
		double cpuTimeToFactorize_;			//!< cpu time to factorize the G matrix

		// Tolerances
		Eigen::VectorXd abs_tolerances_;			//!< absolute tolerances (vector)
		Eigen::VectorXd rel_tolerances_;			//!< relative tolerances (vector)
		double abs_tolerance_;						//!< absolute tolerance (scalar)
		double rel_tolerance_;						//!< relative tolerance (scalar)
		bool abs_tolerances_scalar_;				//!< types of absolute tolerance (scalar or vector)
		bool rel_tolerances_scalar_;				//!< types of relative tolerance (scalar or vector)

		// Strict bound on the independent variable
		bool calculation_with_strict_tmax_;		//!< flag to request a strict bound on the maximum value of the independent variable
		double tStrictMax_;						//!< maximum value of the independent variable (which cannot be overcome)

		// Stop the integration for small first order derivatives
		bool stopIntegrationForSmallYPrimeNorm1_;	//!< flag to request to stop the integration when the norm of the first order derivatives is sufficiently small
		double YPrimeNorm1_;						//!< maximum value of the norm of the first order derivatives

		// Vectors
		Eigen::VectorXd y0_;					//!< vector containing the initial values of the dependent variables
		Eigen::VectorXd one_over_epsilon_;		//!< error vector: one_over_epsilon = (tolAbs + tolRel*abs(y))^(-1)

		// User-defined requests
		double tOut_;					//!< final value of the dependent variable
		unsigned int max_number_steps_;	//!< maximum number of steps allowed
		bool printResults_;				//!< flag to enable the call to the PrintResults function

		// Internal variables (cannot be chosen by the user)
		bool direction_plus_;				//!< integration direction
		OdeErrorStatus errorStatus_;		//!< flag to recognize if the error was ok or too large

		bool fixedOrderMin_;				//!< flag to force the algorithm to use the minimum order

		// Local counters
		unsigned int iterConstSteps_;								//!< number of steps performed without changing the size
		unsigned int iterMinimumOrder_;								//!< number of steps performed under the minimum order constraint
		unsigned int iterOfConstStepsForLargeFactorizationTime_;

		// Default consts
		static const double DEFAULT_TOL_ABS;							//!< default absolute tolerance
		static const double DEFAULT_TOL_REL;							//!< relative absolute tolerance
		static const double DEFAULT_HSCALE_MAX1;						//!< maximum scale factor during the very first steps
		static const double DEFAULT_HSCALE_MAX2;						//!< maximum scale factor during the first steps
		static const double DEFAULT_HSCALE_MAX3;						//!< maximum scale factor after the initial transient
		static const unsigned int DEFAULT_MAXIMUM_NUMBER_OF_STEPS;		//!< maximum step size
		static const unsigned int MAX_ERROR_FAILURE;					//!< maximum number of error failures in the same step
	};

}

#include "MultiValueSolver.hpp"

#endif