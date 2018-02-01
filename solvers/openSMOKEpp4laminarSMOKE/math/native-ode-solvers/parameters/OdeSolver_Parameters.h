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

#ifndef OpenSMOKE_OdeSolver_Parameters_H
#define	OpenSMOKE_OdeSolver_Parameters_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <math/OpenSMOKEFunctions.h>
#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OdeSMOKE
{
	class OdeSolver_Parameters
	{
	public:

		enum ODE_INTEGRATOR { ODE_INTEGRATOR_OPENSMOKEPP, ODE_INTEGRATOR_BZZODE, ODE_INTEGRATOR_CVODE, ODE_INTEGRATOR_DASPK };

	public:
	
		/**
		*@brief Default constructor
		*/
		OdeSolver_Parameters();


		/**
		*@brief Reads the options from a file
		*@param dictionary dictionary from which the oprions are extracted
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		template <typename OdeSolver>
		void TransferDataFromOdeSolver(const OdeSolver& dae_solver, const double cpuTime);

		/**
		*@brief Prints the current status on a file
		*@param fOut file where to print
		*/
		void Status(std::ostream& fOut);

		/**
		*@brief Sets the type of integrator
		*@param type integrator type: ODE_INTEGRATOR_OPENSMOKEPP, ODE_INTEGRATOR_BZZODE, ODE_INTEGRATOR_CVODE, ODE_INTEGRATOR_DASPK
		*/
		void SetType(const ODE_INTEGRATOR type) { type_ = type;}

		/**
		*@brief Sets the relative tolerance
		*@param relative_tolerance the relative tolerance (same for all the variables)
		*/
		void SetRelativeTolerance(const double relative_tolerance) { relative_tolerance_ = relative_tolerance; }

		/**
		*@brief Sets the absolute tolerance
		*@param absolute_tolerance the absolute tolerance (same for all the variables)
		*/
		void SetAbsoluteTolerance(const double absolute_tolerance) { absolute_tolerance_ = absolute_tolerance; }

		/**
		*@brief Sets the minimum step which can be used by the integrator
		*@param minimum_step the minimum step which can be used by the integrator
		*/
		void SetMinimumStep(const double minimum_step) { minimum_step_ = minimum_step; }

		/**
		*@brief Sets the maximum step which can be used by the integrator
		*@param maximum_step the maximum step which can be used by the integrator
		*/
		void SetMaximumStep(const double maximum_step) { maximum_step_ = maximum_step; }

		/**
		*@brief Sets the initial step which can has to be used by the integrator
		*@param intitial_step the initial step which can has to be used by the integrator
		*/
		void SetInitialStep(const double initial_step) { initial_step_ = initial_step; }

		/**
		*@brief Sets the maximum number of steps which can be used by the integrator
		*@param maximum_number_of_steps the maximum number of steps which can be used by the integrator
		*/
		void SetMaximumNumberOfSteps(const int maximum_number_of_steps) { maximum_number_of_steps_ = maximum_number_of_steps; }

		/**
		*@brief Sets the maximum number of error test failures which is allowed
		*@param maximum_err_test_fail the maximum number of error test failures which is allowed
		*/
		void SetMaximumErrorTestFailures(const int maximum_err_test_fails) { maximum_err_test_fails_ = maximum_err_test_fails; }

		/**
		*@brief Sets the maximum number of convergence failures which is allowed
		*@param maximum_conv_fails the maximum number of convergence failures which is allowed
		*/
		void SetMaximumConvergenceFailures(const int maximum_conv_fails) { maximum_conv_fails_ = maximum_conv_fails; }

		/**
		*@brief Sets the maximum number of non linear iterations which is allowed
		*@param maximum_nl_iter the maximum number of non linear iterations which is allowed
		*/
		void SetMaximumNLIterations(const int maximum_nl_iter) { maximum_nl_iter_ = maximum_nl_iter; }

		/**
		*@brief Sets the maximum order adopted by the integrator
		*@param maximum_order the maximum order adopted by the integrator
		*/
		void SetMaximumOrder(const int maximum_order) { maximum_order_ = maximum_order; }

		/**
		*@brief Sets the coefficient to be used in non linear convergence tests
		*@param coefficient_nonlinear_convergence_test the coefficient to be used in non linear convergence tests
		*/
		void SetCoefficientNonLinearConvergenceTest(const double coefficient_nonlinear_convergence_test) { coefficient_nonlinear_convergence_test_ = coefficient_nonlinear_convergence_test; }
		
		/**
		*@brief Sets the maximum number of Jacobian matrices which can be calculated during the integration
		*@param maximum_number_jacobiansthe maximum number of Jacobian matrices which can be calculated during the integration
		*/
		void SetMaximumNumberJacobians(const int maximum_number_jacobians) { maximum_number_jacobians_ = maximum_number_jacobians; }

		/**
		*@brief Sets the minimum residuals for considering the system of equations in steady state conditions
		*@param minimum_ypthe minimum residuals for considering the system of equations in steady state conditions
		*/
		void SetMinimumMeanThreshold(const double minimum_yp) { minimum_yp_ = minimum_yp; }

		/**
		*@brief Turns on/off the minimum constraints on the variables/unknowns
		*@param flag true means the minimum constraints are turned on (by default turned on)
		*/
		void SetMinimumConstraints(const bool flag) { minimum_constraints_ = flag; }

		/**
		*@brief Turns on/off the maximum constraints on the variables/unknowns
		*@param flag true means the maximum constraints are turned on (by default turned off)
		*/
		void SetMaximumConstraints(const bool flag) { maximum_constraints_ = flag; }

		/**
		*@brief True if all the unknowns must be non negative
		*@param flag true means all the unknowns must be non negative (by default turned off)
		*/
		void SetNonNegativeUnknowns(const bool flag) { non_negative_unknowns_ = flag; }

		/**
		*@brief Sets the verbosity level
		*@param verbosity_level the verbosity level (0 means no output)
		*/
		void SetVerbosityLevel(const int verbosity_level) { verbosity_level_ = verbosity_level; }
		
		/**
		*@brief Returns the integrator type: ODE_INTEGRATOR_OPENSMOKEPP, ODE_INTEGRATOR_BZZODE, ODE_INTEGRATOR_IDA, ODE_INTEGRATOR_DASPK
		*/
		ODE_INTEGRATOR type() const { return type_; }

		/**
		*@brief Returns the solver to be used for solving the linear systems associated to the Jacobian matrix
		*/
		OpenSMOKE::SparseSolverType jacobian_solver() const { return jacobian_solver_; }

		/**
		*@brief Returns the preconditioner to be used for solving the linear systems associated to the Jacobian matrix
		*/
		OpenSMOKE::SparsePreconditionerType preconditioner() const { return preconditioner_; }

		/**
		*@brief Returns true is a sparse linear algebra solver is enabled
		*/
		bool sparse_linear_algebra() const { return sparse_linear_algebra_; }

		/**
		*@brief Returns the relative tolerance
		*/
		double relative_tolerance() const { return relative_tolerance_; }

		/**
		*@brief Returns the absolute tolerance
		*/
		double absolute_tolerance() const { return absolute_tolerance_; }

		/**
		*@brief Returns the minimum step (default or defined by the user)
		*/
		double minimum_step() const { return minimum_step_; }

		/**
		*@brief Returns the maximum step (default or defined by the user)
		*/
		double maximum_step() const { return maximum_step_; }

		/**
		*@brief Returns the initial step (default or defined by the user)
		*/
		double initial_step() const { return initial_step_; }

		/**
		*@brief Returns the maximum order (default or defined by the user)
		*/
		int maximum_order() const { return maximum_order_; }

		/**
		*@brief Returns the maximum number of steps (default or defined by the user)
		*/
		int maximum_number_of_steps() const { return maximum_number_of_steps_; }

		/**
		*@brief Returns the maximum number of error test failures (default or defined by the user)
		*/
		int maximum_err_test_fails() const { return maximum_err_test_fails_; }

		/**
		*@brief Returns the maximum number of convergence failures (default or defined by the user)
		*/
		int maximum_conv_fails() const { return maximum_conv_fails_; }

		/**
		*@brief Returns the maximum number of non linear iterations (default or defined by the user)
		*/
		int  maximum_nl_iter() const { return  maximum_nl_iter_; }

		/**
		*@brief Returns the coefficient to be used in non linear convergence tests (default or defined by the user)
		*/
		double coefficient_nonlinear_convergence_test() const { return coefficient_nonlinear_convergence_test_; }

		/**
		*@brief Returns the maximum number of Jacobian matrix evaluations (default or defined by the user)
		*/
		int maximum_number_jacobians() const { return maximum_number_jacobians_; }

		/**
		*@brief Returns the minimum residuals for considering the system of equations in steady state conditions (default or defined by the user)
		*/
		double minimum_yp() const { return minimum_yp_; }

		/**
		*@brief Returns true is the minimum constraints are turned on
		*/
		bool minimum_constraints() const { return minimum_constraints_; }

		/**
		*@brief Returns true is the maximum constraints are turned on
		*/
		bool maximum_constraints() const { return maximum_constraints_; }

		/**
		*@brief Returns true is all the unknowns are non negative
		*/
		bool non_negative_unknowns() const { return non_negative_unknowns_; }

		/**
		*@brief Returns the verbosity level
		*/
		int verbosity_level() const { return verbosity_level_; }

	private:

		ODE_INTEGRATOR type_;									//!< ODE integrator type
		double relative_tolerance_;								//!< relative tolerance
		double absolute_tolerance_;								//!< absolute tolerance
		double minimum_step_;									//!< minimum step to be adopted during the integration
		double maximum_step_;									//!< maximum step to be adopted during the integration
		double initial_step_;									//!< initial step to be adopted 
		int maximum_order_;										//!< maximum order
		int maximum_number_of_steps_;							//!< maximum number of steps
		int maximum_err_test_fails_;							//!< maximum number of error test failures
		int maximum_conv_fails_;								//!< maximum number of convergence failures
		int maximum_nl_iter_;									//!< maximum number of non linear iterations
		OpenSMOKE::SparseSolverType jacobian_solver_;			//!< solver to be used for solving linear systems associated to the Jacobian matrix
		OpenSMOKE::SparsePreconditionerType preconditioner_;	//!< preconditioner to be used for solving linear systems associated to the Jacobian matrix
		bool sparse_linear_algebra_;							//!< true if sparse linear solvers have to be used

		double coefficient_nonlinear_convergence_test_;			//!< coefficient to be used in non linear convergence tests
		int maximum_number_jacobians_;							//!< maximum number of Jacobian matrices evaluations
		double minimum_yp_;										//!< minimum residuals for considering the system of equations in steady state conditions
		bool minimum_constraints_;								//!< true if the minimum constraints are turned on (by default turned on)
		bool maximum_constraints_;								//!< true if the maximum constraints are turned on (by default turned off)
		bool non_negative_unknowns_;							//!< true if the all the variables must be non negative (by default turned off)
		int verbosity_level_;									//!< verbosity level
	};
}

#include "OdeSolver_Parameters.hpp"

#endif	/* OpenSMOKE_OdeSolver_Parameters_H */

