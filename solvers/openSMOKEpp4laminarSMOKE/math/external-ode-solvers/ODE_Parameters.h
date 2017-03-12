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

#ifndef OpenSMOKE_ODE_Parameters_H
#define	OpenSMOKE_ODE_Parameters_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "math/OpenSMOKEFunctions.h"

namespace OpenSMOKE
{
	class ODE_Parameters
	{
	public:

		enum ODE_INTEGRATOR {	ODE_INTEGRATOR_OPENSMOKE, 
								ODE_INTEGRATOR_BZZODE, ODE_INTEGRATOR_CVODE, ODE_INTEGRATOR_DVODE, ODE_INTEGRATOR_DASPK,
								ODE_INTEGRATOR_DLSODE, ODE_INTEGRATOR_DLSODA, ODE_INTEGRATOR_RADAU5, ODE_INTEGRATOR_MEBDF,
								ODE_INTEGRATOR_CHEMEQ2	};

		ODE_Parameters();

		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		template <typename OdeSolver>
		void TransferDataFromOdeSolver(const OdeSolver&ode_solver, const double cpuTime);

		void Status(std::ostream& fOut);

		void SetType(const ODE_INTEGRATOR type) { type_ = type;}
		void SetRelativeTolerance(const double relative_tolerance) { relative_tolerance_ = relative_tolerance; }
		void SetAbsoluteTolerance(const double absolute_tolerance) { absolute_tolerance_ = absolute_tolerance; }
		void SetMinimumStep(const double minimum_step) { minimum_step_ = minimum_step; }
		void SetMaximumStep(const double maximum_step) { maximum_step_ = maximum_step; }
		void SetInitialStep(const double initial_step) { initial_step_ = initial_step; }
		void SetMaximumNumberOfSteps(const int maximum_number_of_steps) { maximum_number_of_steps_ = maximum_number_of_steps; }
		void SetMaximumOrder(const int maximum_order) { maximum_order_ = maximum_order; }
		void SetFullPivoting(const bool flag) { full_pivoting_ = flag; }
		
		void SetCPUTime(const double cpu_time) { cpu_time_ = cpu_time; }
		void SetNumberOfFunctionCalls(const int number_of_function_calls) { number_of_function_calls_ = number_of_function_calls; }
		void SetNumberOfFunctionCallsForJacobian(const int number_of_function_calls_jacobian) { number_of_function_calls_jacobian_ = number_of_function_calls_jacobian; }
		void SetNumberOfJacobians(const int number_of_jacobians) { number_of_jacobians_ = number_of_jacobians; }
		void SetNumberOfFactorizations(const int number_of_factorizations) {number_of_factorizations_ = number_of_factorizations; }
		void SetNumberOfSteps(const int number_of_steps) { number_of_steps_ = number_of_steps; }
		void SetLastOrderUsed(const int last_order_used) { last_order_used_ = last_order_used; }
		void SetLastStepUsed(const double last_step_used) { last_step_used_ = last_step_used; }
		void SetMaxOrderUsed(const int max_order_used) { max_order_used_ = max_order_used; }
		void SetMaxStepUsed(const double max_step_used) { max_step_used_ = max_step_used; }
		void SetMinStepUsed(const double min_step_used) { min_step_used_ = min_step_used; }
		void SetNumberOfNonLinearIterations (const int number_of_nonlinear_iterations)	{ number_of_nonlinear_iterations_ = number_of_nonlinear_iterations; }
		void SetNumberOfConvergenceFailures (const int number_of_convergence_failures)	{ number_of_convergence_failures_ = number_of_convergence_failures; }
		void SetNumberOfErrorTestFailures (const int number_of_error_test_failures)		{ number_of_error_test_failures_ = number_of_error_test_failures; }

		void SetTimeSpentToFactorize (const double time_spent_to_factorize)		{ time_spent_to_factorize_ = time_spent_to_factorize; }
		void SetTimeSpentToEvaluateJacobian (const double time_spent_to_evaluate_jacobian)		{ time_spent_to_evaluate_jacobian_ = time_spent_to_evaluate_jacobian; }
		void SetTimeSpentToSolveLinearSystem (const double time_spent_to_solve_linear_system)		{ time_spent_to_solve_linear_system_ = time_spent_to_solve_linear_system; }

		ODE_INTEGRATOR type() const { return type_; }
		std::string tag() const { return list_names_[type_]; }
		std::string linear_algebra() const { return linear_algebra_; }
		std::string sparse_solver() const { return sparse_solver_; }
		std::string preconditioner() const { return preconditioner_; }
		bool full_pivoting() const { return full_pivoting_; }
		double cpu_time() const { return cpu_time_; }
		double relative_tolerance() const { return relative_tolerance_; }
		double absolute_tolerance() const { return absolute_tolerance_; }
		double minimum_step() const { return minimum_step_; }
		double maximum_step() const { return maximum_step_; }
		double initial_step() const { return initial_step_; }
		int maximum_order() const { return maximum_order_; }
		int number_of_steps() const { return number_of_steps_; }
		int number_of_function_calls() const { return number_of_function_calls_; }
		int number_of_function_calls_jacobian() const { return number_of_function_calls_jacobian_; }
		int number_of_jacobians() const { return number_of_jacobians_; }
		int maximum_number_of_steps() const { return maximum_number_of_steps_; }
		int number_of_nonlinear_iterations() const { return number_of_nonlinear_iterations_; }
		int number_of_convergence_failures() const { return number_of_convergence_failures_; }
		int number_of_error_test_failures() const { return number_of_error_test_failures_; }
		double drop_tolerance() const { return drop_tolerance_; };
		int fill_factor() const { return fill_factor_; };
		

		double time_spent_to_factorize() const { return time_spent_to_factorize_; }
		double time_spent_to_evaluate_jacobian() const { return time_spent_to_evaluate_jacobian_; }
		double time_spent_to_solve_linear_system() const { return time_spent_to_solve_linear_system_; }

	private:

		ODE_INTEGRATOR type_;
		double relative_tolerance_;
		double absolute_tolerance_;
		double minimum_step_;
		double maximum_step_;
		double initial_step_;
		int maximum_order_;
		int maximum_number_of_steps_;
		std::string linear_algebra_;
		std::string sparse_solver_;
		std::string preconditioner_;
		bool full_pivoting_;
		double drop_tolerance_;
		int fill_factor_;

		double cpu_time_; 
		int number_of_function_calls_;
		int number_of_function_calls_jacobian_;
		int number_of_jacobians_;
		int number_of_factorizations_;
		int number_of_steps_;
		int last_order_used_;
		double last_step_used_;
		int max_order_used_;
		double max_step_used_;
		double min_step_used_;
		int number_of_nonlinear_iterations_;
		int number_of_convergence_failures_;
		int number_of_error_test_failures_;

		std::vector<std::string> list_names_;

		double time_spent_to_factorize_ ;
		double time_spent_to_evaluate_jacobian_;
		double time_spent_to_solve_linear_system_;
	};

	class Grammar_ODEParameters_Options : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeSolver", 
																OpenSMOKE::SINGLE_STRING, 
																"ODE Solver: OpenSMOKE | BzzOde | CVODE | DASPK | DVODE | DLSODA | DLSODE | MEBDF | RADAU5", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@LinearAlgebra", 
																OpenSMOKE::SINGLE_STRING, 
																"Linear Algebra: Eigen | BzzMath | Plasma | Flame", 
																false) );	

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparseSolver",
																OpenSMOKE::SINGLE_STRING,
																"Sparse solver: none | EigenBiCGSTAB | EigenGMRES | EigenDGMRES |EigenSparseLU | Pardiso | SuperLUSerial | UMFPack | LIS",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparsePreconditioner",
																OpenSMOKE::SINGLE_STRING,
																"Preconditioner: diagonal | ILUT",
																false));
			
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RelativeTolerance", 
																OpenSMOKE::SINGLE_DOUBLE, 
																"Releative tolerance (default 1.2e-5)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@AbsoluteTolerance", 
																OpenSMOKE::SINGLE_DOUBLE, 
																"Releative tolerance (default 1e-10)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumNumberOfSteps", 
																OpenSMOKE::SINGLE_INT, 
																"Maximum number of steps (default 500000)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumStep", 
																OpenSMOKE::SINGLE_DOUBLE, 
																"Maximum step (default: automatically chosen by the solver)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinimumStep", 
																OpenSMOKE::SINGLE_DOUBLE, 
																"Minimum step (default: automatically chosen by the solver)", 
																false) );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialStep", 
																OpenSMOKE::SINGLE_DOUBLE, 
																"Initial step (default: automatically chosen by the solver)", 
																false) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FullPivoting",
															   OpenSMOKE::SINGLE_BOOL,
															  "Full pivoting during the LU decomposition (default: false)",
															   false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumOrder",
															   OpenSMOKE::SINGLE_INT,
															   "Maximum order to be used during the ODE integration",
															   false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparsePreconditionerDropTolerance",
																OpenSMOKE::SINGLE_DOUBLE,
																"Preconditioner drop tolerance (default 1e-6)",
																false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparsePreconditionerFillFactor",
																OpenSMOKE::SINGLE_INT,
																"Preconditioner fill factor (default 10)",
																false));
		}
	};
}

#include "ODE_Parameters.hpp"

#endif	/* OpenSMOKE_ODE_Parameters_H */

