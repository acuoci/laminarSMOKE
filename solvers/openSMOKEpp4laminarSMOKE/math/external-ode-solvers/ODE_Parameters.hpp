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

namespace OpenSMOKE
{
	void ODE_Parameters::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_ODEParameters_Options grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@OdeSolver") == true)
		{
			std::string name;
			dictionary.ReadString("@OdeSolver", name);
			
			if (name == "OpenSMOKE") 
			{
				type_ = ODE_INTEGRATOR_OPENSMOKE;
			}
			else if (name == "CHEMEQ2")
			{
				type_ = ODE_INTEGRATOR_CHEMEQ2;
			}
			else if (name == "BzzOde")
			{
				#if OPENSMOKE_USE_BZZMATH == 0
					FatalErrorMessage("OpenSMOKE++ was built without the BzzMath support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_BZZODE;
			}
			else if (name == "CVODE")
			{
				#if OPENSMOKE_USE_SUNDIALS == 0
					FatalErrorMessage("OpenSMOKE++ was built without the Sundials (CVODE) support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_CVODE;
			}
			else if (name == "DASPK")
			{
				#if OPENSMOKE_USE_DASPK == 0
					FatalErrorMessage("OpenSMOKE++ was built without the DASPK support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_DASPK;
			}
			else if (name == "DVODE")
			{
				#if OPENSMOKE_USE_DVODE == 0
					FatalErrorMessage("OpenSMOKE++ was built without the DVODE support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_DVODE;
			}
			else if (name == "DLSODA") 
			{
				#if OPENSMOKE_USE_ODEPACK == 0
					FatalErrorMessage("OpenSMOKE++ was built without the ODEPACK (DLSODE and DLSODA) support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_DLSODA;
			}
			else if (name == "DLSODE") 
			{
				#if OPENSMOKE_USE_ODEPACK == 0
					FatalErrorMessage("OpenSMOKE++ was built without the ODEPACK (DLSODE and DLSODA) support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_DLSODE;
			}
			else if (name == "MEBDF") 
			{
				#if OPENSMOKE_USE_MEBDF == 0
					FatalErrorMessage("OpenSMOKE++ was built without the MEBDF support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_MEBDF;
			}
			else if (name == "RADAU5") 
			{
				#if OPENSMOKE_USE_RADAU == 0
					FatalErrorMessage("OpenSMOKE++ was built without the RADAU support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_RADAU5;
			}

			else FatalErrorMessage("Unknown ODE Solver: " + name);
		}

		if (dictionary.CheckOption("@LinearAlgebra") == true)
		{
			std::string name;
			dictionary.ReadString("@LinearAlgebra", name);

			if (name == "Eigen")
			{
				linear_algebra_ = "Eigen";
			}
			else if (name == "BzzMath")
			{
				#if OPENSMOKE_USE_BZZMATH == 0
					FatalErrorMessage("OpenSMOKE++ was built without the BzzMath support. Please select a different @LinearAlgebra option");
				#endif
				linear_algebra_ = "BzzMath";
			}
			else if (name == "Plasma")
			{
				#if OPENSMOKE_USE_PLASMA == 0
					FatalErrorMessage("OpenSMOKE++ was built without the PLASMA support. Please select a different @LinearAlgebra option");
				#endif
					linear_algebra_ = "Plasma";
			}
			else if (name == "Flame")
			{
				#if OPENSMOKE_USE_FLAME == 0
					FatalErrorMessage("OpenSMOKE++ was built without the FLAME support. Please select a different @LinearAlgebra option");
				#endif
				linear_algebra_ = "Flame";
			}
			else
			{
				FatalErrorMessage("Unknown Linear Algebra: " + name);
			}
		}

		if (dictionary.CheckOption("@SparseSolver") == true)
		{
			std::string name;
			dictionary.ReadString("@SparseSolver", name);

			if (name == "none")					sparse_solver_ = "none";
			else if (name == "EigenBiCGSTAB")	sparse_solver_ = "EigenBiCGSTAB";
			else if (name == "EigenGMRES")		sparse_solver_ = "EigenGMRES";
			else if (name == "EigenDGMRES")		sparse_solver_ = "EigenDGMRES";
			else if (name == "EigenSparseLU")	sparse_solver_ = "EigenSparseLU";
			else if (name == "Pardiso")
			{
				#if OPENSMOKE_USE_MKL == 0
					FatalErrorMessage("OpenSMOKE++ was built without the MKL support. Please select a different @SparseSolver option");
				#endif
				sparse_solver_ = "Pardiso";
			}
			else if (name == "SuperLUSerial")
			{
				#if OPENSMOKE_USE_SUPERLU_SERIAL == 0
					FatalErrorMessage("OpenSMOKE++ was built without the SuperLU support. Please select a different @SparseSolver option");
				#endif
				sparse_solver_ = "SuperLUSerial";
			}
			else if (name == "UMFPack")
			{
				#if OPENSMOKE_USE_UMFPACK == 0
					FatalErrorMessage("OpenSMOKE++ was built without the UMFPACK support. Please select a different @SparseSolver option");
				#endif
				sparse_solver_ = "UMFPack";
			}
			else if (name == "LIS")
			{
				#if OPENSMOKE_USE_LIS == 0
					FatalErrorMessage("OpenSMOKE++ was built without the LIS support. Please select a different @SparseSolver option");
				#endif
				sparse_solver_ = "LIS";
			}
			else
			{
				FatalErrorMessage("Unknown Sparse Solver: " + name);
			}
		}

		if (dictionary.CheckOption("@SparsePreconditioner") == true)
		{
			std::string name;
			dictionary.ReadString("@SparsePreconditioner", name);

			if (name == "diagonal")		preconditioner_ = "diagonal";
			else if (name == "ILUT")	preconditioner_ = "ILUT";

			else FatalErrorMessage("Unknown Preconditioner: " + name);
		}

		if (dictionary.CheckOption("@SparsePreconditionerDropTolerance") == true)
			dictionary.ReadDouble("@SparsePreconditionerDropTolerance", drop_tolerance_);

		if (dictionary.CheckOption("@SparsePreconditionerFillFactor") == true)
			dictionary.ReadInt("@SparsePreconditionerFillFactor",fill_factor_);

		if (dictionary.CheckOption("@RelativeTolerance") == true)
			dictionary.ReadDouble("@RelativeTolerance", relative_tolerance_);

		if (dictionary.CheckOption("@AbsoluteTolerance") == true)
			dictionary.ReadDouble("@AbsoluteTolerance", absolute_tolerance_);

		if (dictionary.CheckOption("@MaximumOrder") == true)
			dictionary.ReadInt("@MaximumOrder", maximum_order_);

		if (dictionary.CheckOption("@MaximumStep") == true)
			dictionary.ReadDouble("@MaximumStep", maximum_step_);

		if (dictionary.CheckOption("@MinimumStep") == true)
			dictionary.ReadDouble("@MinimumStep", minimum_step_);

		if (dictionary.CheckOption("@InitialStep") == true)
			dictionary.ReadDouble("@InitialStep", initial_step_);

		if (dictionary.CheckOption("@MaximumNumberOfSteps") == true)
			dictionary.ReadInt("@MaximumNumberOfSteps", maximum_number_of_steps_);

		if (dictionary.CheckOption("@FullPivoting") == true)
			dictionary.ReadBool("@FullPivoting", full_pivoting_);
	}

	ODE_Parameters::ODE_Parameters()
	{
		// Default values
		type_ = ODE_INTEGRATOR_OPENSMOKE;
		linear_algebra_ = "Eigen";
		sparse_solver_ = "none";
		preconditioner_ = "diagonal";
		drop_tolerance_ = 1e-6;
		fill_factor_ = 10;
		relative_tolerance_ = 100.*MachEpsFloat();
		absolute_tolerance_ = 1.e-10;
		maximum_number_of_steps_ = -1;
		minimum_step_ = -1.;
		maximum_step_ = -1.;
		initial_step_ = -1.;
		maximum_order_ = -1;
		full_pivoting_ = false;
		
		// Reset the counters
		time_spent_to_factorize_ = 0.;
		time_spent_to_evaluate_jacobian_ = 0.;
		time_spent_to_solve_linear_system_ = 0.;
		cpu_time_ = 0.;
		number_of_function_calls_ = 0;
		number_of_function_calls_jacobian_ = 0;
		number_of_jacobians_ = 0;
		number_of_factorizations_ = 0;
		number_of_steps_ = 0;
		last_order_used_ = 0;
		last_step_used_ = 0.;
		max_order_used_ = 0;
		max_step_used_ = 0.;
		min_step_used_ = 0.;
		number_of_nonlinear_iterations_ = 0;
		number_of_convergence_failures_ = 0;
		number_of_error_test_failures_ = 0;

		// Names of ODE solvers
		list_names_.push_back("OpenSMOKE");
		list_names_.push_back("BzzOde");
		list_names_.push_back("CVODE");
		list_names_.push_back("DVODE");
		list_names_.push_back("DASPK");
		list_names_.push_back("DLSODE");
		list_names_.push_back("DLSODA");
		list_names_.push_back("RADAU5");
		list_names_.push_back("MEBDF");
		list_names_.push_back("CHEMEQ2");
	}

	void ODE_Parameters::Status(std::ostream& fOut)
	{
		double other_cpu_time = cpu_time_ - (time_spent_to_evaluate_jacobian_ + time_spent_to_factorize_ + time_spent_to_solve_linear_system_);

		fOut << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << "   ODE Integration                                                           " << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << " * Integrator:                          " << list_names_[type_] << std::endl;
		fOut << " * Sparse solver:                       " << sparse_solver_ << std::endl;
		fOut << " * Preconditioner:                      " << preconditioner_ << std::endl;
		fOut << " * Linear algebra:                      " << linear_algebra_ << std::endl;
		fOut << " * CPU Time (s):                        " << std::scientific << cpu_time_ << std::endl;
		fOut << " *   assembling Jacobian:               " << std::scientific << time_spent_to_evaluate_jacobian_;
		fOut                                               << std::fixed << " (" << time_spent_to_evaluate_jacobian_ / cpu_time_*100. << "%)" << std::endl;
		fOut << " *   LU decomposition:                  " << std::scientific << time_spent_to_factorize_;
		fOut                                               << std::fixed << " (" << time_spent_to_factorize_ / cpu_time_*100. << "%)" << std::endl;
		fOut << " *   linear system solution:            " << std::scientific << time_spent_to_solve_linear_system_;
		fOut                                               << std::fixed << " (" << time_spent_to_solve_linear_system_ / cpu_time_*100. << "%)" << std::endl;
		fOut << " *   other:                             " << std::scientific << other_cpu_time;
		fOut                                               << std::fixed << " (" << other_cpu_time / cpu_time_*100. << "%)" << std::endl;
		fOut << " * Relative tolerance:                  " << std::scientific << relative_tolerance_ << std::endl;
		fOut << " * Absolute tolerance:                  " << std::scientific << absolute_tolerance_ << std::endl;
		fOut << " * Number of steps:                     " << std::fixed << number_of_steps_ << std::endl;
		fOut << " * Number of function calls:            " << std::fixed << number_of_function_calls_ + number_of_function_calls_jacobian_ << std::endl;
		fOut << " * Number of function calls (Jacobian): " << std::fixed << number_of_function_calls_jacobian_ << " (" << double(number_of_function_calls_jacobian_) / double(number_of_function_calls_ + number_of_function_calls_jacobian_)*100. << "%)" << std::endl;
		fOut << " * Number of Jacobians:                 " << std::fixed << number_of_jacobians_ << std::endl;
		fOut << " * Number of factorizations:            " << std::fixed << number_of_factorizations_ << std::endl;
		fOut << " * Number of convergence failures:      " << std::fixed << number_of_convergence_failures_ << std::endl;
		fOut << " * Number of error test failures:       " << std::fixed << number_of_error_test_failures_ << std::endl;
		fOut << " * Last order used:                     " << std::fixed << last_order_used_ << std::endl;
		fOut << " * Last step used:                      " << std::scientific << last_step_used_ << std::endl;
		fOut << " * Max order used:                      " << std::fixed << max_order_used_ << std::endl;
		fOut << " * Max step used:                       " << std::scientific << max_step_used_ << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
	}

	template <typename OdeSolver>
	void ODE_Parameters::TransferDataFromOdeSolver(const OdeSolver&ode_solver, const double cpuTime)
	{
		SetCPUTime(cpuTime);
		SetNumberOfFunctionCalls(ode_solver.numberOfFunctionCalls());
		SetNumberOfJacobians(ode_solver.numberOfJacobianEvaluations());
		SetNumberOfFactorizations(ode_solver.numberOfMatrixFactorizations());
		SetNumberOfSteps(ode_solver.numberOfSteps());
		SetLastOrderUsed(ode_solver.lastOrderUsed());
		SetLastStepUsed(ode_solver.lastStepUsed());
		SetMaxOrderUsed(ode_solver.maxOrderUsed());
		SetMinStepUsed(ode_solver.minimumStepUsed());
		SetMaxStepUsed(ode_solver.maximumStepUsed());
		SetNumberOfFunctionCallsForJacobian(ode_solver.numberOfFunctionCallsForJacobian());
		SetTimeSpentToFactorize(ode_solver.cpuTimeToFactorize());
		SetTimeSpentToEvaluateJacobian(ode_solver.cpuTimeToAssembleJacobian());
		SetTimeSpentToSolveLinearSystem(ode_solver.cpuTimeToSolveLinearSystem());
		SetNumberOfConvergenceFailures(ode_solver.numberOfConvergenceFailure());
		SetNumberOfErrorTestFailures(ode_solver.numberOfErrorCheckFailure());
	}
}
