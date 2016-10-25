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

#include "OdeSolver_Parameters_Grammar.h"

namespace OdeSMOKE
{
	OdeSolver_Parameters::OdeSolver_Parameters()
	{
		// Default values
		type_ = ODE_INTEGRATOR_OPENSMOKEPP;
		sparse_linear_algebra_ = false;
		
        jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
		preconditioner_ = OpenSMOKE::PRECONDITIONER_SPARSE_ILUT;

		// Contraints
		minimum_constraints_ = true;
		maximum_constraints_ = false;
		non_negative_unknowns_ = false;
		
		// Default values are adopted for tolerances
		relative_tolerance_ = 100.*OpenSMOKE::MachEpsFloat();
		absolute_tolerance_ = 1.e-10;
		
		// Default values are adopted for initial and maximu steps
		minimum_step_ = -1.;
		maximum_step_ = -1.;
		initial_step_ = -1.;

		// Numerical details
		maximum_number_of_steps_ = 5000000;
		maximum_err_test_fails_ = 50;
		maximum_conv_fails_ = 50;
		maximum_nl_iter_ = 50;
		maximum_order_ = 5;
		coefficient_nonlinear_convergence_test_ = 0.33;
		maximum_number_jacobians_ = -1;
		minimum_yp_ = 1.e-7;
		verbosity_level_ = 1;
	}

	void OdeSolver_Parameters::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		OdeSolver_Parameters_Grammar grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@OdeSolver") == true)
		{
			std::string name;
			dictionary.ReadString("@OdeSolver", name);
			
			if (name == "OpenSMOKE++")
			{
				type_ = ODE_INTEGRATOR_OPENSMOKEPP;
			}
			else if (name == "BzzOde")
			{
				#if OPENSMOKE_USE_BZZMATH == 0
					OpenSMOKE::FatalErrorMessage("OpenSMOKE++ was built without the BzzMath support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_BZZODE;
			}
			else if (name == "Cvode")
			{
				#if OPENSMOKE_USE_SUNDIALS == 0
					OpenSMOKE::FatalErrorMessage("OpenSMOKE++ was built without the Sundials (CVODE) support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_CVODE;
			}
			else if (name == "DASPK")
			{
				#if OPENSMOKE_USE_DASPK == 0
					OpenSMOKE::FatalErrorMessage("OpenSMOKE++ was built without the DASPK support. Please select a different ODE solver");
				#endif
				type_ = ODE_INTEGRATOR_DASPK;
			}

			else OpenSMOKE::FatalErrorMessage("Unknown ODE Solver: " + name);
		}

		// If the OpenSMOKE++ integrator is selected, the user can choose sparse linear algebra
		if (type_ == ODE_INTEGRATOR_OPENSMOKEPP)
		{
			if (dictionary.CheckOption("@Jacobian") == true)
			{
				std::string name;
				dictionary.ReadString("@Jacobian", name);

				if (name == "Band")				sparse_linear_algebra_ = false;
				else if (name == "Sparse")		sparse_linear_algebra_ = true;

				else OpenSMOKE::FatalErrorMessage("Unknown @Jacobian. Available options: Band (default) | Sparse");
			}

			if (dictionary.CheckOption("@SparseSolver") == true)
			{
				std::string name;
				dictionary.ReadString("@SparseSolver", name);

				#ifdef __APPLE__
				if (name == "EigenSparseLU")		jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
				else if (name == "EigenBiCGSTAB")   jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB;
				else if (name == "EigenGMRES")		jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES;
				else if (name == "EigenDGMRES")		jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES;
				else if (name == "Pardiso")			jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO;
				else if (name == "SuperLUSerial")   jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL;
				else if (name == "UMFPack")			jacobian_solver_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK;
                #else
				if (name == "EigenSparseLU")		jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_SPARSE_LU;
				else if (name == "EigenBiCGSTAB")   jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_BICGSTAB;
				else if (name == "EigenGMRES")		jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_GMRES;
				else if (name == "EigenDGMRES")		jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_DGMRES;
				else if (name == "Pardiso")			jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_PARDISO;
				else if (name == "SuperLUSerial")   jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL;
				else if (name == "UMFPack")			jacobian_solver_ = OpenSMOKE::SparseSolverType::SOLVER_SPARSE_EIGEN_UMFPACK;
                #endif

				else OpenSMOKE::FatalErrorMessage("Unknown @SparseSolver. Available options: EigenSparseLU (default) | EigenBiCGSTAB | EigenGMRES | EigenDGMRES | Pardiso | SuperLUSerial | UMFPack");
			}

			if (dictionary.CheckOption("@Preconditioner") == true)
			{
				std::string name;
				dictionary.ReadString("@Preconditioner", name);

				#ifdef __APPLE__
                if (name == "ILUT")			   preconditioner_ = OpenSMOKE::PRECONDITIONER_SPARSE_ILUT;
				else if (name == "diagonal")   preconditioner_ = OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL;
                #else
				if (name == "ILUT")			   preconditioner_ = OpenSMOKE::SparsePreconditionerType::PRECONDITIONER_SPARSE_ILUT;
				else if (name == "diagonal")   preconditioner_ = OpenSMOKE::SparsePreconditionerType::PRECONDITIONER_SPARSE_DIAGONAL;
                #endif

				else OpenSMOKE::FatalErrorMessage("Unknown @Preconditioner. Available options: ILUT (default) | diagonal");
			}
		}
		

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

		if (dictionary.CheckOption("@MeanResidualThreshold") == true)
			dictionary.ReadDouble("@MeanResidualThreshold", minimum_yp_);

		if (dictionary.CheckOption("@VerbosityLevel") == true)
			dictionary.ReadInt("@VerbosityLevel", verbosity_level_);

		if (dictionary.CheckOption("@MinimumConstraints") == true)
			dictionary.ReadBool("@MinimumConstraints", minimum_constraints_);

		if (dictionary.CheckOption("@MaximumConstraints") == true)
			dictionary.ReadBool("@MaximumConstraints", maximum_constraints_);

		if (dictionary.CheckOption("@NonNegativeVariables") == true)
			dictionary.ReadBool("@NonNegativeVariables", non_negative_unknowns_);
	}
}
