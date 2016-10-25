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

#ifndef OpenSMOKE_SensitivityMap_H
#define OpenSMOKE_SensitivityMap_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/src/IterativeSolvers/DGMRES.h>

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include "maps/KineticsMap_CHEMKIN.h"

#if OPENSMOKE_USE_MKL == 1
#include <Eigen/PardisoSupport>
#endif

#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
#include <Eigen/SuperLUSupport>
#endif

#if OPENSMOKE_USE_UMFPACK == 1
#include <Eigen/UmfPackSupport>
#endif

namespace OpenSMOKE
{
	//!  A class to perform sensitivity analysis
	/*!
			A class to perform sensitivity analysis
	*/

	template<typename map> 
	class SensitivityMap
	{
	
	public:

		/**
		*@brief Default constructor
		*@param kineticsMap map containing the kinetic mechanism
		*@param number_of_equations total number of equations of the ODE (or NLS) associated to the sensitivity analysis
		*/
		SensitivityMap(KineticsMap_CHEMKIN<map>& kineticMap, const unsigned int number_of_equations);

		/**
		*@brief Returns the dense solver type
		*/
		DenseSolverType dense_solver_type() const { return dense_solver_type_; }

		/**
		*@brief Returns the sparse solver type
		*/
		SparseSolverType sparse_solver_type() const { return sparse_solver_type_;  }

		/**
		*@brief Set the sparsity pattern 
		*@param row indices of rows of non-zero elements
		*@param col indices of columns of non-zero elements
		*/
		void SetSparsityPattern(const std::vector<unsigned int>& row, const std::vector<unsigned int>& col);
	
		/**
		*@brief Calculates the current sensitivity coefficients (dense version)
		*@param t current time [s]
		*@param T current temperature [K]
		*@param P_Pa current pressure [Pa]
		*@param c current concentrations [kmol/m3]
		*@param J current Jacobian matrix [NExNE]
		*@param scaling_Jp current scaling vector[NE]
		*/
		void Calculate(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J, const OpenSMOKE::OpenSMOKEVectorDouble& scaling_Jp);

		/**
		*@brief Calculates the current sensitivity coefficients (sparse version)
		*@param t current time [s]
		*@param T current temperature [K]
		*@param P_Pa current pressure [Pa]
		*@param c current concentrations [kmol/m3]
		*@param J current Jacobian matrix [NExNE]
		*@param scaling_Jp current scaling vector[NE]
		*/
		void Calculate(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const Eigen::SparseMatrix<double>& J, const OpenSMOKE::OpenSMOKEVectorDouble& scaling_Jp);

		/**
		*@brief Calculates the steady state sensitivity coefficients (dense version)
		*@param T current temperature [K]
		*@param P_Pa current pressure [Pa]
		*@param c current concentrations [kmol/m3]
		*@param J current Jacobian matrix [NExNE]
		*@param scaling_Jp current scaling vector[NE]
		*/
		void Calculate(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J, const OpenSMOKE::OpenSMOKEVectorDouble& scaling_Jp);

		/**
		*@brief Calculates the steady state sensitivity coefficients (sparse version)
		*@param T current temperature [K]
		*@param P_Pa current pressure [Pa]
		*@param c current concentrations [kmol/m3]
		*@param J current Jacobian matrix [NExNE]
		*@param scaling_Jp current scaling vector[NE]
		*/
		void Calculate(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const Eigen::SparseMatrix<double>& J, const OpenSMOKE::OpenSMOKEVectorDouble& scaling_Jp);

		/**
		*@brief Sets the number of substeps to perform the sensitivity analysis (default 2)
		*/
		void SetNumberOfSubSteps(const unsigned int number_of_substeps);
		
		/**
		*@brief Sets the sensitivity type (default SENSITIVITY_KINETIC_CONSTANT)
		*/
		void SetSensitivityType(const PhysicalConstants::sensitivity_type type);

		/**
		*@brief Sets the linear algebra package for solving the linear systems (default SOLVER_DENSE_EIGEN)
		*/
		void SetDenseSolverType(const OpenSMOKE::DenseSolverType type);

		/**
		*@brief Sets the linear algebra package for solving the linear systems (default SOLVER_SPARSE_NONE)
		*/
		void SetSparseSolverType(const OpenSMOKE::SparseSolverType type);

		/**
		*@brief Sets the LU factorization type (default DENSE_DECOMPOSITION_FULL_PIVOTING_LU)
		*/
		void SetDenseDecompositionType(const DenseDecompositionType type);

		/**
		*@brief Sets the sparse preconditioner type (default PRECONDITIONER_SPARSE_ILUT)
		*/
		void SetSparsePreconditionerType(const OpenSMOKE::SparsePreconditionerType type);

		/**
		*@brief Sets the sparse preconditioner drop tolerance
		*/
		void SetSparsePreconditionerDropTolerance(const double droptol_);

		/**
		*@brief Sets the sparse preconditioner fill factor
		*/
		void SetSparsePreconditionerFillFactor(const int fillfactor_);

		/**
		*@brief Sets the type of energy equation (default CONSTANT_VOLUME_SYSTEM)
		*/
		void SetEnergyEquationType(const EnergyEquationType type);

		/**
		*@brief Returns the sensitivity coefficient matrix
		*/
		const Eigen::MatrixXd& sensitivity_coefficients() const { return sensitivity_coeffs_; }

		/**
		*@brief Returns the list of sensitivity parameters, according to the sensitivity_type
		*/
		const OpenSMOKE::OpenSMOKEVectorDouble& parameters() const { return parameters_; }

		/**
		*@brief Returns the total number of parameters
		*/
		unsigned int number_of_parameters() const { return number_of_parameters_; }

		/**
		*@brief Sets the index of temperature equation (1-index based)
		*/
		void SetIndexOfTemperature(const unsigned int index);

		/**
		*@brief Sets the index of density equation (1-index based)
		*/
		void SetIndexOfDensity(const unsigned int index);

		/**
		*@brief Returns the vector containing the derivatives with respect to the parameters for the species
		*/
		inline const OpenSMOKE::OpenSMOKEVectorDouble& Jp_species()		 { return Jp_species_; };

		/**
		*@brief Returns the cpu time to perform a single factorization
		*/
		double cpuTimeSingleFactorization() const { return cpuTimeSingleFactorization_; }

		/**
		*@brief Returns the cpu time to perform a all the single factorizations (cumulative)
		*/
		double cpuTimeFactorization() const { return cpuTimeFactorization_; }

		/**
		*@brief Returns the cpu time to perform a single solution (NP linear systems)
		*/
		double cpuTimeSingleSolution() const { return cpuTimeSingleSolution_; }

		/**
		*@brief Returns the cpu time to perform the single solutions (NP linear systems) (cumulative)
		*/
		double cpuTimeSolution() const { return cpuTimeSolution_; }

		/**
		*@brief Returns the cpu time to perform a single assembling
		*/
		double cpuTimeSingleAssembling() const { return cpuTimeSingleAssembling_; }

		/**
		*@brief Returns the cpu time to perform the single assembling (cumulative)
		*/
		double cpuTimeAssembling() const { return cpuTimeAssembling_; }

	protected:

		KineticsMap_CHEMKIN<map>& kinetics_;	//!< reference to the kinetic map

		unsigned int index_of_temperature_;		//!< index of temperature equation (if available)
		unsigned int index_of_density_;		//!< index of temperature equation (if available)
		unsigned int index_of_species_;			//!< index of first species equation
		unsigned int number_of_species_;		//!< total number of species
		unsigned int number_of_reactions_;		//!< total number of reactions
		unsigned int number_of_parameters_;		//!< total number of parameters
		unsigned int number_of_equations_;		//!< number of equations (for the jacobian construction)

		Eigen::MatrixXd*				A_eigen_;				//!< linear system matrix (dense formulation)
		Eigen::SparseMatrix<double>*	A_eigen_sparse_;		//!< linear system matrix (sparse formulation)
		Eigen::SparseMatrix<double>		ones_;					//!< identity matrix (sparse formulation)
		Eigen::VectorXd*				b_eigen_;				//!< linear system right side vector

		Eigen::MatrixXd sensitivity_coeffs_;			//!< matrix containing the sensitivity coefficients
		Eigen::MatrixXd B_eigen_;						//!< linear system right side matrix
		Eigen::MatrixXd Alfa_;							//!< auxiliary vector

		OpenSMOKE::OpenSMOKEVectorDouble parameters_;					//!< vector containing the sensitivity parameters
		OpenSMOKE::OpenSMOKEVectorDouble Jp_species_;					//!< vector containing the derivatives with respect to the parameters for the species

		double t_;		//!< current time [s]
		double tOld_;	//!< previous time [s]

		EnergyEquationType energy_type_;						//!< type of energy equation (default CONSTANT_VOLUME_SYSTEM)
		PhysicalConstants::sensitivity_type sensitivity_type_;	//!< type of parameters
		unsigned int number_of_substeps_;						//!< number of substeps (default 2)

		// Dense mode
		DenseSolverType dense_solver_type_;						//!< linear algebra package for solving the linear systems (default LINEAR_ALGEBRA_PACKAGE_EIGEN)
		DenseDecompositionType dense_decomposition_type_;		//!< LU factorization type (default LINEAR_DECOMPOSITION_FULL_PIVOTING_LU)

		Eigen::PartialPivLU<Eigen::MatrixXd> partial_LU;
		Eigen::FullPivLU<Eigen::MatrixXd> full_LU;

		// Sparse Mode
		Eigen::SparseLU<Eigen::SparseMatrix<double> > sparse_LU_;														//!< LU solver

		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_bicgstab_diagonal_;	//!< BiCGSTAB solver (diagonal)
		Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_gmres_diagonal_;		//!< GMRES solver (diagonal) 
		Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_dgmres_diagonal_;		//!< DGMRES solver (diagonal) 

		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_bicgstab_ilut_;				//!< BiCGSTAB solver (ILUT)
		Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_gmres_ilut_;					//!< GMRES solver (ILUT) 
		Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_dgmres_ilut_;					//!< DGMRES solver (ILUT) 

		#if OPENSMOKE_USE_MKL == 1
				Eigen::PardisoLU<Eigen::SparseMatrix<double> > sparse_pardiso_;							//!< MKL PARDISO solver 
		#endif

		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
				Eigen::SuperLU<Eigen::SparseMatrix<double> > sparse_superlu_serial_;					//!< SuperLUSerial solver 
		#endif

		#if OPENSMOKE_USE_UMFPACK == 1
				Eigen::UmfPackLU<Eigen::SparseMatrix<double> > sparse_umfpack_;							//!< UMFPack solver 
		#endif

		#if OPENSMOKE_USE_LIS == 1
				// TODO
		#endif

		SparseSolverType			sparse_solver_type_;			//!< sparse solver type (linear algebra)
		SparsePreconditionerType	sparsePreconditionerType_;		//!< preconditioner type (iterative solvers)

		double sparse_preconditioner_droptol_;
		int sparse_preconditioner_fillfactor_;

		double cpuTimeSingleFactorization_;			//!< cpu time to perform a single factorization
		double cpuTimeFactorization_;				//!< cpu time to perform a all the single factorizations (cumulative)
		double cpuTimeSingleSolution_;				//!< cpu time to perform a single solution (NP linear systems)
		double cpuTimeSolution_;					//!< cpu time to perform the single solutions (NP linear systems) (cumulative)
		double cpuTimeSingleAssembling_;			//!< cpu time to perform a single overall operations (factorization, solution assembling) (NP linear systems)
		double cpuTimeAssembling_;					//!< cpu time to perform the single overall operations (factorization, solution assembling) (cumulative)
	};
}

#include "SensitivityMap.hpp"

#endif // OpenSMOKE_SensitivityMap_H
