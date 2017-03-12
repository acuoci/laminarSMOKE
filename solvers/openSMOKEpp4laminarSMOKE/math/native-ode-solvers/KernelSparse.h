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

#ifndef OdeKernelSparse_H
#define OdeKernelSparse_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/src/IterativeSolvers/DGMRES.h>

#if OPENSMOKE_USE_MKL == 1
#include <Eigen/PardisoSupport>
#endif
#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
#include <Eigen/SuperLUSupport>
#endif
#if OPENSMOKE_USE_UMFPACK == 1
#include <Eigen/UmfPackSupport>
#endif


namespace OdeSMOKE
{
	//!  A class to manage dense Jacobian matrices associated to ODE systems
	/*!
	A class to manage dense Jacobian matrices associated to ODE systems
	*/

	template <typename ODESystemObject>
	class KernelSparse : public ODESystemObject
	{
	public:

		/**
		*@brief Default constructor
		*/
		KernelSparse();

		/**
		*@brief Function to set the sparsity pattern
		*@param i vector containing the indices of rows with non-zero elements
		*@param j vector containing the indices of columns with non-zero elements
		*/
		void SetSparsityPattern(const std::vector<unsigned int>& i, const std::vector<unsigned int>& j);

		/**
		*@brief Function to set the solver for the solution of the dense linear system
		*@param linear_algebra_solver EigenSparseLU, EigenBiCGSTAB
		*/
		void SetLinearAlgebraSolver(const std::string linear_algebra_solver);

		/**
		*@brief Function to set the preconditioner for the iterative solvers
		*@param preconditioner: diagonal, ILUT
		*/
		void SetPreconditioner(const std::string preconditioner);

		/**
		*@brief Function to set the preconditioner drop tolerance (ILUT)
		*@param drop the drop tolerance (default 1e-6)
		*/
		void SetPreconditionerDropTolerance(const double drop);

		/**
		*@brief Function to set the preconditioner fill factor (ILUT)
		*@param fillfactor the fille factor (default 10)
		*/
		void SetPreconditionerFillFactor(const int fillfactor);

		/**
		*@brief Returns the number of calls to the system of equations for assembling the Jacobian matrix
		Since the Jacobian matrix is full, the number of calls is equal to the number of Jacobian calls times the number
		of equations
		*/
		unsigned int numberOfFunctionCallsForJacobian() const { return numberOfFunctionCallsForJacobian_; }

		/**
		*@brief Returns the cumulative CPU time to assemble the JAcobian matrix
		*/
		double cpuTimeToAssembleJacobian() const { return cpuTimeToAssembleJacobian_; }

		/**
		*@brief Returns the cumulative CPU time to factorized the G matrix
		*/
		double cpuTimeToFactorize() const { return cpuTimeToFactorize_; }

		/**
		*@brief Returns the cumulative CPU time to solve the linear system (the matrix is supposed already factorized)
		*/
		double cpuTimeToSolveLinearSystem() const { return cpuTimeToSolveLinearSystem_; }

		/**
		*@brief Returns the CPU time for a single assembling of the Jacobian matrix
		*/
		double cpuTimeSingleJacobianAssembling() const { return cpuTimeSingleJacobianAssembling_; }

		/**
		*@brief Returns the CPU time for a single factorization of the G matrix
		*/
		double cpuTimeSingleFactorization() const { return cpuTimeSingleFactorization_; }

		/**
		*@brief Returns the CPU time for a single solution of the linear system (the matrix is supposed already factorized)
		*/
		double cpuTimeSingleLinearSystemSolution() const { return cpuTimeSingleLinearSystemSolution_; }

		/**
		*@brief Summary
		*/
		void OdeSolverKernelSummary(std::ostream& out);

	protected:

		/**
		*@brief Prepares the object (default option, memory allocation, etc.)
		*/
		void MemoryAllocationKernel();

		/**
		*@brief Prepares the object (default option, memory allocation, etc.)
		*/
		void ResetKernel();

		/**
		*@brief Calculates the Jacobian numerically, using the usual differentiation approach
		*@param y the current vector of dependent variables
		*@param t the current value of independent variable
		*@param f the current vector of rhs
		*@param e the current vector of error weights: e(j) = 1. / ( tolAbs + tolRel*y(j) )
		*@param max_constraints true if constraints on the maximum values of unknowns are imposed
		*@param yMax maximum values (if any) imposed on the unknowns
		*/
		void NumericalJacobian(Eigen::VectorXd& y, const double t, const Eigen::VectorXd& f, const double h, const Eigen::VectorXd& e, const bool max_constraints, const Eigen::VectorXd& yMax);

		/**
		*@brief Calculates the Jacobian analytically
		*@param y the current vector of dependent variables
		*@param t the current value of independent variable
		*/
		void UserDefinedJacobian(const Eigen::VectorXd& y, const double t);

		/**
		*@brief Build and factorizes the G matrix to obtain the correction db
		The matrix G is defined as: G = I - hr0*J, where hr0 is a scalar, J is the Jacobian matrix and I the identity matrix
		*@param hr0 the scalar defined above
		*/
		void BuildAndFactorizeMatrixG(const double hr0);

		/**
		*@brief Solve the linear system associated to the G matrix
		*@param db vector of known terms (on input) and vector containing the solution (on output)
		*/
		void SolveLinearSystem(Eigen::VectorXd& db);

		/**
		*@brief Returns the product of the Jacobian matrix times a vector: v_out = J*v_in
		*@param v_in vector to be multiplied
		*@param v_out resulting vector
		*/
		void JacobianTimesVector(const Eigen::VectorXd& v_in, Eigen::VectorXd* v_out);


	private:

		std::vector<unsigned int> i_;
		std::vector<unsigned int> j_;

		Eigen::SparseMatrix<double> J_;		    //!< Jacobian matrix
		Eigen::MatrixXd				Jaux_;		//!< Jacobian matrix
		Eigen::SparseMatrix<double> G_;			//!< matrix to be factorized
		Eigen::SparseMatrix<double> ones_;
		Eigen::VectorXd aux_;					//!< auxiliary vector (dimension equal to the number of equations)

		Eigen::SparseLU<Eigen::SparseMatrix<double> > sparse_LU_;													//!< LU solver

		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_bicgstab_diagonal_;	//!< BiCGSTAB solver (diagonal)
		Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_gmres_diagonal_;		//!< GMRES solver (diagonal) 
		Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_dgmres_diagonal_;		//!< DGMRES solver (diagonal) 

		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_bicgstab_ilut_;			//!< BiCGSTAB solver (ILUT)
		Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_gmres_ilut_;				//!< GMRES solver (ILUT) 
		Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_dgmres_ilut_;				//!< DGMRES solver (ILUT) 

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
		LIS_INT *lis_ptr_;
		LIS_INT *lis_index_;
		LIS_SCALAR *lis_value_;
		LIS_MATRIX lis_G_;
		LIS_VECTOR lis_x_;
		LIS_VECTOR lis_b_;
		LIS_SOLVER lis_solver_;
		#endif
		
		OpenSMOKE::SparseSolverType			solverType_;			//!< sparse solver type (linear algebra)
		OpenSMOKE::SparsePreconditionerType	preconditionerType_;	//!< preconditioner type (iterative solvers)

		double preconditioner_droptol_;
		int preconditioner_fillfactor_;
		
		unsigned int numberOfFunctionCallsForJacobian_;	//!< number of calls to the system of equation for assembling the Jacobian matrix

		double cpuTimeToAssembleJacobian_;				//!< cumulative CPU time for assembling the Jacobian
		double cpuTimeToFactorize_;						//!< cumulative CPU time for factorizing the G matrix
		double cpuTimeToSolveLinearSystem_;				//!< cumulative CPU time for solving the linear system
		double cpuTimeSingleJacobianAssembling_;		//!< CPU time for assembling a single Jacobian matrix
		double cpuTimeSingleFactorization_;				//!< CPU time for factorizing a single G matrix
		double cpuTimeSingleLinearSystemSolution_;		//!< CPU time for solving a single linear system
	};
}

#include "KernelSparse.hpp"

#endif