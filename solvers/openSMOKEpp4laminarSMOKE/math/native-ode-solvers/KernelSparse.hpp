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

namespace OdeSMOKE
{
	template <typename ODESystemObject>
	KernelSparse<ODESystemObject>::KernelSparse()
	{
		solverType_			= OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
		preconditionerType_ = OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL;
		preconditioner_droptol_ = 1e-6;
		preconditioner_fillfactor_ = 10;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::ResetKernel()
	{
		numberOfFunctionCallsForJacobian_ = 0;

		cpuTimeToAssembleJacobian_ = 0.;
		cpuTimeToFactorize_ = 0.;
		cpuTimeToSolveLinearSystem_ = 0.;

		cpuTimeSingleJacobianAssembling_ = 0.;
		cpuTimeSingleFactorization_ = 0.;
		cpuTimeSingleLinearSystemSolution_ = 0.;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::MemoryAllocationKernel()
	{
		// Allocate memory (if needed) for the OdeSystem object
		this->MemoryAllocation();

		// Internal variables
		aux_.resize(this->ne_);
		J_.resize(this->ne_, this->ne_);
		Jaux_.resize(this->ne_, this->ne_);
		G_.resize(this->ne_, this->ne_);
		ones_.resize(this->ne_, this->ne_);

		// Set sparsity pattern Jacobian
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList;
			tripletList.reserve(i_.size());
			for (unsigned int k = 0; k < i_.size(); k++)
				tripletList.push_back(T(i_[k], j_[k], 1.));

			J_.setFromTriplets(tripletList.begin(), tripletList.end());
			J_.makeCompressed();
		}

		// Set sparsity pattern G Matrix (we have to be sure that diagonal elements are present)
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList;
			tripletList.reserve(i_.size() + this->ne_);
			for (unsigned int k = 0; k < i_.size(); k++)
				tripletList.push_back(T(i_[k], j_[k], 1.));
			for (unsigned int k = 0; k < this->ne_; k++)
				tripletList.push_back(T(k, k, 1.));

			G_.setFromTriplets(tripletList.begin(), tripletList.end());
			G_.makeCompressed();
		}

		// Set identity matrix
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList;
			tripletList.reserve(this->ne_);
			for (unsigned int k = 0; k < this->ne_; k++)
				tripletList.push_back(T(k, k, 1.));

			ones_.setFromTriplets(tripletList.begin(), tripletList.end());

			ones_.makeCompressed();
		}

		// Analyze sparsity 
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			sparse_LU_.analyzePattern(G_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_bicgstab_diagonal_.analyzePattern(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_bicgstab_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_bicgstab_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);

				sparse_bicgstab_ilut_.analyzePattern(G_);
			}
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_gmres_diagonal_.analyzePattern(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_gmres_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_gmres_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);

				sparse_gmres_ilut_.analyzePattern(G_);
			}
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_dgmres_diagonal_.analyzePattern(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_dgmres_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_dgmres_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);
				sparse_dgmres_ilut_.analyzePattern(G_);
			}
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			sparse_pardiso_.analyzePattern(G_);
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			sparse_superlu_serial_.analyzePattern(G_);
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			sparse_umfpack_.analyzePattern(G_);
		}
		#endif
		#if OPENSMOKE_USE_LIS == 1
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			lis_value_ = (LIS_SCALAR *)malloc(G_.nonZeros()*sizeof(LIS_SCALAR));
			lis_ptr_ = (LIS_INT *)malloc((this->ne_ + 1)*sizeof(LIS_INT));
			lis_index_ = (LIS_INT *)malloc(G_.nonZeros()*sizeof(LIS_INT));

			for (unsigned int i = 0; i <= this->ne_; i++)
				lis_ptr_[i] = 0;

			unsigned int count = 0;
			for (int k = 0; k < G_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(G_, k); it; ++it)
			{
				lis_index_[count] = it.row();
				lis_ptr_[it.col() + 1]++;
				count++;
			}

			for (unsigned int i = 1; i <= this->ne_; i++)
				lis_ptr_[i] += lis_ptr_[i - 1];

			// Memory allocation
			LIS_INT err;
			err = lis_matrix_create(LIS_COMM_WORLD, &lis_G_);
			CHKERR(err);
			err = lis_matrix_set_size(lis_G_, 0, this->ne_);
			CHKERR(err);


			lis_vector_create(LIS_COMM_WORLD, &lis_x_);
			lis_vector_set_size(lis_x_, this->ne_, 0);

			lis_vector_create(LIS_COMM_WORLD, &lis_b_);
			lis_vector_set_size(lis_b_, this->ne_, 0);


			lis_solver_create(&lis_solver_);
			lis_solver_set_option("-i gmres -p ilut -initx_zeros 0", lis_solver_);
			lis_solver_set_option("-ilut_drop 1e-8 -ilut_rate 2", lis_solver_);
		}
		#endif
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::SetSparsityPattern(const std::vector<unsigned int>& i, const std::vector<unsigned int>& j)
	{
		i_ = i;
		j_ = j;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::SetLinearAlgebraSolver(const std::string linear_algebra_solver)
	{
		if (linear_algebra_solver == "EigenSparseLU")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
		}
		else if (linear_algebra_solver == "EigenBiCGSTAB")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB;
		}
		else if (linear_algebra_solver == "EigenGMRES")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES;
		}
		else if (linear_algebra_solver == "EigenDGMRES")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES;
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (linear_algebra_solver == "Pardiso")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO;
		}
		#else
		else if (linear_algebra_solver == "Pardiso")
		{
			OpenSMOKE::ErrorMessage("KernelSparse<ODESystemObject>", "Requested Pardiso linear algebra solver is not supported because the code was compiled without MKL support!");
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (linear_algebra_solver == "SuperLUSerial")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL;
		}
		#else
		else if (linear_algebra_solver == "SuperLUSerial")
		{
			OpenSMOKE::ErrorMessage("KernelSparse<ODESystemObject>", "Requested SuperLUSerial linear algebra solver is not supported because the code was compiled without SuperLU support!");
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (linear_algebra_solver == "UMFPack")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK;
		}
		#else
		else if (linear_algebra_solver == "UMFPack")
		{
			OpenSMOKE::ErrorMessage("KernelSparse<ODESystemObject>", "Requested UMFPack linear algebra solver is not supported because the code was compiled without UMFPACK support!");
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (linear_algebra_solver == "LIS")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS;
		}
		#else
		else if (linear_algebra_solver == "LIS")
		{
			OpenSMOKE::ErrorMessage("KernelSparse<ODESystemObject>", "Requested LIS linear algebra solver is not supported because the code was compiled without LIS support!");
		}
		#endif
		else
		{
			OpenSMOKE::ErrorMessage("KernelSparse<ODESystemObject>", "Requested linear algebra is not supported! Available options are: EigenSparseLU || EigenBiCGSTAB || EigenGMRES || EigenDGMRES || Pardiso || SuperLUSerial || UMFPack || LIS");
		}
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::SetPreconditioner(const std::string preconditioner)
	{
		if (preconditioner == "diagonal")
		{
			preconditionerType_ = OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL;
		}
		else if (preconditioner == "ILUT")
		{
			preconditionerType_ = OpenSMOKE::PRECONDITIONER_SPARSE_ILUT;
		}
		else
		{
			OpenSMOKE::ErrorMessage("KernelSparse<ODESystemObject>", "Requested preconditioner is not supported! Available options are: diagonal || ILUT");
		}
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::SetPreconditionerDropTolerance(const double drop)
	{
		preconditioner_droptol_ = drop;
	}
	
	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::SetPreconditionerFillFactor(const int fillfactor)
	{
		preconditioner_fillfactor_ = fillfactor;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::JacobianTimesVector(const Eigen::VectorXd& v_in, Eigen::VectorXd* v_out)
	{
		*v_out = J_*v_in;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::UserDefinedJacobian(const Eigen::VectorXd& y, const double t)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		this->Jacobian(y, t, J_);

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleJacobianAssembling_ = tend - tstart;
		cpuTimeToAssembleJacobian_ += cpuTimeSingleJacobianAssembling_;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::NumericalJacobian(Eigen::VectorXd& y, const double t, const Eigen::VectorXd& f, const double h, const Eigen::VectorXd& e,
															const bool max_constraints, const Eigen::VectorXd& yMax)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double BETA = 1.e+3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;


		double hf = BETA * std::fabs(h) * OpenSMOKE::ErrorControl(f, e) * double(this->ne_);
		if (hf < 1.e-10)
			hf = 1.;


		if (max_constraints == false)
		{
			for (unsigned int j = 0; j < this->ne_; j++)
			{
				const double yh = y(j);

				double hJ = ETA2*std::fabs(std::max(yh, 1. / e(j)));

				const double hJf = hf / e(j);

				hJ = std::max(hJ, hJf);

				hJ = std::max(hJ, ZERO_DER);

				hJ = std::min(hJ, 0.001 + 0.001*std::fabs(yh));

				y(j) += hJ;

				this->Equations(y, t, aux_);

				y(j) = yh;

				aux_ -= f;
				const double hInv = 1. / hJ;

				for (unsigned int i = 0; i < this->ne_; i++)
					Jaux_(i, j) = hInv*aux_(i);
			}
		}
		else
		{
			for (unsigned int j = 0; j < this->ne_; j++)
			{
				const double yh = y(j);

				double hJ = ETA2*OpenSMOKE::MaxAbs(yh, 1. / e(j));
				
				const double hJf = hf / e(j);
				
				hJ = std::max(hJ, hJf);
				
				hJ = std::max(hJ, ZERO_DER);
				
				if (yh + hJ > yMax(j))
					hJ = -hJ;
				
				y(j) += hJ;

				this->Equations(y, t, aux_);

				y(j) = yh;

				aux_ -= f;
				const double hInv = 1. / hJ;

				for (unsigned int i = 0; i < this->ne_; i++)
					Jaux_(i, j) = hInv*aux_(i);
			}
		}

		for (int k = 0; k<J_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(J_, k); it; ++it)
			{
				it.valueRef() = Jaux_(it.row(), it.col());
			}
			
		numberOfFunctionCallsForJacobian_ += this->ne_;

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleJacobianAssembling_ = tend - tstart;
		cpuTimeToAssembleJacobian_ += cpuTimeSingleJacobianAssembling_;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::BuildAndFactorizeMatrixG(const double hr0)
	{
		G_ = J_;
		G_ *= -hr0;
		G_ += ones_;

		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			sparse_LU_.compute(G_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_bicgstab_diagonal_.compute(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_bicgstab_ilut_.compute(G_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_gmres_diagonal_.compute(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_gmres_ilut_.compute(G_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_dgmres_diagonal_.compute(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_dgmres_ilut_.compute(G_);
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			sparse_pardiso_.compute(G_);
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			sparse_superlu_serial_.compute(G_);
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			sparse_umfpack_.compute(G_);
		}
		#endif
		#if OPENSMOKE_USE_LIS == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
		}
		#endif

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();
		cpuTimeSingleFactorization_ = tend - tstart;
		cpuTimeToFactorize_ += cpuTimeSingleFactorization_;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::SolveLinearSystem(Eigen::VectorXd& db)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		Eigen::VectorXd v = db;

		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			db = sparse_LU_.solve(v);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				db = sparse_bicgstab_diagonal_.solve(v);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				db = sparse_bicgstab_ilut_.solve(v);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				db = sparse_gmres_diagonal_.solve(v);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				db = sparse_gmres_ilut_.solve(v);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				db = sparse_dgmres_diagonal_.solve(v);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				db = sparse_dgmres_ilut_.solve(v);
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			db = sparse_pardiso_.solve(v);
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			db = sparse_superlu_serial_.solve(v);
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			db = sparse_umfpack_.solve(v);
		}
		#endif
		#if OPENSMOKE_USE_LIS == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			unsigned int count = 0;
			for (int k = 0; k < G_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(G_, k); it; ++it)
			{
				lis_value_[count] = it.value();
				count++;
			}

			lis_matrix_set_csc(G_.nonZeros(), lis_ptr_, lis_index_, lis_value_, lis_G_);
			lis_matrix_assemble(lis_G_);

			for (unsigned int i = 0; i < this->ne_; i++)
				lis_vector_set_value(LIS_INS_VALUE, i, v(i), lis_b_);
			for (unsigned int i = 0; i < this->ne_; i++)
				lis_vector_set_value(LIS_INS_VALUE, i, aux_(i), lis_x_);

			lis_solve(lis_G_, lis_b_, lis_x_, lis_solver_);

			for (unsigned int i = 0; i < this->ne_; i++)
			{	
				LIS_SCALAR dummy;
				lis_vector_get_value(lis_x_, i, &dummy);
				db(i) = dummy;
			}
		}
		#endif

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleLinearSystemSolution_ = tend - tstart;
		cpuTimeToSolveLinearSystem_ += cpuTimeSingleLinearSystemSolution_;
	}

	template <typename ODESystemObject>
	void KernelSparse<ODESystemObject>::OdeSolverKernelSummary(std::ostream& out)
	{
		double totalCpu = cpuTimeToAssembleJacobian_ + cpuTimeToFactorize_ + cpuTimeToSolveLinearSystem_;
		double totalSingleCpu = cpuTimeSingleFactorization_ + cpuTimeSingleJacobianAssembling_ + cpuTimeSingleLinearSystemSolution_;

		out << std::endl;
		out << "Data for the Dense ODE solver Kernel" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << "Number of function calls (only to assemble Jacobian): " << numberOfFunctionCallsForJacobian_ << " (" << numberOfFunctionCallsForJacobian_ / this->ne_ << ")" << std::endl;
		out << "Cumulative CPU time for assembling Jacobian:          " << cpuTimeToAssembleJacobian_ << " (" << cpuTimeToAssembleJacobian_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU time for LU decomposition:             " << cpuTimeToFactorize_ << " (" << cpuTimeToFactorize_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU time for solving the linear system:    " << cpuTimeToSolveLinearSystem_ << " (" << cpuTimeToSolveLinearSystem_ / totalCpu *100. << "%)" << std::endl;
		out << "CPU time for assembling Jacobian:                     " << cpuTimeSingleJacobianAssembling_ << " (" << cpuTimeSingleFactorization_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU time for LU decomposition:                        " << cpuTimeSingleFactorization_ << " (" << cpuTimeSingleJacobianAssembling_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU time for solving the linear system:               " << cpuTimeSingleLinearSystemSolution_ << " (" << cpuTimeSingleLinearSystemSolution_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << std::endl;
	}
}
