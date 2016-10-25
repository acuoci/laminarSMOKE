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

#include "math/OpenSMOKEUtilities.h"

namespace OpenSMOKE
{
	template<typename map> 
	SensitivityMap<map>::SensitivityMap(KineticsMap_CHEMKIN<map>& kinetics, const unsigned int number_of_equations) :
	kinetics_(kinetics), 
	number_of_equations_(number_of_equations)
	{
		// Default values
		sensitivity_type_			= PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR;
		energy_type_				= CONSTANT_VOLUME_SYSTEM;
		number_of_substeps_			= 2;

		dense_decomposition_type_	= DENSE_DECOMPOSITION_FULL_PIVOTING_LU;
		dense_solver_type_			= SOLVER_DENSE_EIGEN;

		number_of_species_ = boost::lexical_cast<unsigned int>(kinetics.NamesOfSpecies().size());
		number_of_reactions_ = boost::lexical_cast<unsigned int>(kinetics.NumberOfReactions());
		index_of_species_ = 1;			

		index_of_temperature_ = 0;			// means that there is no temperature equation
		index_of_density_ = 0;				// means that there is no density equation

		if (sensitivity_type_ == PhysicalConstants::SENSITIVITY_KINETIC_CONSTANT)
			number_of_parameters_ = kinetics.NumberOfReactions();
		else if (sensitivity_type_ == PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR)
			number_of_parameters_ = kinetics.NumberOfReactions() + kinetics.NumberOfFallOffReactions() + kinetics.NumberOfCABRReactions();

		ChangeDimensions(number_of_parameters_, &parameters_, true);
		ChangeDimensions(number_of_species_, &Jp_species_, true);
		
		sensitivity_coeffs_.resize(number_of_equations_, number_of_parameters_);
		sensitivity_coeffs_.setConstant(0.);

		Alfa_.resize(number_of_equations_, number_of_parameters_);
		Alfa_.setConstant(0.);

		B_eigen_.resize(number_of_equations_, number_of_parameters_);
		B_eigen_.setConstant(0.);
		
		// Dense formulation
		A_eigen_ = new Eigen::MatrixXd();

		// Sparse formulation
		A_eigen_sparse_ = new Eigen::SparseMatrix<double>();
		sparse_solver_type_ = SOLVER_SPARSE_EIGEN_SPARSE_LU;
		sparsePreconditionerType_ = PRECONDITIONER_SPARSE_DIAGONAL;
		sparse_preconditioner_droptol_ = 1e-6;
		sparse_preconditioner_fillfactor_ = 10;

		// Vectors
		b_eigen_ = new Eigen::VectorXd(number_of_equations_);
		b_eigen_->setConstant(0.);

		
		tOld_ = 0.;

		cpuTimeSingleFactorization_ = 0.;
		cpuTimeFactorization_		= 0.;
		cpuTimeSingleAssembling_    = 0.;
		cpuTimeAssembling_    		= 0.;
		cpuTimeSingleSolution_		= 0.;
		cpuTimeSolution_			= 0.;
	}

	template<typename map> 
	void SensitivityMap<map>::SetSensitivityType(const PhysicalConstants::sensitivity_type type)
	{
		sensitivity_type_ = type;
		if (sensitivity_type_ == PhysicalConstants::SENSITIVITY_KINETIC_CONSTANT)
			number_of_parameters_ = kinetics_.NumberOfReactions();
		else if (sensitivity_type_ == PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR)
			number_of_parameters_ = kinetics_.NumberOfReactions() + kinetics_.NumberOfFallOffReactions() + kinetics_.NumberOfCABRReactions();

		ChangeDimensions(number_of_parameters_, &parameters_, true);
		sensitivity_coeffs_.resize(number_of_equations_, number_of_parameters_);
		sensitivity_coeffs_.setConstant(0.);
	}

	template<typename map> 
	void SensitivityMap<map>::SetIndexOfTemperature(const unsigned int index)
	{
		index_of_temperature_ = index;
	}

	template<typename map> 
	void SensitivityMap<map>::SetIndexOfDensity(const unsigned int index)
	{
		index_of_density_ = index;
	}

	template<typename map> 
	void SensitivityMap<map>::SetDenseSolverType(const DenseSolverType type)
	{
		dense_solver_type_ = type;
	}

	template<typename map>
	void SensitivityMap<map>::SetSparseSolverType(const SparseSolverType type)
	{
		sparse_solver_type_ = type;
	}

	template<typename map> 
	void SensitivityMap<map>::SetDenseDecompositionType(const DenseDecompositionType type)
	{
		dense_decomposition_type_ = type;
	}

	template<typename map>
	void SensitivityMap<map>::SetSparsePreconditionerType(const SparsePreconditionerType type)
	{
		sparsePreconditionerType_ = type;
	}

	template<typename map>
	void SensitivityMap<map>::SetSparsePreconditionerDropTolerance(const double droptol)
	{
		sparse_preconditioner_droptol_ = droptol;
	}

	template<typename map>
	void SensitivityMap<map>::SetSparsePreconditionerFillFactor(const int fillfactor)
	{
		sparse_preconditioner_fillfactor_ = fillfactor;
	}

	template<typename map> 
	void SensitivityMap<map>::SetEnergyEquationType(const EnergyEquationType type)
	{
		energy_type_ = type;
	}

	template<typename map> 
	void SensitivityMap<map>::SetNumberOfSubSteps(const unsigned int number_of_substeps)
	{
		number_of_substeps_ = number_of_substeps;
	}

	template<typename map>
	void SensitivityMap<map>::SetSparsityPattern(const std::vector<unsigned int>& row, const std::vector<unsigned int>& col)
	{
		// Set size of sparse matrix
		A_eigen_sparse_->resize(number_of_equations_, number_of_equations_);

		// Main matrix
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList;
			tripletList.reserve(row.size() + number_of_equations_);
			for (unsigned int k = 0; k < row.size(); k++)
				tripletList.push_back(T(row[k], col[k], 1.));
			for (unsigned int k = 0; k < number_of_equations_; k++)
				tripletList.push_back(T(k, k, 1.));

			A_eigen_sparse_->setFromTriplets(tripletList.begin(), tripletList.end());
			A_eigen_sparse_->makeCompressed();
		}
		
		// Identity matrix
		{
			ones_.resize(number_of_equations_, number_of_equations_);
			ones_.setIdentity();
			ones_.makeCompressed();
		}

		// Analyze sparsity pattern
		if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			sparse_LU_.analyzePattern(*A_eigen_sparse_);
		}
		else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
			{
				sparse_bicgstab_diagonal_.analyzePattern(*A_eigen_sparse_);
			}
			else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_bicgstab_ilut_.preconditioner().setDroptol(sparse_preconditioner_droptol_);
				sparse_bicgstab_ilut_.preconditioner().setFillfactor(sparse_preconditioner_fillfactor_);
				sparse_bicgstab_ilut_.analyzePattern(*A_eigen_sparse_);
			}
		}
		else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
			{
				sparse_gmres_diagonal_.analyzePattern(*A_eigen_sparse_);
			}
			else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_gmres_ilut_.preconditioner().setDroptol(sparse_preconditioner_droptol_);
				sparse_gmres_ilut_.preconditioner().setFillfactor(sparse_preconditioner_fillfactor_);
				sparse_gmres_ilut_.analyzePattern(*A_eigen_sparse_);
			}
		}
		else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
			{
				sparse_dgmres_diagonal_.analyzePattern(*A_eigen_sparse_);
			}
			else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_dgmres_ilut_.preconditioner().setDroptol(sparse_preconditioner_droptol_);
				sparse_dgmres_ilut_.preconditioner().setFillfactor(sparse_preconditioner_fillfactor_);
				sparse_dgmres_ilut_.analyzePattern(*A_eigen_sparse_);
			}
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_PARDISO)
		{
			sparse_pardiso_.analyzePattern(*A_eigen_sparse_);
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			sparse_superlu_serial_.analyzePattern(*A_eigen_sparse_);
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			sparse_umfpack_.analyzePattern(*A_eigen_sparse_);
		}
		#endif
		#if OPENSMOKE_USE_LIS == 1
		if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_LIS)
		{
			// TODO
		}
		#endif
	}

	template<typename map> 
	void SensitivityMap<map>::Calculate(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J, const OpenSMOKE::OpenSMOKEVectorDouble& scaling_Jp)
	{
		t_ = t;

		const double u_number_of_substeps_ = 1./double(number_of_substeps_);
		const double deltat_ = t_-tOld_;
		const double u_number_of_substeps_deltat_ = u_number_of_substeps_ * deltat_;

		OpenSMOKE::OpenSMOKEVectorDouble mole_fractions(number_of_species_);
		const double cTot = P_Pa / T / PhysicalConstants::R_J_kmol;
		for (unsigned int i=1;i<=number_of_species_;i++)
			mole_fractions[i] = c[i]/cTot;

		if (dense_solver_type_ == SOLVER_DENSE_EIGEN)
		{
			// Factorize the matrix
			{
				const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

				if (A_eigen_->rows() == 0)
				A_eigen_->resize(number_of_equations_, number_of_equations_);

				const double* ptJ = J.GetHandle();
				for (unsigned int i = 0; i < number_of_equations_; i++)
				for (unsigned int j = 0; j < number_of_equations_; j++)
					(*A_eigen_)(i, j) = *ptJ++ * (-u_number_of_substeps_deltat_);
				for (unsigned int i = 0; i < number_of_equations_; i++)
					(*A_eigen_)(i, i) += 1.;

				if (dense_decomposition_type_ == DENSE_DECOMPOSITION_FULL_PIVOTING_LU)
					full_LU.compute(*A_eigen_);
				else if (dense_decomposition_type_ == DENSE_DECOMPOSITION_PARTIAL_PIVOTING_LU)
					partial_LU.compute(*A_eigen_);

				const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeSingleFactorization_ = tend - tstart;
				cpuTimeFactorization_ += cpuTimeSingleFactorization_;
			}

			// Assembling
			{
				const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

				for (unsigned int j = 1; j <= number_of_parameters_; j++)
				{
					double Jp_T_;
					double Jp_rho_;
					if (index_of_temperature_ == 0 && index_of_density_ == 0)
						kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, j, c, &Jp_species_, parameters_[j]);
					else if (index_of_temperature_ != 0 && index_of_density_ == 0)
						kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, energy_type_, j, c, mole_fractions, &Jp_species_, Jp_T_, parameters_[j]);
					else if (index_of_temperature_ != 0 && index_of_density_ != 0)
						kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, energy_type_, j, c, mole_fractions, &Jp_species_, Jp_T_, Jp_rho_, parameters_[j]);

					for (unsigned int i = 1; i <= number_of_species_; i++)
						(*b_eigen_)(i - 1) = Jp_species_[i] * scaling_Jp[i] * u_number_of_substeps_deltat_;

					if (index_of_temperature_ != 0)
						(*b_eigen_)(index_of_temperature_ - 1) = Jp_T_*scaling_Jp[index_of_temperature_] * u_number_of_substeps_deltat_;

					if (index_of_density_ != 0)
						(*b_eigen_)(index_of_density_ - 1) = Jp_T_*scaling_Jp[index_of_density_] * u_number_of_substeps_deltat_;

					Alfa_.col(j - 1) = (*b_eigen_);
				}

				const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeSingleAssembling_ = tend - tstart;
				cpuTimeAssembling_ += cpuTimeSingleAssembling_;
			}

			// Solution
			{
				const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

				for (unsigned int k = 1; k <= number_of_substeps_; k++)
				{
					B_eigen_ = Alfa_ + sensitivity_coeffs_;

					if (dense_decomposition_type_ == DENSE_DECOMPOSITION_FULL_PIVOTING_LU)
						sensitivity_coeffs_ = full_LU.solve(B_eigen_);
					else if (dense_decomposition_type_ == DENSE_DECOMPOSITION_PARTIAL_PIVOTING_LU)
						sensitivity_coeffs_ = partial_LU.solve(B_eigen_);
				}

				const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeSingleSolution_ = tend - tstart;
				cpuTimeSolution_ += cpuTimeSingleSolution_;
			}
		}
		else
		{
			OpenSMOKE::FatalErrorMessage("Sensitivity Analysis: The current version of OpenSMOKE was compiled only with Eigen support");
		}
		
		tOld_ = t_;
	}

	template<typename map>
	void SensitivityMap<map>::Calculate(const double t, const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const Eigen::SparseMatrix<double>& J, const OpenSMOKE::OpenSMOKEVectorDouble& scaling_Jp)
	{
		t_ = t;

		const double u_number_of_substeps_ = 1. / double(number_of_substeps_);
		const double deltat_ = t_ - tOld_;
		const double u_number_of_substeps_deltat_ = u_number_of_substeps_ * deltat_;

		OpenSMOKE::OpenSMOKEVectorDouble mole_fractions(number_of_species_);
		const double cTot = P_Pa / T / PhysicalConstants::R_J_kmol;
		for (unsigned int i = 1; i <= number_of_species_; i++)
			mole_fractions[i] = c[i] / cTot;

		// Factorize the matrix
		{
			const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

			*A_eigen_sparse_ = ones_ - u_number_of_substeps_deltat_*J;
			A_eigen_sparse_->makeCompressed();

			if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SPARSE_LU)
			{
				sparse_LU_.compute(*A_eigen_sparse_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_BICGSTAB)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sparse_bicgstab_diagonal_.compute(*A_eigen_sparse_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sparse_bicgstab_ilut_.compute(*A_eigen_sparse_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_GMRES)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sparse_gmres_diagonal_.compute(*A_eigen_sparse_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sparse_gmres_ilut_.compute(*A_eigen_sparse_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_DGMRES)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sparse_dgmres_diagonal_.compute(*A_eigen_sparse_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sparse_dgmres_ilut_.compute(*A_eigen_sparse_);
			}
			#if OPENSMOKE_USE_MKL == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_PARDISO)
			{
				sparse_pardiso_.compute(*A_eigen_sparse_);
			}
			#endif
			#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
			{
				sparse_superlu_serial_.compute(*A_eigen_sparse_);
			}
			#endif
			#if OPENSMOKE_USE_UMFPACK == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_UMFPACK)
			{
				sparse_umfpack_.compute(*A_eigen_sparse_);
			}
			#endif
			#if OPENSMOKE_USE_LIS == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_LIS)
			{
				// TODO
			}
			#endif

			const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeSingleFactorization_ = tend - tstart;
			cpuTimeFactorization_ += cpuTimeSingleFactorization_;
		}

		// Assembling
		{
			const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

			for (unsigned int j = 1; j <= number_of_parameters_; j++)
			{
				double Jp_T_;
				double Jp_rho_;
				if (index_of_temperature_ == 0 && index_of_density_ == 0)
					kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, j, c, &Jp_species_, parameters_[j]);
				else if (index_of_temperature_ != 0 && index_of_density_ == 0)
					kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, energy_type_, j, c, mole_fractions, &Jp_species_, Jp_T_, parameters_[j]);
				else if (index_of_temperature_ != 0 && index_of_density_ != 0)
					kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, energy_type_, j, c, mole_fractions, &Jp_species_, Jp_T_, Jp_rho_, parameters_[j]);

				for (unsigned int i = 1; i <= number_of_species_; i++)
					(*b_eigen_)(i - 1) = Jp_species_[i] * scaling_Jp[i] * u_number_of_substeps_deltat_;

				if (index_of_temperature_ != 0)
					(*b_eigen_)(index_of_temperature_ - 1) = Jp_T_*scaling_Jp[index_of_temperature_] * u_number_of_substeps_deltat_;

				if (index_of_density_ != 0)
					(*b_eigen_)(index_of_density_ - 1) = Jp_T_*scaling_Jp[index_of_density_] * u_number_of_substeps_deltat_;

				Alfa_.col(j - 1) = (*b_eigen_);
			}

			const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeSingleAssembling_ = tend - tstart;
			cpuTimeAssembling_ += cpuTimeSingleAssembling_;
		}
			
		// Solution
		{
			const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

			for (unsigned int k = 1; k <= number_of_substeps_; k++)
			{
				B_eigen_ = Alfa_ + sensitivity_coeffs_;

				if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SPARSE_LU )
				{
					sensitivity_coeffs_ = sparse_LU_.solve(B_eigen_);
				}
				else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_BICGSTAB)
				{
					if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
						sensitivity_coeffs_ = sparse_bicgstab_diagonal_.solve(B_eigen_);
					else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
						sensitivity_coeffs_ = sparse_bicgstab_ilut_.solve(B_eigen_);
				}
				else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_GMRES)
				{
					if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
						sensitivity_coeffs_ = sparse_gmres_diagonal_.solve(B_eigen_);
					else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
						sensitivity_coeffs_ = sparse_gmres_ilut_.solve(B_eigen_);
				}
				else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_DGMRES)
				{
					if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
						sensitivity_coeffs_ = sparse_dgmres_diagonal_.solve(B_eigen_);
					else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
						sensitivity_coeffs_ = sparse_dgmres_ilut_.solve(B_eigen_);
				}
				#if OPENSMOKE_USE_MKL == 1
				else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_PARDISO)
				{
					sensitivity_coeffs_ = sparse_pardiso_.solve(B_eigen_);
				}
				#endif
				#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
				else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
				{
					sensitivity_coeffs_ = sparse_superlu_serial_.solve(B_eigen_);
				}
				#endif
				#if OPENSMOKE_USE_UMFPACK == 1
				else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_UMFPACK)
				{
					sensitivity_coeffs_ = sparse_umfpack_.solve(B_eigen_);
				}
				#endif
				#if OPENSMOKE_USE_LIS == 1
				else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_LIS)
				{
					// TODO
				}
				#endif
			}

			const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeSingleSolution_ = tend - tstart;
			cpuTimeSolution_ += cpuTimeSingleSolution_;
		}

		tOld_ = t_;
	}

	template<typename map> 
	void SensitivityMap<map>::Calculate(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const OpenSMOKE::OpenSMOKEMatrixDouble& J, const OpenSMOKE::OpenSMOKEVectorDouble& scaling_Jp)
	{
		OpenSMOKE::OpenSMOKEVectorDouble mole_fractions(number_of_species_);
		const double cTot = P_Pa / T / PhysicalConstants::R_J_kmol;
		for (unsigned int i=1;i<=number_of_species_;i++)
			mole_fractions[i] = c[i]/cTot;

		if (dense_solver_type_ == SOLVER_DENSE_EIGEN)
		{
			// Factorize the matrix
			{
				const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

				if (A_eigen_->rows() == 0)
					A_eigen_->resize(number_of_equations_, number_of_equations_);

				const double* ptJ = J.GetHandle();
				for (unsigned int i = 0; i < number_of_equations_; i++)
				for (unsigned int j = 0; j < number_of_equations_; j++)
					(*A_eigen_)(i, j) = *ptJ++;

				if (dense_decomposition_type_ == DENSE_DECOMPOSITION_FULL_PIVOTING_LU)
					full_LU.compute(*A_eigen_);
				else if (dense_decomposition_type_ == DENSE_DECOMPOSITION_PARTIAL_PIVOTING_LU)
					partial_LU.compute(*A_eigen_);

				const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeSingleFactorization_ = tend - tstart;
				cpuTimeFactorization_ += cpuTimeSingleFactorization_;
			}

			// Assembling
			{
				const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

				for (unsigned int j = 1; j <= number_of_parameters_; j++)
				{
					double Jp_T_;
					if (index_of_temperature_ == 0 && index_of_density_ == 0)
						kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, j, c, &Jp_species_, parameters_[j]);
					else if (index_of_temperature_ != 0 && index_of_density_ == 0)
						kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, energy_type_, j, c, mole_fractions, &Jp_species_, Jp_T_, parameters_[j]);

					for (unsigned int i = 1; i <= number_of_species_; i++)
						(*b_eigen_)(i - 1) = -Jp_species_[i] * scaling_Jp[i];

					if (index_of_temperature_ != 0)
						(*b_eigen_)(index_of_temperature_ - 1) = -Jp_T_*scaling_Jp[index_of_temperature_];

					Alfa_.col(j - 1) = (*b_eigen_);
				}

				const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeSingleAssembling_ = tend - tstart;
				cpuTimeAssembling_ += cpuTimeSingleAssembling_;
			}

			// Solution
			{
				const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

				if (dense_decomposition_type_ == DENSE_DECOMPOSITION_FULL_PIVOTING_LU)
					sensitivity_coeffs_ = full_LU.solve(Alfa_);
				else if (dense_decomposition_type_ == DENSE_DECOMPOSITION_PARTIAL_PIVOTING_LU)
					sensitivity_coeffs_ = partial_LU.solve(Alfa_);

				const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuTimeSingleSolution_ = tend - tstart;
				cpuTimeSolution_ += cpuTimeSingleSolution_;
			}

			/*
			const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

			for (unsigned int j = 1; j <= number_of_parameters_; j++)
			{
				double Jp_T_;
				if (index_of_temperature_ == 0 && index_of_density_ == 0)
					kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, j, c, &Jp_species_, parameters_[j]);
				else if (index_of_temperature_ != 0 && index_of_density_ == 0)
					kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, energy_type_, j, c, mole_fractions, &Jp_species_, Jp_T_, parameters_[j]);

				for (unsigned int i = 1; i <= number_of_species_; i++)
					(*b_eigen_)(i - 1) = -Jp_species_[i] * scaling_Jp[i];

				if (index_of_temperature_ != 0)
					(*b_eigen_)(index_of_temperature_ - 1) = -Jp_T_*scaling_Jp[index_of_temperature_];

				if (dense_decomposition_type_ == DENSE_DECOMPOSITION_FULL_PIVOTING_LU)
					*x_eigen_ = full_LU.solve(*b_eigen_);
				else if (dense_decomposition_type_ == DENSE_DECOMPOSITION_PARTIAL_PIVOTING_LU)
					*x_eigen_ = partial_LU.solve(*b_eigen_);

				for (unsigned int i = 1; i <= number_of_equations_; i++)
					sensitivity_coefficients_[i][j] = (*x_eigen_)(i - 1);
			}

			const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeSingleSolution_ = tend - tstart;
			cpuTimeSolution_ += cpuTimeSingleSolution_;
			*/
		}
		else
		{
			OpenSMOKE::FatalErrorMessage("Sensitivity Analysis: The current version of OpenSMOKE was compiled with Eigen support");
		}
	}

	template<typename map>
	void SensitivityMap<map>::Calculate(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c, const Eigen::SparseMatrix<double>& J, const OpenSMOKE::OpenSMOKEVectorDouble& scaling_Jp)
	{
		OpenSMOKE::OpenSMOKEVectorDouble mole_fractions(number_of_species_);
		const double cTot = P_Pa / T / PhysicalConstants::R_J_kmol;
		for (unsigned int i = 1; i <= number_of_species_; i++)
			mole_fractions[i] = c[i] / cTot;

		// Factorize the matrix
		{
			const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

			*A_eigen_sparse_ = J;
			
			if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SPARSE_LU)
			{
				sparse_LU_.compute(*A_eigen_sparse_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_BICGSTAB)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sparse_bicgstab_diagonal_.compute(*A_eigen_sparse_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sparse_bicgstab_ilut_.compute(*A_eigen_sparse_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_GMRES)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sparse_gmres_diagonal_.compute(*A_eigen_sparse_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sparse_gmres_ilut_.compute(*A_eigen_sparse_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_DGMRES)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sparse_dgmres_diagonal_.compute(*A_eigen_sparse_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sparse_dgmres_ilut_.compute(*A_eigen_sparse_);
			}
			#if OPENSMOKE_USE_MKL == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_PARDISO)
			{
				sparse_pardiso_.compute(*A_eigen_sparse_);
			}
			#endif
			#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
			{
				sparse_superlu_serial_.compute(*A_eigen_sparse_);
			}
			#endif
			#if OPENSMOKE_USE_UMFPACK == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_UMFPACK)
			{
				sparse_umfpack_.compute(*A_eigen_sparse_);
			}
			#endif
			#if OPENSMOKE_USE_LIS == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_LIS)
			{
				// TODO
			}
			#endif

			const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeSingleFactorization_ = tend - tstart;
			cpuTimeFactorization_ += cpuTimeSingleFactorization_;
		}

		// Assembling
		{
			const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

			for (unsigned int j = 1; j <= number_of_parameters_; j++)
			{
				double Jp_T_;
				if (index_of_temperature_ == 0 && index_of_density_ == 0)
					kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, j, c, &Jp_species_, parameters_[j]);
				else if (index_of_temperature_ != 0 && index_of_density_ == 0)
					kinetics_.SensitivityWithRespectKineticParameter(sensitivity_type_, energy_type_, j, c, mole_fractions, &Jp_species_, Jp_T_, parameters_[j]);

				for (unsigned int i = 1; i <= number_of_species_; i++)
					(*b_eigen_)(i - 1) = -Jp_species_[i] * scaling_Jp[i];

				if (index_of_temperature_ != 0)
					(*b_eigen_)(index_of_temperature_ - 1) = -Jp_T_*scaling_Jp[index_of_temperature_];

				Alfa_.col(j - 1) = (*b_eigen_);
			}

			const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeSingleAssembling_ = tend - tstart;
			cpuTimeAssembling_ += cpuTimeSingleAssembling_;
		}

		// Solution
		{
			const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

			if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SPARSE_LU )
			{
				sensitivity_coeffs_ = sparse_LU_.solve(Alfa_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_BICGSTAB)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sensitivity_coeffs_ = sparse_bicgstab_diagonal_.solve(Alfa_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sensitivity_coeffs_ = sparse_bicgstab_ilut_.solve(Alfa_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_GMRES)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sensitivity_coeffs_ = sparse_gmres_diagonal_.solve(Alfa_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sensitivity_coeffs_ = sparse_gmres_ilut_.solve(Alfa_);
			}
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_DGMRES)
			{
				if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_DIAGONAL)
					sensitivity_coeffs_ = sparse_dgmres_diagonal_.solve(Alfa_);
				else if (sparsePreconditionerType_ == PRECONDITIONER_SPARSE_ILUT)
					sensitivity_coeffs_ = sparse_dgmres_ilut_.solve(Alfa_);
			}
			#if OPENSMOKE_USE_MKL == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_PARDISO)
			{
				sensitivity_coeffs_ = sparse_pardiso_.solve(Alfa_);
			}
			#endif
			#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
			{
				sensitivity_coeffs_ = sparse_superlu_serial_.solve(Alfa_);
			}
			#endif
			#if OPENSMOKE_USE_UMFPACK == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_UMFPACK)
			{
				sensitivity_coeffs_ = sparse_umfpack_.solve(Alfa_);
			}
			#endif
			#if OPENSMOKE_USE_LIS == 1
			else if (sparse_solver_type_ == SOLVER_SPARSE_EIGEN_LIS)
			{
				// TODO
			}
			#endif
				

			const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

			cpuTimeSingleSolution_ = tend - tstart;
			cpuTimeSolution_ += cpuTimeSingleSolution_;
		}
	}
}
