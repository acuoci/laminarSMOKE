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
	KernelDense<ODESystemObject>::KernelDense()
	{
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::ResetKernel()
	{
		solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		full_pivoting_ = false;

		numberOfFunctionCallsForJacobian_ = 0;

		cpuTimeToAssembleJacobian_ = 0.;
		cpuTimeToFactorize_ = 0.;
		cpuTimeToSolveLinearSystem_ = 0.;

		cpuTimeSingleJacobianAssembling_ = 0.;
		cpuTimeSingleFactorization_ = 0.;
		cpuTimeSingleLinearSystemSolution_ = 0.;
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::MemoryAllocationKernel()
	{
		// Allocate memory (if needed) for the OdeSystem object
		this->MemoryAllocation();

		// Internal variables
		aux_.resize(this->ne_);
		J_.resize(this->ne_, this->ne_);
		G_.resize(this->ne_, this->ne_);
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::SetLinearAlgebraSolver(const std::string linear_algebra_solver)
	{
		if (linear_algebra_solver == "Eigen")
		{
			solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		}
		else
			OpenSMOKE::ErrorMessage("KernelDense<ODESystemObject>", "Requested linear algebra is not supported!");
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::SetFullPivoting(const bool flag)
	{
		full_pivoting_ = flag;
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::JacobianTimesVector(const Eigen::VectorXd& v_in, Eigen::VectorXd* v_out)
	{
		*v_out = J_*v_in;
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::UserDefinedJacobian(const Eigen::VectorXd& y, const double t)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		// TODO
		// this->Jacobian(y, t, J_);
		OpenSMOKE::ErrorMessage("KernelDense<ODESystemObject>", "User defined Jacobian is still not supported!");

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleJacobianAssembling_ = tend - tstart;
		cpuTimeToAssembleJacobian_ += cpuTimeSingleJacobianAssembling_;
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::NumericalJacobian(Eigen::VectorXd& y, const double t, const Eigen::VectorXd& f, const double h, const Eigen::VectorXd& e,
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
					J_(i, j) = hInv*aux_(i);
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
					J_(i, j) = hInv*aux_(i);
			}
		}

		numberOfFunctionCallsForJacobian_ += this->ne_;

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleJacobianAssembling_ = tend - tstart;
		cpuTimeToAssembleJacobian_ += cpuTimeSingleJacobianAssembling_;
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::BuildAndFactorizeMatrixG(const double hr0)
	{
		G_ = J_;
		G_ *= -hr0;
		for (unsigned int i = 0; i < this->ne_; i++)
			G_(i, i) += 1.;

		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (solverType_ == OpenSMOKE::SOLVER_DENSE_EIGEN)
		{
			if (full_pivoting_ == true)
			{
				full_LU_.compute(G_);
			}
			else
			{
				partial_LU_.compute(G_);
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();
		cpuTimeSingleFactorization_ = tend - tstart;
		cpuTimeToFactorize_ += cpuTimeSingleFactorization_;
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::SolveLinearSystem(Eigen::VectorXd& db)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (solverType_ == OpenSMOKE::SOLVER_DENSE_EIGEN)
		{
			Eigen::VectorXd v = db;

			if (full_pivoting_ == true)
			{
				db = full_LU_.solve(v);
			}
			else
			{
				db = partial_LU_.solve(v);
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleLinearSystemSolution_ = tend - tstart;
		cpuTimeToSolveLinearSystem_ += cpuTimeSingleLinearSystemSolution_;
	}

	template <typename ODESystemObject>
	void KernelDense<ODESystemObject>::OdeSolverKernelSummary(std::ostream& out)
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
		out << "CPU time for assembling Jacobian:                     " << cpuTimeSingleFactorization_ << " (" << cpuTimeSingleFactorization_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU time for LU decomposition:                        " << cpuTimeSingleJacobianAssembling_ << " (" << cpuTimeSingleJacobianAssembling_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU time for solving the linear system:               " << cpuTimeSingleLinearSystemSolution_ << " (" << cpuTimeSingleLinearSystemSolution_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << std::endl;
	}
}
