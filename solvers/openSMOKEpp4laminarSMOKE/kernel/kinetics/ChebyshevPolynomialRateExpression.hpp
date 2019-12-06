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

	ChebyshevPolynomialRateExpression::ChebyshevPolynomialRateExpression()
	{
		is_violation_allowed_ = false;
	}
	
	ChebyshevPolynomialRateExpression::ChebyshevPolynomialRateExpression(const ChebyshevPolynomialRateExpression& orig)
	{
	}
	
	//ChebyshevPolynomialRateExpression::~ChebyshevPolynomialRateExpression()
	//{
	//}

	void ChebyshevPolynomialRateExpression::ErrorMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ChebyshevPolynomialRateExpression"	<< std::endl;
		std::cout << "Error:  " << message							<< std::endl;
		std::cout << "Press a key to continue... "					<< std::endl;
		getchar();
		exit(-1);
	}

	void ChebyshevPolynomialRateExpression::WarningMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ChebyshevPolynomialRateExpression"	<< std::endl;
		std::cout << "Warning:  "	<< message						<< std::endl;
		std::cout << "Press a key to continue... "					<< std::endl;
		getchar();
	}

	void ChebyshevPolynomialRateExpression::Setup(std::vector<double> coefficients_, std::vector<double> pressures_, std::vector<double> temperatures_)
	{
		N = (unsigned int)(coefficients_[0]);
		M = (unsigned int)(coefficients_[1]);
		a.resize(N,M);
		phi_n.resize(N);
		phi_m.resize(M);

		unsigned int j=2;
		for(unsigned int n=0;n<N;n++)
			for(unsigned int m=0;m<M;m++)
				a(n,m) = coefficients_[j++];

		conversion = coefficients_[coefficients_.size()-1];

		Tmin = temperatures_[0];
		Tmax = temperatures_[1];
		Pmin = pressures_[0]*101325.;
		Pmax = pressures_[1]*101325.;

		log10_Pmin = std::log10(Pmin);
		log10_Pmax = std::log10(Pmax);
	}

	void ChebyshevPolynomialRateExpression::SetViolationAllowed(const bool flag)
	{
		is_violation_allowed_ = flag;
	}

	double ChebyshevPolynomialRateExpression::KineticConstant(const double T, const double P)
	{
		double Tc = T;
		double Pc = P;

		if (is_violation_allowed_ == false)
		{
			if (T > Tmax) 
			{ 
				std::cout << "Current T[K]: " << T << " - Maximum T[K]: " << Tmax << std::endl; 
				ErrorMessage("Temperature is larger than Tmax in Chebyshev Polynomial expression"); 
			}

			if (T < Tmin)
			{
				std::cout << "Current T[K]: " << T << " - Minimum T[K]: " << Tmin << std::endl;
				ErrorMessage("Temperature is lower than Tmin in Chebyshev Polynomial expression");
			}

			if (P > Pmax)
			{
				std::cout << "Current P[Pa]: " << P << " - Maximum P[Pa]: " << Pmax << std::endl;
				ErrorMessage("Pressure is larger than Pmax in Chebyshev Polynomial expression");
			}

			if (P < Pmin)
			{
				std::cout << "Current P[Pa]: " << P << " - Minimum P[Pa]: " << Pmin << std::endl;
				ErrorMessage("Pressure is lower than Pmin in Chebyshev Polynomial expression");
			}
		}
		else
		{
			if (T > Tmax) Tc = Tmax; 
			if (T < Tmin) Tc = Tmin;
			if (P > Pmax) Pc = Pmax;
			if (P < Pmin) Pc = Pmin;
		}

		double Ttilde = (2. / Tc - 1. / Tmin - 1. / Tmax) / (1. / Tmax - 1. / Tmin);
		double Ptilde = (2.*std::log10(Pc) - log10_Pmin - log10_Pmax) / (log10_Pmax - log10_Pmin);

		for (unsigned int n = 1; n <= N; n++)	phi_n(n - 1) = Phi(n, Ttilde);
		for (unsigned int m = 1; m <= M; m++)	phi_m(m - 1) = Phi(m, Ptilde);

		double sum = 0.;
		for (unsigned int n = 0; n < N; n++)
			for (unsigned int m = 0; m < M; m++)
				sum += a(n, m)*phi_n(n)*phi_m(m);

		return std::pow(10., sum) * conversion;
	}

	void ChebyshevPolynomialRateExpression::WriteStatus(std::ostream& fOutput) const
	{
		fOutput << " Chebyshev Polynomial: " << N << " x " << M << std::endl; 
		for(unsigned int n=0;n<N;n++)
		{
			for(unsigned int m=0;m<M;m++)
				fOutput << "   " << std::setw(12) << std::left << std::setprecision(6) << a(n,m);
			fOutput << std::endl;
		}

		fOutput << " Temperature limits [K]:  " << Tmin << " " << Tmax << std::endl;
		fOutput << " Pressure limits [atm]:   " << Pmin/101325. << " " << Pmax/101325. << std::endl;
	}

	void ChebyshevPolynomialRateExpression::ReadFromASCIIFile(std::istream& fInput)
	{
		double n;
		std::vector<double> coefficients;
		std::vector<double> pressures(2);
		std::vector<double> temperatures(2);

		fInput >> n;
		
		coefficients.resize(int(n));
		for(unsigned int i=0;i<(unsigned int)(n);i++)
			fInput >> coefficients[i];

		fInput >> pressures[0];
		fInput >> pressures[1];
		fInput >> temperatures[0];
		fInput >> temperatures[1];

		Setup(coefficients, pressures, temperatures);
	}

	void ChebyshevPolynomialRateExpression::WriteOnASCIIFileOldStyle(std::ostream& fOutput) const
	{
		fOutput << N << std::endl;
		fOutput << M << std::endl;
		for(unsigned int n=0;n<N;n++)
			for(unsigned int m=0;m<M;m++)
				fOutput << a(n,m) << std::endl;

		fOutput << 1./conversion << std::endl;	// TODO (units)

		fOutput << Tmin << std::endl << Tmax << std::endl;
		fOutput << Pmin/101325. << std::endl << Pmax/101325. << std::endl;
	}

	void ChebyshevPolynomialRateExpression::WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const
	{
		fOutput << std::setw(9) << " "; 
		fOutput << "Chebyshev Polynomial Reaction Rate" << std::endl;
		
		fOutput << std::setw(9) << " "; 
		fOutput << "Matrix of coefficients: " << N << " x " << M << std::endl;
		for(unsigned int n=0;n<N;n++)
		{
			fOutput << std::setw(9) << " "; 
			for(unsigned int m=0;m<M;m++)
				fOutput << std::setw(12) << std::setprecision(6) << std::left << a(n,m);
			fOutput << std::endl;
		}

		fOutput << std::setw(9) << " "; fOutput << "Conversion: " << 1./conversion	<< std::endl;	// TODO (units)
		fOutput << std::setw(9) << " "; fOutput << "Tmin: " << Tmin << " K" << std::endl;
		fOutput << std::setw(9) << " "; fOutput << "Tmax: " << Tmax << " K" << std::endl;
		fOutput << std::setw(9) << " "; fOutput << "Pmin: " << Pmin/101325. << " atm" << std::endl;	
		fOutput << std::setw(9) << " "; fOutput << "Pmax: " << Pmax/101325. << " atm" << std::endl;	
	}

	double ChebyshevPolynomialRateExpression::Phi(const int n, const double x)
	{
		return std::cos( double(n-1)*std::acos(x) );
	}

	void ChebyshevPolynomialRateExpression::FitttingArrheniusLaw(std::ostream& fOutput)
	{
		// Write on file
		{
			fOutput << "Temperature limits [K]:  " << Tmin << " " << Tmax << std::endl;
			fOutput << "Pressure limits [atm]:   " << Pmin / 101325. << " " << Pmax / 101325. << std::endl;
			fOutput << std::endl;

			fOutput << "Matrix of coefficients: " << N << " x " << M << std::endl;
			for (unsigned int n = 0; n < N; n++)
			{
				for (unsigned int m = 0; m < M; m++)
					fOutput << std::setw(16) << std::scientific << std::setprecision(6) << std::right << a(n, m);
				fOutput << std::endl;
			}
			fOutput << std::endl;
		}

		const unsigned int nT = 20;
		const unsigned int nP = 10;
		std::vector<double> list_T(nT);
		std::vector<double> P(nP);
		std::vector<double> logP(nP);

		// List of temperatures
		const double delta_T = (Tmax - Tmin - 40) / static_cast<double>(nT - 1);
		for (unsigned int i = 0; i < nT; i++)
			list_T[i] = (Tmin + 20) + i * delta_T;

		// List of pressures
		const double logPmin = std::log(Pmin);
		const double logPmax = std::log(Pmax);
		const double delta_logP = (logPmax - logPmin) / static_cast<double>(nP - 1);
		for (unsigned int i = 0; i < nP; i++)
			logP[i] = logPmin + i * delta_logP;

		P[0] = Pmin * 1.01;
		for (unsigned int i = 1; i < nP - 1; i++)
			P[i] = std::exp(logP[i]);
		P[nP - 1] = Pmax * 0.999;

		auto const it = std::lower_bound(logP.begin(), logP.end(), std::log(101325.));
		if (it != logP.end())
		{
			logP[std::distance(logP.begin(), it)] = std::log(101325.);
			P[std::distance(logP.begin(), it)] = 101325.;
		}

		// Two-parameter Arrhenius law
		{
			const unsigned int nparameters = 2;
			Eigen::Matrix2d XTX;
			Eigen::MatrixXd X(nT, nparameters);
			Eigen::MatrixXd XT(nparameters, nT);

			for (unsigned int i = 0; i < nT; i++)
			{
				X(i, 0) = 1.;
				X(i, 1) = 1. / list_T[i];
			}
			XT = X.transpose();
			XTX = XT * X;

			std::vector<double> A(nP);
			std::vector<double> E(nP);
			Eigen::VectorXd Y(nT);
			Eigen::VectorXd parameters(nparameters);
			for (unsigned int i = 0; i < nP; i++)
			{
				std::vector<double> list_k(nT);
				Eigen::VectorXd y(nT);
				for (unsigned int j = 0; j < nT; j++)
				{
					list_k[j] = KineticConstant(list_T[j], P[i]) / conversion;	// [mol, cm, s]
					y(j) = std::log(list_k[j]);
				}

				Y = XT * y;
				parameters = XTX.fullPivLu().solve(Y);

				A[i] = std::exp(parameters(0));						// [mol, cm, s]
				E[i] = -parameters(1)*PhysicalConstants::R_cal_mol;	// [cal/mol]
			}

			// Write in CHEMKIN format
			for (unsigned int i = 0; i < P.size(); i++)
			{
				fOutput << "PLOG / ";
				fOutput << std::setw(12) << std::scientific << std::setprecision(4) << std::right << std::exp(logP[i]) / 101325.;
				fOutput << std::setw(12) << std::scientific << std::setprecision(4) << std::right << A[i];
				fOutput << std::setw(7) << std::right << std::left << std::setprecision(2) << std::fixed << std::right << 0;
				fOutput << std::setw(12) << std::right << std::left << std::setprecision(2) << std::fixed << std::right << E[i];
				fOutput << " /" << std::endl;
			}
			fOutput << std::endl;
		}

		// Three-parameter Arrhenius law
		{
			const unsigned int nparameters = 3;
			Eigen::Matrix3d XTX;
			Eigen::MatrixXd X(nT, nparameters);
			Eigen::MatrixXd XT(nparameters, nT);

			for (unsigned int i = 0; i < nT; i++)
			{
				X(i, 0) = 1.;
				X(i, 1) = 1. / list_T[i];
				X(i, 2) = std::log(list_T[i]);
			}
			XT = X.transpose();
			XTX = XT * X;

			std::vector<double> A(nP);
			std::vector<double> E(nP);
			std::vector<double> Beta(nP);
			Eigen::VectorXd Y(nT);
			Eigen::VectorXd parameters(nparameters);
			for (unsigned int i = 0; i < nP; i++)
			{
				std::vector<double> list_k(nT);
				Eigen::VectorXd y(nT);
				for (unsigned int j = 0; j < nT; j++)
				{
					list_k[j] = KineticConstant(list_T[j], P[i]) / conversion;	// [mol, cm, s]
					y(j) = std::log(list_k[j]);
				}

				Y = XT * y;
				parameters = XTX.fullPivLu().solve(Y);

				A[i] = std::exp(parameters(0));						// [mol, cm, s]
				E[i] = -parameters(1)*PhysicalConstants::R_cal_mol;	// [cal/mol]
				Beta[i] = parameters(2);							// [-]
			}

			// Write in CHEMKIN format
			for (unsigned int i = 0; i < nP; i++)
			{
				fOutput << "PLOG / ";
				fOutput << std::setw(12) << std::scientific << std::setprecision(4) << std::right << std::exp(logP[i]) / 101325.;
				fOutput << std::setw(12) << std::scientific << std::setprecision(4) << std::right << A[i];
				fOutput << std::setw(7) << std::right << std::left << std::setprecision(2) << std::fixed << std::right << Beta[i];
				fOutput << std::setw(12) << std::right << std::left << std::setprecision(2) << std::fixed << std::right << E[i];
				fOutput << " /" << std::endl;
			}
			fOutput << std::endl;
		}
	}
}
