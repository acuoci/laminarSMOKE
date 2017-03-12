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

	double ChebyshevPolynomialRateExpression::KineticConstant(const double T, const double P)
	{
	//	if (T>Tmax)	{ cout << T << " " << Tmax << endl; ErrorMessage("Temperature is larger than Tmax");}
	//	if (T<Tmin)	{ cout << T << " " << Tmin << endl; ErrorMessage("Temperature is smaller than Tmin");}
	//	if (P>Pmax)	{ cout << P << " " << Pmax << endl; ErrorMessage("Pressure is larger than Pmax");}
	//	if (P<Pmin)	{ cout << P << " " << Pmin << endl; ErrorMessage("Pressure is smaller than Pmin");}

		double Ttilde = (2./T-1./Tmin-1./Tmax) / (1./Tmax-1./Tmin);
		double Ptilde = (2.*std::log10(P)-log10_Pmin-log10_Pmax) / (log10_Pmax-log10_Pmin);

		for(unsigned int n=1;n<=N;n++)	phi_n(n-1) = Phi(n, Ttilde);
		for(unsigned int m=1;m<=M;m++)	phi_m(m-1) = Phi(m, Ptilde);

		double sum = 0.;
		for(unsigned int n=0;n<N;n++)
			for(unsigned int m=0;m<M;m++)	
				sum += a(n,m)*phi_n(n)*phi_m(m);

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

}
