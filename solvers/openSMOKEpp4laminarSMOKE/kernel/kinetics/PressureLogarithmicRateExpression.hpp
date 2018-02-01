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

	PressureLogarithmicRateExpression::PressureLogarithmicRateExpression()
	{
	}
	
	PressureLogarithmicRateExpression::PressureLogarithmicRateExpression(const PressureLogarithmicRateExpression& orig)
	{
	}
	
	//PressureLogarithmicRateExpression::~PressureLogarithmicRateExpression()
	//{
	//}

	void PressureLogarithmicRateExpression::ErrorMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  PressureLogarithmicRateExpression"	<< std::endl;
		std::cout << "Error:  " << message							<< std::endl;
		std::cout << "Press a key to continue... "					<< std::endl;
		getchar();
		exit(-1);
	}

	void PressureLogarithmicRateExpression::WarningMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  PressureLogarithmicRateExpression"	<< std::endl;
		std::cout << "Warning:  "	<< message						<< std::endl;
		std::cout << "Press a key to continue... "					<< std::endl;
		getchar();
	}

	void PressureLogarithmicRateExpression::Setup(std::vector<double> coefficients_)
	{
		const double threshold = 1.e-32;

		N = (unsigned int)((coefficients_.size() - 2) / 4);
		lnA_.resize(N);
		Beta_.resize(N);
		E_over_R_.resize(N);
		p_.resize(N);
		lnp_.resize(N);

		double conversion_A_ = coefficients_[coefficients_.size()-2];
		double conversion_E_ = coefficients_[coefficients_.size()-1];

		unsigned count = 0;
		for(unsigned int i=0;i<N;i++)
		{
			if (coefficients_[count] < 0.)
				ErrorMessage("The PLOG option can be used only with non-negative frequency factors!");

			p_[i] = coefficients_[count++]*101325.;
			lnp_[i] = std::log(p_[i]);
			lnA_[i] = std::log(std::max(coefficients_[count++]*conversion_A_, threshold)) ;
			Beta_[i] = coefficients_[count++];
			E_over_R_[i] = (coefficients_[count++] * conversion_E_) / PhysicalConstants::R_J_kmol;
		}

		// Checking if the values are provided in the correct order
		for(unsigned int i=1;i<N;i++)
			if (p_[i] <= p_[i-1])
				ErrorMessage("The points on the pressure axis (PLOG) are not provided in the correct order!");
	}

	double PressureLogarithmicRateExpression::KineticConstant(const double T, const double P)
	{
		if	(P<=p_[0])	return std::exp(lnA_[0]+Beta_[0]*std::log(T)-E_over_R_[0]/T);
		else if (P>=p_[N-1])	return std::exp(lnA_[N-1]+Beta_[N-1]*std::log(T)-E_over_R_[N-1]/T);
		else
		{
			unsigned int i=0;
			for(i=0;i<N-1;i++)
				if (P<p_[i+1])
					break;
			
			double ln_kA = lnA_[i]+Beta_[i]*std::log(T)-E_over_R_[i]/T;
			double ln_kB = lnA_[i+1]+Beta_[i+1]*std::log(T)-E_over_R_[i+1]/T;
			return	std::exp( ln_kA+(ln_kB-ln_kA)*(std::log(P)-lnp_[i])/(lnp_[i+1]-lnp_[i]) );
		}
	}

	void PressureLogarithmicRateExpression::ReadFromASCIIFile(std::istream& fInput)
	{
		std::vector<double> coefficients;
		double n;
		fInput >> n;
		
		coefficients.resize(int(n));
		for(unsigned int i=0;i<(unsigned int)(n);i++)
			fInput >> coefficients[i];

		Setup(coefficients);
	}

	void PressureLogarithmicRateExpression::WriteOnASCIIFileOldStyle(std::ostream& fOutput) const
	{
		fOutput << N << std::endl;
		for (unsigned int j=0;j<N;j++)
			fOutput << p_[j]/101325. << std::endl;
		for (unsigned int j=0;j<N;j++)
			fOutput << std::exp(lnA_[j]) << std::endl;
		for (unsigned int j=0;j<N;j++)
			fOutput << Beta_[j] << std::endl;
		for (unsigned int j=0;j<N;j++)
			fOutput << E_over_R_[j] * PhysicalConstants::R_J_kmol/Conversions::J_from_kcal << std::endl;
	}

	void PressureLogarithmicRateExpression::WriteShortSummaryOnASCIIFile(std::ostream& fOutput, const double conversion_factor_A) const
	{
		fOutput << std::setw(9) << " "; 
		fOutput << "Pressure Logarithmic Interpolation" << std::endl;
		for (unsigned int j=0;j<N;j++)
		{
			fOutput << std::setw(9)		<< " "; 
			fOutput << j+1				<< "\t";
			fOutput << std::scientific << std::setprecision(6) << std::right << p_[j] / 101325. << "\t";
			fOutput << std::scientific << std::setprecision(6) << std::right << std::exp(lnA_[j])/conversion_factor_A << "\t";
			fOutput << std::setw(8) << std::setprecision(2) << std::fixed << std::right << Beta_[j];
			fOutput << std::setw(14) << std::setprecision(2) << std::fixed << std::right << E_over_R_[j] * PhysicalConstants::R_J_kmol / Conversions::J_from_kcal << std::endl;
		}
	}

}
