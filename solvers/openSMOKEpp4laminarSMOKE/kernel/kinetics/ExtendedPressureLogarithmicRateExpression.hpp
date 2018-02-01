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

	ExtendedPressureLogarithmicRateExpression::ExtendedPressureLogarithmicRateExpression()
	{
	}
	
	ExtendedPressureLogarithmicRateExpression::ExtendedPressureLogarithmicRateExpression(const ExtendedPressureLogarithmicRateExpression& orig)
	{
	}
	
	//ExtendedPressureLogarithmicRateExpression::~ExtendedPressureLogarithmicRateExpression()
	//{
	//}

	void ExtendedPressureLogarithmicRateExpression::ErrorMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ExtendedPressureLogarithmicRateExpression"	<< std::endl;
		std::cout << "Error:  " << message									<< std::endl;
		std::cout << "Press a key to continue... "							<< std::endl;
		getchar();
		exit(-1);
	}

	void ExtendedPressureLogarithmicRateExpression::WarningMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ExtendedPressureLogarithmicRateExpression"	<< std::endl;
		std::cout << "Warning:  "	<< message								<< std::endl;
		std::cout << "Press a key to continue... "							<< std::endl;
		getchar();
	}

	void ExtendedPressureLogarithmicRateExpression::Setup(std::vector<std::string> coefficients)
	{
		const unsigned int n = (unsigned int)((coefficients.size() - 2) / 6);

		             conversion_A_ = boost::lexical_cast<double>(coefficients[coefficients.size() - 2]);
		const double conversion_E  = boost::lexical_cast<double>(coefficients[coefficients.size() - 1]);

		// Resize species
		species_.resize(0);

		// Create list of species
		for (unsigned int i = 0; i < coefficients.size()-2; i += 6)
		{
			bool found = false;
			for (unsigned int k = 0; k < species_.size(); k++)
				if (coefficients[i] == species_[k])
				{
					found = true;
					break;
				}
			
			if (found == false)
			{
				species_.push_back(coefficients[i]);
				species_indices_.push_back(boost::lexical_cast<int>(coefficients[i + 1]));
			}
		}

		// Populate local coefficients
		std::vector< std::vector<double> > local_coefficients(species_.size());
		for (unsigned int i = 0; i < coefficients.size()-2; i += 6)
		{
			for (unsigned int k = 0; k < species_.size(); k++)
				if (coefficients[i] == species_[k])
				{
					local_coefficients[k].push_back(boost::lexical_cast<double>(coefficients[i + 2]));
					local_coefficients[k].push_back(boost::lexical_cast<double>(coefficients[i + 3]));
					local_coefficients[k].push_back(boost::lexical_cast<double>(coefficients[i + 4]));
					local_coefficients[k].push_back(boost::lexical_cast<double>(coefficients[i + 5]));
				}
		}

		// Check for data consistency
		{
			CheckForMixtureParameters();
		}

		// Memory allocation
		N_.resize(species_.size());
		lnA_.resize(species_.size());
		Beta_.resize(species_.size());
		E_over_R_.resize(species_.size());
		p_.resize(species_.size());
		lnp_.resize(species_.size());

		// Setup for single species
		for (unsigned int k = 0; k < species_.size(); k++)
			Setup(k, local_coefficients[k], conversion_E);
	}

	void ExtendedPressureLogarithmicRateExpression::Setup(unsigned int index, std::vector<double> coefficients, const double conversion_E)
	{
		N_[index] = (unsigned int)(coefficients.size()/4);
		
		lnA_[index].resize(N_[index]);
		Beta_[index].resize(N_[index]);
		E_over_R_[index].resize(N_[index]);
		p_[index].resize(N_[index]);
		lnp_[index].resize(N_[index]);

		const double threshold = 1.e-32;
		
		unsigned count = 0;
		for(int i=0;i<N_[index];i++)
		{
			if (coefficients[count] < 0.)
				ErrorMessage("The PLOGMX/PLOGSP option can be used only with non-negative frequency factors!");

			p_[index][i] = coefficients[count++]*101325.;
			lnp_[index][i] = std::log(p_[index][i]);
			lnA_[index][i] = std::log(std::max(coefficients[count++]*conversion_A_, threshold)) ;
			Beta_[index][i] = coefficients[count++];
			E_over_R_[index][i] = (coefficients[count++] * conversion_E) / PhysicalConstants::R_J_kmol;
		}

		// Checking if the values are provided in the correct order
		for(int i=1;i<N_[index];i++)
			if (p_[index][i] <= p_[index][i-1])
				ErrorMessage("The points on the pressure axis (PLOGMX/PLOGSP) are not provided in the correct order!");
	}

	double ExtendedPressureLogarithmicRateExpression::KineticConstant(const double T, const double P, const double cTot, const double* c)
	{
		double ctot_minus_cspecies = cTot;
		for (unsigned int k = 0; k < species_.size(); k++)
			if (species_indices_[k] != -1)	ctot_minus_cspecies -= c[species_indices_[k]];

		double kinetic_constant = 0.;
		for (unsigned int k = 0; k < species_.size(); k++)
		{
			if (species_indices_[k] != -1)
				kinetic_constant += c[species_indices_[k]] * KineticConstant(k, T, P);
			else
				kinetic_constant += ctot_minus_cspecies * KineticConstant(k, T, P);
		}
		kinetic_constant /= cTot;

		return kinetic_constant;
	}

	double ExtendedPressureLogarithmicRateExpression::KineticConstant(unsigned int index, const double T, const double P)
	{
		if (P <= p_[index][0])
		{
			return std::exp(lnA_[index][0] + Beta_[index][0] * std::log(T) - E_over_R_[index][0] / T);
		}
		else if (P >= p_[index][N_[index] - 1])
		{
			return std::exp(lnA_[index][N_[index] - 1] + Beta_[index][N_[index] - 1] * std::log(T) - E_over_R_[index][N_[index] - 1] / T);
		}
		else
		{
			int i=0;
			for(i=0;i<N_[index]-1;i++)
				if (P<p_[index][i+1])
					break;
			
			double ln_kA = lnA_[index][i]+Beta_[index][i]*std::log(T)-E_over_R_[index][i]/T;
			double ln_kB = lnA_[index][i+1]+Beta_[index][i+1]*std::log(T)-E_over_R_[index][i+1]/T;
			return	std::exp( ln_kA+(ln_kB-ln_kA)*(std::log(P)-lnp_[index][i])/(lnp_[index][i+1]-lnp_[index][i]) );
		}
	}

	void ExtendedPressureLogarithmicRateExpression::ReadFromASCIIFile(std::istream& fInput)
	{
		std::vector<std::string> coefficients;
		double n;
		fInput >> n;

		coefficients.resize(int(n));
		for (unsigned int i = 0; i<(unsigned int)(n); i++)
			fInput >> coefficients[i];

		Setup(coefficients);
	}

	void ExtendedPressureLogarithmicRateExpression::WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const
	{
		fOutput << std::setw(9) << " ";
		fOutput << "Extended Pressure Logarithmic Interpolation" << std::endl;
		for (unsigned int k = 0; k < species_.size(); k++)
			WriteShortSummaryOnASCIIFile(k, fOutput);
	}

	void ExtendedPressureLogarithmicRateExpression::WriteShortSummaryOnASCIIFile(const unsigned int index, std::ostream& fOutput) const
	{
		fOutput << std::setw(9) << " "; 
		fOutput << "Bath species: " << species_[index] << std::endl;
		for (int j=0;j<N_[index];j++)
		{
			fOutput << std::setw(9)									<< " "; 
			fOutput << std::scientific << j+1						<< "\t";
			fOutput << std::scientific << std::setprecision(6) << std::right << p_[index][j]/101325.     << "\t";
			fOutput << std::scientific << std::setprecision(6) << std::right << std::exp(lnA_[index][j])/conversion_A_ << "\t";
			fOutput << std::setw(8) << std::setprecision(2) << std::fixed << std::right << Beta_[index][j];
			fOutput << std::setw(14) << std::setprecision(2) << std::fixed << std::right << E_over_R_[index][j] * PhysicalConstants::R_J_kmol / Conversions::J_from_kcal << std::endl;
		}
	}

	void ExtendedPressureLogarithmicRateExpression::WriteCHEMKINOnASCIIFile(std::ostream& fOutput) const
	{
		fOutput.unsetf(std::ios_base::floatfield);
		fOutput.precision(6);

		for (unsigned int k = 0; k < species_.size(); k++)
		{
			for (int j = 0; j < N_[k]; j++)
			{
				if (species_indices_[k] == -1)
				{
					fOutput << " PLOGMX / ";
				}
				else
				{
					fOutput << " PLOGSP / ";
					fOutput << species_[k] << "  ";
				}

				fOutput << std::showpoint << std::setw(12) << std::left << p_[k][j] / 101325.;
				fOutput << std::showpoint << std::setw(12) << std::left << std::exp(lnA_[k][j]) / conversion_A_;
				fOutput << std::showpoint << std::setw(12) << std::left << Beta_[k][j];
				fOutput << std::showpoint << std::setw(12) << std::left << E_over_R_[k][j] * PhysicalConstants::R_J_kmol / Conversions::J_from_kcal;;

				fOutput << "/" << std::endl;
			}
		}
	}

	// Check if the kinetic parameters for the mixture are provided
	void ExtendedPressureLogarithmicRateExpression::CheckForMixtureParameters()
	{
		for (unsigned int k = 0; k < species_.size(); k++)
			if (species_indices_[k] == -1)
				return;

		ErrorMessage("The kinetic parameters have to be provided for the mixture (through PLOGMX)");
	}

}
