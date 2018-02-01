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

	ExtendedFallOff::ExtendedFallOff()
	{
	}
	
	ExtendedFallOff::ExtendedFallOff(const ExtendedFallOff& orig)
	{
	}
	
	//ExtendedFallOff::~ExtendedFallOff()
	//{
	//}

	void ExtendedFallOff::ErrorMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ExtendedFallOff"	<< std::endl;
		std::cout << "Error:  " << message									<< std::endl;
		std::cout << "Press a key to continue... "							<< std::endl;
		getchar();
		exit(-1);
	}

	void ExtendedFallOff::WarningMessage(const std::string message)
	{
		std::cout << std::endl;
		std::cout << "Class:  ExtendedFallOff"	<< std::endl;
		std::cout << "Warning:  "	<< message								<< std::endl;
		std::cout << "Press a key to continue... "							<< std::endl;
		getchar();
	}

	void ExtendedFallOff::Setup(const std::vector< std::vector<std::string> >& coefficients)
	{
		// Count number of entries
		n_ = 0;
		for (unsigned int i = 0; i < coefficients.size(); i++)
			if (coefficients[i][0] == "low")	n_++;

		// Conversion factors
		             conversion_A0_   = boost::lexical_cast<double>(coefficients[coefficients.size() - 1][0]);
		             conversion_AInf_ = boost::lexical_cast<double>(coefficients[coefficients.size() - 1][1]);
		const double conversion_E     = boost::lexical_cast<double>(coefficients[coefficients.size() - 1][2]);

		// Conversion factors
		AInf_ = boost::lexical_cast<double>(coefficients[coefficients.size() - 2][0]);
		BetaInf_ = boost::lexical_cast<double>(coefficients[coefficients.size() - 2][1]);
		EInf_over_R_ = boost::lexical_cast<double>(coefficients[coefficients.size() - 2][2]);

		// Third bodies
		for (unsigned int i = 0; i < coefficients[coefficients.size() - 3].size(); i++)
			third_body_indices_.push_back(boost::lexical_cast<int>(coefficients[coefficients.size() - 3][i]));
		for (unsigned int i = 0; i < coefficients[coefficients.size() - 4].size(); i++)
			third_body_efficiencies_.push_back(boost::lexical_cast<double>(coefficients[coefficients.size() - 4][i])-1.);

		// Create list of species
		species_.resize(0);
		for (unsigned int i = 0; i < coefficients.size()-4; i++)
		{
			bool found = false;
			for (unsigned int k = 0; k < species_.size(); k++)
				if (coefficients[i][1] == species_[k])
				{
					found = true;
					break;
				}

			if (found == false)
			{
				species_.push_back(coefficients[i][1]);
				species_indices_.push_back(boost::lexical_cast<int>(coefficients[i][2]));
			}
		}

		// Checking consistency of data
		{
			CheckForRedundancy(coefficients);
			CheckForMixtureLowPressureParameters(coefficients);
			CheckForLowPressureKineticParameters(coefficients);
			CheckForUnitaryThirdBodyEfficiencies();
		}

		// Memory allocation
		A0_.resize(n_);
		Beta0_.resize(n_);
		E0_over_R_.resize(n_);
		k_.resize(n_);
		teta_.resize(n_);
		type_.resize(n_);
		for(int i=0;i<n_;i++)
			type_[i] = EXTENDED_FALLOFF_LINDEMANN;

		// Low pressure values
		for (unsigned int i = 0; i < coefficients.size() - 4; i++)
		{
			if (coefficients[i][0] == "low")
			{
				for (unsigned int k = 0; k < species_.size(); k++)
					if (coefficients[i][1] == species_[k])
					{
						A0_[k] = boost::lexical_cast<double>(coefficients[i][3]);
						Beta0_[k] = boost::lexical_cast<double>(coefficients[i][4]);
						E0_over_R_[k] = boost::lexical_cast<double>(coefficients[i][5]);
					}
			}
			else
			{
				for (unsigned int k = 0; k < species_.size(); k++)
					if (coefficients[i][1] == species_[k])
					{
						if (coefficients[i][0] == "troe")
							type_[k] = EXTENDED_FALLOFF_TROE;
						else if (coefficients[i][0] == "sri")
							type_[k] = EXTENDED_FALLOFF_SRI;

						for (unsigned int j = 3; j<coefficients[i].size(); j++)
							teta_[k].push_back(boost::lexical_cast<double>(coefficients[i][j]));
					}
			}
		}

		// Conversions
		AInf_ *= conversion_AInf_;
		EInf_over_R_ *= conversion_E / PhysicalConstants::R_J_kmol;
		for (unsigned int k = 0; k < species_.size(); k++)
		{
			A0_[k] *= conversion_A0_;
			E0_over_R_[k] *= conversion_E / PhysicalConstants::R_J_kmol;
		}
	}

	double ExtendedFallOff::KineticConstant(const double T, const double P, const double cTot, const double* c)
	{
		// High pressure limit [kmol,m3,s]
		const double kInf = AInf_ * std::pow(T, BetaInf_)*std::exp(-EInf_over_R_ / T);

		// Third body concentration [kmol/m3]
		double M = cTot;
		for (unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
			M += third_body_efficiencies_[j] * c[third_body_indices_[j]];

		for (int i = 0; i < n_; i++)
		{
			// Low pressure limit [kmol,m3,s]
			const double k0 = A0_[i] * std::pow(T, Beta0_[i])*std::exp(-E0_over_R_[i] / T);

			// Reduced pressure
			const double Pr = k0 / kInf * M;

			// Function F
			double F = 1.;
			if (type_[i] == EXTENDED_FALLOFF_LINDEMANN)
			{
				// Do nothing (default values): F=1
			}
			else if (type_[i] == EXTENDED_FALLOFF_TROE)
			{
				double Fcent = (1. - teta_[i][0])*std::exp(-T / teta_[i][1]) + teta_[i][0] * std::exp(-T / teta_[i][2]);
				if (teta_[i].size() == 4)
					Fcent += std::exp(-teta_[i][3] / T);

				double logFcent = -300.;
				if (Fcent > 1e-300)
					logFcent = std::log10(Fcent);

				if (Pr > 1.e-32)
				{
					const double c = -0.4 - 0.67*logFcent;
					const double n = 0.75 - 1.27*logFcent;
					const double e = std::log10(Pr) + c;
					F = std::pow(10., logFcent / (1. + boost::math::pow<2>(e / (n - 0.14*e))));
				}
				else
				{
					F = std::pow(10., logFcent / (1. + boost::math::pow<2>(1./0.14)));
				}
			}
			else if (type_[i] == EXTENDED_FALLOFF_SRI)
			{
				const double X = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
				F = std::pow(teta_[i][0] * std::exp(-teta_[i][1] / T) + std::exp(-T / teta_[i][2]), X);
				if (teta_[i].size() == 5)
					F *= teta_[i][3] * std::pow(T, teta_[i][4]);
			}

			// Kinetic constant
			k_[i] = kInf * (Pr / (1. + Pr)) * F;
		}

		// Weight the kinetic constants
		double ctot_minus_cspecies = cTot;
		for (unsigned int i = 0; i < species_.size(); i++)
			if (species_indices_[i] != -1)	ctot_minus_cspecies -= c[species_indices_[i]];
		
		double kinetic_constant = 0.;
		for (unsigned int i = 0; i < species_.size(); i++)
		{
			if (species_indices_[i] != -1)
				kinetic_constant += c[species_indices_[i]] * k_[i];
			else
				kinetic_constant += ctot_minus_cspecies * k_[i];
		}
		kinetic_constant /= cTot;

		return kinetic_constant;
	}

	void ExtendedFallOff::ReadFromASCIIFile(std::istream& fInput)
	{
		double n;
		fInput >> n;

		std::vector< std::vector<std::string> > coefficients(static_cast<int>(n));
		for (unsigned int i = 0; i < static_cast<unsigned int>(n); i++)
		{
			double nn;
			fInput >> nn;

			coefficients[i].resize(static_cast<int>(nn));
			for (unsigned int j = 0; j < static_cast<unsigned int>(nn); j++)
				fInput >> coefficients[i][j];
		}

		Setup(coefficients);
	}

	void ExtendedFallOff::WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const
	{
		fOutput << std::setw(9) << " ";
		fOutput << std::setw(9) << std::left << "kInf:";
		fOutput << std::scientific << std::setprecision(6) << std::right << AInf_/conversion_AInf_ << "\t";
		fOutput << std::setw(8) << std::setprecision(2) << std::fixed << std::right << BetaInf_;
		fOutput << std::setw(14) << std::setprecision(2) << std::fixed << std::right << EInf_over_R_ * PhysicalConstants::R_J_kmol / Conversions::J_from_kcal << std::endl;

		for (unsigned int k = 0; k < species_.size(); k++)
			WriteShortSummaryOnASCIIFile(k, fOutput);
	}

	void ExtendedFallOff::WriteShortSummaryOnASCIIFile(const unsigned int index, std::ostream& fOutput) const
	{
		fOutput << std::setw(9) << " ";
		fOutput << "Bath species: " << species_[index] << std::endl;
		{
			fOutput << std::setw(12) << " ";
			fOutput << std::setw(9) << std::left << "k0:";
			fOutput << std::scientific << std::setprecision(6) << std::right << A0_[index]/conversion_A0_ << "\t";
			fOutput << std::setw(8) << std::setprecision(2) << std::fixed << std::right << Beta0_[index];
			fOutput << std::setw(14) << std::setprecision(2) << std::fixed << std::right << E0_over_R_[index]*PhysicalConstants::R_J_kmol / Conversions::J_from_kcal << std::endl;
		}

		if (type_[index] == EXTENDED_FALLOFF_TROE)
		{
			fOutput << std::setw(12) << " "; fOutput << "Troe Parameters" << std::endl;
			fOutput << std::setw(12) << " "; fOutput << std::scientific << "a    " << teta_[index][0] << std::endl;
			fOutput << std::setw(12) << " "; fOutput << std::scientific << "T*** " << teta_[index][1] << std::endl;
			fOutput << std::setw(12) << " "; fOutput << std::scientific << "T*   " << teta_[index][2] << std::endl;

			if (teta_[index].size() == 3)
			{
				fOutput << std::setw(12) << " "; fOutput << std::scientific << "T**  " << 0. << std::endl;
			}

			if (teta_[index].size() == 4)
			{
				fOutput << std::setw(12) << " "; fOutput << std::scientific << "T**  " << teta_[index][3] << std::endl;
			}
		}

		if (type_[index] == EXTENDED_FALLOFF_SRI)
		{
			fOutput << std::setw(12) << " "; fOutput << "SRI Parameters" << std::endl;
			fOutput << std::setw(12) << " "; fOutput << std::scientific << "a  " << teta_[index][0] << std::endl;
			fOutput << std::setw(12) << " "; fOutput << std::scientific << "b  " << teta_[index][1] << std::endl;
			fOutput << std::setw(12) << " "; fOutput << std::scientific << "c  " << teta_[index][2] << std::endl;

			if (teta_[index].size() == 3)
			{
				fOutput << std::setw(12) << " "; fOutput << std::scientific << "d  " << 1. << std::endl;
				fOutput << std::setw(12) << " "; fOutput << std::scientific << "e  " << 0. << std::endl;
			}

			if (teta_[index].size() == 5)
			{
				fOutput << std::setw(12) << " "; fOutput << std::scientific << "d  " << teta_[index][3] << std::endl;
				fOutput << std::setw(12) << " "; fOutput << std::scientific << "e  " << teta_[index][4] << std::endl;
			}
		}
	}

	void ExtendedFallOff::WriteCHEMKINOnASCIIFile(std::ostream& fOutput) const
	{
		for (unsigned int k = 0; k < species_.size(); k++)
		{
			// Low-pressure kinetic parameters
			{
				if (species_indices_[k] == -1)
					fOutput << "LOWMX/ ";
				else
				{
					fOutput << "LOWSP/ ";
					fOutput << species_[k];
				}

				fOutput.width(11);
				fOutput.unsetf(std::ios_base::floatfield);
				fOutput.precision(6);
				fOutput << std::showpoint << "   " << A0_[k]/conversion_A0_;
				fOutput << std::showpoint << "   " << Beta0_[k];
				fOutput << std::showpoint << "   " << E0_over_R_[k] * PhysicalConstants::R_J_kmol / Conversions::J_from_kcal;

				fOutput << "/" << std::endl;
			}

			// Attitional options
			{
				if (type_[k] == EXTENDED_FALLOFF_TROE)
				{
					if (species_indices_[k] == -1)
						fOutput << "TROEMX/";
					else
					{
						fOutput << "TROESP/";
						fOutput << species_[k];
					}
				}

				else if (type_[k] == EXTENDED_FALLOFF_SRI)
				{
					if (species_indices_[k] == -1)
						fOutput << "SRIMX/ ";
					else
					{
						fOutput << "SRISP/ ";
						fOutput << species_[k];
					}
				}

				if (type_[k] != EXTENDED_FALLOFF_LINDEMANN)
				{
					fOutput.width(11);
					fOutput.unsetf(std::ios_base::floatfield);
					for (unsigned int j = 0; j < teta_[k].size(); j++)
					{
						fOutput.precision(6);
						fOutput << std::showpoint << "   " << teta_[k][j];
					}

					fOutput << "/" << std::endl;
				}
			}
		}
	}

	// Check if redundant data are provided
	void ExtendedFallOff::CheckForRedundancy(const std::vector< std::vector<std::string> >& coefficients)
	{
		std::vector<int> counter(n_);
		std::fill(counter.begin(), counter.end(), 0);
		for (unsigned int i = 0; i < coefficients.size() - 4; i++)
		{
			for (unsigned int k = 0; k < species_.size(); k++)
				if (coefficients[i][1] == species_[k])
				{
					counter[k]++;
					break;
				}
		}

		std::vector<int>::iterator result = std::max_element(counter.begin(), counter.end());
		const int position = std::distance(counter.begin(), result);

		if (counter[position] > 2)
			ErrorMessage("Too many options were specified for species: " + species_[position]);
	}

	// Check if the low-pressure parameters for the mixture are provided
	void ExtendedFallOff::CheckForMixtureLowPressureParameters(const std::vector< std::vector<std::string> >& coefficients)
	{
		bool found = false;
		for (unsigned int i = 0; i < coefficients.size() - 4; i++)
		{
			if (coefficients[i][2] == "-1" && coefficients[i][0] == "low")
			{
				found = true;
				break;
			}
		}

		if (found == false)
			ErrorMessage("The low-pressure kinetic parameters have to be provided for the mixture (through LOWMX)");
	}

	// Check if the low-pressure parameters are provided for every species
	void ExtendedFallOff::CheckForLowPressureKineticParameters(const std::vector< std::vector<std::string> >& coefficients)
	{
		std::vector<bool> found(n_);
		std::fill(found.begin(), found.end(), false);
		for (unsigned int k = 0; k < species_.size(); k++)
			for (unsigned int i = 0; i < coefficients.size() - 4; i++)
			{
				if (coefficients[i][1] == species_[k] && coefficients[i][0] == "low")
				{
					found[k] = true;
					break;
				}
			}

		for (unsigned int k = 0; k < species_.size(); k++)
			if (found[k] == false)
			{
				ErrorMessage("Too low-pressure kinetic parameters were specified for species: " + species_[k]);
			}
	}

	// Check if the third-body efficiencies of specified species are equal to 1
	void ExtendedFallOff::CheckForUnitaryThirdBodyEfficiencies()
	{
		for (unsigned int k = 0; k < species_.size(); k++)
			for (unsigned int i = 0; i<third_body_indices_.size(); i++)
				if (species_indices_[k] == third_body_indices_[i])
				{
					if (third_body_efficiencies_[i] != 0.)
						ErrorMessage("The third-body efficiency for species " + species_[k] + " must be equal to 1");
				}
	}

	void ExtendedFallOff::WriteAdditionalParameters(std::ofstream& fOut, const int i)
	{
		if (type_[i] == EXTENDED_FALLOFF_TROE)
		{
			fOut << " Troe parameters:  " << std::endl;
			fOut << std::setw(6) << std::left << " a:" << teta_[i][0] << std::endl;
			fOut << std::setw(6) << std::left << " T***:" << teta_[i][1] << std::endl;
			fOut << std::setw(6) << std::left << " T*:" << teta_[i][2] << std::endl;
			if (teta_[i].size() == 4)
				fOut << std::setw(6) << std::left << " T**:" << teta_[i][3] << std::endl;
		}
		else if (type_[i] == EXTENDED_FALLOFF_SRI)
		{
			fOut << " SRI parameters:  " << std::endl;
			fOut << std::setw(6) << std::left << " a:" << teta_[i][0] << std::endl;
			fOut << std::setw(6) << std::left << " b:" << teta_[i][1] << std::endl;
			fOut << std::setw(6) << std::left << " c:" << teta_[i][2] << std::endl;
			if (teta_[i].size() == 5)
			{
				fOut << std::setw(6) << std::left << " d:" << teta_[i][3] << std::endl;
				fOut << std::setw(6) << std::left << " e:" << teta_[i][4] << std::endl;
			}
			else
			{
				fOut << std::setw(6) << std::left << " d:" << 1. << std::endl;
				fOut << std::setw(6) << std::left << " e:" << 0. << std::endl;
			}
		}
	}
}
