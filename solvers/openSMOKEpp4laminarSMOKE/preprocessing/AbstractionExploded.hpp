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
|   Copyright(C) 2018  Alberto Cuoci                                      |
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
	AbstractionExploded::AbstractionExploded(const std::vector<std::string>& species_R, const std::string& species_RpH,
		const std::string& species_Rp, const std::vector<std::string>& species_RH,
		const double n_H, const unsigned int type_H)
	{
		species_R_ = species_R;
		species_RpH_ = species_RpH;
		species_Rp_ = species_Rp;
		species_RH_ = species_RH;
		n_H_ = n_H;
		type_H_ = type_H;

		n_ = species_R_.size();
		A_.resize(n_);
		Beta_.resize(n_);
		E_.resize(n_);
		is_dummy_.resize(n_);
		is_existing_.resize(n_);

		std::vector<int> species_R_length(n_);
		std::vector<int> species_RH_length(n_);
		for (unsigned int i = 0; i < n_; i++)
		{
			species_R_length[i] = species_R_[i].length();
			species_RH_length[i] = species_RH_[i].length();
		}

		const std::vector<int>::iterator length_R_ = std::max_element(species_R_length.begin(), species_R_length.end());
		const std::vector<int>::iterator length_RH_ = std::max_element(species_RH_length.begin(), species_RH_length.end());

		for (unsigned int i = 0; i < n_; i++)
		{
			std::stringstream	dummy;
			dummy << std::setw(*length_R_) << std::left << species_R_[i];
			dummy << std::setw(3) << " + ";
			dummy << std::setw(species_RpH_.length()) << std::left << species_RpH_;

			dummy << std::setw(4) << " => ";

			dummy << std::setw(species_Rp_.length()) << std::left << species_Rp_;
			dummy << std::setw(3) << " + ";
			dummy << std::setw(*length_RH_) << std::left << species_RH_[i];

			label_.push_back(dummy.str());
		}

		// Remove dummy reactions
		for (unsigned int i = 0; i < species_R_.size(); i++)
		{
			if (species_R_[i] == species_Rp_ && species_RpH_ == species_RH_[i])
				is_dummy_[i] = true;
			else
				is_dummy_[i] = false;
		}

		// Reactions already available in the kinetic mechanism
		for (unsigned int i = 0; i < species_R_.size(); i++)
			is_existing_[i] = false;
	}

	void AbstractionExploded::PrintExplodedReactions(std::ostream& out) const
	{
		out << std::endl;
		out << "! ABSTRACTION REACTION: R + " << species_RpH_ << " => " << species_Rp_ + " + RH " << std::endl;
		out << "! Number of H abstracted: " << std::setprecision(2) << std::fixed << n_H_ << std::endl;
		out << "! Type of H abstracted:   " << "[" << type_H_ << "] " << std::endl;

		for (unsigned int i = 0; i < n_; i++)
		{
			if (is_dummy_[i] == false && is_existing_[i] == false)
			{
				out << "  ";
				out << label_[i];
				out << std::setw(16) << std::right << std::scientific << std::setprecision(5) << A_[i];
				out << std::setw(12) << std::right << std::fixed << std::setprecision(3) << Beta_[i];
				out << std::setw(12) << std::right << std::fixed << std::setprecision(3) << E_[i];
				out << std::endl;
			}

			else if (is_dummy_[i] == true || is_existing_[i] == true)
			{
				out << "! ";
				out << label_[i];
				out << std::setw(16) << std::right << std::scientific	<< std::setprecision(5) << A_[i];
				out << std::setw(12) << std::right << std::fixed		<< std::setprecision(3) << Beta_[i];
				out << std::setw(12) << std::right << std::fixed		<< std::setprecision(3) << E_[i];

				if (is_existing_[i] == true)	out << "    OVERWRITTEN";
				else if (is_dummy_[i] == true)	out << "    DUMMY";

				out << std::endl;
			}
		}

		out << std::endl;
	}

	unsigned int AbstractionExploded::RemoveExistingReactions(const std::vector<std::string>& reactant_species, const std::vector<std::string>& product_species)
	{
		unsigned int n_existing = 0;
		if (std::find(reactant_species.begin(), reactant_species.end(), species_RpH_) != reactant_species.end())
		{
			if (std::find(product_species.begin(), product_species.end(), species_Rp_) != product_species.end())
			{
				for (unsigned int i = 0; i < n_; i++)
				{
					if (std::find(reactant_species.begin(), reactant_species.end(), species_R_[i]) != reactant_species.end())
					{
						if (std::find(product_species.begin(), product_species.end(), species_RH_[i]) != product_species.end())
						{
							is_existing_[i] = true;
							n_existing++;
						}
					}
				}
			}
		}

		return n_existing;
	}

	void AbstractionExploded::ExplodeReactions(	const std::vector<std::string>& abstractors,
												const std::vector<double>& abstractors_A,
												const std::vector<double>& abstractors_Beta,
												const std::vector<double>& abstractors_E,
												const std::vector<double>& corrections_A,
												const std::vector<double>& corrections_E,
												const std::vector<double>& corrections_Er,
												const double Eref, const double alpha, const double Tref)
	{
		int index_Rprime = -1;
		std::vector<std::string>::const_iterator iRprime = std::find(abstractors.begin(), abstractors.end(), species_Rp_);
		if (iRprime == abstractors.end())
		{
			const std::string label = "R + " + species_RpH_ + " => " + species_Rp_ + " + RH ";
			std::cout << "   WARNING: The following species (R') is not available in the ABSTRACTORS section: " << species_Rp_ << std::endl;
			std::cout << "            No reversibility correction will be applied to reactions: " << label << std::endl;
		}
		else
		{
			index_Rprime = std::distance(abstractors.begin(), iRprime);
		}

		for (unsigned int i = 0; i < n_; i++)
		{
			std::vector<std::string>::const_iterator iR = std::find(abstractors.begin(), abstractors.end(), species_R_[i]);
			if (iR == abstractors.end())
			{
				std::cout << "The following species (R) is not available in the ABSTRACTORS section: " << species_R_[i] << std::endl;
				exit(-1);
			}
			const unsigned int index_R = std::distance(abstractors.begin(), iR);
			
			A_[i] = abstractors_A[index_R];
			Beta_[i] = abstractors_Beta[index_R];
			E_[i] = abstractors_E[index_R];

			// Correction
			if (type_H_ != 0)
			{
				A_[i] *= n_H_ * corrections_A[type_H_ - 1] * std::pow(Tref, -Beta_[i])*std::exp(-Beta_[i]);
				E_[i] += corrections_E[type_H_ - 1] * std::pow(abstractors_E[index_R] / Eref, alpha);
				E_[i] -= Beta_[i] * Tref*PhysicalConstants::R_cal_mol;

				// Reversibility correction
				if (index_Rprime != -1)
					E_[i] -= (1. - std::pow(abstractors_E[index_Rprime] / Eref, alpha))*corrections_Er[i];
			}
			else
			{
				A_[i] *= n_H_ * std::pow(Tref, -Beta_[i])*std::exp(-Beta_[i]);
				E_[i] -= Beta_[i] * Tref*PhysicalConstants::R_cal_mol;
			}
		}
	}

}
