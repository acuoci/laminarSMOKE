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

#include <boost/algorithm/string.hpp>

namespace OpenSMOKE
{
	AbstractionReactions::AbstractionReactions(InputFileCHEMKIN& kinetics)
	{
		is_active_ = false;
		Tref_ = 1000.;
		alpha_ = 1. / 3.;
		Eref_ = 13500.;

		first_line_of_abstractors_  = 0;
		first_line_of_corrections_  = 0;
		first_line_of_abstractions_ = 0;

		for (unsigned int j = 0; j < kinetics.good_lines().size(); j++)
		{
			std::string line = kinetics.good_lines()[j];

			boost::replace_all(line, "abstractors", "ABSTRACTORS    ");
			boost::replace_all(line, "corrections", "CORRECTIONS    ");
			boost::replace_all(line, "abstractions", "ABSTRACTIONS    ");
			boost::replace_all(line, "end", "END    ");

			{
				size_t found = line.find("ABSTRACTORS");
				if (found != std::string::npos)	first_line_of_abstractors_ = j;
			}
			{
				size_t found = line.find("CORRECTIONS");
				if (found != std::string::npos)	first_line_of_corrections_ = j;
			}
			{
				size_t found = line.find("ABSTRACTIONS");
				if (found != std::string::npos)	first_line_of_abstractions_ = j;
			}
			{
				size_t found = line.find("END");
				if (found != std::string::npos)	abstractions_end_lines_.push_back(j);
			}
		}

		if (first_line_of_abstractors_ != 0 && first_line_of_corrections_ == 0)
		{
			std::cout << "Reading kinetic scheme. No CORRECTIONS section is present." << std::endl;
			std::cout << "CORRECTIONS section is strictly needed if ABSTRACTORS section is present" << std::endl;
			OpenSMOKE::FatalErrorMessage("Press enter to exit...");
		}
		if (first_line_of_abstractors_ == 0 && first_line_of_corrections_ != 0)
		{
			std::cout << "Reading kinetic scheme. No ABSTRACTORS section is present." << std::endl;
			std::cout << "ABSTRACTORS section is strictly needed if CORRECTIONS section is present" << std::endl;
			OpenSMOKE::FatalErrorMessage("Press enter to exit...");
		}
		if (first_line_of_abstractions_ != 0 && (first_line_of_abstractors_ == 0 && first_line_of_corrections_ == 0))
		{
			std::cout << "Reading kinetic scheme. No ABSTRACTORS and/or CORRECTIONS section is present." << std::endl;
			std::cout << "ABSTRACTORS and CORRECTIONS sections are strictly needed if ABSTRACTIONS section is present" << std::endl;
			OpenSMOKE::FatalErrorMessage("Press enter to exit...");
		}

		if (first_line_of_abstractors_*first_line_of_corrections_ > 0)
			is_active_ = true;
	}

	bool AbstractionReactions::Checking(InputFileCHEMKIN& kinetics)
	{
		unsigned int last_line_of_abstractors = kinetics.good_lines().size();
		unsigned int last_line_of_corrections = kinetics.good_lines().size();
		unsigned int last_line_of_abstractions = kinetics.good_lines().size();
		for (unsigned int i = 0; i < abstractions_end_lines_.size(); i++)
		{
			if (abstractions_end_lines_[i] > first_line_of_abstractors_)
			{
				if (abstractions_end_lines_[i] < last_line_of_abstractors)
					last_line_of_abstractors = abstractions_end_lines_[i];
			}

			if (abstractions_end_lines_[i] > first_line_of_corrections_)
			{
				if (abstractions_end_lines_[i] < last_line_of_corrections)
					last_line_of_corrections = abstractions_end_lines_[i];
			}

			if (abstractions_end_lines_[i] > first_line_of_abstractions_)
			{
				if (abstractions_end_lines_[i] < last_line_of_abstractions)
					last_line_of_abstractions = abstractions_end_lines_[i];
			}
		}

		for (unsigned int j = first_line_of_abstractors_ + 1; j < last_line_of_abstractors; j++)
			lines_abstractors_.push_back(kinetics.good_lines()[j]);
		for (unsigned int j = first_line_of_corrections_ + 1; j < last_line_of_corrections; j++)
			lines_corrections_.push_back(kinetics.good_lines()[j]);
		for (unsigned int j = first_line_of_abstractions_ + 1; j < last_line_of_abstractions; j++)
			lines_abstractions_.push_back(kinetics.good_lines()[j]);

		std::vector<unsigned int> lines_to_remove;
		for (unsigned int j = first_line_of_abstractors_; j <= last_line_of_abstractors; j++)
			lines_to_remove.push_back(j);
		for (unsigned int j = first_line_of_corrections_; j <= last_line_of_corrections; j++)
			lines_to_remove.push_back(j);
		for (unsigned int j = first_line_of_abstractions_; j <= last_line_of_abstractions; j++)
			lines_to_remove.push_back(j);

		kinetics.ConvertGoodLinesIntoBlankLines(lines_to_remove);

		return true;
	}

	bool AbstractionReactions::ExtractTables()
	{
		ExtractTableAbstractors();
		ExtractTableCorrections();
		ExtractTableAbstractions();

		return true;
	}

	bool AbstractionReactions::ExtractTableAbstractors()
	{
		for (unsigned int j = 0; j < lines_abstractors_.size(); j++)
		{
			typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
			boost::char_separator<char> sep(" ");
			tokenizer tokens(lines_abstractors_[j], sep);

			unsigned int count = 0;
			for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter)
				count++;

			if (count != 5)
			{
				std::cout << "Reading ABSTRACTORS section. Error in the following line:" << std::endl;
				std::cout << lines_abstractors_[j] << std::endl;
				return false;
			}

			count = 0;
			for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter)
			{
				if (count == 0)			abstractors_.push_back(*tok_iter);
				else if (count == 1)	abstracted_.push_back(*tok_iter);
				else if (count == 2)	A_.push_back(boost::lexical_cast<double>(*tok_iter));
				else if (count == 3)	Beta_.push_back(boost::lexical_cast<double>(*tok_iter));
				else if (count == 4)	E_.push_back(boost::lexical_cast<double>(*tok_iter));
				count++;
			}
		}

		return true;
	}

	bool AbstractionReactions::ExtractTableCorrections()
	{
		for (unsigned int j = 0; j < lines_corrections_.size(); j++)
		{
			typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
			boost::char_separator<char> sep(" ");
			tokenizer tokens(lines_corrections_[j], sep);

			unsigned int count = 0;
			for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter)
				count++;

			if (count != 3)
			{
				std::cout << "Reading CORRECTIONS section. Error in the following line:" << std::endl;
				std::cout << lines_corrections_[j] << std::endl;
				return false;
			}

			count = 0;
			for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter)
			{
				if (count == 0)			corrections_types_.push_back(*tok_iter);
				else if (count == 1)	corrections_A_.push_back(boost::lexical_cast<double>(*tok_iter));
				else if (count == 2)	corrections_E_.push_back(boost::lexical_cast<double>(*tok_iter));
				count++;
			}
		}

		return true;
	}


	bool AbstractionReactions::ExtractTableAbstractions()
	{
		for (unsigned int j = 0; j < lines_abstractions_.size(); j++)
		{
			typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
			boost::char_separator<char> sep(" ");
			tokenizer tokens(lines_abstractions_[j], sep);

			unsigned int count = 0;
			for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter)
				count++;

			if (count != 2)
			{
				std::cout << "Reading ABSTRACTIONS section. Error in the following line:" << std::endl;
				std::cout << lines_abstractions_[j] << std::endl;
				return false;
			}

			std::string element1 = "";
			std::string element2 = "";
			count = 0;
			for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter)
			{
				if (count == 0)			element1 = *tok_iter;
				else if (count == 1)	element2 = *tok_iter;
				count++;
			}

			if (element1 == "TREF")
			{
				Tref_ = boost::lexical_cast<double>(element2);
			}
			else if (element1 == "ALPHA")
			{
				alpha_ = boost::lexical_cast<double>(element2);
			}
			else if (element1 == "EREF")
			{
				Eref_ = boost::lexical_cast<double>(element2);
			}
			else
			{
				std::cout << "Reading ABSTRACTIONS section. Error in the following line:" << std::endl;
				std::cout << lines_abstractions_[j] << std::endl;
				return false;
			}

		}

		return true;
	}

	void AbstractionReactions::Summary(std::ostream& out)
	{
		out << std::endl;
		out << "---------------------------------------------------------------------------------" << std::endl;
		out << "                                   ABSTRACTORS                                   " << std::endl;
		out << "---------------------------------------------------------------------------------" << std::endl;
		out << " Abstractor          Abstracted          A[mol,cm,s]     Beta       E[cal/mol]   " << std::endl;
		out << "---------------------------------------------------------------------------------" << std::endl;
		for (unsigned int j = 0; j < A_.size(); j++)
		{
			out << " ";
			out << std::setw(20) << std::left << abstractors_[j];
			out << std::setw(20) << std::left << abstracted_[j];
			out << std::setw(16) << std::left << std::scientific << std::setprecision(3) << A_[j];
			out << std::setw(11) << std::left << std::fixed << std::setprecision(3) << Beta_[j];
			out << std::setw(16) << std::left << std::fixed << std::setprecision(3) << E_[j];
			out << std::endl;
		}
		out << "---------------------------------------------------------------------------------" << std::endl;

		out << std::endl;
		out << "---------------------------------------------------------------------------------" << std::endl;
		out << "                                   CORRECTIONS                                   " << std::endl;
		out << "---------------------------------------------------------------------------------" << std::endl;
		out << " Index  Type                A[-]       E[cal/mol]                                " << std::endl;
		out << "---------------------------------------------------------------------------------" << std::endl;
		for (unsigned int j = 0; j < corrections_A_.size(); j++)
		{
			out << " ";
			out << std::setw(7)  << std::left << j+1;
			out << std::setw(20) << std::left << corrections_types_[j];
			out << std::setw(11) << std::left << std::fixed << std::setprecision(3) << corrections_A_[j];
			out << std::setw(16) << std::left << std::fixed << std::setprecision(3) << corrections_E_[j];
			out << std::endl;
		}
		out << "---------------------------------------------------------------------------------" << std::endl;


	}
	
	int AbstractionReactions::ParseReaction(ReactionPolicy_CHEMKIN& reaction, const std::vector<std::string>& map_of_species)
	{
		std::vector<std::string> reactant_species;
		for (unsigned int i = 0; i<reaction.reactant_nu().size(); i++)
			reactant_species.push_back(map_of_species[reaction.reactant_nu_indices()[i]]);

		std::vector<std::string> product_species;
		for (unsigned int i = 0; i<reaction.product_nu().size(); i++)
			product_species.push_back(map_of_species[reaction.product_nu_indices()[i]]);

		if (std::find(reactant_species.begin(), reactant_species.end(), "R") != reactant_species.end() &&
			std::find(product_species.begin(), product_species.end(), "RH") != product_species.end())
		{
			if (reactant_species.size() == 2 && product_species.size() == 2)
			{
				std::string RpH = "";
				for (unsigned int i = 0; i < reaction.reactant_nu().size(); i++)
					if (reactant_species[i] != "R")
						RpH=reactant_species[i];

				std::string Rp = "";
				for (unsigned int i = 0; i < reaction.product_nu().size(); i++)
					if (product_species[i] != "RH")
						Rp = product_species[i];

				const double n_H = reaction.A()*1e3;
				const int type_H = std::round(reaction.E_over_R()*PhysicalConstants::R_cal_mol);

				exploded_.push_back(AbstractionExploded(abstractors_, RpH, Rp, abstracted_, n_H, type_H));
				
				return 1;
			}
		}

		if (std::find(reactant_species.begin(), reactant_species.end(), "R") == reactant_species.end() &&
			std::find(product_species.begin(), product_species.end(), "RH") != product_species.end())
		{
			std::cout << "Found a reaction with RH species on the product side, without R species on the reactiant side." << std::endl;
			std::cout << "The R and RH names for species are reserved to the ABSTRACTIONS formalism." << std::endl;
			return -1;
		}

		if (std::find(reactant_species.begin(), reactant_species.end(), "R") != reactant_species.end() &&
			std::find(product_species.begin(), product_species.end(), "RH") == product_species.end())
		{
			std::cout << "Found a reaction with R species on the reactant side, without RH species on the product side." << std::endl;
			std::cout << "The R and RH names for species are reserved to the ABSTRACTIONS formalism." << std::endl;
			return -1;
		}

		return 0;
	}

	unsigned int AbstractionReactions::RemoveExistingReactions(ReactionPolicy_CHEMKIN& reaction, const std::vector<std::string>& map_of_species)
	{
		unsigned int n_existing = 0;

		if (reaction.reactant_nu().size() == 2 && reaction.product_nu().size() == 2)
		{
			std::vector<std::string> reactant_species;
			for (unsigned int i = 0; i<reaction.reactant_nu().size(); i++)
				reactant_species.push_back(map_of_species[reaction.reactant_nu_indices()[i]]);

			std::vector<std::string> product_species;
			for (unsigned int i = 0; i<reaction.product_nu().size(); i++)
				product_species.push_back(map_of_species[reaction.product_nu_indices()[i]]);

			for (unsigned int j = 0; j < exploded_.size(); j++)
				n_existing += exploded_[j].RemoveExistingReactions(reactant_species, product_species);
		}

		return n_existing;
	}

	std::string my_trim(const std::string& str, const std::string& whitespace = " \t")
	{
		const auto strBegin = str.find_first_not_of(whitespace);
		if (strBegin == std::string::npos)
			return ""; // no content

		const auto strEnd = str.find_last_not_of(whitespace);
		const auto strRange = strEnd - strBegin + 1;

		return str.substr(strBegin, strRange);
	}

	std::string my_reduce(std::string& str, const std::string& fill = " ", const std::string& whitespace = " \t")
	{
		// trim first
		auto result = my_trim(str, whitespace);
		// replace sub ranges
		auto beginSpace = result.find_first_of(whitespace);
		while (beginSpace != std::string::npos)
		{
			const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
			const auto range = endSpace - beginSpace;

			result.replace(beginSpace, range, fill);

			const auto newStart = beginSpace + fill.length();
			beginSpace = result.find_first_of(whitespace, newStart);
		}

		return result;
	}

	int AbstractionReactions::ReplaceAbstractionReaction(std::string& line)
	{
		if (line.find("HNUMBER") != std::string::npos && line.find("HTYPE") == std::string::npos)
		{
			std::cout << "The HTYPE parameter is miising from the abstraction reaction" << std::endl;
			std::cout << "Please use the following format for abstraction reactions: R+R'H => R'+RH    HNUMBER=X  HTYPE=X" << std::endl;
			return -1;
		}

		if (line.find("HNUMBER") == std::string::npos && line.find("HTYPE") != std::string::npos)
		{
			std::cout << "The HNUMBER parameter is miising from the abstraction reaction" << std::endl;
			std::cout << "Please use the following format for abstraction reactions: R+R'H => R'+RH    HNUMBER=X  HTYPE=X" << std::endl;
			return -1;
		}

		if (line.find("HNUMBER") != std::string::npos && line.find("HTYPE") != std::string::npos)
		{
			const std::string tmp = my_reduce(line);
			line = tmp;

			// Check the order
			if (line.find("HNUMBER") > line.find("HTYPE"))
			{
				std::cout << "The HNUMBER parameter must preceed the HTYPE option" << std::endl;
				std::cout << "Please use the following format for abstraction reactions: R+R'H => R'+RH    HNUMBER=X  HTYPE=X" << std::endl;
				return -1;
			}

			// Replace
			boost::replace_all(line, "HNUMBER=", "");
			boost::replace_all(line, "HTYPE=", " 0 ");

			// Replace
			for (unsigned int i = 0; i < corrections_types_.size(); i++)
			{
				std::stringstream index; index << i + 1;
				boost::replace_all(line, corrections_types_[i], index.str());
			}

			// Separate the line
			std::vector<std::string> tokens;
			boost::algorithm::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);

			std::vector<std::string> elements;
			for (auto token : tokens)
				elements.push_back(token);

			try
			{
				boost::lexical_cast<int>(elements[elements.size()-1]);		// must be int
				boost::lexical_cast<int>(elements[elements.size()-2]);		// must be 0
				boost::lexical_cast<double>(elements[elements.size()-3]);	// must be double	
			}
			catch (...) 
			{
				std::cout << "Error in parameters of abstraction reaction." << std::endl;
				std::cout << "Please use the following format for abstraction reactions: R+R'H => R'+RH    HNUMBER=X  HTYPE=X" << std::endl;
				return -1;
			}

			const double nh = boost::lexical_cast<double>(elements[elements.size() - 3]);
			if (nh<=0.)
			{
				std::cout << "The number of hydrogen atoms in abstraction reactions must be strictly positive." << std::endl;
				std::cout << "Please check the HNUMBER option." << std::endl;
				return -1;
			}

			const int itype = boost::lexical_cast<int>(elements[elements.size() - 1]);
			if (itype < 0 || itype > corrections_types_.size())
			{
				std::cout << "The type of hydrogen atoms must be specified between 0 and " << corrections_types_.size() << std::endl;
				std::cout << "Please check the HTYPE option." << std::endl;
				return -1;
			}

			return 1;
		}

		return 0;
	}

	void AbstractionReactions::ExplodeReactions()
	{
		corrections_Er_.resize(abstractors_.size());
		for (unsigned int i = 0; i < abstractors_.size(); i++)
			corrections_Er_[i] = 0;

		for (unsigned int i = 0; i < abstractors_.size(); i++)
		{
			bool iFound = false;
			for (unsigned int j = 0; j < exploded_.size(); j++)
				if (exploded_[j].species_RpH() == abstracted_[i])
				{
					if (exploded_[j].type_H() != 0)
						corrections_Er_[i] = corrections_E_[exploded_[j].type_H() - 1];

					iFound = true;
					break;
				}

			if (iFound == false)
			{
				const std::string label_target = abstractors_[i] + " + R'H" + " => " + "R' + " + abstracted_[i];
				const std::string label_missing = "R + " + abstracted_[i] + " => " + abstractors_[i] + " + RH";
				std::cout << "   WARNING: No abstraction reaction with the following form was declared: " << label_missing << std::endl;
				std::cout << "            No reversibility correction will be applied to reactions: " << label_target << std::endl;
			}
		}

		for (unsigned int i = 0; i < exploded_.size(); i++)
			exploded_[i].ExplodeReactions(abstractors_, A_, Beta_, E_, corrections_A_, corrections_E_, corrections_Er_, Eref_, alpha_, Tref_);
	}

	bool AbstractionReactions::CheckListOfSpecies(const std::vector<std::string> names_species)
	{
		for (unsigned int j = 0; j < A_.size(); j++)
		{
			if (std::count(names_species.begin(), names_species.end(), abstractors_[j]) == 0)
			{
				std::cout << "The following abstractor species is not included in the main list of species: " << abstractors_[j] << std::endl;
				std::cout << "Please correct the ABSTRACTORS section or the SPECIES section" << std::endl;
				return false;
			}
			if (std::count(names_species.begin(), names_species.end(), abstracted_[j]) == 0)
			{
				std::cout << "The following abstracted species is not included in the main list of species: " << abstracted_[j] << std::endl;
				std::cout << "Please correct the ABSTRACTORS section or the SPECIES section" << std::endl;
				return false;
			}
		}

		return true;
	}
}
