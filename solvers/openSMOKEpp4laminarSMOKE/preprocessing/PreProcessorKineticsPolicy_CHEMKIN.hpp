/*-----------------------------------------------------------------------*\
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

#include "math/OpenSMOKEStdInclude.h"
#include "math/PhysicalConstants.h"
#include "kernel/thermo/AtomicCompositionTable.h"
#include <boost/algorithm/string.hpp>
#include <iterator>

namespace OpenSMOKE
{

	template<typename Reactions>
	PreProcessorKineticsPolicy_CHEMKIN<Reactions>::PreProcessorKineticsPolicy_CHEMKIN() {
	}

	template<typename Reactions>
	PreProcessorKineticsPolicy_CHEMKIN<Reactions>::PreProcessorKineticsPolicy_CHEMKIN(const PreProcessorKineticsPolicy_CHEMKIN& orig) {
	}

	template<typename Reactions>
	PreProcessorKineticsPolicy_CHEMKIN<Reactions>::~PreProcessorKineticsPolicy_CHEMKIN() 
	{
		delete myKinetics;
		reaction_lines.clear();	
		iReactionLines.clear();

		reactions_.clear();
		names_species_.clear();	
		map_of_species.clear();	
	}

	template<typename Reactions>
	bool PreProcessorKineticsPolicy_CHEMKIN<Reactions>::ReadFromASCIIFile(const std::string file_name)
	{
		std::cout << " * Reading kinetic file... " << std::endl;
		myKinetics = new InputFileCHEMKIN(file_name);
		//myKinetics->Status(std::cout);

		{
			// Abstractions (CRECK Modeling standard, not available in original CHEMKIN)
			std::cout << " * Check if Abstraction Reaction module is available... " << std::endl;
			abstractions_ = new AbstractionReactions(*myKinetics);
			if (abstractions_->is_active() == true)
			{
				std::cout << " * Found Abstraction Reaction Module... " << std::endl;
				abstractions_->Checking(*myKinetics);
				abstractions_->ExtractTables();
				abstractions_->Summary(std::cout);
			}

			for (unsigned int j = 0; j<myKinetics->good_lines().size(); j++)
			{
				size_t found = myKinetics->good_lines()[j].find("REACTIONS");
				if (found != std::string::npos)
				{
					iReactionLines.push_back(j);
					break;
				}
				found = myKinetics->good_lines()[j].find("reactions");
				if (found != std::string::npos)
				{
					iReactionLines.push_back(j);
					break;
				}
			}

			if (iReactionLines.size() == 0)
			{
				std::cout << "Reading kinetic scheme. No REACTIONS keyword found!" << std::endl;
				return false;
			}

			std::vector<unsigned int> iElementLines;
			std::vector<unsigned int> iSpeciesLines;
			std::vector<unsigned int> iEndLines;
			std::vector<unsigned int> iThermoLines;
			for (unsigned int j=0;j<iReactionLines[0];j++)
			{
				std::string line = myKinetics->good_lines()[j];
			
				boost::replace_all(line, "ELEMENTS", "ELEM    ");
				boost::replace_all(line, "elements", "ELEM    ");

				boost::replace_all(line, "SPECIES", "SPEC");
				boost::replace_all(line, "species", "SPEC");

				boost::replace_all(line, "end", "END");

				boost::replace_all(line, "thermo", "THERMO    ");

				{
					size_t found=line.find("ELEM");
					if (found!=std::string::npos)	iElementLines.push_back(j);
				}
				{
					size_t found=line.find("SPEC");
					if (found!=std::string::npos)	iSpeciesLines.push_back(j);
				}
				{
					size_t found=line.find("END");
					if (found!=std::string::npos)	iEndLines.push_back(j);
				}
				{
					size_t found=line.find("THERMO");
					if (found!=std::string::npos)	iThermoLines.push_back(j);
				}
			}

			if (iElementLines.size() == 0)
			{
				std::cout << "Reading kinetic scheme. No ELEMENTS section found!" << std::endl;
				return false;
			}

			if (iSpeciesLines.size() == 0)
			{
				std::cout << "Reading kinetic scheme. No SPECIES section found!" << std::endl;
				return false;
			}

			if (iThermoLines.size() > 0)
			{
				if ( (iThermoLines.size() > 1) || (iThermoLines[0] < iSpeciesLines[0]) )
				{
					std::cout << "Reading kinetic scheme. The kinetic interpreter does not support the possibility " << std::endl;
					std::cout << "to write thermodynamic data in the kinetic file." << std::endl;
					std::cout << "Please move the thermodinamyc section in a different file." << std::endl;
					return false;
				}

				unsigned int last_line_thermo = 0;
				for (unsigned int i=0;i<iEndLines.size();i++)
					if (iEndLines[i] >= iThermoLines[0])
					{
						last_line_thermo = iEndLines[i];
						break;
					}

				if (last_line_thermo == 0)
				{
					std::cout << "Reading kinetic scheme. The kinetic interpreter does not support the possibility " << std::endl;
					std::cout << "to write thermodynamic data in the kinetic file." << std::endl;
					std::cout << "Please move the thermodinamyc section in a different file." << std::endl;
					return false;
				}
				else
				{
					for(unsigned int i=iThermoLines[0]+1;i<=last_line_thermo-1;i++)
					{
						for (unsigned int j=0;j<=myKinetics->good_lines()[i].size();j++)
							if(myKinetics->good_lines()[i][j] != ' ')
							{
								std::cout << "Reading kinetic scheme. The kinetic interpreter does not support the possibility " << std::endl;
								std::cout << "to write thermodynamic data in the kinetic file." << std::endl;
								std::cout << "Please move the thermodinamyc section in a different file." << std::endl;
								return false;
							}
					}
				}
			}

			{
				if (iElementLines.back() >= iSpeciesLines.front())
				{
					std::cout << "Reading kinetic scheme. Please specifiy all the elements before the species!" << std::endl;
					return false;
				}
			
				if (iSpeciesLines.back() >= iReactionLines.front())
				{
					std::cout << "Reading kinetic scheme. Please specifiy all the species before the reactions!" << std::endl;
					return false;
				}

				// Element section 
				{
					// Limits
					unsigned int last_line_of_elements = 0;
					for (unsigned int i=0;i<iEndLines.size();i++)
						if (iEndLines[i] < iSpeciesLines.front())
							last_line_of_elements = iEndLines[i];
					if (last_line_of_elements < iElementLines.back())
						last_line_of_elements = iElementLines.back();

					// Element analysis
					std::string line_elements;
					for (unsigned int i=0;i<=last_line_of_elements;i++)
						line_elements += myKinetics->good_lines()[i] + " ";
					boost::replace_all(line_elements, "ELEMENTS", " ");
					boost::replace_all(line_elements, "elements", " ");
					boost::replace_all(line_elements, "ELEM", " ");
					boost::replace_all(line_elements, "elem", " ");
					boost::replace_all(line_elements, "END", " ");
					boost::replace_all(line_elements, "end", " ");

					boost::replace_all(line_elements, "/", " / ");


					std::vector<std::string> name_provisional;
					typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
					boost::char_separator<char> sep_blank(" ");
					tokenizer_blank tokens(line_elements, sep_blank);
					for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
						name_provisional.push_back(*tok_iter);
			
					std::vector<std::string>	name_elements;
					std::vector<double>			weights;
					unsigned int i=0;
					for (;;)
					{
						if (i >= name_provisional.size())
							break;

						if (name_provisional[i] != "/")
						{
							name_elements.push_back(name_provisional[i]);
							weights.push_back(0.);
							i++;
						}
						else
						{
							if ( (i+2+1)<=name_provisional.size() )
							{
								if (name_provisional[i+2] == "/")
								{
									try
									{
										weights[weights.size()-1] = boost::lexical_cast<double>(name_provisional[i+1]);
									}
									catch(boost::bad_lexical_cast &)
									{
										std::cout << "Reading kinetic scheme. The isotopic weight is not written properly!" << std::endl;
										return false;
									}
								}
								else
								{
									std::cout << "Reading kinetic scheme. The isotopic element is not declared properly!" << std::endl;
									return false;
								}
								i+=3;
							}
							else
							{
								std::cout << "Reading kinetic scheme. The isotopic element is not declared properly!" << std::endl;
								return false;
							}
						}
					}

					// Checking
				}
			}

			// Species section 
			{
				// Limits
				unsigned int last_line_of_species = iSpeciesLines.front();
				for (unsigned int i=0;i<iEndLines.size();i++)
				{
					if (iEndLines[i] < iReactionLines.front())
						last_line_of_species = iEndLines[i];
				}
				for (unsigned int i=0;i<iEndLines.size();i++)
				{
					if (iThermoLines.size() == 1)
						if (iEndLines[i] < iThermoLines[0])
							last_line_of_species = iEndLines[i];
				}
				if (last_line_of_species < iSpeciesLines.back())
					last_line_of_species = iSpeciesLines.back();

				// Element analysis
				std::string line_species;
				for (unsigned int i=iSpeciesLines.front();i<=last_line_of_species;i++)
					line_species += myKinetics->good_lines()[i] + " ";
				boost::trim(line_species);

				boost::replace_all(line_species, "SPECIES", " ");
				boost::replace_all(line_species, "species", " ");
				boost::replace_all(line_species, "SPEC", " ");
				boost::replace_all(line_species, "spec", " ");
				boost::replace_all(line_species, "END", " ");
				boost::replace_all(line_species, "end", " ");

				
				typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
				boost::char_separator<char> sep_blank(" ");
				tokenizer_blank tokens(line_species, sep_blank);
				for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
					names_species_.push_back(*tok_iter);

				// Abstractions
				if (abstractions_->is_active() == true)
				{
					names_species_.push_back("R");
					names_species_.push_back("RH");

					// Check the species
					if ( abstractions_->CheckListOfSpecies(names_species_) == false)
						return false;
				}
				
				// Checking if a species appears more tha once
				{
					for (unsigned int i = 0; i < names_species_.size(); i++)
						if (std::count(names_species_.begin(), names_species_.end(), names_species_[i]) > 1)
						{
							std::cout << "Reading kinetic scheme. The following species is declared more than once: " << names_species_[i] << std::endl;
							return false;
						}
				}

				// Create the map of species
				for(unsigned int i=0;i<names_species_.size();i++)
					map_of_species.insert(std::make_pair(names_species_[i],  i));
			}
		}

		return true;
	}

	template<typename Reactions>
	bool PreProcessorKineticsPolicy_CHEMKIN<Reactions>::KineticsFromASCIIFile(AtomicCompositionTable& atomicComposition, std::ostream& flog)
	{

		// Reaction units
		PhysicalConstants::UNITS_REACTION e_units = PhysicalConstants::UNITS_CAL_MOLE;
		PhysicalConstants::UNITS_REACTION a_units = PhysicalConstants::UNITS_MOLES;
		{
			typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
			boost::char_separator<char> sep_blank(" ");
			tokenizer_blank tokens(myKinetics->good_lines()[iReactionLines[0]], sep_blank);
			const std::size_t n = std::distance (tokens.begin(), tokens.end());
			if (n>3)
			{
				std::cout << "Reading kinetic scheme. Wrong number of arguments on the REACTION line (line " << myKinetics->indices_of_good_lines()[iReactionLines[0]] << ")!" << std::endl;
				return false;
			}

			//for (tokenizer_blank::iterator tok_iter = std::next(tokens.begin()); tok_iter != tokens.end(); ++tok_iter)
			for (tokenizer_blank::iterator tok_iter = boost::next(tokens.begin()); tok_iter != tokens.end(); ++tok_iter)
			{
				std::string tag = *tok_iter;
				
				if (tag == "CAL" || tag == "CAL/MOLE" || tag == "cal" || tag == "cal/mole")
					e_units = PhysicalConstants::UNITS_CAL_MOLE;
				else if (tag == "EVOL" || tag == "EVOLTS" || tag == "evol" || tag == "evolts")
					e_units = PhysicalConstants::UNITS_EVOLTS;
				else if (tag == "JOUL" || tag == "JOULES/MOLE" || tag == "joul" || tag == "joules/mole")
					e_units = PhysicalConstants::UNITS_JOULES_MOLE;
				else if (tag == "KCAL" || tag == "KCAL/MOLE" || tag == "kcal" || tag == "kcal/mole")
					e_units = PhysicalConstants::UNITS_KCAL_MOLE;
				else if (tag == "KJOU" || tag == "KJOULES/MOLE" || tag == "kjou" || tag == "kjoules/mole")
					e_units = PhysicalConstants::UNITS_KJOULES_MOLE;
				else if (tag == "KELV" || tag == "KELVINS" || tag == "kelv" || tag == "kelvins")
					e_units = PhysicalConstants::UNITS_KELVINS;
				
				else if (tag == "MOLEC" || tag == "MOLECULES" || tag == "molec" || tag == "molecules")
					a_units = PhysicalConstants::UNITS_MOLECULES;
				else if (tag == "MOLE" || tag == "MOLES" || tag == "mole" || tag == "moles")
					a_units = PhysicalConstants::UNITS_MOLES;
				
				else
				{
					std::cout << "Reading kinetic scheme. Wrong units on the REACTION line (line " << myKinetics->indices_of_good_lines()[iReactionLines[0]] << ")!" << std::endl;
					std::cout << "Available options: CAL || CAL/MOLE || EVOL || EVOLTS || JOUL || JOULES/MOLE || KCAL || KCAL/MOLE || KELV || KELVINS || MOLEC || MOLECULES || MOLE || MOLES" << std::endl;
					return false;
				}
			}
		}

		unsigned int iEndReactionLine = 0;
		for (unsigned int j=iReactionLines[0];j<myKinetics->good_lines().size();j++)
		{
			size_t found_END=myKinetics->good_lines()[j].find("END");
			if (found_END!=std::string::npos)
			{
				iEndReactionLine = j+1;
				break;
			}
			else
			{
				size_t found_end=myKinetics->good_lines()[j].find("end");
				if (found_end!=std::string::npos)
				{
					iEndReactionLine = j+1;
					break;
				}
			}
		}
		if (iEndReactionLine == 0)
		{
			std::cout << "Reading kinetic scheme. Missing END (end) keyword at the end of the list of the reactions!" << std::endl;
			return false;
		}

		unsigned int number_of_reactions = 0;
		for (unsigned int j=0;j<iEndReactionLine;j++)
		{
			size_t found=myKinetics->good_lines()[j].find("=");
			if (found!=std::string::npos)
			{
				reaction_lines.push_back(j);
				number_of_reactions++;
			}
		}

		reaction_lines.push_back(iEndReactionLine-1);

		reactions_.resize(number_of_reactions);
		if (number_of_reactions > 20)
			std::cout << " * Parsing " << number_of_reactions << " reactions: ";

		// Abstractions
		if (abstractions_->is_active() == true)
		{
			for (unsigned int j = 0; j < number_of_reactions; j++)
			{
				std::string line = myKinetics->good_lines()[reaction_lines[j]];
				const int success = abstractions_->ReplaceAbstractionReaction(line);

				if (success == -1)
				{
					const unsigned int index_line = myKinetics->indices_of_good_lines()[reaction_lines[j]];
					std::cout << "Fatal error: please check the abstraction reaction at line: " << index_line << std::endl;
					return false;
				}

				myKinetics->ReplaceGoodLine(reaction_lines[j], line);
			}
		}

		unsigned int large_error_in_stoichiometries = 0;
		unsigned int small_error_in_stoichiometries = 0;
		for (unsigned int j=0;j<number_of_reactions;j++)
		{
			
			if (number_of_reactions > 20)
			{
				if (j%(number_of_reactions/20) == 0) 
					std::cout << "%";
				if (j==number_of_reactions-1)
					std::cout << std::endl;
			}

			std::vector<std::string> list_of_lines;
			for (unsigned int i=reaction_lines[j];i<reaction_lines[j+1];i++)
				list_of_lines.push_back(myKinetics->good_lines()[i]);
			try
			{
				reactions_[j].SetDefaultUnits();
				reactions_[j].SetUnits(a_units, e_units);
				
				bool successReading = reactions_[j].ReadReactionFromStrings(list_of_lines, map_of_species);
				if (successReading == false)
					throw reaction_lines[j];
				
				unsigned int successStoichiometry = atomicComposition.CheckStoichiometry(flog, reactions_[j], 1e-3);
				if (successStoichiometry > 0)
				{
					if (successStoichiometry == 1)
					{
						large_error_in_stoichiometries++;
						std::string reaction_string; reactions_[j].GetReactionString(names_species(),reaction_string);
						boost::erase_all(reaction_string, " ");
						flog << "Error in reaction (line " << reaction_lines[j]+1 << "): " << reaction_string << std::endl;
					}
					if (successStoichiometry == 2)
					{
						small_error_in_stoichiometries++;
						
						// Write uncorrected reaction
						{
							std::string reaction_string; reactions_[j].GetReactionString(names_species(), reaction_string);
							boost::erase_all(reaction_string, " ");
							flog << "Warning in reaction (line " << reaction_lines[j] + 1 << "): " << reaction_string << std::endl;
						}

						// Correct the stoichiometry and write the corrected reaction
						{
							unsigned int flag = atomicComposition.CorrectStoichiometry(reactions_[j]);
							std::string reaction_string; reactions_[j].GetReactionString(names_species(), reaction_string);
							boost::erase_all(reaction_string, " ");
							if (flag == 1)	flog << "Corrected reaction (line " << reaction_lines[j] + 1 << "): " << reaction_string << std::endl;
							else            flog << "Correction failed for reaction (line " << reaction_lines[j] + 1 << "): " << reaction_string << std::endl;
						}
					}					
					flog << std::endl;
				}
			}
			catch(unsigned int k)
			{
				std::cout << "Reading kinetic scheme: error in reaction starting at line " << myKinetics->indices_of_good_lines()[k] << std::endl;
				return false;
			}
		}

		if (abstractions_->is_active() == true)
		{
			// Search for abstraction reactions
			std::vector<unsigned int> list_lines_to_be_exploded;
			{
				std::cout << " * Search for abstraction reactions..." << std::endl;
				
				unsigned int n_abstraction_reactions = 0;
				for (unsigned int j = 0; j < number_of_reactions; j++)
				{
					const int is_abstraction = abstractions_->ParseReaction(reactions_[j], names_species_);

					if (is_abstraction == 1)	// abstraction reaction
					{
						n_abstraction_reactions++;
						std::vector<unsigned int> list_of_lines;
						for (unsigned int i = reaction_lines[j]; i < reaction_lines[j + 1]; i++)
							list_of_lines.push_back(i);

						if (list_of_lines.size() != 1)
						{
							std::cout << "The abstraction reactions must be written in a single line. No multiple lines are allowed." << std::endl;
							std::cout << "Please check the abstraction reaction starting at line: " << reaction_lines[j] << std::endl;
							return false;
						}
						else
						{
							list_lines_to_be_exploded.push_back(list_of_lines[0]);
						}
					}
					else if (is_abstraction == -1) // conventional reaction
					{
						std::cout << "Please check the abstraction reaction starting at line: " << myKinetics->indices_of_good_lines()[reaction_lines[j]] << std::endl;
						return false;
					}
				}

				std::cout << "   Found (to be exploded): " << n_abstraction_reactions << std::endl;
			}

			// Check for exisiting abstraction reactions
			{
				std::cout << " * Checking for abstraction reactions already available in the mechanism..." << std::endl;
				unsigned int n_existing_reactions = 0;
				for (unsigned int j = 0; j < number_of_reactions; j++)
					n_existing_reactions += abstractions_->RemoveExistingReactions(reactions_[j], names_species_);
				std::cout << "   Found: " << n_existing_reactions << std::endl;
			}

			// Exploding abstraction reactions
			{
				std::cout << " * Exploding abstraction reactions..." << std::endl;
				abstractions_->ExplodeReactions();
			}

			// Write abstraction reactions on EXT file
			{
				std::cout << " * Writing abstraction reactions on file..." << std::endl;
				boost::filesystem::path exploded_file_name = myKinetics->folder_path() / (myKinetics->file_name().string() + ".EXT");
				std::ofstream fExploded(exploded_file_name.string(), std::ios::out);
				for (unsigned int j = 0; j < myKinetics->good_lines().size(); j++)
				{
					std::vector<unsigned int>::iterator it = std::find(list_lines_to_be_exploded.begin(), list_lines_to_be_exploded.end(), j);

					if (it != list_lines_to_be_exploded.end())
					{
						const unsigned int k = std::distance(list_lines_to_be_exploded.begin(), it);
						abstractions_->exploded()[k].PrintExplodedReactions(fExploded);
					}
					else
					{
						fExploded << myKinetics->good_lines()[j] << " " << myKinetics->strong_comments()[j] << std::endl;
					}
				}
				fExploded.close();

				// Final message before leaving the pre-processor
				std::cout << "   Abstractions correctly written on the output file: " << exploded_file_name.string() << std::endl;
				exit(-1);
			}
		}

		if (large_error_in_stoichiometries > 0)
		{
			std::cout << std::endl;
			std::cout << " ! ERROR MESSAGE: Large inconsistencies were found in the stoichiometries of " << large_error_in_stoichiometries << " reactions." << std::endl;
			std::cout << "                  Please check the log file for additional details." << std::endl;
			std::cout << std::endl;
			return false;
		}

		if (small_error_in_stoichiometries > 0)
		{
			std::cout << std::endl;
			std::cout << " ! WARNING MESSAGE: Small inconsistencies were found in the stoichiometries of " << small_error_in_stoichiometries << " reactions." << std::endl;
			std::cout << "                    Please check the log file for additional details." << std::endl;
			std::cout << std::endl;
		}

		// Check for duplicate reactions
		{
			std::cout << " * Looking for duplicate reactions... " << std::endl;
			
			for (unsigned int i=0; i<reactions_.size(); i++)
				for (unsigned int j=i+1; j<reactions_.size(); j++)
				{
					// Check direct vs direct
					if ( OpenSMOKE_Utilities::compare_vectors( reactions_[i].reactant_nu_indices(), reactions_[j].reactant_nu_indices() ) == true)
					{
						if ( OpenSMOKE_Utilities::compare_vectors( reactions_[i].product_nu_indices(), reactions_[j].product_nu_indices() ) == true)
						{
							if ( reactions_[i].Tag() == reactions_[j].Tag() )
							{
								if (reactions_[i].IsDuplicate() == false || reactions_[j].IsDuplicate() == false )
								{
									if (reactions_[i].pressureDependentSpeciesIndex() == reactions_[j].pressureDependentSpeciesIndex())
									{
										std::cout << "The following reactions must be declared as DUPLICATE" << std::endl;
										std::cout << "Reaction " << i+1 << " starting at line: " << myKinetics->indices_of_good_lines()[reaction_lines[i]] << std::endl;
										std::cout << "Reaction " << j+1 << " starting at line: " << myKinetics->indices_of_good_lines()[reaction_lines[j]] << std::endl;
										return false;
									}
								}
							}
						}
					}

					// Check direct vs reverse
					if (reactions_[i].IsReversible() || reactions_[j].IsReversible())
					{
						if (OpenSMOKE_Utilities::compare_vectors(reactions_[i].reactant_nu_indices(), reactions_[j].product_nu_indices()) == true)
						{
							if (OpenSMOKE_Utilities::compare_vectors(reactions_[i].product_nu_indices(), reactions_[j].reactant_nu_indices()) == true)
							{
								if (reactions_[i].Tag() == reactions_[j].Tag())
								{
									if (reactions_[i].IsDuplicate() == false || reactions_[j].IsDuplicate() == false)
									{
										if (reactions_[i].pressureDependentSpeciesIndex() == reactions_[j].pressureDependentSpeciesIndex())
										{
											std::cout << "The following reactions must be declared as DUPLICATE" << std::endl;
											std::cout << "Reaction " << i + 1 << " starting at line: " << myKinetics->indices_of_good_lines()[reaction_lines[i]] << std::endl;
											std::cout << "Reaction " << j + 1 << " starting at line: " << myKinetics->indices_of_good_lines()[reaction_lines[j]] << std::endl;
											return false;
										}
									}
								}
							}
						}
					}
				}
		}

		std::cout << " * Reactions correctly imported!" << std::endl;

		return true;
	}

	template<typename Reactions>
	bool PreProcessorKineticsPolicy_CHEMKIN<Reactions>::KineticsToCHEMKINFile(AtomicCompositionTable& atomicComposition, boost::filesystem::path file_path)
	{
		std::ofstream fOut;
		fOut.open(file_path.string().c_str(), std::ios::out);

		// Elements section
		{
			fOut << "ELEMENTS" << std::endl;
			for (unsigned int j = 0; j < atomicComposition.element_names_list().size(); j++)
				fOut << atomicComposition.element_names_list()[j] << std::endl;
			fOut << "END" << std::endl;
			fOut << std::endl;
		}

		// Species section
		{
			fOut << "SPECIES" << std::endl;
			unsigned int count = 1;
			for (;;)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					if (count <= names_species().size())
						fOut << std::setw(30) << std::left << names_species()[count - 1];
					count++;
				}
				fOut << std::endl;
				if (count > names_species().size())
					break;
			}
			fOut << "END" << std::endl;
			fOut << std::endl;
		}

		// Reactions section
		{
			fOut << "REACTIONS" << std::endl;
			fOut << std::endl;

			for (unsigned int j = 0; j < reactions_.size(); j++)
			{
				std::stringstream reaction_data;
				reactions_[j].GetReactionStringCHEMKIN(names_species(), reaction_data, myKinetics->strong_comments()[reaction_lines[j]]);
				fOut << reaction_data.str();
				fOut << std::endl;
			}
			
			fOut << "END" << std::endl;
			fOut << std::endl;
		}

		fOut.close();

		return true;
	}

	template<typename Reactions>
	template<typename Thermodynamics>
	void PreProcessorKineticsPolicy_CHEMKIN<Reactions>::WriteShortSummaryOnASCIIFile(const std::string file_name, Thermodynamics& thermodynamics) const
	{
		std::cout << " * Writing the summary of gas-phase kinetic mechanism..." << std::endl;

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		thermodynamics.WriteElementTableOnASCIIFile(fOutput);

		fOutput << "---------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                                  CHEMICAL REACTIONS                                   " << std::endl;
		fOutput << std::endl;
		fOutput << "                          Units: [mol, cm3, s] and [cal/mol]                           " << std::endl;
		fOutput << "---------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
		fOutput << std::endl;

		unsigned int count = 1;
		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
		{
			fOutput << std::setw(7) << std::right << count++;
			fOutput << ". ";
			(*it).WriteShortSummary(fOutput, names_species_);
			fOutput << std::endl;
			fOutput << std::endl;
		}

		fOutput.close();
	}

	template<typename T>
	void WriteObjectASCIIFileOldStyle(const T& v, std::ostream& fOut)
	{
		fOut << v.size() << std::endl;
		for(int i=0;i<v.size();i++)
			fOut << 	v(i) << " ";
		fOut << std::endl;
	}

	template<typename Reactions>
	bool PreProcessorKineticsPolicy_CHEMKIN<Reactions>::WriteASCIIFileOldStyle(const std::string file_name) const
	{
		std::cout << " * Writing the interpreted kinetic file in ASCII format for old versions of OpenSMOKE..." << std::endl;

		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			if ( (*it).IsExplicitlyReversible()==true)
				ErrorMessage("PreProcessorKineticsPolicy_CHEMKIN<Reactions>::WriteASCIIFileOldStyle(const std::string file_name)",
							 "The kinetic mechanism contains one or more REV reactions. Please write them as separate reactions to continue...");

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		// 1. Number of species
		fOutput << names_species_.size() << std::endl;

		// 2. Number of reactions
		fOutput << reactions_.size() << std::endl;

		// 3. Kinetic data for each reaction
		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			(*it).WriteKineticsDataOnASCIIFileOldStyle(fOutput);
	
		// 4. Stoichiometric data
		WriteStoichiometricDataOnASCIIFile(fOutput);

		// 6. Writing reaction strings
		{
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				std::string line_reaction;
				(*it).GetReactionString(names_species_, line_reaction);
				boost::erase_all(line_reaction, " ");
				fOutput << line_reaction << std::endl; 
			}
		}

		// 7. Reaction orders
		{
			Eigen::VectorXd forwardOrders(reactions_.size()); 
			Eigen::VectorXd backwardOrders(reactions_.size()); 
			forwardOrders.setConstant(0.);
			backwardOrders.setConstant(0.);
			for (unsigned int k=0; k<reactions_.size(); k++)
			{
				forwardOrders(k) = reactions_[k].sumLambdaReactants();
				backwardOrders(k) = reactions_[k].sumLambdaProducts();
			}

			WriteObjectASCIIFileOldStyle(forwardOrders, fOutput);
			WriteObjectASCIIFileOldStyle(backwardOrders, fOutput);
		}

		// 8. Old data (TODO reaction orders different than stoichiometric coefficients)
		{
			fOutput << 0 << std::endl;
			fOutput << 0 << std::endl;
		}

		fOutput.close();




		return true;
	}

	template<typename Reactions>
	bool PreProcessorKineticsPolicy_CHEMKIN<Reactions>::WriteStoichiometricDataOnASCIIFile(std::ostream& fOutput) const	
	{
		{
			Eigen::VectorXi numDir1(names_species_.size()); numDir1.setConstant(0);
			Eigen::VectorXi numDir2(names_species_.size()); numDir2.setConstant(0);
			Eigen::VectorXi numDir3(names_species_.size()); numDir3.setConstant(0);
			Eigen::VectorXi numDir4(names_species_.size()); numDir4.setConstant(0);
			Eigen::VectorXi numDir5(names_species_.size()); numDir5.setConstant(0);

			Eigen::VectorXi jDir1;
			Eigen::VectorXi jDir2;
			Eigen::VectorXi jDir3;
			Eigen::VectorXi jDir4;
			Eigen::VectorXi jDir5;
			Eigen::VectorXd vDir5;

			Eigen::VectorXi numInvTot1(names_species_.size()); numInvTot1.setConstant(0);
			Eigen::VectorXi numInvTot2(names_species_.size()); numInvTot2.setConstant(0);
			Eigen::VectorXi numInvTot3(names_species_.size()); numInvTot3.setConstant(0);
			Eigen::VectorXi numInvTot4(names_species_.size()); numInvTot4.setConstant(0);
			Eigen::VectorXi numInvTot5(names_species_.size()); numInvTot5.setConstant(0);

			Eigen::VectorXi jInvTot1;
			Eigen::VectorXi jInvTot2;
			Eigen::VectorXi jInvTot3;
			Eigen::VectorXi jInvTot4;
			Eigen::VectorXi jInvTot5;
			Eigen::VectorXd vInvTot5;

			Eigen::VectorXi numInvEq1(names_species_.size()); numInvEq1.setConstant(0);
			Eigen::VectorXi numInvEq2(names_species_.size()); numInvEq2.setConstant(0);
			Eigen::VectorXi numInvEq3(names_species_.size()); numInvEq3.setConstant(0);
			Eigen::VectorXi numInvEq4(names_species_.size()); numInvEq4.setConstant(0);
			Eigen::VectorXi numInvEq5(names_species_.size()); numInvEq5.setConstant(0);

			Eigen::VectorXi jInvEq1;
			Eigen::VectorXi jInvEq2;
			Eigen::VectorXi jInvEq3;
			Eigen::VectorXi jInvEq4;
			Eigen::VectorXi jInvEq5;
			Eigen::VectorXd vInvEq5;

			{
				for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				{
					for (unsigned int i=0;i<(*it).reactant_nu_indices().size();i++)
					{
						if ( (*it).reactant_nu()[i] == 1. )			numDir1((*it).reactant_nu_indices()[i])++;
						else if ( (*it).reactant_nu()[i] == 2.  )	numDir2((*it).reactant_nu_indices()[i])++;
						else if ( (*it).reactant_nu()[i] == 3.  )	numDir3((*it).reactant_nu_indices()[i])++;
						else if ( (*it).reactant_nu()[i] == 0.5 )	numDir4((*it).reactant_nu_indices()[i])++;
						else 										numDir5((*it).reactant_nu_indices()[i])++;
					}
				}

				jDir1.resize(numDir1.sum()); jDir1.setConstant(0);
				jDir2.resize(numDir2.sum()); jDir2.setConstant(0);
				jDir3.resize(numDir3.sum()); jDir3.setConstant(0);
				jDir4.resize(numDir4.sum()); jDir4.setConstant(0);
				jDir5.resize(numDir5.sum()); jDir5.setConstant(0);
				vDir5.resize(numDir5.sum()); vDir5.setConstant(0.);

				const unsigned int sumDirTot = numDir1.sum()+numDir2.sum()+numDir3.sum()+numDir4.sum()+numDir5.sum();

				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					for (unsigned int k=0; k<reactions_.size(); k++)
					{
						for (unsigned int i=0;i<reactions_[k].reactant_nu_indices().size();i++)
						{
							if ( reactions_[k].reactant_nu_indices()[i] == j )
							{
								if (reactions_[k].reactant_nu()[i] == 1.)
								{
									jDir1(countGlobal1++) = k+1;
								}
								else if (reactions_[k].reactant_nu()[i] == 2.)
								{
									jDir2(countGlobal2++) = k+1;
								}
								else if (reactions_[k].reactant_nu()[i] == 3.)
								{
									jDir3(countGlobal3++) = k+1;
								}
								else if (reactions_[k].reactant_nu()[i] == 0.5)
								{
									jDir4(countGlobal4++) = k+1;
								}
								else
								{
									jDir5(countGlobal5) = k+1;
									vDir5(countGlobal5) = reactions_[k].reactant_nu()[i];
									countGlobal5++;
								}
								countGlobal++;
							}
						}
						if (countGlobal == sumDirTot)
							break;
					}
				}
			}
			
			{
				for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				{
					for (unsigned int i=0;i<(*it).product_nu_indices().size();i++)
					{
						if ( (*it).product_nu()[i] == 1. )			numInvTot1((*it).product_nu_indices()[i])++;
						else if ( (*it).product_nu()[i] == 2.  )	numInvTot2((*it).product_nu_indices()[i])++;
						else if ( (*it).product_nu()[i] == 3.  )	numInvTot3((*it).product_nu_indices()[i])++;
						else if ( (*it).product_nu()[i] == 0.5 )	numInvTot4((*it).product_nu_indices()[i])++;
						else 										numInvTot5((*it).product_nu_indices()[i])++;
					}
				}

				jInvTot1.resize(numInvTot1.sum()); jInvTot1.setConstant(0);
				jInvTot2.resize(numInvTot2.sum()); jInvTot2.setConstant(0);
				jInvTot3.resize(numInvTot3.sum()); jInvTot3.setConstant(0);
				jInvTot4.resize(numInvTot4.sum()); jInvTot4.setConstant(0);
				jInvTot5.resize(numInvTot5.sum()); jInvTot5.setConstant(0);
				vInvTot5.resize(numInvTot5.sum()); vInvTot5.setConstant(0.);

				const unsigned int sumInvTot = numInvTot1.sum()+numInvTot2.sum()+numInvTot3.sum()+numInvTot4.sum()+numInvTot5.sum();
		
				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					for (unsigned int k=0; k<reactions_.size(); k++)
					{
						for (unsigned int i=0;i<reactions_[k].product_nu_indices().size();i++)
						{
							if ( reactions_[k].product_nu_indices()[i] == j )
							{
								if (reactions_[k].product_nu()[i] == 1.)
								{
									jInvTot1(countGlobal1++) = k+1;
								}
								else if (reactions_[k].product_nu()[i] == 2.)
								{
									jInvTot2(countGlobal2++) = k+1;
								}
								else if (reactions_[k].product_nu()[i] == 3.)
								{
									jInvTot3(countGlobal3++) = k+1;
								}
								else if (reactions_[k].product_nu()[i] == 0.5)
								{
									jInvTot4(countGlobal4++) = k+1;
								}
								else
								{
									jInvTot5(countGlobal5) = k+1;
									vInvTot5(countGlobal5) = reactions_[k].product_nu()[i];
									countGlobal5++;
								}
								countGlobal++;
							}
						}
						if (countGlobal == sumInvTot)
							break;
					}
				}
			}
			
			{
				for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				{
					if ((*it).IsReversible() == true)
					{
						for (unsigned int i=0;i<(*it).product_nu_indices().size();i++)
						{
							if ( (*it).product_nu()[i] == 1. )			numInvEq1((*it).product_nu_indices()[i])++;
							else if ( (*it).product_nu()[i] == 2.  )	numInvEq2((*it).product_nu_indices()[i])++;
							else if ( (*it).product_nu()[i] == 3.  )	numInvEq3((*it).product_nu_indices()[i])++;
							else if ( (*it).product_nu()[i] == 0.5 )	numInvEq4((*it).product_nu_indices()[i])++;
							else 										numInvEq5((*it).product_nu_indices()[i])++;
						}
					}
				}

				jInvEq1.resize(numInvEq1.sum()); jInvEq1.setConstant(0);
				jInvEq2.resize(numInvEq2.sum()); jInvEq2.setConstant(0);
				jInvEq3.resize(numInvEq3.sum()); jInvEq3.setConstant(0);
				jInvEq4.resize(numInvEq4.sum()); jInvEq4.setConstant(0);
				jInvEq5.resize(numInvEq5.sum()); jInvEq5.setConstant(0);
				vInvEq5.resize(numInvEq5.sum()); vInvEq5.setConstant(0.);

				const unsigned int sumInvEq = numInvEq1.sum()+numInvEq2.sum()+numInvEq3.sum()+numInvEq4.sum()+numInvEq5.sum();
		
				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					for (unsigned int k=0; k<reactions_.size(); k++)
					{
						if (reactions_[k].IsReversible() == true)
						{
							for (unsigned int i=0;i<reactions_[k].product_nu_indices().size();i++)
							{
								if ( reactions_[k].product_nu_indices()[i] == j )
								{
									if (reactions_[k].product_nu()[i] == 1.)
									{
										jInvEq1(countGlobal1++) = k+1;
									}
									else if (reactions_[k].product_nu()[i] == 2.)
									{
										jInvEq2(countGlobal2++) = k+1;
									}
									else if (reactions_[k].product_nu()[i] == 3.)
									{
										jInvEq3(countGlobal3++) = k+1;
									}
									else if (reactions_[k].product_nu()[i] == 0.5)
									{
										jInvEq4(countGlobal4++) = k+1;
									}
									else
									{
										jInvEq5(countGlobal5) = k+1;
										vInvEq5(countGlobal5) = reactions_[k].product_nu()[i];
										countGlobal5++;
									}
									countGlobal++;
								}
							}
						}
						if (countGlobal == sumInvEq)
							break;
					}
				}
			}
			
			WriteObjectASCIIFileOldStyle(numDir1, fOutput);
			WriteObjectASCIIFileOldStyle(numDir2, fOutput);
			WriteObjectASCIIFileOldStyle(numDir3, fOutput);
			WriteObjectASCIIFileOldStyle(numDir4, fOutput);
			WriteObjectASCIIFileOldStyle(numDir5, fOutput);
			
			WriteObjectASCIIFileOldStyle(numInvTot1, fOutput);
			WriteObjectASCIIFileOldStyle(numInvTot2, fOutput);
			WriteObjectASCIIFileOldStyle(numInvTot3, fOutput);
			WriteObjectASCIIFileOldStyle(numInvTot4, fOutput);
			WriteObjectASCIIFileOldStyle(numInvTot5, fOutput);

			WriteObjectASCIIFileOldStyle(numInvEq1, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq2, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq3, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq4, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq5, fOutput);

			WriteObjectASCIIFileOldStyle(jDir1, fOutput);
			WriteObjectASCIIFileOldStyle(jDir2, fOutput);
			WriteObjectASCIIFileOldStyle(jDir3, fOutput);
			WriteObjectASCIIFileOldStyle(jDir4, fOutput);
			WriteObjectASCIIFileOldStyle(jDir5, fOutput);
			WriteObjectASCIIFileOldStyle(vDir5, fOutput);

			WriteObjectASCIIFileOldStyle(jInvTot1, fOutput);
			WriteObjectASCIIFileOldStyle(jInvTot2, fOutput);
			WriteObjectASCIIFileOldStyle(jInvTot3, fOutput);
			WriteObjectASCIIFileOldStyle(jInvTot4, fOutput);
			WriteObjectASCIIFileOldStyle(jInvTot5, fOutput);
			WriteObjectASCIIFileOldStyle(vInvTot5, fOutput);

			WriteObjectASCIIFileOldStyle(jInvEq1, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq2, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq3, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq4, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq5, fOutput);
			WriteObjectASCIIFileOldStyle(vInvEq5, fOutput);
		}
		
		// 5. Sum of stoichiometric coefficients
		{
			Eigen::VectorXd sumStoichiometricCoefficients(reactions_.size()); 
			sumStoichiometricCoefficients.setConstant(0.);
			for (unsigned int k=0; k<reactions_.size(); k++)
				sumStoichiometricCoefficients(k) = reactions_[k].sumNuProducts()-reactions_[k].sumNuReactants();

			WriteObjectASCIIFileOldStyle(sumStoichiometricCoefficients, fOutput);
		}
		
		return true;
	}

	template<typename Reactions>
	bool PreProcessorKineticsPolicy_CHEMKIN<Reactions>::WriteReactionOrdersOnASCIIFile(std::ostream& fOutput) const	
	{
		bool write_reaction_orders = false;
		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			if ((*it).IsFORD() == true || (*it).IsRORD() == true)
			{
				write_reaction_orders = true;
				break;
			}

		fOutput << write_reaction_orders << std::endl;

		if (write_reaction_orders == true)
		{
		{
			Eigen::VectorXi numDir1(names_species_.size()); numDir1.setConstant(0);
			Eigen::VectorXi numDir2(names_species_.size()); numDir2.setConstant(0);
			Eigen::VectorXi numDir3(names_species_.size()); numDir3.setConstant(0);
			Eigen::VectorXi numDir4(names_species_.size()); numDir4.setConstant(0);
			Eigen::VectorXi numDir5(names_species_.size()); numDir5.setConstant(0);

			Eigen::VectorXi jDir1;
			Eigen::VectorXi jDir2;
			Eigen::VectorXi jDir3;
			Eigen::VectorXi jDir4;
			Eigen::VectorXi jDir5;
			Eigen::VectorXd vDir5;

			Eigen::VectorXi numInvEq1(names_species_.size()); numInvEq1.setConstant(0);
			Eigen::VectorXi numInvEq2(names_species_.size()); numInvEq2.setConstant(0);
			Eigen::VectorXi numInvEq3(names_species_.size()); numInvEq3.setConstant(0);
			Eigen::VectorXi numInvEq4(names_species_.size()); numInvEq4.setConstant(0);
			Eigen::VectorXi numInvEq5(names_species_.size()); numInvEq5.setConstant(0);

			Eigen::VectorXi jInvEq1;
			Eigen::VectorXi jInvEq2;
			Eigen::VectorXi jInvEq3;
			Eigen::VectorXi jInvEq4;
			Eigen::VectorXi jInvEq5;
			Eigen::VectorXd vInvEq5;

			{
				for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				{
					for (unsigned int i=0;i<(*it).reactant_lambda_indices().size();i++)
					{
						if ( (*it).reactant_lambda()[i] == 1. )			numDir1((*it).reactant_lambda_indices()[i])++;
						else if ( (*it).reactant_lambda()[i] == 2.  )	numDir2((*it).reactant_lambda_indices()[i])++;
						else if ( (*it).reactant_lambda()[i] == 3.  )	numDir3((*it).reactant_lambda_indices()[i])++;
						else if ( (*it).reactant_lambda()[i] == 0.5 )	numDir4((*it).reactant_lambda_indices()[i])++;
						else 											numDir5((*it).reactant_lambda_indices()[i])++;
					}
				}

				jDir1.resize(numDir1.sum()); jDir1.setConstant(0);
				jDir2.resize(numDir2.sum()); jDir2.setConstant(0);
				jDir3.resize(numDir3.sum()); jDir3.setConstant(0);
				jDir4.resize(numDir4.sum()); jDir4.setConstant(0);
				jDir5.resize(numDir5.sum()); jDir5.setConstant(0);
				vDir5.resize(numDir5.sum()); vDir5.setConstant(0.);

				const unsigned int sumDirTot = numDir1.sum()+numDir2.sum()+numDir3.sum()+numDir4.sum()+numDir5.sum();

				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					for (unsigned int k=0; k<reactions_.size(); k++)
					{
						for (unsigned int i=0;i<reactions_[k].reactant_lambda_indices().size();i++)
						{
							if ( reactions_[k].reactant_lambda_indices()[i] == j )
							{
								if (reactions_[k].reactant_lambda()[i] == 1.)
								{
									jDir1(countGlobal1++) = k+1;
								}
								else if (reactions_[k].reactant_lambda()[i] == 2.)
								{
									jDir2(countGlobal2++) = k+1;
								}
								else if (reactions_[k].reactant_lambda()[i] == 3.)
								{
									jDir3(countGlobal3++) = k+1;
								}
								else if (reactions_[k].reactant_lambda()[i] == 0.5)
								{
									jDir4(countGlobal4++) = k+1;
								}
								else
								{
									jDir5(countGlobal5) = k+1;
									vDir5(countGlobal5) = reactions_[k].reactant_lambda()[i];
									countGlobal5++;
								}
								countGlobal++;
							}
						}
						if (countGlobal == sumDirTot)
							break;
					}
				}
			}
			
			{
				for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				{
					if ((*it).IsReversible() == true)
					{
						for (unsigned int i=0;i<(*it).product_lambda_indices().size();i++)
						{
							if ( (*it).product_lambda()[i] == 1. )			numInvEq1((*it).product_lambda_indices()[i])++;
							else if ( (*it).product_lambda()[i] == 2.  )	numInvEq2((*it).product_lambda_indices()[i])++;
							else if ( (*it).product_lambda()[i] == 3.  )	numInvEq3((*it).product_lambda_indices()[i])++;
							else if ( (*it).product_lambda()[i] == 0.5 )	numInvEq4((*it).product_lambda_indices()[i])++;
							else 											numInvEq5((*it).product_lambda_indices()[i])++;
						}
					}
				}

				jInvEq1.resize(numInvEq1.sum()); jInvEq1.setConstant(0);
				jInvEq2.resize(numInvEq2.sum()); jInvEq2.setConstant(0);
				jInvEq3.resize(numInvEq3.sum()); jInvEq3.setConstant(0);
				jInvEq4.resize(numInvEq4.sum()); jInvEq4.setConstant(0);
				jInvEq5.resize(numInvEq5.sum()); jInvEq5.setConstant(0);
				vInvEq5.resize(numInvEq5.sum()); vInvEq5.setConstant(0.);

				const unsigned int sumInvEq = numInvEq1.sum()+numInvEq2.sum()+numInvEq3.sum()+numInvEq4.sum()+numInvEq5.sum();
		
				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					for (unsigned int k=0; k<reactions_.size(); k++)
					{
						if (reactions_[k].IsReversible() == true)
						{
							for (unsigned int i=0;i<reactions_[k].product_lambda_indices().size();i++)
							{
								if ( reactions_[k].product_lambda_indices()[i] == j )
								{
									if (reactions_[k].product_lambda()[i] == 1.)
									{
										jInvEq1(countGlobal1++) = k+1;
									}
									else if (reactions_[k].product_lambda()[i] == 2.)
									{
										jInvEq2(countGlobal2++) = k+1;
									}
									else if (reactions_[k].product_lambda()[i] == 3.)
									{
										jInvEq3(countGlobal3++) = k+1;
									}
									else if (reactions_[k].product_lambda()[i] == 0.5)
									{
										jInvEq4(countGlobal4++) = k+1;
									}
									else
									{
										jInvEq5(countGlobal5) = k+1;
										vInvEq5(countGlobal5) = reactions_[k].product_lambda()[i];
										countGlobal5++;
									}
									countGlobal++;
								}
							}
						}
						if (countGlobal == sumInvEq)
							break;
					}
				}
			}
			
			WriteObjectASCIIFileOldStyle(numDir1, fOutput);
			WriteObjectASCIIFileOldStyle(numDir2, fOutput);
			WriteObjectASCIIFileOldStyle(numDir3, fOutput);
			WriteObjectASCIIFileOldStyle(numDir4, fOutput);
			WriteObjectASCIIFileOldStyle(numDir5, fOutput);

			WriteObjectASCIIFileOldStyle(numInvEq1, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq2, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq3, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq4, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq5, fOutput);

			WriteObjectASCIIFileOldStyle(jDir1, fOutput);
			WriteObjectASCIIFileOldStyle(jDir2, fOutput);
			WriteObjectASCIIFileOldStyle(jDir3, fOutput);
			WriteObjectASCIIFileOldStyle(jDir4, fOutput);
			WriteObjectASCIIFileOldStyle(jDir5, fOutput);
			WriteObjectASCIIFileOldStyle(vDir5, fOutput);

			WriteObjectASCIIFileOldStyle(jInvEq1, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq2, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq3, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq4, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq5, fOutput);
			WriteObjectASCIIFileOldStyle(vInvEq5, fOutput);
		}
		}// write_reaction_orders
		
		return true;
	}


	template<typename Reactions>
	bool PreProcessorKineticsPolicy_CHEMKIN<Reactions>::WriteASCIIFile(const std::string file_name) const
	{
		std::cout << " * Writing the interpreted kinetic file in ASCII format..." << std::endl;

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out | std::ios::binary);
		fOutput.setf(std::ios::scientific);

		// 1. Number of species
		fOutput << "number-of-species" << std::endl;
		fOutput << names_species_.size() << std::endl;

		// 2. Number of reactions
		fOutput << "number-of-reactions" << std::endl;
		fOutput << reactions_.size() << std::endl;

		// 3. Kinetic data for each reaction
		unsigned int number_of_irreversible_reactions = 0;
		unsigned int number_of_reversible_reactions = 0;
		unsigned int number_of_explicitly_reversible_reactions = 0;
		unsigned int number_of_thermodynamic_reversible_reactions = 0;
		unsigned int number_of_thirdbody_reactions = 0;
		unsigned int number_of_falloff_reactions = 0;
		unsigned int number_of_cabr_reactions = 0;
		unsigned int number_of_chebyshev_reactions = 0;
		
		unsigned int number_of_pressurelog_reactions = 0;
		unsigned int number_of_extendedpressurelog_reactions = 0;
		unsigned int number_of_extendedfalloff_reactions = 0;
		unsigned int number_of_fit1_reactions = 0;
		unsigned int number_of_janevlanger_reactions = 0;
		unsigned int number_of_landauteller_reactions = 0;

		unsigned int number_of_falloff_lindemann_reactions = 0;
		unsigned int number_of_falloff_troe_reactions = 0;
		unsigned int number_of_falloff_sri_reactions = 0;
		unsigned int number_of_cabr_lindemann_reactions = 0;
		unsigned int number_of_cabr_troe_reactions = 0;
		unsigned int number_of_cabr_sri_reactions = 0;

		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
		{
			if ( (*it).IsReversible()==true)
			{
				number_of_reversible_reactions++;
				if ( (*it).IsExplicitlyReversible()==true)	number_of_explicitly_reversible_reactions++;
				else										number_of_thermodynamic_reversible_reactions++;
			}
			else
				number_of_irreversible_reactions++;
			
			if ( (*it).Tag() == PhysicalConstants::REACTION_THIRDBODY) number_of_thirdbody_reactions++;
			else if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF) { number_of_falloff_reactions++; number_of_falloff_lindemann_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_FALLOFF) {number_of_falloff_reactions++; number_of_falloff_troe_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_FALLOFF) {number_of_falloff_reactions++; number_of_falloff_sri_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR) {number_of_cabr_reactions++; number_of_cabr_lindemann_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_CABR) {number_of_cabr_reactions++; number_of_cabr_troe_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_CABR) {number_of_cabr_reactions++; number_of_cabr_sri_reactions++; }
			else if ((*it).Tag() == PhysicalConstants::REACTION_CHEBYSHEV) number_of_chebyshev_reactions++;
			else if ((*it).Tag() == PhysicalConstants::REACTION_EXTENDEDFALLOFF) number_of_extendedfalloff_reactions++;
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE)
			{
				if ( (*it).IsPressureLog() == true )				number_of_pressurelog_reactions++;
				else if ((*it).IsExtendedPressureLog() == true)		number_of_extendedpressurelog_reactions++;
				else if ( (*it).IsJanevLanger() == true )			number_of_janevlanger_reactions++;
				else if ( (*it).IsFit1() == true )					number_of_fit1_reactions++;
				else if ( (*it).IsLandauTeller() == true )			number_of_landauteller_reactions++;
			}
		}

		std::vector<unsigned int> indices_of_irreversible_reactions(number_of_irreversible_reactions);
		std::vector<unsigned int> indices_of_reversible_reactions(number_of_reversible_reactions);
		std::vector<unsigned int> indices_of_explicitly_reversible_reactions(number_of_explicitly_reversible_reactions);
		std::vector<unsigned int> indices_of_thermodynamic_reversible_reactions(number_of_thermodynamic_reversible_reactions);
		std::vector<unsigned int> indices_of_thirdbody_reactions(number_of_thirdbody_reactions);
		std::vector<unsigned int> indices_of_falloff_reactions(number_of_falloff_reactions);
		std::vector<unsigned int> indices_of_cabr_reactions(number_of_cabr_reactions);
		std::vector<unsigned int> indices_of_chebyshev_reactions(number_of_chebyshev_reactions);
		std::vector<unsigned int> indices_of_pressurelog_reactions(number_of_pressurelog_reactions);
		std::vector<unsigned int> indices_of_extendedpressurelog_reactions(number_of_extendedpressurelog_reactions);
		std::vector<unsigned int> indices_of_extendedfalloff_reactions(number_of_extendedfalloff_reactions);
		std::vector<unsigned int> indices_of_janevlanger_reactions(number_of_janevlanger_reactions);
		std::vector<unsigned int> indices_of_fit1_reactions(number_of_fit1_reactions);
		std::vector<unsigned int> indices_of_landauteller_reactions(number_of_landauteller_reactions);

		std::vector<unsigned int> indices_of_falloff_lindemann_reactions(number_of_falloff_lindemann_reactions);
		std::vector<unsigned int> indices_of_cabr_lindemann_reactions(number_of_cabr_lindemann_reactions);
		std::vector<unsigned int> indices_of_falloff_troe_reactions(number_of_falloff_troe_reactions);
		std::vector<unsigned int> indices_of_cabr_troe_reactions(number_of_cabr_troe_reactions);
		std::vector<unsigned int> indices_of_falloff_sri_reactions(number_of_falloff_sri_reactions);
		std::vector<unsigned int> indices_of_cabr_sri_reactions(number_of_cabr_sri_reactions);

		{		
			unsigned int count_irreversible=0;
			unsigned int count_reversible=0;
			unsigned int count_explicit=0;
			unsigned int count_thermodynamic=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).IsReversible()==true) 
				{
						indices_of_reversible_reactions[count_reversible++] = it - reactions_.begin() + 1;
						if ( (*it).IsExplicitlyReversible() == true)  indices_of_explicitly_reversible_reactions[count_explicit++] = it - reactions_.begin() + 1;
						else										  indices_of_thermodynamic_reversible_reactions[count_thermodynamic++] = it - reactions_.begin() + 1;
				}
				else
					indices_of_irreversible_reactions[count_irreversible++] = it - reactions_.begin() + 1;
			}

			fOutput << "irreversible-reactions" << std::endl;
			fOutput << reactions_.size() - number_of_reversible_reactions << std::endl;
			for(unsigned int i=0;i<reactions_.size() - number_of_reversible_reactions;i++)
				fOutput << indices_of_irreversible_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "reversible-reactions" << std::endl;
			fOutput << number_of_reversible_reactions << std::endl;
			for(unsigned int i=0;i<number_of_reversible_reactions;i++)
				fOutput << indices_of_reversible_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "thermodynamic-reversible-reactions" << std::endl;
			fOutput << number_of_thermodynamic_reversible_reactions << std::endl;
			for(unsigned int i=0;i<number_of_thermodynamic_reversible_reactions;i++)
				fOutput << indices_of_thermodynamic_reversible_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "explicitly-reversible-reactions" << std::endl;
			fOutput << number_of_explicitly_reversible_reactions << std::endl;
			for(unsigned int i=0;i<number_of_explicitly_reversible_reactions;i++)
				fOutput << indices_of_explicitly_reversible_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_THIRDBODY) indices_of_thirdbody_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "thirdbody-reactions" << std::endl;
			fOutput << number_of_thirdbody_reactions << std::endl;
			for(unsigned int i=0;i<number_of_thirdbody_reactions;i++)
				fOutput << indices_of_thirdbody_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			unsigned int count_lindemann=0;
			unsigned int count_troe=0;
			unsigned int count_sri=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_falloff_lindemann_reactions[count_lindemann++] = it - reactions_.begin() + 1;
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_falloff_troe_reactions[count_troe++] = it - reactions_.begin() + 1;
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_falloff_sri_reactions[count_sri++] = it - reactions_.begin() + 1;
				}
			}

			fOutput << "falloff-reactions" << std::endl;
			fOutput << number_of_falloff_reactions << std::endl;
			for(unsigned int i=0;i<number_of_falloff_reactions;i++)
				fOutput << indices_of_falloff_reactions[i] << " ";
			fOutput << std::endl;
			
			fOutput << "falloff-lindemann-reactions" << std::endl;
			fOutput << number_of_falloff_lindemann_reactions << std::endl;
			for(unsigned int i=0;i<number_of_falloff_lindemann_reactions;i++)
				fOutput << indices_of_falloff_lindemann_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "falloff-troe-reactions" << std::endl;
			fOutput << number_of_falloff_troe_reactions << std::endl;
			for(unsigned int i=0;i<number_of_falloff_troe_reactions;i++)
				fOutput << indices_of_falloff_troe_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "falloff-sri-reactions" << std::endl;
			fOutput << number_of_falloff_sri_reactions << std::endl;
			for(unsigned int i=0;i<number_of_falloff_sri_reactions;i++)
				fOutput << indices_of_falloff_sri_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			unsigned int count_lindemann=0;
			unsigned int count_troe=0;
			unsigned int count_sri=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR ) 
				{
					indices_of_cabr_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_cabr_lindemann_reactions[count_lindemann++] = it - reactions_.begin() + 1;
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_CABR ) 
				{
					indices_of_cabr_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_cabr_troe_reactions[count_troe++] = it - reactions_.begin() + 1;
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_CABR ) 
				{
					indices_of_cabr_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_cabr_sri_reactions[count_sri++] = it - reactions_.begin() + 1;
				}
			}

			fOutput << "cabr-reactions" << std::endl;
			fOutput << number_of_cabr_reactions << std::endl;
			for(unsigned int i=0;i<number_of_cabr_reactions;i++)
				fOutput << indices_of_cabr_reactions[i] << " ";
			fOutput << std::endl;
			
			fOutput << "cabr-lindemann-reactions" << std::endl;
			fOutput << number_of_cabr_lindemann_reactions << std::endl;
			for(unsigned int i=0;i<number_of_cabr_lindemann_reactions;i++)
				fOutput << indices_of_cabr_lindemann_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "cabr-troe-reactions" << std::endl;
			fOutput << number_of_cabr_troe_reactions << std::endl;
			for(unsigned int i=0;i<number_of_cabr_troe_reactions;i++)
				fOutput << indices_of_cabr_troe_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "cabr-sri-reactions" << std::endl;
			fOutput << number_of_cabr_sri_reactions << std::endl;
			for(unsigned int i=0;i<number_of_cabr_sri_reactions;i++)
				fOutput << indices_of_cabr_sri_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_CHEBYSHEV ) indices_of_chebyshev_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "chebyshev-reactions" << std::endl;
			fOutput << number_of_chebyshev_reactions << std::endl;
			for(unsigned int i=0;i<number_of_chebyshev_reactions;i++)
				fOutput << indices_of_chebyshev_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsPressureLog() == true) indices_of_pressurelog_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "pressurelog-reactions" << std::endl;
			fOutput << number_of_pressurelog_reactions << std::endl;
			for(unsigned int i=0;i<number_of_pressurelog_reactions;i++)
				fOutput << indices_of_pressurelog_reactions[i] << " ";
			fOutput << std::endl;
		}

		{
			unsigned int count = 0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsExtendedPressureLog() == true) indices_of_extendedpressurelog_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "extended-pressurelog-reactions" << std::endl;
			fOutput << number_of_extendedpressurelog_reactions << std::endl;
			for (unsigned int i = 0; i<number_of_extendedpressurelog_reactions; i++)
				fOutput << indices_of_extendedpressurelog_reactions[i] << " ";
			fOutput << std::endl;
		}

		{
			unsigned int count = 0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_EXTENDEDFALLOFF) indices_of_extendedfalloff_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "extended-falloff-reactions" << std::endl;
			fOutput << number_of_extendedfalloff_reactions << std::endl;
			for (unsigned int i = 0; i<number_of_extendedfalloff_reactions; i++)
				fOutput << indices_of_extendedfalloff_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsFit1() == true) indices_of_fit1_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "fit1-reactions" << std::endl;
			fOutput << number_of_fit1_reactions << std::endl;
			for(unsigned int i=0;i<number_of_fit1_reactions;i++)
				fOutput << indices_of_fit1_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsJanevLanger() == true) indices_of_janevlanger_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "janevlanger-reactions" << std::endl;
			fOutput << number_of_janevlanger_reactions << std::endl;
			for(unsigned int i=0;i<number_of_janevlanger_reactions;i++)
				fOutput << indices_of_janevlanger_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsLandauTeller() == true) indices_of_landauteller_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "landauteller-reactions" << std::endl;
			fOutput << number_of_landauteller_reactions << std::endl;
			for(unsigned int i=0;i<number_of_landauteller_reactions;i++)
				fOutput << indices_of_landauteller_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			OpenSMOKEVectorDouble lnA(reactions_.size());
			OpenSMOKEVectorDouble Beta(reactions_.size());
			OpenSMOKEVectorDouble E_over_R(reactions_.size());
			OpenSMOKEVectorInt    negative_lnA;
			
			unsigned int j=1;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				lnA[j] = ((*it).A() == 0.) ? log(1.e-100) : log(std::fabs((*it).A()));

				// Negative frequency factors: list of reactions (1-index based)
				if ((*it).A() < 0.)
					negative_lnA.Append(j);
				
				Beta[j] = (*it).Beta();
				E_over_R[j] = (*it).E_over_R();
				j++;
			}

			fOutput << "lnA" << std::endl;
			lnA.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "negative-lnA" << std::endl;
			negative_lnA.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "Beta" << std::endl;
			Beta.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "E_over_R" << std::endl;
			E_over_R.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
		}

		{		
			OpenSMOKEVectorDouble lnA_reversible(number_of_explicitly_reversible_reactions);
			OpenSMOKEVectorDouble Beta_reversible(number_of_explicitly_reversible_reactions);
			OpenSMOKEVectorDouble E_over_R_reversible(number_of_explicitly_reversible_reactions);

			unsigned int j=1;
			for (unsigned int i=0;i<number_of_explicitly_reversible_reactions;i++)
			{
				double A = reactions_[indices_of_explicitly_reversible_reactions[i]-1].A_reversible();
				lnA_reversible[j] = (A == 0.) ? log(1.e-100) : log(A);
				Beta_reversible[j] = reactions_[indices_of_explicitly_reversible_reactions[i]-1].Beta_reversible();
				E_over_R_reversible[j] = reactions_[indices_of_explicitly_reversible_reactions[i]-1].E_over_R_reversible();
				j++;
			}

			fOutput << "reversible-lnA" << std::endl;
			lnA_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "reversible-Beta" << std::endl;
			Beta_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "reversible-E_over_R" << std::endl;
			E_over_R_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
		}

		{		
			fOutput << "thirdbody-parameters" << std::endl;
			for (unsigned int i=0;i<number_of_thirdbody_reactions;i++)
				reactions_[indices_of_thirdbody_reactions[i]-1].WriteThirdBodyParametersOnASCIIFile(fOutput);
		}

		{		
			fOutput << "falloff-kinetics" << std::endl;

			fOutput << "lnA-falloff-inf" << std::endl;
			fOutput << number_of_falloff_reactions << std::endl;
			for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				fOutput << log( reactions_[indices_of_falloff_reactions[i]-1].A_inf() )<< " ";
			fOutput << std::endl;

			fOutput << "Beta-falloff-inf" << std::endl;
			fOutput << number_of_falloff_reactions << std::endl;
			for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				fOutput << reactions_[indices_of_falloff_reactions[i]-1].Beta_inf() << " ";
			fOutput << std::endl;

			fOutput << "E_over_R-falloff-inf" << std::endl;
			fOutput << number_of_falloff_reactions << std::endl;
			for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				fOutput << reactions_[indices_of_falloff_reactions[i]-1].E_over_R_inf() << " ";
			fOutput << std::endl;
		}

		{		
			fOutput << "falloff-parameters" << std::endl;
			for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				reactions_[indices_of_falloff_reactions[i]-1].WritePressureDependentParametersOnASCIIFile(fOutput);
		}

		{		
			fOutput << "cabr-kinetics" << std::endl;
			fOutput << "lnA-cabr-inf" << std::endl;
			fOutput << number_of_cabr_reactions << std::endl;
			for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				fOutput << log( reactions_[indices_of_cabr_reactions[i]-1].A_inf() )<< " ";
			fOutput << std::endl;

			fOutput << "Beta-cabr-inf" << std::endl;
			fOutput << number_of_cabr_reactions << std::endl;
			for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				fOutput << reactions_[indices_of_cabr_reactions[i]-1].Beta_inf()<< " ";
			fOutput << std::endl;

			fOutput << "E_over_R-cabr-inf" << std::endl;
			fOutput << number_of_cabr_reactions << std::endl;
			for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				fOutput << reactions_[indices_of_cabr_reactions[i]-1].E_over_R_inf() << " ";
			fOutput << std::endl;
		}

		{		
			fOutput << "cabr-parameters" << std::endl;
			for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				reactions_[indices_of_cabr_reactions[i]-1].WritePressureDependentParametersOnASCIIFile(fOutput);
		}

		{
			fOutput << "chebyshev-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_chebyshev_reactions;j++)
				reactions_[indices_of_chebyshev_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "pressurelog-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_pressurelog_reactions;j++)
				reactions_[indices_of_pressurelog_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "extended-pressurelog-parameters " << std::endl;
			for (unsigned int j = 0; j<number_of_extendedpressurelog_reactions; j++)
				reactions_[indices_of_extendedpressurelog_reactions[j] - 1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "extended-falloff-parameters " << std::endl;
			for (unsigned int j = 0; j<number_of_extendedfalloff_reactions; j++)
				reactions_[indices_of_extendedfalloff_reactions[j] - 1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "fit1-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_fit1_reactions;j++)
				reactions_[indices_of_fit1_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "janevlanger-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_janevlanger_reactions;j++)
				reactions_[indices_of_janevlanger_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "landauteller-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_landauteller_reactions;j++)
				reactions_[indices_of_landauteller_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}

		{		
			fOutput << "stoichiometric-coefficients" << std::endl;
			WriteStoichiometricDataOnASCIIFile(fOutput);
			WriteReactionOrdersOnASCIIFile(fOutput);
		}

		return true;

	}

	template<typename T>
	void FormatXML(std::stringstream& xml_string, const std::string tag, std::vector<T>& v, bool close=true)
	{
			xml_string<< "<" << tag << ">" << std::endl;
			xml_string << v.size() << std::endl;
			if (v.size()>0)
			{
				for(unsigned int i=0;i<v.size();i++)
				{
					xml_string << v[i] << " ";
					if ((i+1)%30==0) xml_string << std::endl;
				}
				xml_string << std::endl;
			}
			if (close == true) xml_string << "</" << tag << ">" << std::endl;
	}


	template<typename Reactions>
	bool PreProcessorKineticsPolicy_CHEMKIN<Reactions>::WriteXMLFile(std::stringstream& fOutput) const
	{
		std::cout << " * Writing the interpreted kinetic file in XML format..." << std::endl;

		fOutput << "<Kinetics type=\"OpenSMOKE\" version=\"04-22-2013\">" << std::endl;
		
		// 2. Number of reactions
		fOutput << "<NumberOfReactions>" << std::endl;
		fOutput << reactions_.size() << std::endl;
		fOutput << "</NumberOfReactions>" << std::endl;

		// 3. Kinetic data for each reaction
		unsigned int number_of_irreversible_reactions = 0;
		unsigned int number_of_reversible_reactions = 0;
		unsigned int number_of_explicitly_reversible_reactions = 0;
		unsigned int number_of_thermodynamic_reversible_reactions = 0;
		unsigned int number_of_thirdbody_reactions = 0;
		unsigned int number_of_falloff_reactions = 0;
		unsigned int number_of_cabr_reactions = 0;
		unsigned int number_of_chebyshev_reactions = 0;
		
		unsigned int number_of_pressurelog_reactions = 0;
		unsigned int number_of_extendedpressurelog_reactions = 0;
		unsigned int number_of_extendedfalloff_reactions = 0;
		unsigned int number_of_fit1_reactions = 0;
		unsigned int number_of_janevlanger_reactions = 0;
		unsigned int number_of_landauteller_reactions = 0;

		unsigned int number_of_falloff_lindemann_reactions = 0;
		unsigned int number_of_falloff_troe_reactions = 0;
		unsigned int number_of_falloff_sri_reactions = 0;
		unsigned int number_of_cabr_lindemann_reactions = 0;
		unsigned int number_of_cabr_troe_reactions = 0;
		unsigned int number_of_cabr_sri_reactions = 0;

		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
		{
			if ( (*it).IsReversible()==true)
			{
				number_of_reversible_reactions++;
				if ( (*it).IsExplicitlyReversible()==true)	number_of_explicitly_reversible_reactions++;
				else										number_of_thermodynamic_reversible_reactions++;
			}
			else
				number_of_irreversible_reactions++;
			
			if ( (*it).Tag() == PhysicalConstants::REACTION_THIRDBODY) number_of_thirdbody_reactions++;
			else if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF) { number_of_falloff_reactions++; number_of_falloff_lindemann_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_FALLOFF) {number_of_falloff_reactions++; number_of_falloff_troe_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_FALLOFF) {number_of_falloff_reactions++; number_of_falloff_sri_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR) {number_of_cabr_reactions++; number_of_cabr_lindemann_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_CABR) {number_of_cabr_reactions++; number_of_cabr_troe_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_CABR) {number_of_cabr_reactions++; number_of_cabr_sri_reactions++; }
			else if ((*it).Tag() == PhysicalConstants::REACTION_CHEBYSHEV) number_of_chebyshev_reactions++;
			else if ((*it).Tag() == PhysicalConstants::REACTION_EXTENDEDFALLOFF) number_of_extendedfalloff_reactions++;
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE)
			{
				if ( (*it).IsPressureLog() == true )		number_of_pressurelog_reactions++;
				if ((*it).IsExtendedPressureLog() == true)	number_of_extendedpressurelog_reactions++;
				else if ( (*it).IsJanevLanger() == true )	number_of_janevlanger_reactions++;
				else if ( (*it).IsFit1() == true )			number_of_fit1_reactions++;
				else if ( (*it).IsLandauTeller() == true )	number_of_landauteller_reactions++;
			}
		}

		std::vector<unsigned int> indices_of_irreversible_reactions(number_of_irreversible_reactions);
		std::vector<unsigned int> indices_of_reversible_reactions(number_of_reversible_reactions);
		std::vector<unsigned int> indices_of_explicitly_reversible_reactions(number_of_explicitly_reversible_reactions);
		std::vector<unsigned int> indices_of_thermodynamic_reversible_reactions(number_of_thermodynamic_reversible_reactions);
		std::vector<unsigned int> indices_of_thirdbody_reactions(number_of_thirdbody_reactions);
		std::vector<unsigned int> indices_of_falloff_reactions(number_of_falloff_reactions);
		std::vector<unsigned int> indices_of_cabr_reactions(number_of_cabr_reactions);
		std::vector<unsigned int> indices_of_chebyshev_reactions(number_of_chebyshev_reactions);
		std::vector<unsigned int> indices_of_pressurelog_reactions(number_of_pressurelog_reactions);
		std::vector<unsigned int> indices_of_extendedpressurelog_reactions(number_of_extendedpressurelog_reactions);
		std::vector<unsigned int> indices_of_extendedfalloff_reactions(number_of_extendedfalloff_reactions);
		std::vector<unsigned int> indices_of_janevlanger_reactions(number_of_janevlanger_reactions);
		std::vector<unsigned int> indices_of_fit1_reactions(number_of_fit1_reactions);
		std::vector<unsigned int> indices_of_landauteller_reactions(number_of_landauteller_reactions);

		std::vector<unsigned int> indices_of_falloff_lindemann_reactions(number_of_falloff_lindemann_reactions);
		std::vector<unsigned int> indices_of_cabr_lindemann_reactions(number_of_cabr_lindemann_reactions);
		std::vector<unsigned int> indices_of_falloff_troe_reactions(number_of_falloff_troe_reactions);
		std::vector<unsigned int> indices_of_cabr_troe_reactions(number_of_cabr_troe_reactions);
		std::vector<unsigned int> indices_of_falloff_sri_reactions(number_of_falloff_sri_reactions);
		std::vector<unsigned int> indices_of_cabr_sri_reactions(number_of_cabr_sri_reactions);

		{		
			unsigned int count_irreversible=0;
			unsigned int count_reversible=0;
			unsigned int count_explicit=0;
			unsigned int count_thermodynamic=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).IsReversible()==true) 
				{
						indices_of_reversible_reactions[count_reversible++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
						if ((*it).IsExplicitlyReversible() == true)  indices_of_explicitly_reversible_reactions[count_explicit++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
						else										  indices_of_thermodynamic_reversible_reactions[count_thermodynamic++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
				}
				else
					indices_of_irreversible_reactions[count_irreversible++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
			}

			FormatXML(fOutput, "Irreversible", indices_of_irreversible_reactions);
			FormatXML(fOutput, "Reversible", indices_of_reversible_reactions);
			FormatXML(fOutput, "Reversible-Thermodynamics", indices_of_thermodynamic_reversible_reactions);
			FormatXML(fOutput, "Reversible-Explicit", indices_of_explicitly_reversible_reactions);

		}
	
		// Three-body reactions
		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_THIRDBODY) indices_of_thirdbody_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);

			FormatXML(fOutput, "ThreeBody", indices_of_thirdbody_reactions);
		}

		{		
			unsigned int count=0;
			unsigned int count_lindemann=0;
			unsigned int count_troe=0;
			unsigned int count_sri=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
					indices_of_falloff_lindemann_reactions[count_lindemann++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
					indices_of_falloff_troe_reactions[count_troe++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
					indices_of_falloff_sri_reactions[count_sri++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
				}
			}
			
			FormatXML(fOutput, "FallOff", indices_of_falloff_reactions, false);
				FormatXML(fOutput, "Lindemann", indices_of_falloff_lindemann_reactions);
				FormatXML(fOutput, "Troe", indices_of_falloff_troe_reactions);
				FormatXML(fOutput, "SRI", indices_of_falloff_sri_reactions);
			fOutput << "</FallOff>" << std::endl;
		}

		{		
			unsigned int count=0;
			unsigned int count_lindemann=0;
			unsigned int count_troe=0;
			unsigned int count_sri=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR ) 
				{
					indices_of_cabr_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
					indices_of_cabr_lindemann_reactions[count_lindemann++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_CABR ) 
				{
					indices_of_cabr_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
					indices_of_cabr_troe_reactions[count_troe++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_CABR ) 
				{
					indices_of_cabr_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
					indices_of_cabr_sri_reactions[count_sri++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);
				}
			}

			FormatXML(fOutput, "CABR", indices_of_cabr_reactions, false);
				FormatXML(fOutput, "Lindemann", indices_of_cabr_lindemann_reactions);
				FormatXML(fOutput, "Troe", indices_of_cabr_troe_reactions);
				FormatXML(fOutput, "SRI", indices_of_cabr_sri_reactions);
			fOutput << "</CABR>" << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_CHEBYSHEV) indices_of_chebyshev_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);

			FormatXML(fOutput, "Chebyshev", indices_of_chebyshev_reactions);
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsPressureLog() == true) indices_of_pressurelog_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);

			FormatXML(fOutput, "PressureLog", indices_of_pressurelog_reactions);
		}

		{
			unsigned int count = 0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsExtendedPressureLog() == true) indices_of_extendedpressurelog_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);

			FormatXML(fOutput, "ExtendedPressureLog", indices_of_extendedpressurelog_reactions);
		}

		{
			unsigned int count = 0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_EXTENDEDFALLOFF) indices_of_extendedfalloff_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);

			FormatXML(fOutput, "ExtendedFallOff", indices_of_extendedfalloff_reactions);
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsFit1() == true) indices_of_fit1_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);

			FormatXML(fOutput, "FIT1", indices_of_fit1_reactions);
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsJanevLanger() == true) indices_of_janevlanger_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);

			FormatXML(fOutput, "JanevLanger", indices_of_janevlanger_reactions);
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ((*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsLandauTeller() == true) indices_of_landauteller_reactions[count++] = boost::lexical_cast<int>(it - reactions_.begin() + 1);

			FormatXML(fOutput, "LandauTeller", indices_of_landauteller_reactions);
		}

		{		
			OpenSMOKEVectorDouble lnA(boost::lexical_cast<int>(reactions_.size()));
			OpenSMOKEVectorDouble Beta(boost::lexical_cast<int>(reactions_.size()));
			OpenSMOKEVectorDouble E_over_R(boost::lexical_cast<int>(reactions_.size()));
			OpenSMOKEVectorInt    negative_lnA;
			
			unsigned int j=1;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				lnA[j] = ((*it).A() == 0.) ? log(1.e-100) : log(std::fabs((*it).A()));
				
				// Negative frequency factors: list of reactions (1-index based)
				if ((*it).A() < 0.)
					negative_lnA.Append(j);

				Beta[j] = (*it).Beta();
				E_over_R[j] = (*it).E_over_R();
				j++;
			}

			fOutput << "<KineticParameters>" << std::endl;

			fOutput << "<Direct>" << std::endl;

			fOutput << "<lnA>" << std::endl;
			fOutput.precision(6);
			lnA.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
			fOutput << "</lnA>" << std::endl;

			fOutput << "<negative-lnA>" << std::endl;
			negative_lnA.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
			fOutput << "</negative-lnA>" << std::endl;

			fOutput << "<Beta>" << std::endl;
			Beta.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
			fOutput << "</Beta>" << std::endl;

			fOutput << "<E_over_R>" << std::endl;
			E_over_R.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
			fOutput << "</E_over_R>" << std::endl;
			
			fOutput << "</Direct>" << std::endl;
		}

		if (number_of_explicitly_reversible_reactions != 0)
		{		
			OpenSMOKEVectorDouble lnA_reversible(number_of_explicitly_reversible_reactions);
			OpenSMOKEVectorDouble Beta_reversible(number_of_explicitly_reversible_reactions);
			OpenSMOKEVectorDouble E_over_R_reversible(number_of_explicitly_reversible_reactions);

			unsigned int j=1;
			for (unsigned int i=0;i<number_of_explicitly_reversible_reactions;i++)
			{
				double A = reactions_[indices_of_explicitly_reversible_reactions[i]-1].A_reversible();
				lnA_reversible[j] = (A == 0.) ? log(1.e-100) : log(A);
				Beta_reversible[j] = reactions_[indices_of_explicitly_reversible_reactions[i]-1].Beta_reversible();
				E_over_R_reversible[j] = reactions_[indices_of_explicitly_reversible_reactions[i]-1].E_over_R_reversible();
				j++;
			}

			fOutput << "<Reverse>" << std::endl;

			fOutput << "<lnA>" << std::endl;
			lnA_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
			fOutput << "</lnA>" << std::endl;

			fOutput << "<Beta>" << std::endl;
			Beta_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
			fOutput << "</Beta>" << std::endl;

			fOutput << "<E_over_R>" << std::endl;
			E_over_R_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
			fOutput << "</E_over_R>" << std::endl;

			fOutput << "</Reverse>" << std::endl;
		}

		
		if (number_of_thirdbody_reactions != 0)
		{		
			fOutput << "<ThreeBody>" << std::endl;
			for (unsigned int i=0;i<number_of_thirdbody_reactions;i++)
			{
				fOutput << "<reaction i=\"" << indices_of_thirdbody_reactions[i] << "\">" << std::endl;
				reactions_[indices_of_thirdbody_reactions[i]-1].WriteThirdBodyParametersOnASCIIFile(fOutput);
				fOutput << "</reaction>" << std::endl;
			}
			fOutput << "</ThreeBody>" << std::endl;
		}

		if (number_of_falloff_reactions != 0)
		{
			fOutput << "<FallOff>" << std::endl;
			{		
				fOutput << "<lnA>" << std::endl;
				fOutput << number_of_falloff_reactions << std::endl;
				for (unsigned int i=0;i<number_of_falloff_reactions;i++)
					fOutput << log( reactions_[indices_of_falloff_reactions[i]-1].A_inf() )<< " ";
				fOutput << std::endl;
				fOutput << "</lnA>" << std::endl;

				fOutput << "<Beta>" << std::endl;
				fOutput << number_of_falloff_reactions << std::endl;
				for (unsigned int i=0;i<number_of_falloff_reactions;i++)
					fOutput << reactions_[indices_of_falloff_reactions[i]-1].Beta_inf() << " ";
				fOutput << std::endl;
				fOutput << "</Beta>" << std::endl;

				fOutput << "<E_over_R>" << std::endl;
				fOutput << number_of_falloff_reactions << std::endl;
				for (unsigned int i=0;i<number_of_falloff_reactions;i++)
					fOutput << reactions_[indices_of_falloff_reactions[i]-1].E_over_R_inf() << " ";
				fOutput << std::endl;
				fOutput << "</E_over_R>" << std::endl;
			}
		
			{		
				fOutput << "<Parameters>" << std::endl;
				for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				{
					fOutput << "<reaction i=\"" << indices_of_falloff_reactions[i] << "\">" << std::endl;
					reactions_[indices_of_falloff_reactions[i]-1].WritePressureDependentParametersOnASCIIFile(fOutput);
					fOutput << "</reaction>" << std::endl;
				}
				fOutput << "</Parameters>" << std::endl;
			}
			fOutput << "</FallOff>" << std::endl;
		}
		
		if (number_of_cabr_reactions != 0)
		{
			fOutput << "<CABR>" << std::endl;
			{		
				fOutput << "<lnA>" << std::endl;
				fOutput << number_of_cabr_reactions << std::endl;
				for (unsigned int i=0;i<number_of_cabr_reactions;i++)
					fOutput << log( reactions_[indices_of_cabr_reactions[i]-1].A_inf() )<< " ";
				fOutput << std::endl;
				fOutput << "</lnA>" << std::endl;

				fOutput << "<Beta>" << std::endl;
				fOutput << number_of_cabr_reactions << std::endl;
				for (unsigned int i=0;i<number_of_cabr_reactions;i++)
					fOutput << reactions_[indices_of_cabr_reactions[i]-1].Beta_inf()<< " ";
				fOutput << std::endl;
				fOutput << "</Beta>" << std::endl;

				fOutput << "<E_over_R>" << std::endl;
				fOutput << number_of_cabr_reactions << std::endl;
				for (unsigned int i=0;i<number_of_cabr_reactions;i++)
					fOutput << reactions_[indices_of_cabr_reactions[i]-1].E_over_R_inf() << " ";
				fOutput << std::endl;
				fOutput << "</E_over_R>" << std::endl;
			}

			{		
				fOutput << "<Parameters>" << std::endl;
				for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				{
					fOutput << "<reaction i=\"" << indices_of_cabr_reactions[i] << "\">" << std::endl;
					reactions_[indices_of_cabr_reactions[i]-1].WritePressureDependentParametersOnASCIIFile(fOutput);
					fOutput << "</reaction>" << std::endl;
				}
				fOutput << "</Parameters>" << std::endl;
			}
			fOutput << "</CABR>" << std::endl;
		}
		
		{
			if (number_of_chebyshev_reactions != 0)
			{
				fOutput << "<Chebyshev>" << std::endl;
				for(unsigned int j=0;j<number_of_chebyshev_reactions;j++)
					reactions_[indices_of_chebyshev_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
				fOutput << "</Chebyshev>" << std::endl;
			}

			if (number_of_pressurelog_reactions != 0)
			{
				fOutput << "<PressureLog>" << std::endl;
				for(unsigned int j=0;j<number_of_pressurelog_reactions;j++)
					reactions_[indices_of_pressurelog_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
				fOutput << "</PressureLog>" << std::endl;
			}

			if (number_of_extendedpressurelog_reactions != 0)
			{
				fOutput << "<ExtendedPressureLog>" << std::endl;
				for (unsigned int j = 0; j<number_of_extendedpressurelog_reactions; j++)
					reactions_[indices_of_extendedpressurelog_reactions[j] - 1].WriteAdditionalDataOnASCIIFile(fOutput);
				fOutput << "</ExtendedPressureLog>" << std::endl;
			}

			if (number_of_extendedfalloff_reactions != 0)
			{
				fOutput << "<ExtendedFallOff>" << std::endl;
				for (unsigned int j = 0; j<number_of_extendedfalloff_reactions; j++)
					reactions_[indices_of_extendedfalloff_reactions[j] - 1].WriteAdditionalDataOnASCIIFile(fOutput);
				fOutput << "</ExtendedFallOff>" << std::endl;
			}

			if (number_of_fit1_reactions != 0)
			{
				fOutput << "<FIT1>" << std::endl;
				for(unsigned int j=0;j<number_of_fit1_reactions;j++)
					reactions_[indices_of_fit1_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
				fOutput << "</FIT1>" << std::endl;
			}

			if (number_of_janevlanger_reactions != 0)
			{
				fOutput << "<JanevLanger>" << std::endl;
				for(unsigned int j=0;j<number_of_janevlanger_reactions;j++)
					reactions_[indices_of_janevlanger_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
				fOutput << "</JanevLanger>" << std::endl;
			}

			if (number_of_landauteller_reactions != 0)
			{
				fOutput << "<LandauTeller>" << std::endl;
				for(unsigned int j=0;j<number_of_landauteller_reactions;j++)
					reactions_[indices_of_landauteller_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
				fOutput << "</LandauTeller>" << std::endl;
			}
		}

		fOutput << "</KineticParameters>" << std::endl;
		
		// Write stoichiometric data
		{		
			fOutput << "<Stoichiometry type=\"OpenSMOKE\" version=\"04-22-2013\">" << std::endl;
			WriteStoichiometricDataOnASCIIFile(fOutput);
			WriteReactionOrdersOnASCIIFile(fOutput);
			fOutput << false << std::endl; // TOREMOVE

			fOutput << "</Stoichiometry>" << std::endl;
		}

		// Globar reactions (TODO)
		{
			fOutput << "<GlobalReactions type=\"none\">" << std::endl;
			fOutput << "</GlobalReactions>" << std::endl;
		}
		
		fOutput << "</Kinetics>" << std::endl;

		return true;

	}

}
