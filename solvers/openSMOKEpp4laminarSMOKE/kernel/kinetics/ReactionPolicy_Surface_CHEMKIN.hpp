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

#include <iomanip>
#include "math/PhysicalConstants.h"
#include "KineticsUtilityFunctions.h"
#include "ChebyshevPolynomialRateExpression.h"
#include "PressureLogarithmicRateExpression.h"

namespace OpenSMOKE
{
	ReactionPolicy_Surface_CHEMKIN::ReactionPolicy_Surface_CHEMKIN() 
	{
		SetDefaultUnits();
	}

	void ReactionPolicy_Surface_CHEMKIN::SetDefaultUnits()
	{
		iReversible_ = false;
		iExplicitlyReversible_ = false;
		iFord_ = false;
		iRord_ = false;
		iDuplicate_ = false;
		iStick_ = false;
		iMotzWiseCorrection_ = false;
		iCoverageDependent_ = false;
		iLangmuir_ = false;
		iLumped_ = false;
		name_of_lumped_function_ = "";
		iUBIQEP_ = false;
		ubiqep_reaction_class_ = 0;
		ubiqep_reaction_type_ = PhysicalConstants::UBIQEP_TYPE_DUMMY;
		surface_reaction_membership_ = 0;
		delta_nu_gas_ = 0.;
		delta_occupancy_sites_ = 0.;
	}

	ReactionPolicy_Surface_CHEMKIN::ReactionPolicy_Surface_CHEMKIN(const ReactionPolicy_Surface_CHEMKIN& orig) 
	{
	}

	ReactionPolicy_Surface_CHEMKIN::~ReactionPolicy_Surface_CHEMKIN() {
	}

	PhysicalConstants::TAG_REACTION_SURFACE ReactionPolicy_Surface_CHEMKIN::Tag() const
	{
		return tag_reaction_;
	}

	std::string ReactionPolicy_Surface_CHEMKIN::TagASCII() const
	{
		if (tag_reaction_ == PhysicalConstants::REACTION_SURFACE_SIMPLE)
			return "Simple";
		// TODO
		else
			return "Unknown reaction";
	}

	bool ReactionPolicy_Surface_CHEMKIN::SetUnits(const PhysicalConstants::UNITS_REACTION a_units, const PhysicalConstants::UNITS_REACTION e_units, const PhysicalConstants::UNITS_REACTION_COMPOSITION composition_units)
	{
		a_units_ = a_units;
		e_units_ = e_units;
		composition_units_ = composition_units;
		
		return true;
	}

	bool ReactionPolicy_Surface_CHEMKIN::SetMotzWiseCorrection(const bool iMotzWiseCorrection)
	{
		iMotzWiseCorrection_ = iMotzWiseCorrection;
		return true;
	}

	bool ReactionPolicy_Surface_CHEMKIN::SetGlobalNonConservationOfSites(const bool global_non_conservation_of_sites)
	{
		global_non_conservation_of_sites_ = global_non_conservation_of_sites;
		return true;
	}

	bool ReactionPolicy_Surface_CHEMKIN::ReadReactionFromStrings(	const std::vector<std::string>& lines, const std::map<std::string, unsigned int>& map_of_species,
																	const unsigned int number_of_gas_species, const unsigned int number_of_site_species, const unsigned int number_of_bulk_species,
																	const std::vector<double>& occupancy_site_species, const std::vector<unsigned int>& site_phase_membership, const std::vector<unsigned int>& bulk_phase_membership)
	{
		number_of_gas_species_ = number_of_gas_species;
		number_of_site_species_ = number_of_site_species;
		number_of_bulk_species_ = number_of_bulk_species;

		const std::string firstline = lines[0];
		
		std::string reconstructedline;
		std::string line_reactants;
		std::string line_products;
		
		std::vector<std::string> reactant_species;
		std::vector<std::string> product_species;

		if ( OpenSMOKE_Utilities::AnalyzeReactionLine(firstline, A_, beta_, E_, reconstructedline) == false) return false;
		if ( OpenSMOKE_Utilities::SeparateReactantSideFromProductSide(reconstructedline, line_reactants, line_products, iReversible_) == false) return false;

		if ( OpenSMOKE_Utilities::SeparateSpeciesAndStoichiometricCoefficients(line_reactants, reactant_species, reactant_nu_) == false) return false;
		if ( OpenSMOKE_Utilities::SeparateSpeciesAndStoichiometricCoefficients(line_products, product_species, product_nu_) == false) return false;

		// Reactants indices
		reactant_nu_indices_.resize(reactant_nu_.size());
		for(unsigned int i=0;i<reactant_nu_.size();i++)
		{
			std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(reactant_species[i]);
			if( it == map_of_species.end())
			{
				std::cout << "The following species is not available: " << reactant_species[i] << std::endl;
				return false;
			}
			else
			{
				reactant_nu_indices_[i] = (*it).second;
			}
		}

		// Product indices
		product_nu_indices_.resize(product_nu_.size());
		for(unsigned int i=0;i<product_nu_.size();i++)
		{
			std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(product_species[i]);
			if( it == map_of_species.end())
			{
				std::cout << "The following species is not available: " << product_species[i] << std::endl;
				return false;
			}
			else
			{
				product_nu_indices_[i] = (*it).second;
			}
		}

		// Save original stoichiometry (useful for UBI-QEP kinetic mechanisms)
		{
			original_reactant_nu_indices_ = reactant_nu_indices_;
			original_product_nu_indices_ = product_nu_indices_;
			original_reactant_nu_ = reactant_nu_;
			original_product_nu_ = product_nu_;
		}

		// Reorder indices
		{
			OpenSMOKE_Utilities::ReorderPairsOfVectors(reactant_nu_indices_, reactant_nu_);
			OpenSMOKE_Utilities::ReorderPairsOfVectors(product_nu_indices_, product_nu_);
		}

		// Clean indices (only after they are reordered)
		{
			OpenSMOKE_Utilities::CleanPairsOfVectors(reactant_nu_indices_, reactant_nu_);
			OpenSMOKE_Utilities::CleanPairsOfVectors(product_nu_indices_, product_nu_);
		}

		// Kinetic orders
		reactant_lambda_indices_ = reactant_nu_indices_;
		reactant_lambda_ = reactant_nu_;

		if (iReversible_ == true)
		{
			product_lambda_indices_  = product_nu_indices_;
			product_lambda_  = product_nu_;	
		}
		
		// Next lines
		{
			for (unsigned int l=1;l<lines.size();l++)
			{
				std::string line = lines[l];
				std::string keyword;
				OpenSMOKE_Utilities::LookForKeyWord(line, keyword);

				std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(keyword);
				if( it != map_of_species.end())
				{
					std::cout << "Third body efficiencies: The following species is not available: " << std::endl;
					return false;
				}
				else
				{
					if (boost::iequals(keyword, "COV"))
					{
						if (iLangmuir_ == true)
						{
							std::cout << "LANG and COV reactions are mutually exclusive!" << std::endl;
							return false;
						}
						
						if (iLumped_ == true)
						{
							std::cout << "LUMPED and COV reactions are mutually exclusive!" << std::endl;
							return false;
						}

						std::string species;
						std::vector<double> coefficients;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Coverage dependent reaction (COV). ", line, 4, species, coefficients);
						if (tag == false) return false;

						std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(species);
						if( it == map_of_species.end())
						{
							std::cout << "The following species is not available: " << species << std::endl;
							return false;
						}
						else
						{
							if ((*it).second < number_of_gas_species_)
							{
								std::cout << "COV species must bi site- or bulk- species: " << species << std::endl;
								return false;
							}
							else
							{
								if ( (*it).second < number_of_gas_species_+number_of_site_species_)
								{
									coverage_dependent_species_site_type_.push_back(true);
									coverage_dependent_species_index_.push_back((*it).second+1-number_of_gas_species_);
								}
								else
								{
									coverage_dependent_species_site_type_.push_back(false);
									coverage_dependent_species_index_.push_back((*it).second-number_of_gas_species_+1-number_of_site_species_);
								}

								coverage_dependent_eta_.push_back(coefficients[0]);
								coverage_dependent_mu_.push_back(coefficients[1]);
								coverage_dependent_epsilon_.push_back(coefficients[2]);
								iCoverageDependent_ = true;
							}
						}
					}

					else if (boost::iequals(keyword, "LANG"))
					{
						if (iStick_ == true)
						{
							std::cout << "LANG and STICK reactions are mutually exclusive!" << std::endl;
							return false;
						}

						if (iCoverageDependent_ == true)
						{
							std::cout << "LANG and COV reactions are mutually exclusive!" << std::endl;
							return false;
						}

						std::string species;
						std::vector<double> coefficients;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Langmuir-Hinshelwood reaction. ", line, 5, species, coefficients);
						if (tag == false) return false;

						std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(species);
						if( it == map_of_species.end())
						{
							std::cout << "The following species is not available: " << species << std::endl;
							return false;
						}
						else
						{
							if ((*it).second >= number_of_gas_species_)
							{
								std::cout << "LANG species must be gas-species: " << species << std::endl;
								return false;
							}
							else
							{
									langmuir_species_index_.push_back((*it).second+1);
									langmuir_A_.push_back( coefficients[0] );
									langmuir_Beta_.push_back(coefficients[1]);
									langmuir_H_over_R_.push_back(coefficients[2]/PhysicalConstants::R_cal_mol); // TODO
									langmuir_order_.push_back(coefficients[3]);
									langmuir_denominator_order_ = 2.0;
									langmuir_units_ = PhysicalConstants::UNITS_STD;
									iLangmuir_ = true;
							}
						}
					}
					
					else if (boost::iequals(keyword, "LUMPED"))
					{
						if (iStick_ == true)
						{
							std::cout << "LUMPED and STICK reactions are mutually exclusive!" << std::endl;
							return false;
						}

						if (iCoverageDependent_ == true)
						{
							std::cout << "LUMPED and COV reactions are mutually exclusive!" << std::endl;
							return false;
						}

						std::vector<std::string> values(1);
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusWords("Lumped reaction. ", line, 1, values);
						name_of_lumped_function_ = values[0];
						iLumped_ = true;
						if (tag == false) return false;
					}
					
					else if (boost::iequals(keyword, "LHDE"))
					{
						if (iLangmuir_ == false)
						{
							std::cout << "LHDE option must be preceeded by the LANG option!" << std::endl;
							return false;
						}
						
						if (iLumped_ == false)
						{
							std::cout << "LUMPED option cannot be used with the LHDE option!" << std::endl;
							return false;
						}

						std::vector<double> coefficients;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Langmuir-Hinshelwood denominator exponent parameter", line, 1, coefficients);
						langmuir_denominator_order_ = coefficients[0];
					}
					else if (boost::iequals(keyword, "LHNU"))
					{
						if (iLangmuir_ == false)
						{
							std::cout << "LHNU option must be preceeded by the LANG option!" << std::endl;
							return false;
						}
						
						if (iLumped_ == false)
						{
							std::cout << "LUMPED option cannot be used with the LHNU option!" << std::endl;
							return false;
						}						

						std::vector<std::string> words;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusWords("Langmuir-Hinshelwood denominator exponent parameter", line, 1, words);
						
						std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(words[0]);
						if( it == map_of_species.end())
						{
							std::cout << "The following species is not available: " << words[0] << std::endl;
							return false;
						}
						else
						{
							if ((*it).second >= number_of_gas_species_)
							{
								std::cout << "LHNU species must be gas-species: " << words[0] << std::endl;
								return false;
							}
							else
							{
								langmuir_numerator_species_index_provisional_.push_back((*it).second+1);
							}
						}
					}
					
					else if (boost::iequals(keyword, "LHPR"))
					{
						if (iLangmuir_ == false)
						{
							std::cout << "LHPR option must be preceeded by the LANG option!" << std::endl;
							return false;
						}
						
						if (iLumped_ == false)
						{
							std::cout << "LUMPED option cannot be used with the LHPR option!" << std::endl;
							return false;
						}

						std::vector<std::string> words;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusWords("Langmuir-Hinshelwood units", line, 1, words);
						
						if (words[0] != "atm" && words[0] != "ATM" &&
							words[0] != "bar" && words[0] != "BAR" &&
							words[0] != "torr" && words[0] != "TORR" &&
							words[0] != "pasc" && words[0] != "PASC" &&
							words[0] != "dyne" && words[0] != "DYNE"  )
						{
							std::cout << "LHPR option must be equal to one of the following: atm || bar || torr || pasc || dyne!" << std::endl;
							return false;
						}

						if (words[0] == "atm" || words[0] == "ATM")			langmuir_units_ = PhysicalConstants::UNITS_ATM;
						else if (words[0] == "bar" || words[0] == "BAR")	langmuir_units_ = PhysicalConstants::UNITS_BAR;
						else if (words[0] == "torr" || words[0] == "TORR")	langmuir_units_ = PhysicalConstants::UNITS_TORR;
						else if (words[0] == "pasc" || words[0] == "PASC")	langmuir_units_ = PhysicalConstants::UNITS_PASCALS;
						else if (words[0] == "dyne" || words[0] == "DYNE")	langmuir_units_ = PhysicalConstants::UNITS_DYNES;					
					}

					else if (boost::iequals(keyword, "STICK"))
					{
						if (iLangmuir_ == true)
						{
							std::cout << "LANG and STICK reactions are mutually exclusive!" << std::endl;
							return false;
						}

						iStick_ = true;
					}
					else if (boost::iequals(keyword, "DUP") || boost::iequals(keyword, "DUPLICATE"))
					{
						if (iDuplicate_ == true)
						{
							std::cout << "The DUPLICATE (or DUP) option is used more than once!" << std::endl;
							return false;
						}

						// Check if additional options are specified
						{
							typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
							boost::char_separator<char> sep_blank(" ");
							tokenizer_blank tokens(line, sep_blank);
							const std::size_t n = std::distance(tokens.begin(), tokens.end());
							if (n != 1)
							{
								std::cout << "The DUPLICATE (or DUP) keyword is used together with an additional option on the same line. Please, split the line." << std::endl;
								return false;
							}
						}
						iDuplicate_ = true;
					}
					else if (boost::iequals(keyword, "FORD"))
					{
						std::string species;
						std::vector<double> coefficients;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Forward Reaction Order Parameter (FORD). ", line, 2, species, coefficients);
						if (tag == false) return false;
						else
						{
							iFord_ = true;
							std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(species);
							if( it == map_of_species.end())
							{
								std::cout << "The following species is not available: " << species << std::endl;
								return false;
							}
							else
							{
								bool iFound = false;
								for(unsigned int i=0;i<reactant_lambda_.size();i++)
								{
									if (reactant_lambda_indices_[i] == (*it).second)
									{
										reactant_lambda_[i] = coefficients[0];
										iFound = true;
										break;
									}
								}
								if (iFound==false)
								{
									reactant_lambda_indices_.push_back((*it).second);
									reactant_lambda_.push_back(coefficients[0]);
								}
							}
						}			
					}
					else if (boost::iequals(keyword, "MWON"))
					{
						iMotzWiseCorrection_ = true;
					}
					else if (boost::iequals(keyword, "MWOFF"))
					{
						iMotzWiseCorrection_ = false;
					}
					else if (boost::iequals(keyword, "REV"))
					{
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}

						if (iExplicitlyReversible_ == true)
						{
							std::cout << "The REV option is used more than once!" << std::endl;
							return false;
						}

						std::vector<double> coefficients;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Explicit reversible kinetics (REV). ", line, 3, coefficients);
						if (tag == false) return false;
						if (iReversible_ == false)
						{
							std::cout << "The REV tag is used, but the reaction is not reversible!" << std::endl;
							return false;
						}
						ARev_ = coefficients[0];
						betaRev_ = coefficients[1];
						ERev_ = coefficients[2];
						iExplicitlyReversible_ = true;
					}
					else if (boost::iequals(keyword, "RORD"))
					{
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}

						if (iReversible_ == false)
						{
							std::cout << "Backward Reaction Order Parameter (RORD). " << std::endl;
							std::cout << "The reaction is not reversible." << std::endl;
							return false;
						}

						std::string species;
						std::vector<double> coefficients;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Backward Reaction Order Parameter (FORD). ", line, 2, species, coefficients);
						if (tag == false) return false;
						else
						{
							iRord_ = true;
							std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(species);
							if( it == map_of_species.end())
							{
								std::cout << "Backward Reaction Order Parameter (RORD). " << std::endl;
								std::cout << "The following species is not available: " << species << std::endl;
								return false;
							}
							else
							{
								bool iFound = false;
								for(unsigned int i=0;i<product_lambda_.size();i++)
								{
									if (product_lambda_indices_[i] == (*it).second)
									{
										product_lambda_[i] = coefficients[0];
										iFound = true;
										break;
									}
								}
								if (iFound==false)
								{
									product_lambda_indices_.push_back((*it).second);
									product_lambda_.push_back(coefficients[0]);
								}
							}
						}
					}
					else if (boost::iequals(keyword, "UNITS"))
					{
						std::vector<std::string> words;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusWords("Units. ", line, 1, words);
						if (tag == false) return false;
						else
						{
							if (words[0] == "CAL" || words[0] == "CAL/MOLE" || words[0] == "cal" || words[0] == "cal/mole")
								e_units_ = PhysicalConstants::UNITS_CAL_MOLE;
							else if (words[0] == "EVOL" || words[0] == "EVOLTS" || words[0] == "evol" || words[0] == "evolts")
								e_units_ = PhysicalConstants::UNITS_EVOLTS;
							else if (words[0] == "JOUL" || words[0] == "JOULES/MOLE" || words[0] == "joul" || words[0] == "joules/mole")
								e_units_ = PhysicalConstants::UNITS_JOULES_MOLE;
							else if (words[0] == "KJOU" || words[0] == "KJOULES/MOLE" || words[0] == "kjou" || words[0] == "kjoules/mole")
								e_units_ = PhysicalConstants::UNITS_KJOULES_MOLE;
							else if (words[0] == "KCAL" || words[0] == "KCAL/MOLE" || words[0] == "kcal" || words[0] == "kcal/mole")
								e_units_ = PhysicalConstants::UNITS_KCAL_MOLE;
							else if (words[0] == "KELV" || words[0] == "KELVINS" || words[0] == "kelv" || words[0] == "kelvins")
								e_units_ = PhysicalConstants::UNITS_KELVINS;
							
							else if (words[0] == "MOLEC" || words[0] == "MOLECULES" || words[0] == "molec" || words[0] == "molecules")
								a_units_ = PhysicalConstants::UNITS_MOLECULES;
							else if (words[0] == "MOLE" || words[0] == "MOLES" || words[0] == "mole" || words[0] == "moles")
								a_units_ = PhysicalConstants::UNITS_MOLES;


							else if (words[0] == "ATM" || words[0] == "atm")
								composition_units_ = PhysicalConstants::UNITS_ATM;
							else if (words[0] == "BAR" || words[0] == "bar")
								composition_units_ = PhysicalConstants::UNITS_BAR;
							else if (words[0] == "DYN" || words[0] == "dyn" || words[0] == "DYNES" || words[0] == "dynes")
								composition_units_ = PhysicalConstants::UNITS_DYNES;
							else if (words[0] == "PAS" || words[0] == "pas" || words[0] == "PASCALS" || words[0] == "pascals")
								composition_units_ = PhysicalConstants::UNITS_PASCALS;
							else if (words[0] == "TOR" || words[0] == "tor" || words[0] == "TORR" || words[0] == "torr")
								composition_units_ = PhysicalConstants::UNITS_TORR;
							
							else
							{
								std::cout << "Wrong units" << std::endl;
								std::cout << "Available units: CAL || CAL/MOLE || EVOL || EVOLTS || JOUL || JOULES/MOLE || KCAL || KCAL/MOLE || KELV || KELVINS || MOLEC || MOLECULES || MOLE || MOLES || ATM || BAR || DYN || DYNES || PAS || PASCALS || TOR || TORR" << std::endl;
								return false;
							}
						}
					}
					else if (boost::iequals(keyword, "USRPROG"))
					{
						std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
						return false;
					}

					else if (boost::iequals(keyword, "UBIQEP"))
					{
						iUBIQEP_ = true;

						if (iStick_ == true)
						{
							std::cout << "UBIQEP and STICK reactions are mutually exclusive!" << std::endl;
							return false;
						}

						if (iCoverageDependent_ == true)
						{
							std::cout << "UBIQEP and COV reactions are mutually exclusive!" << std::endl;
							return false;
						}

						if (iLangmuir_ == true)
						{
							std::cout << "UBIQEP and LANG reactions are mutually exclusive!" << std::endl;
							return false;
						}
						
						if (iLumped_ == true)
						{
							std::cout << "UBIQEP and LUMPED reactions are mutually exclusive!" << std::endl;
							return false;
						}						

						std::vector<std::string> coefficients;
						bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusWords("UBI-QEP reaction. ", line, coefficients);
						if (tag == false) return false;

						ubiqep_reaction_class_ = boost::lexical_cast<unsigned int>(coefficients[0]);
						if (ubiqep_reaction_class_ < 1 || ubiqep_reaction_class_ >8)
						{
							std::cout << "The type of UBI-QEP reaction must be a integer value between 1 and 8" << std::endl;
							return false;
						}

						if (coefficients[1] == "ADS" || coefficients[1] == "ads")
							ubiqep_reaction_type_ = PhysicalConstants::UBIQEP_TYPE_ADSORPTION;
						else if (coefficients[1] == "DES" || coefficients[1] == "des")
							ubiqep_reaction_type_ = PhysicalConstants::UBIQEP_TYPE_DESORPTION;
						else if (coefficients[1] == "SUP" || coefficients[1] == "sup")
							ubiqep_reaction_type_ = PhysicalConstants::UBIQEP_TYPE_SURFACE;
						else
						{
							std::cout << "Only the following types are allowed for UBI-QEP reactions: ADS || DES || SUP" << std::endl;
							return false;
						}

						if (ubiqep_reaction_type_ == PhysicalConstants::UBIQEP_TYPE_ADSORPTION && ubiqep_reaction_class_>4)
						{
							std::cout << "The ADS type is allowed only for UBI-QEP reactions of class 1, 2, 3 or 4" << std::endl;
							return false;
						}

						if (coefficients[2] == "DIR" || coefficients[2] == "dir")
							ubiqep_direct_ = true;
						else if (coefficients[2] == "REV" || coefficients[2] == "rev")
							ubiqep_direct_ = false;
						else
						{
							std::cout << "Only the following types are allowed for UBI-QEP reactions: DIR || REV" << std::endl;
							return false;
						}

						const std::size_t n = coefficients.size() - 3;

						if (ubiqep_reaction_class_ == 1 && n != 0) 
						{
							std::cout << "UBI-QEP reactions of type 1 require the user specifies 0 gas phase species in the option line" << std::endl;
							return false;
						}
						
						if (ubiqep_reaction_class_ == 2 && n != 1) 
						{
							std::cout << "UBI-QEP reactions of type 2 require the user specifies 1 gas phase species in the option line" << std::endl;
							return false;
						}
						if (ubiqep_reaction_class_ == 2 && n == 1) 
						{
							std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[3]);
							if( it == map_of_species.end())
							{
								std::cout << "The following species is not available: " << coefficients[3] << std::endl;
								return false;
							}
							else
							{
								ubiqep_index_A_ = (*it).second + 1;
							}
						}
						if (ubiqep_reaction_class_ == 3 && n != 1) 
						{
							std::cout << "UBI-QEP reactions of type 3 require the user specifies 1 gas phase species in the option line" << std::endl;
							return false;
						}
						if (ubiqep_reaction_class_ == 3 && n == 1) 
						{
							std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[3]);
							if( it == map_of_species.end())
							{
								std::cout << "The following species is not available: " << coefficients[3] << std::endl;
								return false;
							}
							else
							{
								ubiqep_index_A_ = (*it).second + 1;
							}
						}
						if (ubiqep_reaction_class_ == 4 && n != 2) 
						{
							std::cout << "UBI-QEP reactions of type 4 require the user specifies 2 gas phase species in the option line" << std::endl;
							return false;
						}
						if (ubiqep_reaction_class_ == 4 && n == 2) 
						{
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[3]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[3] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_A_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[4]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[4] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_B_ = (*it).second + 1;
								}
							}
						}
						if (ubiqep_reaction_class_ == 5 && n != 3) 
						{
							std::cout << "UBI-QEP reactions of type 5 require the user specifies 3 gas phase species in the option line" << std::endl;
							return false;
						}
						if (ubiqep_reaction_class_ == 5 && n == 3) 
						{
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[3]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[3] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_A_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[4]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[4] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_B_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[5]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[5] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_AB_ = (*it).second + 1;
								}
							}
						}
						if (ubiqep_reaction_class_ == 6 && n != 4) 
						{
							std::cout << "UBI-QEP reactions of type 6 require the user specifies 4 gas phase species in the option line" << std::endl;
							return false;
						}
						if (ubiqep_reaction_class_ == 6 && n == 4) 
						{
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[3]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[3] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_A_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[4]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[4] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_B_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[5]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[5] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_C_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[6]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[6] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_D_ = (*it).second + 1;
								}
							}
						}
						if (ubiqep_reaction_class_ == 7 && n != 3) 
						{
							std::cout << "UBI-QEP reactions of type 7 require the user specifies 3 gas phase species in the option line" << std::endl;
							return false;
						}
						if (ubiqep_reaction_class_ == 7 && n == 3) 
						{
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[3]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[3] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_A_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[4]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[4] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_B_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[5]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[5] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_C_ = (*it).second + 1;
								}
							}
						}
						if (ubiqep_reaction_class_ == 8 && n != 3) 
						{
							std::cout << "UBI-QEP reactions of type 8 require the user specifies 3 gas phase species in the option line" << std::endl;
							return false;
						}
						if (ubiqep_reaction_class_ == 8 && n == 3) 
						{
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[3]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[3] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_A_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[4]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[4] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_C_ = (*it).second + 1;
								}
							}
							{
								std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(coefficients[5]);
								if( it == map_of_species.end())
								{
									std::cout << "The following species is not available: " << coefficients[5] << std::endl;
									return false;
								}
								else
								{
									ubiqep_index_D_ = (*it).second + 1;
								}
							}
						}

						ubiqep_A_			= A_;
						ubiqep_Beta_		= beta_;
						ubiqep_BondIndex_	= E_;
					}
					else
					{
						std::cout << "Sorry! The " << keyword << " option is not available." << std::endl;
						return false;
					}
				}
			}
		}

		if (iReversible_ == true && iExplicitlyReversible_ == false)
		{
			if (iUBIQEP_ == true)
			{
				std::cout << "Sorry! Reversible reactions are allowed only if the kinetic parameters of the inverse reaction are explicitly provided through the REV keyword!" << std::endl;
				return false;
			}
		}

		if (iReversible_ == true)
		{
			if (iStick_ == true)
			{
				std::cout << "The STICK option can be applied only to irreversible reactions" << std::endl;
				return false;
			}
			if (iCoverageDependent_ == true)
			{
				std::cout << "The COV option can be applied only to irreversible reactions" << std::endl;
				return false;
			}
			if (iLangmuir_ == true)
			{
				std::cout << "The LANG option can be applied only to irreversible reactions" << std::endl;
				return false;
			}
		}

		// Reorder and clean indices
		if (iFord_ == true)
		{
			OpenSMOKE_Utilities::ReorderPairsOfVectors(reactant_lambda_indices_, reactant_lambda_);
			OpenSMOKE_Utilities::CleanPairsOfVectors(reactant_lambda_indices_, reactant_lambda_);
			
		}
		if (iRord_ == true)
		{
			OpenSMOKE_Utilities::ReorderPairsOfVectors(product_lambda_indices_, product_lambda_);
			OpenSMOKE_Utilities::CleanPairsOfVectors(product_lambda_indices_, product_lambda_);
		}

		if (iStick_ == true)
		{
			double stick_stoichiometric_coefficient = 0.;
			unsigned int number_of_gas_species = 0;
			for(unsigned int i=0;i<reactant_lambda_.size();i++)
				if (reactant_lambda_indices_[i] < number_of_gas_species_)	
				{
					number_of_gas_species++;
					stick_gas_species_ = reactant_lambda_indices_[i] + 1;
					stick_stoichiometric_coefficient = reactant_nu_[i];
				}
			if (number_of_gas_species != 1)
			{
				std::cout << "Sticking coefficients can be adopted only for reactions involving exactly one species in gas phase!" << std::endl;
				return false;
			}
			if (stick_stoichiometric_coefficient != 1)
			{
				std::cout << "The stoichiometric coefficient of the gas-phase species in a STICK reaction must be equal to 1!" << std::endl;
				return false;
			}

			// Sum of all the stoichiometric coefficients of reactants that are surface species
			stick_power_ = 0.;
			for(unsigned int i=0;i<reactant_lambda_.size();i++)
			{
				if (reactant_lambda_indices_[i] < number_of_gas_species_)	continue;
				else if (reactant_lambda_indices_[i] < number_of_gas_species_+number_of_site_species_)
				{
					stick_power_+=reactant_lambda_[i];
					stick_indices_site_species_.push_back(reactant_lambda_indices_[i]+1-number_of_gas_species_);
					stick_exponents_site_species_.push_back(reactant_lambda_[i]);
				}
			}
		}

		if (iLangmuir_ == true)
		{
			// Checking units
			if (composition_units_ != langmuir_units_)
			{
				if (composition_units_ != PhysicalConstants::UNITS_STD)
				{
					std::cout << "LHPR units do not match the current units for composition." << std::endl;
					std::cout << "Please use the UNITS option to force the agreement between the two." << std::endl;
					return false;
				}
					
				composition_units_ = langmuir_units_;
			}

			// Check if a species is specified two or more times
			{
				for(unsigned int j=0;j<langmuir_species_index_.size();j++)
					for(unsigned int k=0;k<langmuir_species_index_.size();k++)
						if (langmuir_species_index_[j] == langmuir_species_index_[k] && j!=k)
						{
							std::cout << "The LANG option was applied two times to the same species: " << std::endl;
							return false;
						}
			}

			// Finalising setup...
			langmuir_numerator_species_.resize(langmuir_species_index_.size());
			for(unsigned int j=0;j<langmuir_numerator_species_.size();j++)
				langmuir_numerator_species_[j] = false;

			for(unsigned int j=0;j<langmuir_species_index_.size();j++)
				if (langmuir_species_index_[j] >= number_of_gas_species_)	
				{
					std::cout << "LANG species must be gas-species " << std::endl;
					return false;
				}
				
			for(unsigned int j=0;j<langmuir_numerator_species_index_provisional_.size();j++)
			{
				bool iFound = false;
				for(unsigned int k=0;k<langmuir_species_index_.size();k++)
					if (langmuir_species_index_[k] == langmuir_numerator_species_index_provisional_[j])
					{
						iFound = true;
						langmuir_numerator_species_[k] = true;
						break;
					}
				
					if (iFound == false)
					{
						std::cout << "The LHNU option can be used only if the LANG option was previously used for the same species!" << std::endl;
						return false;
					}
			}

			// Conversions
			if (langmuir_units_ == PhysicalConstants::UNITS_STD)
			{
				for(unsigned int k=0;k<langmuir_A_.size();k++)
					langmuir_A_[k] /= std::pow(1.e3, langmuir_order_[k]);
			}
			else
			{
				double conversion_factor = 1.;

				if (composition_units_ == PhysicalConstants::UNITS_ATM)				conversion_factor = 101325.;
				else if (composition_units_ == PhysicalConstants::UNITS_BAR)		conversion_factor = 100000.;
				else if (composition_units_ == PhysicalConstants::UNITS_PASCALS)	conversion_factor = 1.;
				else if (composition_units_ == PhysicalConstants::UNITS_TORR)		conversion_factor = 101325./760.;
				else if (composition_units_ == PhysicalConstants::UNITS_DYNES)		conversion_factor = 0.1;

				for(unsigned int k=0;k<langmuir_A_.size();k++)
					langmuir_A_[k] /= std::pow(conversion_factor, langmuir_order_[k]);
			}
		}

		// Conservation of sites
		{
			double number_of_reactant_sites = 0.;
			for(unsigned int i=0;i<reactant_nu_.size();i++)
				if (reactant_nu_indices_[i] >= number_of_gas_species_ && reactant_nu_indices_[i]<number_of_gas_species_+number_of_site_species_)	
				{
					unsigned int index = reactant_nu_indices_[i] - number_of_gas_species_;
					number_of_reactant_sites += occupancy_site_species[index]*reactant_nu_[i];
				}

			double number_of_product_sites = 0.;
			for(unsigned int i=0;i<product_nu_.size();i++)
				if (product_nu_indices_[i] >= number_of_gas_species_ && product_nu_indices_[i]<number_of_gas_species_+number_of_site_species_)	
				{
					unsigned int index = product_nu_indices_[i] - number_of_gas_species_;
					number_of_product_sites += occupancy_site_species[index]*product_nu_[i];
				}

			delta_occupancy_sites_ = number_of_product_sites - number_of_reactant_sites;

			
			if (std::fabs(delta_occupancy_sites_) > 1.e-16)
			{
				if (global_non_conservation_of_sites_ == false)
				{
					std::cout << "The reaction does not conserve the number of surface species!" << std::endl;
					std::cout << "The non conservation of surface species can be enabled using the NONCON keyword." << std::endl;
					return false;
				}
				else
				{
					std::cout << "The reaction does not conserve the number of surface species!" << std::endl;
					std::cout << "This was explicitly required by the user through the NONCON keyword." << std::endl;
				}
			}
		}

		// Consistency of phases
		{
			std::vector<unsigned int> memberships;
			for(unsigned int i=0;i<reactant_nu_.size();i++)
				if (reactant_nu_indices_[i] >= number_of_gas_species_ && reactant_nu_indices_[i]<number_of_gas_species_+number_of_site_species_)	
				{
					unsigned int index = reactant_nu_indices_[i] - number_of_gas_species_;
					memberships.push_back(site_phase_membership[index]);
				}
			for(unsigned int i=0;i<product_nu_.size();i++)
				if (product_nu_indices_[i] >= number_of_gas_species_ && product_nu_indices_[i]<number_of_gas_species_+number_of_site_species_)	
				{
					unsigned int index = product_nu_indices_[i] - number_of_gas_species_;
					memberships.push_back(site_phase_membership[index]);
				}

			for(unsigned int i=0;i<memberships.size();i++)
				for(unsigned int j=i+1;j<memberships.size();j++)
					if (memberships[i] != memberships[j])
					{
						std::cout << "The reaction mixes surface species belonging to different phases!" << std::endl;
						return false;
					}

			if (memberships.size() != 0)
				surface_reaction_membership_ = memberships[0];
		}

		if (iUBIQEP_ == true)
		{
			if (ubiqep_direct_ == true)
			{
				double number_of_gas_species_reactant_side	= 0;
				double number_of_site_species_reactant_side	= 0;
				double number_of_bulk_species_reactant_side	= 0;
				double number_of_gas_species_product_side	= 0;
				double number_of_site_species_product_side	= 0;
				double number_of_bulk_species_product_side	= 0;

				for(unsigned int i=0;i<original_reactant_nu_.size();i++)
				{
					if (original_reactant_nu_indices_[i] < number_of_gas_species_)	
						number_of_gas_species_reactant_side += original_reactant_nu_[i];
					else if (original_reactant_nu_indices_[i] < number_of_gas_species_+number_of_site_species_)
						number_of_site_species_reactant_side += original_reactant_nu_[i];
					else number_of_bulk_species_reactant_side += original_reactant_nu_[i];;
				}

				for(unsigned int i=0;i<original_product_nu_.size();i++)
				{
					if (original_product_nu_indices_[i] < number_of_gas_species_)	
						number_of_gas_species_product_side += original_product_nu_[i];
					else if (original_product_nu_indices_[i] < number_of_gas_species_+number_of_site_species_)
						number_of_site_species_product_side += original_product_nu_[i];
					else number_of_bulk_species_product_side += original_product_nu_[i];
				}

				if (number_of_bulk_species_reactant_side !=0 || number_of_bulk_species_product_side != 0)
				{
					std::cout << "UBI-QEP reactions must involve only gas-phase and surface-phase species!" << std::endl;
					return false;
				}
				

				// 1. Non-activated atomic or non-dissociative molecular adsorption (A + * = A*)
				if (ubiqep_reaction_class_ == 1)
				{
					std::string message = "UBI-QEP reaction of type 1 must have the following form: [ A + * = A* ]";
					if (number_of_gas_species_reactant_side  != 1)	return FatalErrorMessage(message);
					if (number_of_gas_species_product_side   != 0)	return FatalErrorMessage(message);
					if (number_of_site_species_reactant_side != 1)	return FatalErrorMessage(message);
					if (number_of_site_species_product_side  != 1)	return FatalErrorMessage(message);

					// Checking reactant side
					if (original_reactant_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0]>=number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[1]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[1] != 1.)								return FatalErrorMessage(message);

					// Checking product side
					if (original_product_nu_.size() != 1)							return FatalErrorMessage(message);
					if (original_product_nu_indices_[0]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[0] != 1.)								return FatalErrorMessage(message);
					

					ubiqep_index_A_     = original_reactant_nu_indices_[0] +1;
					ubiqep_index_Star_  = original_reactant_nu_indices_[1] -number_of_gas_species_+1;
					ubiqep_index_AStar_ = original_product_nu_indices_[0]  -number_of_gas_species_+1;
				}

				// 2. Non-activated homonuclear dissociative adsorption (A2 + 2* = 2A*)
				// 3. Activated homonuclear dissociative adsorption (A2 + 2* = 2A*)
				if (ubiqep_reaction_class_ == 2 || ubiqep_reaction_class_ == 3)
				{
					std::string message = "UBI-QEP reaction of types 2 and 3 must have the following form: [ A2 + 2* = 2A* ]";
					if (number_of_gas_species_reactant_side  != 1)	return FatalErrorMessage(message);
					if (number_of_gas_species_product_side   != 0)	return FatalErrorMessage(message);
					if (number_of_site_species_reactant_side != 2)	return FatalErrorMessage(message);
					if (number_of_site_species_product_side  != 2)	return FatalErrorMessage(message);

					// Checking reactant side
					if (original_reactant_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0]>=number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[1]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[1] != 2.)								return FatalErrorMessage(message);

					// Checking product side
					if (original_product_nu_.size() != 1)							return FatalErrorMessage(message);
					if (original_product_nu_indices_[0]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[0] != 2.)								return FatalErrorMessage(message);

					ubiqep_index_A2_	= original_reactant_nu_indices_[0]+1;
					ubiqep_index_Star_	= original_reactant_nu_indices_[1]-number_of_gas_species_+1;
					ubiqep_index_AStar_ = original_product_nu_indices_[0]-number_of_gas_species_+1;
				}

				// 4. Activated homonuclear dissociative adsorption (AB + 2* = A* + B*)
				if (ubiqep_reaction_class_ == 4)
				{
					std::string message = "UBI-QEP reaction of type 4 must have the following form: [ AB + 2* = A* + B* ]";
					if (number_of_gas_species_reactant_side  != 1)	return FatalErrorMessage(message);
					if (number_of_gas_species_product_side   != 0)	return FatalErrorMessage(message);
					if (number_of_site_species_reactant_side != 2)	return FatalErrorMessage(message);
					if (number_of_site_species_product_side  != 2)	return FatalErrorMessage(message);

					// Checking reactant side
					if (original_reactant_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0]>=number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[1]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[1] != 2.)								return FatalErrorMessage(message);

					// Checking product side
					if (original_product_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_product_nu_indices_[0]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_product_nu_indices_[1]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[1] != 1.)								return FatalErrorMessage(message);
					if (original_product_nu_indices_[0] == 
						original_product_nu_indices_[1] )							return FatalErrorMessage(message);

					ubiqep_index_AB_	= original_reactant_nu_indices_[0]+1;
					ubiqep_index_Star_	= original_reactant_nu_indices_[1]-number_of_gas_species_+1;
					ubiqep_index_AStar_ = original_product_nu_indices_[0]-number_of_gas_species_+1;
					ubiqep_index_BStar_ = original_product_nu_indices_[1]-number_of_gas_species_+1;
				}

				// 5. Heteronuclear surface dissociation (AB + * = A* + B*)
				if (ubiqep_reaction_class_ == 5)
				{
					std::string message = "UBI-QEP reaction of type 5 must have the following form: [ AB* + * = A* + B* ]";
					if (number_of_gas_species_reactant_side  != 0)	return FatalErrorMessage(message);
					if (number_of_gas_species_product_side   != 0)	return FatalErrorMessage(message);
					if (number_of_site_species_reactant_side != 2)	return FatalErrorMessage(message);
					if (number_of_site_species_product_side  != 2)	return FatalErrorMessage(message);

					// Checking reactant side
					if (original_reactant_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[1]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[1] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0] == 
						original_reactant_nu_indices_[1] )							return FatalErrorMessage(message);

					// Checking product side
					if (original_product_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_product_nu_indices_[0]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_product_nu_indices_[1]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[1] != 1.)								return FatalErrorMessage(message);
					if (original_product_nu_indices_[0] == 
						original_product_nu_indices_[1] )							return FatalErrorMessage(message);
						
					ubiqep_index_ABStar_ = original_reactant_nu_indices_[0]-number_of_gas_species_+1;
					ubiqep_index_Star_   = original_reactant_nu_indices_[1]-number_of_gas_species_+1;

					ubiqep_index_AStar_ = original_product_nu_indices_[0]-number_of_gas_species_+1;
					ubiqep_index_BStar_ = original_product_nu_indices_[1]-number_of_gas_species_+1;
				}

				// 6. Surface disproportion (A* + B* = C* + D*)
				if (ubiqep_reaction_class_ == 6)
				{
					std::string message = "UBI-QEP reaction of type 6 must have the following form: [ A* + B* = C* + D* ]";
					if (number_of_gas_species_reactant_side  != 0)	return FatalErrorMessage(message);
					if (number_of_gas_species_product_side   != 0)	return FatalErrorMessage(message);
					if (number_of_site_species_reactant_side != 2)	return FatalErrorMessage(message);
					if (number_of_site_species_product_side  != 2)	return FatalErrorMessage(message);

					// Checking reactant side
					if (original_reactant_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[1]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[1] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0] == 
						original_reactant_nu_indices_[1] )							return FatalErrorMessage(message);

					// Checking product side
					if (original_product_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_product_nu_indices_[0]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_product_nu_indices_[1]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[1] != 1.)								return FatalErrorMessage(message);
					if (original_product_nu_indices_[0] == 
						original_product_nu_indices_[1] )							return FatalErrorMessage(message);

					ubiqep_index_AStar_ = original_reactant_nu_indices_[0]-number_of_gas_species_+1;
					ubiqep_index_BStar_   = original_reactant_nu_indices_[1]-number_of_gas_species_+1;

					ubiqep_index_CStar_ = original_product_nu_indices_[0]-number_of_gas_species_+1;
					ubiqep_index_DStar_ = original_product_nu_indices_[1]-number_of_gas_species_+1;
				}

				// 7. Homonuclear surface disproportion (2C* = A* + B*)
				if (ubiqep_reaction_class_ == 7)
				{
					std::string message = "UBI-QEP reaction of type 7 must have the following form: [ 2C* = A* + B* ]";
					if (number_of_gas_species_reactant_side  != 0)	return FatalErrorMessage(message);
					if (number_of_gas_species_product_side   != 0)	return FatalErrorMessage(message);
					if (number_of_site_species_reactant_side != 2)	return FatalErrorMessage(message);
					if (number_of_site_species_product_side  != 2)	return FatalErrorMessage(message);

					// Checking reactant side
					if (original_reactant_nu_.size() != 1)							return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[0] != 2.)								return FatalErrorMessage(message);

					// Checking product side
					if (original_product_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_product_nu_indices_[0]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_product_nu_indices_[1]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[1] != 1.)								return FatalErrorMessage(message);
					if (original_product_nu_indices_[0] == 
						original_product_nu_indices_[1] )							return FatalErrorMessage(message);

					ubiqep_index_CStar_ = original_reactant_nu_indices_[0]-number_of_gas_species_+1;

					ubiqep_index_AStar_ = original_product_nu_indices_[0]-number_of_gas_species_+1;
					ubiqep_index_BStar_ = original_product_nu_indices_[1]-number_of_gas_species_+1;
				}

				// 8. Surface disproportionation to homonuclear product (C* + D* = 2A)
				if (ubiqep_reaction_class_ == 8)
				{
					std::string message = "UBI-QEP reaction of type 8 must have the following form: [ C* + D* = 2A* ]";
					if (number_of_gas_species_reactant_side  != 0)	return FatalErrorMessage(message);
					if (number_of_gas_species_product_side   != 0)	return FatalErrorMessage(message);
					if (number_of_site_species_reactant_side != 2)	return FatalErrorMessage(message);
					if (number_of_site_species_product_side  != 2)	return FatalErrorMessage(message);

					// Checking reactant side
					if (original_reactant_nu_.size() != 2)							return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[0] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[1]<number_of_gas_species_)	return FatalErrorMessage(message);
					if (original_reactant_nu_[1] != 1.)								return FatalErrorMessage(message);
					if (original_reactant_nu_indices_[0] == 
						original_reactant_nu_indices_[1] )							return FatalErrorMessage(message);

					// Checking product side
					if (original_product_nu_.size() != 1)							return FatalErrorMessage(message);
					if (original_product_nu_indices_[0]<number_of_gas_species_)		return FatalErrorMessage(message);
					if (original_product_nu_[0] != 2.)								return FatalErrorMessage(message);

					ubiqep_index_CStar_ = original_reactant_nu_indices_[0]-number_of_gas_species_+1;
					ubiqep_index_DStar_ = original_reactant_nu_indices_[1]-number_of_gas_species_+1;
					ubiqep_index_AStar_ = original_product_nu_indices_[0]-number_of_gas_species_+1;
				}
			}
		}

		if (iReversible_ == true && iExplicitlyReversible_ == false)
		{
			double number_of_gas_species_reactant_side	= 0;
			double number_of_gas_species_product_side	= 0;

			for(unsigned int i=0;i<original_reactant_nu_.size();i++)
			{
				if (original_reactant_nu_indices_[i] < number_of_gas_species_)	
					number_of_gas_species_reactant_side += original_reactant_nu_[i];
			}

			for(unsigned int i=0;i<original_product_nu_.size();i++)
			{
				if (original_product_nu_indices_[i] < number_of_gas_species_)	
					number_of_gas_species_product_side += original_product_nu_[i];
			}

			delta_nu_gas_ = number_of_gas_species_product_side - number_of_gas_species_reactant_side;
		}

		tag_reaction_ = PhysicalConstants::REACTION_SURFACE_SIMPLE;
		if (iStick_ == true && iCoverageDependent_ == false) 
			tag_reaction_ = PhysicalConstants::REACTION_SURFACE_STICK;
		if (iStick_ == false && iCoverageDependent_ == true) 
			tag_reaction_ = PhysicalConstants::REACTION_SURFACE_COVERAGE_DEPENDENT;
		if (iStick_ == true && iCoverageDependent_ == true) 
			tag_reaction_ = PhysicalConstants::REACTION_SURFACE_STICK_COVERAGE_DEPENDENT;
		if (iLangmuir_ == true) 
			tag_reaction_ = PhysicalConstants::REACTION_SURFACE_LANGMUIR;
		if (iLumped_ == true) 
			tag_reaction_ = PhysicalConstants::REACTION_SURFACE_LUMPED;

		ConvertUnits();

		return true;
	}

	void ReactionPolicy_Surface_CHEMKIN::ReactionOrders() 
	{

		sumLambdaGasReactants_  = 0.;
		sumLambdaSiteReactants_ = 0.;
		sumLambdaBulkReactants_ = 0.;
		for(unsigned int i=0;i<reactant_lambda_.size();i++)
		{
			if (reactant_lambda_indices_[i] < number_of_gas_species_)								sumLambdaGasReactants_ += reactant_lambda_[i];
			else if (reactant_lambda_indices_[i] < number_of_gas_species_+number_of_site_species_)	sumLambdaSiteReactants_ += reactant_lambda_[i];
			else																					sumLambdaBulkReactants_ += reactant_lambda_[i];
		}

		sumLambdaGasProducts_  = 0.;
		sumLambdaSiteProducts_ = 0.;
		sumLambdaBulkProducts_ = 0.;
		for(unsigned int i=0;i<product_lambda_.size();i++)
		{
			if (product_lambda_indices_[i] < number_of_gas_species_)								sumLambdaGasProducts_  += product_lambda_[i];
			else if (product_lambda_indices_[i] < number_of_gas_species_+number_of_site_species_)	sumLambdaSiteProducts_ += product_lambda_[i];
			else 																					sumLambdaBulkProducts_ += product_lambda_[i];
		}
	}

	void ReactionPolicy_Surface_CHEMKIN::ConvertUnits() 
	{
		ConversionFactors();
		ReactionOrders();

		// Conversion of activation energy
		{
			E_ *= e_conversion_;																	// J/kmol
			
			if (iExplicitlyReversible_ == true)
				ERev_ *= e_conversion_;																// J/kmol
		}

		// Conversion of pre-exponential factors
		if (iStick_ == false)
		{			
			if (composition_units_ == PhysicalConstants::UNITS_STD)
			{
				A_ *=	pow(1.e-3, 1.-sumLambdaGasReactants_-sumLambdaSiteReactants_) / 
						pow(1.e-4, 1.-1.5*sumLambdaGasReactants_-sumLambdaSiteReactants_)   / 
						pow(a_conversion_, 1.-sumLambdaGasReactants_-sumLambdaSiteReactants_);				// [m, kmol, s]
			}
			else 
			{
				double conversion_factor = 1.;

				if (composition_units_ == PhysicalConstants::UNITS_ATM)				conversion_factor = 101325.;
				else if (composition_units_ == PhysicalConstants::UNITS_BAR)		conversion_factor = 100000.;
				else if (composition_units_ == PhysicalConstants::UNITS_PASCALS)	conversion_factor = 1.;
				else if (composition_units_ == PhysicalConstants::UNITS_TORR)		conversion_factor = 101325./760.;
				else if (composition_units_ == PhysicalConstants::UNITS_DYNES)		conversion_factor = 0.1;

				A_ *=	pow(conversion_factor, -sumLambdaGasReactants_)*std::pow(10., 1.-sumLambdaSiteReactants_) /
						pow(a_conversion_, 1.-sumLambdaSiteReactants_);									// [m, kmol, s];
			}

			if (iLangmuir_ == true)
			{
				double sum = 0.;
				for(unsigned int j=0;j<langmuir_species_index_.size();j++)
					if (langmuir_numerator_species_[j] == true)
						sum += langmuir_order_[j];

					A_ *= std::pow(1.e3, sum);
			}

			if (iExplicitlyReversible_ == true)
				ARev_ *=	std::pow(1.e-3, 1. - sumLambdaGasProducts_ - sumLambdaSiteProducts_) /
							std::pow(1.e-4, 1. - 1.5*sumLambdaGasProducts_ - sumLambdaSiteProducts_) /
							std::pow(a_conversion_, 1. - sumLambdaGasProducts_ - sumLambdaSiteProducts_);		// [m, kmol, s]
		}

		if (iCoverageDependent_ == true)
		{
			for (unsigned int j = 0; j < coverage_dependent_species_site_type_.size(); j++)
				coverage_dependent_epsilon_[j] *= e_conversion_;
		}
	}

	void ReactionPolicy_Surface_CHEMKIN::ConversionFactors() 
	{
		if		(e_units_ == PhysicalConstants::UNITS_CAL_MOLE)		e_conversion_ = 1.e3*Conversions::J_from_cal;
		else if (e_units_ == PhysicalConstants::UNITS_KCAL_MOLE)	e_conversion_ = 1.e6*Conversions::J_from_cal;
		else if (e_units_ == PhysicalConstants::UNITS_JOULES_MOLE)	e_conversion_ = 1.e3;
		else if (e_units_ == PhysicalConstants::UNITS_KJOULES_MOLE)	e_conversion_ = 1.e6;
		else if (e_units_ == PhysicalConstants::UNITS_EVOLTS)		e_conversion_ = 1.e3*PhysicalConstants::Nav_mol*Conversions::J_from_eV;
		else if (e_units_ == PhysicalConstants::UNITS_KELVINS)		e_conversion_ = 1.e3*PhysicalConstants::R_J_mol;

		if		(a_units_ == PhysicalConstants::UNITS_MOLES)		a_conversion_ = 1.;
		else if	(a_units_ == PhysicalConstants::UNITS_MOLECULES)	a_conversion_ = PhysicalConstants::Nav_mol;

		// Calculate number of moles
		sumNuReactants_ = 0.;
		for(unsigned int i=0;i<reactant_nu_.size();i++)
			sumNuReactants_ += reactant_nu_[i];
		sumNuProducts_ = 0.;
		for(unsigned int i=0;i<product_nu_.size();i++)
			sumNuProducts_ += product_nu_[i];
	}

	void ReactionPolicy_Surface_CHEMKIN::GetReactionString(const std::vector<std::string>& list_species, std::string& line_reaction) const
	{
		// Reactants
		for(unsigned int i=0;i<reactant_nu_.size();i++)
		{
			if (reactant_nu_[i] != 1)
			{
				std::stringstream nu; nu << reactant_nu_[i];
				line_reaction += nu.str();
			}
			line_reaction += list_species[reactant_nu_indices_[i]];
			
			{
				if (i == reactant_nu_.size()-1) line_reaction += " ";
				else line_reaction += " + ";
			}
		}

		// Rversible vs Non reversible
		if (iReversible_ == true) line_reaction += " = ";
		else                      line_reaction += " => ";

		// Products
		for(unsigned int i=0;i<product_nu_.size();i++)
		{
			if (product_nu_[i] != 1)
			{
				std::stringstream nu; nu << product_nu_[i];
				line_reaction += nu.str();
			}
			line_reaction += list_species[product_nu_indices_[i]];

			{
				if (i == product_nu_.size()-1) line_reaction += " ";
				else line_reaction += " + ";
			}
		}
	}

	void ReactionPolicy_Surface_CHEMKIN::GetReactionStringCHEMKIN(	const std::vector<std::string>& list_species, std::stringstream& reaction_data) const
	{
		std::vector<bool> isReducedSpecies(list_species.size());
		for(unsigned int i=0;i<list_species.size();i++)
			isReducedSpecies[i] = true;

		GetReactionStringCHEMKIN( list_species, reaction_data, isReducedSpecies);
	}

	void ReactionPolicy_Surface_CHEMKIN::GetReactionStringCHEMKIN(const std::vector<std::string>& list_species,
                std::stringstream& reaction_data, const std::vector<bool>& isReducedSpecies) const
	{
        std::string reaction_string;
        GetReactionString(list_species, reaction_string);
        boost::erase_all(reaction_string, " ");
        reaction_data.precision(4);
        
        {
            reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << A_ / A_conversion();
            reaction_data.precision(3);
            reaction_data.width(9);
			reaction_data << std::fixed << std::right << beta_;
            reaction_data.precision(2);
            reaction_data << std::setw(13) << std::right <<E_over_R() * PhysicalConstants::R_cal_mol << std::endl;
            
            if(IsDuplicate() == true)
                reaction_data << " DUPLICATE" << std::endl;
            
            if(IsExplicitlyReversible() == true)
            {
                reaction_data << " REV /  ";
                reaction_data.precision(4);
                reaction_data << std::scientific << ARev_ / Arev_conversion() << "  ";
                reaction_data.precision(3);
		reaction_data <<std::fixed <<  betaRev_ << "  ";
                reaction_data.precision(2);
                reaction_data << E_over_R_reversible() * PhysicalConstants::R_cal_mol;
                reaction_data << "  /" << std::endl;
            }
        }
        
        if( IsFORD() == true)
        {
            reaction_data.unsetf(std::ios_base::floatfield);
            reaction_data.precision(4);
            for(unsigned int k = 0; k < reactant_lambda_.size(); k++)
            {
                reaction_data << " FORD /  ";
                
                int index = reactant_lambda_indices_[k];
                
                reaction_data << list_species[index] << "  ";
                
                reaction_data << std::showpoint <<std::fixed << reactant_lambda_[k];
                reaction_data << "/" << std::endl;
            }
        }
	}

	double ReactionPolicy_Surface_CHEMKIN::A_conversion() const
	{
		double conversion_factor = 0;
		// simple
		{
			conversion_factor = std::pow(1.e-3, 1. - sumLambdaGasReactants_ - sumLambdaSiteReactants_) /
								std::pow(1.e-4, 1. - 1.5*sumLambdaGasReactants_ - sumLambdaSiteReactants_) /
								std::pow(a_conversion_, 1. - sumLambdaGasReactants_ - sumLambdaSiteReactants_);
		}
            
		return conversion_factor;
	}
        
	double ReactionPolicy_Surface_CHEMKIN::Arev_conversion() const
    {
        double conversion_factor = 0;
        if(	tag_reaction_ == PhysicalConstants::REACTION_SURFACE_SIMPLE )
        {
            if(iExplicitlyReversible_ == true)
				conversion_factor = std::pow(1.e-3, 1. - sumLambdaGasReactants_ - sumLambdaSiteReactants_) /
									std::pow(1.e-4, 1. - 1.5*sumLambdaGasReactants_ - sumLambdaSiteReactants_) /
									std::pow(a_conversion_, 1. - sumLambdaGasReactants_ - sumLambdaSiteReactants_);
            else
                ErrorMessage("double Arev_conversion() const", "The reaction is not reversible!");
        }
        else
        {
            ErrorMessage("double Arev_conversion() const", "Conversion of A reversible allowed only for "
                    "simple and third body reactions");
        }
            
        return conversion_factor;
    }

	void ReactionPolicy_Surface_CHEMKIN::WriteAdditionalDataOnASCIIFile(std::ostream& fOut) const
	{
		if (tag_reaction_ == PhysicalConstants::REACTION_SURFACE_SIMPLE)
		{
			
		}
	}

	void ReactionPolicy_Surface_CHEMKIN::WriteShortSummary(std::ostream& fOut, const std::vector<std::string>& list_species) const
	{
		std::string line_reaction;
		GetReactionString(list_species, line_reaction);

		fOut << line_reaction << std::endl;

		// Conventional reactions
		{
			fOut << std::setw(9)  << " ";
			fOut << std::setw(9)  << std::left << "k:";
			fOut << std::scientific	<< std::setprecision(6) << std::right << A_ << "\t";
			fOut << std::setw(8) << std::setprecision(2) << std::fixed << std::right << beta_;
			fOut << std::setw(14) << std::setprecision(2) << std::fixed << std::right << E_/Conversions::J_from_kcal << std::endl;
		}

		// Units
		if (composition_units_ != PhysicalConstants::UNITS_STD)
		{
			fOut << std::setw(9) << " ";
			if (composition_units_ == PhysicalConstants::UNITS_ATM)
				fOut << "Gas-phase species units are atmospheres (ATM)" << std::endl;
			else if (composition_units_ == PhysicalConstants::UNITS_BAR)
				fOut << "Gas-phase species units are bars (BAR)" << std::endl;
			else if (composition_units_ == PhysicalConstants::UNITS_TORR)
				fOut << "Gas-phase species units are torrs (TORR)" << std::endl;
			else if (composition_units_ == PhysicalConstants::UNITS_PASCALS)
				fOut << "Gas-phase species units are Pascals (PASCALS)" << std::endl;
			else if (composition_units_ == PhysicalConstants::UNITS_DYNES)
				fOut << "Gas-phase species units are dynes per cm2 (DYNE)" << std::endl;
		}

		// Reverse Rate
		if (iExplicitlyReversible_ == true)
		{
			fOut << std::setw(9) << " ";
			fOut << "This reaction has the explicit reverse reaction: ";
			fOut << std::setw(9)  << " ";
			fOut << std::setw(9)  << std::left << "krev:";
			fOut << std::scientific	<< std::setprecision(6) << std::right << ARev_ << "\t";
			fOut << std::setw(8) << std::setprecision(2) << std::fixed << std::right << betaRev_;
			fOut << std::setw(14) << std::setprecision(2) << std::fixed << std::right << ERev_/Conversions::J_from_kcal << std::endl;
		}

		// Ford reactions
		if (iFord_ == true)
		{
			for(unsigned int j=0;j<reactant_lambda_indices_.size();j++)
			{
				fOut << std::setw(9) << " ";
				fOut << "Reaction order of species ";
				fOut << std::setw(18) << std::left << list_species[reactant_lambda_indices_[j]] << " equal to ";	
				fOut << std::setprecision(4) << std::fixed << reactant_lambda_[j] << std::endl;	
			}
		}

		// Stick reaction
		if (iStick_ == true)
		{
			fOut << std::setw(9) << " ";
			fOut << "Coefficients are sticking parameters...";
			
			fOut << std::endl;
			fOut << std::setw(9) << " ";
			if (iMotzWiseCorrection_ == false)	fOut << std::setw(9)  << std::left << "The Motz-Wise correction is turned off" << std::endl;
			else                                fOut << std::setw(9)  << std::left << "The Motz-Wise correction is turned on" << std::endl;
			
			for(unsigned int j=0;j<stick_indices_site_species_.size();j++)
			{
				fOut << std::setw(9) << " ";
				fOut << std::setw(9) << std::left << "site species: ";
				fOut << std::scientific	<< std::setprecision(6) << std::right << list_species[number_of_gas_species_+stick_indices_site_species_[j]-1] << "\t" << stick_exponents_site_species_[j] << std::endl;
			}
		}

		// Surface dependent reactions
		if (iCoverageDependent_ == true)
		{
			for(unsigned int j=0;j<coverage_dependent_species_site_type_.size();j++)
			{
				fOut << std::setw(9) << " ";
				if (coverage_dependent_species_site_type_[j] == true)
					fOut << "Coverage dependent parameters for species " << list_species[number_of_gas_species_+coverage_dependent_species_index_[j]-1] << ": " << std::endl;	
				else
					fOut << "Coverage dependent parameters for species " << list_species[number_of_gas_species_+number_of_bulk_species_+coverage_dependent_species_index_[j]-1] << ": " << std::endl;	
				fOut << std::setw(18) << " ";
				fOut << coverage_dependent_eta_[j] << "\t";	
				fOut << coverage_dependent_mu_[j] << "\t";
				fOut << coverage_dependent_epsilon_[j] / Conversions::J_from_kcal << "\t";
				fOut << std::endl;
			}
		}

		// Langmuir-Hinshelwood reactions
		if (iLangmuir_ == true)
		{
			fOut << std::setw(9) << " ";
			fOut << "Langmuir-Hinshelwood reaction with denominator exponent parameter of " << langmuir_denominator_order_ << std::endl;

			if (std::find(langmuir_numerator_species_.begin(), langmuir_numerator_species_.end(), true)!=langmuir_numerator_species_.end())
			{
				fOut << std::setw(9) << " ";
				fOut << "Langmuir-Hinshelwood with explicit inclusion of equilibrium constants for: ";
				for(unsigned int j=0;j<langmuir_species_index_.size();j++)
					if (langmuir_numerator_species_[j] == true)	fOut << list_species[langmuir_species_index_[j]-1] << "\t";
				fOut << std::endl;
			}

			for(unsigned int j=0;j<langmuir_species_index_.size();j++)
			{
				fOut << std::setw(9) << " ";
				fOut << "Langmuir-Hinshelwood parameters for species " << list_species[langmuir_species_index_[j]-1] << ": " << std::endl;	
				fOut << std::setw(18) << " ";
				fOut << langmuir_A_[j] << "\t";	
				fOut << langmuir_Beta_[j] << "\t";
				fOut << langmuir_H_over_R_[j]*PhysicalConstants::R_cal_mol << "\t";	
				fOut << langmuir_order_[j] << "\t";	
				fOut << std::endl;
			}
		}
		
		// Lumped reactions
		if (iLumped_ == true)
		{
			fOut << std::setw(9) << " ";
			fOut << "Lumped reaction: " << name_of_lumped_function_ << std::endl;
		}
	}

	double ReactionPolicy_Surface_CHEMKIN::GetForwardConversionFactor() const
	{
		return 1. / (	std::pow(1.e-3, 1. - sumLambdaGasReactants_ - sumLambdaSiteReactants_) /
						std::pow(1.e-4, 1. - 1.5*sumLambdaGasReactants_ - sumLambdaSiteReactants_));
	}

	double ReactionPolicy_Surface_CHEMKIN::GetBackwardConversionFactor() const
	{
		return 1. / (	std::pow(1.e-3, 1. - sumLambdaGasReactants_ - sumLambdaSiteReactants_) /
						std::pow(1.e-4, 1. - 1.5*sumLambdaGasReactants_ - sumLambdaSiteReactants_));
	}

	void ReactionPolicy_Surface_CHEMKIN::WriteSummary(std::ofstream& fOut, const std::vector<std::string>& list_species, const unsigned int index) const
	{
		std::string line_reaction;
		GetReactionString(list_species, line_reaction);
	
		fOut << "================================================================================================================================" << std::endl;
		fOut << " KINETIC DATA - REACTION  " << index << " " << TagASCII() << std::endl;
		fOut << "  " << line_reaction << std::endl; 
		fOut << "================================================================================================================================" << std::endl;
		fOut << " Change in moles in the reaction = " << sumNuProducts_ - sumNuReactants_ << std::endl;
	
		fOut << " Reaction order (forward, gas)   = " << std::setprecision(3) << sumLambdaGasReactants_ << std::endl;
		fOut << " Reaction order (forward, site)  = " << std::setprecision(3) << sumLambdaSiteReactants_ << std::endl;
		fOut << " Reaction order (forward, bulk)  = " << std::setprecision(3) << sumLambdaBulkReactants_ << std::endl;
		if (iReversible_ == true)
		{
			fOut << " Reaction order (backward, gas)  = " << std::setprecision(3) << sumLambdaGasProducts_ << std::endl;
			fOut << " Reaction order (backward, site) = " << std::setprecision(3) << sumLambdaSiteProducts_ << std::endl;
			fOut << " Reaction order (backward, bulk) = " << std::setprecision(3) << sumLambdaBulkProducts_ << std::endl;
		}

		{
			std::string line_kForwardStringSI;
			std::string line_kForwardStringCGS;
		//	OpenSMOKE_Utilities::GetKineticConstantString(A_, beta_, E_, sumLambdaReactants_, line_kForwardStringSI, line_kForwardStringCGS);
		//	fOut << " " << line_kForwardStringSI  << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsSI(sumLambdaReactants_)  << " and [J/kmol]" << std::endl;
		//	fOut << " " << line_kForwardStringCGS << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsCGS(sumLambdaReactants_) << " and [cal/mol]" << std::endl;
		}
	
	//	if (iReversible_ == true)
	//		fOut << " Reverse reaction units: " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsSI(sumLambdaProducts_) << " or " 
	//											<< OpenSMOKE_Utilities::GetUnitsOfKineticConstantsCGS(sumLambdaProducts_) << std::endl;
	}

	bool ReactionPolicy_Surface_CHEMKIN::WriteUBIParametersOnFile(std::ostream& fOutput) const
	{
		fOutput << ubiqep_A_ << " " << ubiqep_Beta_ << " " << ubiqep_BondIndex_ << " " << ubiqep_reaction_class_ << " " << ubiqep_reaction_type_ << std::endl;
		
		if (ubiqep_reaction_type_ == PhysicalConstants::UBIQEP_TYPE_ADSORPTION)
			fOutput << original_reactant_nu_indices_[0]+1 << std::endl;

		if (ubiqep_direct_ == true)
		{
			if (ubiqep_reaction_class_ == 1)
				fOutput << ubiqep_index_A_ << " " << ubiqep_index_Star_ << " " << ubiqep_index_AStar_ << std::endl;
			else if (ubiqep_reaction_class_ == 2)
				fOutput << ubiqep_index_A2_ << " " << ubiqep_index_Star_ << " " << ubiqep_index_AStar_ << " " 
						<< ubiqep_index_A_ << std::endl;
			else if (ubiqep_reaction_class_ == 3)
				fOutput << ubiqep_index_A2_ << " " << ubiqep_index_Star_ << " " << ubiqep_index_AStar_ << " " 
						<< ubiqep_index_A_ << std::endl;
			else if (ubiqep_reaction_class_ == 4)
				fOutput << ubiqep_index_AB_ << " " << ubiqep_index_Star_ << " " << ubiqep_index_AStar_ << " " << ubiqep_index_BStar_ << " " 
						<< ubiqep_index_A_ << " " << ubiqep_index_B_ << std::endl;
			else if (ubiqep_reaction_class_ == 5)
				fOutput << ubiqep_index_ABStar_ << " " << ubiqep_index_Star_ << " " << ubiqep_index_AStar_ << " " << ubiqep_index_BStar_ << " " 
						<< ubiqep_index_A_ << " " << ubiqep_index_B_ << " " << ubiqep_index_AB_ << std::endl;
			else if (ubiqep_reaction_class_ == 6)
				fOutput << ubiqep_index_AStar_ << " " << ubiqep_index_BStar_ << " " << ubiqep_index_CStar_ << " " << ubiqep_index_DStar_ << " " 
						<< ubiqep_index_A_ << " " << ubiqep_index_B_ << " " << ubiqep_index_C_ << " " << ubiqep_index_D_<< std::endl;
			else if (ubiqep_reaction_class_ == 7)
				fOutput << ubiqep_index_CStar_ << " " << ubiqep_index_AStar_ << " " << ubiqep_index_BStar_ << " " 
						<< ubiqep_index_A_ << " " << ubiqep_index_B_ << " " << ubiqep_index_C_ << std::endl;
			else if (ubiqep_reaction_class_ == 8)
				fOutput << ubiqep_index_CStar_ << " " << ubiqep_index_DStar_ << " " << ubiqep_index_AStar_ << " " 
						<< ubiqep_index_A_ << " " << ubiqep_index_C_ << " " << ubiqep_index_D_ << std::endl;
		}
		return true;
	}

	bool ReactionPolicy_Surface_CHEMKIN::FatalErrorMessage(const std::string message)
	{
		std::cout << message << std::endl;
		return false;
	}
}



