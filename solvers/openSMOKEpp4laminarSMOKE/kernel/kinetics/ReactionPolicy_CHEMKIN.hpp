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
#include "ExtendedPressureLogarithmicRateExpression.h"
#include "ExtendedFallOff.h"
#include "ReactionPolicy_CHEMKIN.h"

namespace OpenSMOKE
{
	ReactionPolicy_CHEMKIN::ReactionPolicy_CHEMKIN() 
	{
		SetDefaultUnits();
	}

	void ReactionPolicy_CHEMKIN::SetDefaultUnits()
	{
		iThirdBody_ = false;
		iHigh_ = false;
		iLow_ = false;
		iSRI_ = false;
		iLindemann_ = false;
		iTroe_ = false;
		iReversible_ = false;
		iExplicitlyReversible_ = false;
		iFord_ = false;
		iRord_ = false;
		iDuplicate_ = false;
		iChebyshev_ = false;
		iFit1_ = false;
		iJan_ = false;
		iLandauTeller_ = false;
		iPlog_ = false;
		iExtPlog_ = false;
		iExtLow_ = false;
		pressureDependentSpeciesIndex_ = -1;
	}

	ReactionPolicy_CHEMKIN::ReactionPolicy_CHEMKIN(const ReactionPolicy_CHEMKIN& orig) 
	{
	}

	//ReactionPolicy_CHEMKIN::~ReactionPolicy_CHEMKIN() {
	//}

	PhysicalConstants::TAG_REACTION ReactionPolicy_CHEMKIN::Tag() const
	{
		return tag_reaction_;
	}

	std::string ReactionPolicy_CHEMKIN::TagASCII() const
	{
		if (tag_reaction_ == PhysicalConstants::REACTION_SIMPLE)
			return "Simple";
		else if (tag_reaction_ == PhysicalConstants::REACTION_THIRDBODY)
			return "Third-body";
		else if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF)
			return "Fall-off (Lindemann)";
		else if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR)
			return "CABR (Lindemann)";
		else if (tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF)
			return "Fall-off (Troe)";
		else if (tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR)
			return "CABR (Troe)";
		else if (tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF)
			return "Fall-off (SRI)";
		else if (tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)
			return "CABR (SRI)";
		else if (tag_reaction_ == PhysicalConstants::REACTION_CHEBYSHEV)
			return "Chebyshev Polynomial Rate Expression";
		else if (tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
			return "Extended Fall-off";
		else
			return "Unknown reaction";
	}

	bool ReactionPolicy_CHEMKIN::SetUnits(const PhysicalConstants::UNITS_REACTION a_units, const PhysicalConstants::UNITS_REACTION e_units)
	{
		a_units_ = a_units;
		e_units_ = e_units;
		
		return true;
	}

	bool ReactionPolicy_CHEMKIN::ReadReactionFromStrings(const std::vector<std::string>& lines, const std::map<std::string, unsigned int>& map_of_species)
	{
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

			std::string pressureDependentReaction;
		
			// Look for Pressure dependent reaction
			// At the end of this part of the code the pressureDependentReaction tag is filled
			{
				pressureDependentReaction = "none";
				std::vector<unsigned int> indices_to_erase;
			
				for(unsigned int i=0;i<reactant_species.size();i++)
				{
					if (reactant_species[i] == "$PDR$")
						indices_to_erase.push_back(i);
				}

				if(indices_to_erase.size() > 1)
				{
					std::cout << "The pressure dependent reaction tag is specified more than once" << std::endl;
					return false;
				}
				if(indices_to_erase.size() == 1)
				{
					pressureDependentReaction = "all";
					reactant_species.erase(reactant_species.begin()+indices_to_erase[0]);
					reactant_nu_.erase(reactant_nu_.begin()+indices_to_erase[0]);
					indices_to_erase.resize(0);
				}

				for(unsigned int i=0;i<product_species.size();i++)
				{
					if (product_species[i] == "$PDR$")
						indices_to_erase.push_back(i);
				}
				if(indices_to_erase.size() > 1)
				{
					std::cout << "The pressure dependent reaction tag is specified more than once" << std::endl;
					return false;
				}
				if(indices_to_erase.size() == 1)
				{
					pressureDependentReaction = "all";
					product_species.erase(product_species.begin()+indices_to_erase[0]);
					product_nu_.erase(product_nu_.begin()+indices_to_erase[0]);
					indices_to_erase.clear();
				}
			}

			// Look for Pressure dependent reaction species
			{
				std::vector<unsigned int> indices_to_erase;
			
				for(unsigned int i=0;i<reactant_species.size();i++)
				{
					std::string substring=reactant_species[i].substr(0,5);
					if (substring == "$PDS$")
					{
						boost::replace_all(reactant_species[i], "$PDS$", "");
						indices_to_erase.push_back(i);
					}
				}

				if(indices_to_erase.size() > 1)
				{
					std::cout << "The pressure dependent reaction tag is specified more than once" << std::endl;
					return false;
				}
				if(indices_to_erase.size() == 1)
				{
					if (pressureDependentReaction == "all")
					{
						std::cout << "The pressure dependent reaction is not specified correctly." << std::endl;
						std::cout << "Two different tags are used for the pressure dependence (reactant side)." << std::endl;
						return false;
					}
					pressureDependentReaction = "partial";
					std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(reactant_species[indices_to_erase[0]]);
					if( it == map_of_species.end())
					{
						std::cout << "The following species is not available (pressure dependent reaction, reactant side): " << reactant_species[indices_to_erase[0]] << std::endl;
						return false;
					}
					else
					{
						pressureDependentSpeciesIndex_ = (*it).second;
					}
				
					reactant_species.erase(reactant_species.begin()+indices_to_erase[0]);
					reactant_nu_.erase(reactant_nu_.begin()+indices_to_erase[0]);
					indices_to_erase.resize(0);
				}
			
				for(unsigned int i=0;i<product_species.size();i++)
				{
					std::string substring=product_species[i].substr(0,5);
					if (substring == "$PDS$")
					{
						boost::replace_all(product_species[i], "$PDS$", "");
						indices_to_erase.push_back(i);
					}
				}
				if(indices_to_erase.size() > 1)
				{
					std::cout << "The pressure dependent reaction tag is specified more than once" << std::endl;
					return false;
				}
				if(indices_to_erase.size() == 1)
				{
					if (pressureDependentReaction == "all")
					{
						std::cout << "The pressure dependent reaction is not specified correctly." << std::endl;
						std::cout << "Two different tags are used for the pressure dependence (product side)." << std::endl;
						return false;
					}
					pressureDependentReaction = "partial";
					std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(product_species[indices_to_erase[0]]);
					if(it == map_of_species.end())
					{
						std::cout << "The following species is not available (pressure dependent reaction, product side): " << product_species[indices_to_erase[0]] << std::endl;
						return false;
					}
					else
					{
						if (pressureDependentSpeciesIndex_ != -1)
							if ( pressureDependentSpeciesIndex_ != (*it).second )
							{
								std::cout << "The pressure dependent reaction is not specified correctly." << std::endl;
								std::cout << "Please check the species acting as third/body." << std::endl;
								return false;
							}
					}

					product_species.erase(product_species.begin()+indices_to_erase[0]);
					product_nu_.erase(product_nu_.begin()+indices_to_erase[0]);
					indices_to_erase.resize(0);
				}
			}

			// Look for Third-body reaction
			{
				std::vector<unsigned int> indices_to_erase;
			
				for(unsigned int i=0;i<reactant_species.size();i++)
				{
					if (reactant_species[i] == "$TBR$")
						indices_to_erase.push_back(i);
				}

				if(indices_to_erase.size() > 1)
				{
					std::cout << "The third body reaction tag is specified more than once" << std::endl;
					return false;
				}
				if(indices_to_erase.size() == 1)
				{
					if (pressureDependentReaction != "none")
					{
						std::cout << "The third body reaction is not specified correctly" << std::endl;
						return false;
					}

					reactant_species.erase(reactant_species.begin()+indices_to_erase[0]);
					reactant_nu_.erase(reactant_nu_.begin()+indices_to_erase[0]);
					indices_to_erase.resize(0);

					iThirdBody_ = true;
				}

				for(unsigned int i=0;i<product_species.size();i++)
				{
					if (product_species[i] == "$TBR$")
						indices_to_erase.push_back(i);
				}
				if(indices_to_erase.size() > 1)
				{
					std::cout << "The third body reaction tag is specified more than once" << std::endl;
					return false;
				}
				if(indices_to_erase.size() == 1)
				{
					if (pressureDependentReaction != "none")
					{
						std::cout << "The third body reaction is not specified correctly" << std::endl;
						return false;
					}

					product_species.erase(product_species.begin()+indices_to_erase[0]);
					product_nu_.erase(product_nu_.begin()+indices_to_erase[0]);
					indices_to_erase.clear();

					iThirdBody_ = true;
				}
			}
		
			reactant_nu_indices_.resize(reactant_nu_.size());
			for(unsigned int i=0;i<reactant_nu_.size();i++)
			{
				std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(reactant_species[i]);
				if( it == map_of_species.end())
				{
                                    if (reactant_species[i] == "M" || reactant_species[i] == "m")
                                        std::cout << "The third body species (M or m) must be always the last species, both on the reactant and the product side (e.g. A+B+M = C+D+M)." << std::endl;
                                    else
					std::cout << "The following species is not available: " << reactant_species[i] << std::endl;
                                    return false;
				}
				else
				{
					reactant_nu_indices_[i] = (*it).second;
				}
			}

			product_nu_indices_.resize(product_nu_.size());
			for(unsigned int i=0;i<product_nu_.size();i++)
			{
				std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(product_species[i]);
				if( it == map_of_species.end())
				{
                                    if (product_species[i] == "M" || product_species[i] == "m")
                                        std::cout << "The third body species (M or m) must be always the last species, both on the reactant and the product side (e.g. A+B+M = C+D+M)." << std::endl;
                                    else
					std::cout << "The following species is not available: " << product_species[i] << std::endl;
				    return false;
				}
				else
				{
					product_nu_indices_[i] = (*it).second;
				}
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

			// Provisional data
		
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
						// Third body efficiencies
						boost::algorithm::trim(line);
						typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
						boost::char_separator<char> sep_slash("/");
						tokenizer_slash tokens(line, sep_slash);
						const std::size_t n = std::distance(tokens.begin(), tokens.end());
				
						//if (n%2 != 0 || line.back() != '/')
						if (n%2 != 0 || line.at(line.size()-1) != '/')
						{
								std::cout << "Third body efficiencies: syntax error!" << std::endl;
								return false;
						}

						tokenizer_slash::iterator tok_slash = tokens.begin();
						for (std::size_t i = 1; i <= n / 2; i++)
						{
							std::string name = *tok_slash;
							boost::algorithm::trim(name);
							++tok_slash;
							std::string number = *(tok_slash);
							boost::algorithm::trim(number);
							++tok_slash;

							std::map<std::string, unsigned int>::const_iterator it=map_of_species.find(name);
							if( it == map_of_species.end())
							{
								std::cout << "Third body efficiencies: The following species is not available: " << name << std::endl;
								return false;
							}
							else
							{
								third_body_indices_.push_back( (*it).second);
								try
								{
									double efficiency = boost::lexical_cast<double>(number);
									third_body_efficiencies_.push_back(efficiency);
								}
								catch(boost::bad_lexical_cast &)
								{
									std::cout << "Third body efficiencies: The efficiencies are not written properly (they are not numbers)" << std::endl;
									return false;
								}  
							}
						}

						if (pressureDependentSpeciesIndex_ != -1)
						{
							std::cout << "Third body reactions require (+M)" << std::endl;
							return false;
						}
					}
					else
					{
						if (boost::iequals(keyword, "CHEB"))
						{
							if (pressureDependentReaction != "all")
							{
								std::cout << "The CHEB option can be used only for pressure-dependent reactions with the (+M) option!" << std::endl;
								return false;
							}
							if (iLow_ == true || iHigh_ == true || troe_.size()!=0 || sri_.size()!=0 || iFit1_ == true || iJan_ == true || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The CHEB option cannot be used together with the following options: LOW | HIGH | TROE | SRI | FIT1 | JAN | LT | PLOG | PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}

							iChebyshev_ = true;
							std::vector<double> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadCoefficients("Chebishev Polynomial Rate Expressions (CHEB). ", line, coefficients);
							if (tag == false) return false;
							for(unsigned int i=0;i<coefficients.size();i++)
								chebyshev_coefficients_.push_back(coefficients[i]);
						}
						else if (boost::iequals(keyword, "COLLEFF"))
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}
						else if (boost::iequals(keyword, "DUP") || boost::iequals(keyword, "DUPLICATE"))
						{
							if (iDuplicate_ == true)
							{
								std::cout << "The DUPLICATE (or DUP) option is used more than once!" << std::endl;
								return false;
							}
							iDuplicate_ = true;
						}
						else if (boost::iequals(keyword, "EXCI"))
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}
						else if (boost::iequals(keyword, "FIT1"))
						{
							if (iFit1_ == true)
							{
								std::cout << "The FIT1 option is used more than once!" << std::endl;
								return false;
							}
							if (iChebyshev_ == true || chebyshev_pressure_limits_.size()!=0 || chebyshev_temperature_limits_.size()!=0 || iJan_ == true || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The FIT1 option cannot be used together with the following options: CHEB | TCHEB | PCHEB | JAN | LT | PLOG | PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Power series modified Arrhenus law (FIT1). ", line, 4, fit1_coefficients_);
							if (tag == false) return false;
							iFit1_ = true;
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
						else if (boost::iequals(keyword, "HIGH"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The HIGH option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}
							if (iHigh_ == true)
							{
								std::cout << "The HIGH option is used more than once!" << std::endl;
								return false;
							}
							if (iLow_ == true || iChebyshev_ == true || chebyshev_pressure_limits_.size()!=0 || chebyshev_temperature_limits_.size()!=0 ||
								iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The HIGH option cannot be used together with the following options: LOW | CHEB | TCHEB | PCHEB | LT | PLOG | PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}

							std::vector<double> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("High-pressure limit (HIGH). ", line, 3, coefficients);
							if (tag == false) return false;
							else
							{
								iHigh_ = true;
								if (iTroe_ == false && iSRI_ == false)	iLindemann_ = true;
								AInf_ = coefficients[0];
								betaInf_ = coefficients[1];
								EInf_ = coefficients[2];

								if (AInf_ < 0.)
								{
									std::cout << "Negative frequency factors are not allowed for HIGH type reactions!" << std::endl;
									CheckForFatalError(false);
								}
							}
						}
						else if (boost::iequals(keyword, "JAN"))
						{
							if (iJan_ == true)
							{
								std::cout << "The JAN option is used more than once!" << std::endl;
								return false;
							}
							if (iChebyshev_ == true || chebyshev_pressure_limits_.size()!=0 || chebyshev_temperature_limits_.size()!=0 || iFit1_ == true || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The JAN option cannot be used together with the following options: CHEB | TCHEB | PCHEB | FIT1 | LT | PLOG || PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Janev-Langer reaction rate (JAN). ", line, 9, janev_langer_coefficients_);
							if (tag == false) return false;
							iJan_ = true;
						}
						else if (boost::iequals(keyword, "LOW"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The LOW option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}
							if (iLow_ == true)
							{
								std::cout << "The LOW option is used more than once!" << std::endl;
								return false;
							}
							if (iHigh_ == true || iChebyshev_ == true || chebyshev_pressure_limits_.size()!=0 || chebyshev_temperature_limits_.size()!=0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true)
							{
								std::cout << "The LOW option cannot be used together with the following options: HIGH | CHEB | TCHEB | PCHEB | LT | PLOG || PLOGMX | PLOGSP" << std::endl;
								return false;
							}

							std::vector<double> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Low-pressure limit (LOW). ", line, 3, coefficients);
							if (tag == false) return false;
							else
							{
								iLow_ = true;
								if (iTroe_ == false && iSRI_ == false)	iLindemann_ = true;
								AInf_ = A_;
								betaInf_ = beta_;
								EInf_ = E_;

								A_ = coefficients[0];
								beta_ = coefficients[1];
								E_ = coefficients[2];

								if (AInf_ < 0.)
								{
									std::cout << "Negative frequency factors are not allowed for LOW type reactions!" << std::endl;
									CheckForFatalError(false);
								}
							}
						}
						else if (boost::iequals(keyword, "LOWMX"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The LOWMX option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}
						
							if (iLow_ == true || iHigh_ == true || iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true)
							{
								std::cout << "The LOWMX option cannot be used together with the following options: LOW | HIGH | CHEB | TCHEB | PCHEB | LT | PLOG || PLOGMX | PLOGSP" << std::endl;
								return false;
							}

							std::vector<std::string> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Low-pressure limit (LOWMX). ", line, 3, coefficients);
							if (tag == false) return false;
							else
							{
								for (unsigned int i = 0; i < coefficients.size(); i++)
								{
									// Check if they are really numbers
									try
									{
										const double dummy = boost::lexical_cast<double>(coefficients[i]);
									}
									catch (boost::bad_lexical_cast &)
									{
										std::cout << "LOWMX: The coefficients are not written properly (they are not numbers)" << std::endl;
										return false;
									}
								}

								coefficients.insert(coefficients.begin(), "-1");
								coefficients.insert(coefficients.begin(), "Mixture");
								coefficients.insert(coefficients.begin(), "low");
								extendedfalloff_coefficients_.push_back(coefficients);
								iExtLow_ = true;
							}
						}
						else if (boost::iequals(keyword, "LOWSP"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The LOWSP option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}

							if (iLow_ == true || iHigh_ == true || iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true)
							{
								std::cout << "The LOWSP option cannot be used together with the following options: LOW | HIGH | CHEB | TCHEB | PCHEB | LT | PLOG || PLOGMX | PLOGSP" << std::endl;
								return false;
							}

							std::vector<std::string> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Low-pressure limit (LOWSP). ", line, 4, coefficients);
							if (tag == false) return false;
							else
							{
								// Bath species
								{
									std::map<std::string, unsigned int>::const_iterator it = map_of_species.find(coefficients[0]);
									if (it == map_of_species.end())
									{
										std::cout << "LOWSP: The following bath species is not available: " << coefficients[0] << std::endl;
										return false;
									}

									coefficients.insert(coefficients.begin() + 1, boost::lexical_cast<std::string>((*it).second));
								}

								for (unsigned int i = 2; i < coefficients.size(); i++)
								{
									// Check if they are really numbers
									try
									{
										const double dummy = boost::lexical_cast<double>(coefficients[i]);
									}
									catch (boost::bad_lexical_cast &)
									{
										std::cout << "LOWSP: The coefficients are not written properly (they are not numbers)" << std::endl;
										return false;
									}
								}

								coefficients.insert(coefficients.begin(), "low");
								extendedfalloff_coefficients_.push_back(coefficients);
								iExtLow_ = true;
							}
						}
						else if (boost::iequals(keyword, "LT"))
						{
							if (iLandauTeller_ == true)
							{
								std::cout << "The LT option is used more than once!" << std::endl;
								return false;
							}
							if (iChebyshev_ == true || chebyshev_pressure_limits_.size()!=0 || chebyshev_temperature_limits_.size()!=0 || iFit1_ == true || iJan_ == true ||
								iLow_ == true || iHigh_ == true || iTroe_ == true || iSRI_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The LT option cannot be used together with the following options: CHEB | TCHEB | PCHEB | FIT1 | JAN | LOW | HIGH | TROE | SRI | PLOG | PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Landau Teller reaction rate (LT). ", line, 2, landau_teller_coefficients_);
							if (tag == false) return false;
							iLandauTeller_ = true;
						}
						else if (boost::iequals(keyword, "MOME"))
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}
						else if (boost::iequals(keyword, "PCHEB"))
						{
							if (pressureDependentReaction != "all")
							{
								std::cout << "The PCHEB option can be used only for pressure-dependent reactions with the (+M) option!" << std::endl;
								return false;
							}
							if (chebyshev_pressure_limits_.size() != 0)
							{
								std::cout << "The PCHEB option is used more than once!" << std::endl;
								return false;
							}
							if (iLow_ == true || iHigh_ == true || troe_.size()!=0 || sri_.size()!=0 || iFit1_==true || iJan_ == true || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The PCHEB option cannot be used together with the following options: LOW | HIGH | TROE | SRI | FIT1 | JAN | LT | PLOG | PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}

							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure limit for Chebyshev Polynomial Rate Expression (PCHEB). ", line, 2, chebyshev_pressure_limits_);
							if (tag == false) return false;
						}
						else if (boost::iequals(keyword, "PLOG"))
						{
							if (iThirdBody_ == true)
							{
								std::cout << "The PLOG option cannot be used for third-body reactions!" << std::endl;
								return false;
							}
							if (iLow_ == true || iHigh_ == true || troe_.size()!=0 || sri_.size()!=0 || iFit1_ == true || iJan_ == true || iLandauTeller_ == true ||
								iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The PLOG option cannot be used together with the following options: LOW | HIGH | TROE | SRI | FIT1 | JAN | LT | CHEB | PCHEB | TCHEB | PLOG | LOWMX | LOWSP" << std::endl;
								return false;
							}

							iPlog_ = true;
							std::vector<double> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure Dependence through Logarithmic Interpolation Rate Expressions (PLOG). ", line, 4, coefficients);
							if (tag == false) return false;
							for(unsigned int i=0;i<coefficients.size();i++)
								plog_coefficients_.push_back(coefficients[i]);
						}
						else if (boost::iequals(keyword, "PLOGMX"))
						{
							if (iThirdBody_ == true)
							{
								std::cout << "The PLOGMX option cannot be used for third-body reactions!" << std::endl;
								return false;
							}
							if (iLow_ == true || iHigh_ == true || troe_.size() != 0 || sri_.size() != 0 || iFit1_ == true || iJan_ == true || iLandauTeller_ == true ||
								iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The PLOGMX option cannot be used together with the following options: LOW | HIGH | TROE | SRI | FIT1 | JAN | LT | CHEB | PCHEB | TCHEB | PLOG | LOWMX | LOWSP" << std::endl;
								return false;
							}

							iExtPlog_ = true;
							std::vector<std::string> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusWords("Extended Pressure Dependence through Logarithmic Interpolation Rate Expressions (PLOGMX). ", line, 4, coefficients);
							if (tag == false) return false;
							
							// Mixture (same structure for specific species)
							extendedplog_coefficients_.push_back("Mixture");
							extendedplog_coefficients_.push_back("-1");

							// Coefficients
							for (unsigned int i = 0; i < coefficients.size(); i++)
							{
								// Check if they are really numbers
								try
								{
									const double dummy = boost::lexical_cast<double>(coefficients[i]);
								}
								catch (boost::bad_lexical_cast &)
								{
									std::cout << "PLOGMX: The coefficients are not written properly (they are not numbers)" << std::endl;
									return false;
								}

								extendedplog_coefficients_.push_back(coefficients[i]);
							}
						}
						else if (boost::iequals(keyword, "PLOGSP"))
						{
							if (iThirdBody_ == true)
							{
								std::cout << "The PLOGSP option cannot be used for third-body reactions!" << std::endl;
								return false;
							}
							if (iLow_ == true || iHigh_ == true || troe_.size() != 0 || sri_.size() != 0 || iFit1_ == true || iJan_ == true || iLandauTeller_ == true ||
								iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The PLOGSP option cannot be used together with the following options: LOW | HIGH | TROE | SRI | FIT1 | JAN | LT | CHEB | PCHEB | TCHEB | PLOG | LOWMX | LOWSP" << std::endl;
								return false;
							}

							iExtPlog_ = true;
							std::vector<std::string> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusWords("Extended Pressure Dependence through Logarithmic Interpolation Rate Expressions (PLOGSP). ", line, 5, coefficients);
							if (tag == false) return false;

							// Bath species
							{
								std::map<std::string, unsigned int>::const_iterator it = map_of_species.find(coefficients[0]);
								if (it == map_of_species.end())
								{
									std::cout << "PLOGSP: The following bath species is not available: " << coefficients[0] << std::endl;
									return false;
								}

								extendedplog_coefficients_.push_back(coefficients[0]);
								extendedplog_coefficients_.push_back(boost::lexical_cast<std::string>((*it).second));
							}

							// Coefficients
							for (unsigned int i = 1; i < coefficients.size(); i++)
							{
								// Check if they are really numbers
								try
								{
									const double dummy = boost::lexical_cast<double>(coefficients[i]);
								}
								catch (boost::bad_lexical_cast &)
								{
									std::cout << "PLOGSP: The coefficients are not written properly (they are not numbers)" << std::endl;
									return false;
								}

								extendedplog_coefficients_.push_back(coefficients[i]);
							}
						}
						else if (boost::iequals(keyword, "REV"))
						{
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
						else if (boost::iequals(keyword, "RLT"))
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
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
						else if (boost::iequals(keyword, "SRI"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The SRI option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}
							if (iSRI_ == true)
							{
								std::cout << "The SRI option is used more than once!" << std::endl;
								return false;
							}
							if (troe_.size()!=0 || iChebyshev_ == true || chebyshev_pressure_limits_.size()!=0 || chebyshev_temperature_limits_.size()!=0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The SRI option cannot be used together with the following options: TROE | CHEB | TCHEB | PCHEB | LT | PLOG | PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}

							std::vector<double> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure-dependent reaction (SRI). ", line, 3, 5, coefficients);
							if (tag == false) return false;
							else
							{
								iSRI_ = true;
								iLindemann_ = false;
								sri_ = coefficients;
							}
						}
						else if (boost::iequals(keyword, "SRIMX"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The SRIMX option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}
							if (iLow_ == true || iHigh_ == true || troe_.size() != 0 || sri_.size() != 0 || iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true)
							{
								std::cout << "The SRIMX option cannot be used together with the following options: LOW | HIGH | TROE | SRI | CHEB | TCHEB | PCHEB | LT | PLOG | PLOGMX | PLOGSP" << std::endl;
								return false;
							}

							std::vector<std::string> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure-dependent reaction (SRIMX). ", line, 3, 5, coefficients);
							if (tag == false) return false;
							else
							{
								for (unsigned int i = 0; i < coefficients.size(); i++)
								{
									// Check if they are really numbers
									try
									{
										const double dummy = boost::lexical_cast<double>(coefficients[i]);
									}
									catch (boost::bad_lexical_cast &)
									{
										std::cout << "SRIMX: The coefficients are not written properly (they are not numbers)" << std::endl;
										return false;
									}
								}

								coefficients.insert(coefficients.begin(), "-1");
								coefficients.insert(coefficients.begin(), "Mixture");
								coefficients.insert(coefficients.begin(), "sri");
								extendedfalloff_coefficients_.push_back(coefficients);
								iExtLow_ = true;
							}
						}
						else if (boost::iequals(keyword, "SRISP"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The SRISP option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}
							if (iLow_ == true || iHigh_ == true || troe_.size() != 0 || sri_.size() != 0 || iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true)
							{
								std::cout << "The SRISP option cannot be used together with the following options: LOW | HIGH | TROE | SRI | CHEB | TCHEB | PCHEB | LT | PLOG | PLOGMX | PLOGSP" << std::endl;
								return false;
							}

							std::vector<std::string> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure-dependent reaction (SRISP). ", line, 4, 6, coefficients);
							if (tag == false) return false;
							else
							{
								// Bath species
								{
									std::map<std::string, unsigned int>::const_iterator it = map_of_species.find(coefficients[0]);
									if (it == map_of_species.end())
									{
										std::cout << "SRISP: The following bath species is not available: " << coefficients[0] << std::endl;
										return false;
									}

									coefficients.insert(coefficients.begin() + 1, boost::lexical_cast<std::string>((*it).second));
								}

								for (unsigned int i = 2; i < coefficients.size(); i++)
								{
									// Check if they are really numbers
									try
									{
										const double dummy = boost::lexical_cast<double>(coefficients[i]);
									}
									catch (boost::bad_lexical_cast &)
									{
										std::cout << "SRISP: The coefficients are not written properly (they are not numbers)" << std::endl;
										return false;
									}
								}

								coefficients.insert(coefficients.begin(), "sri");
								extendedfalloff_coefficients_.push_back(coefficients);
								iExtLow_ = true;
							}
						}
						else if (boost::iequals(keyword, "TCHEB"))
						{
							if (pressureDependentReaction != "all")
							{
								std::cout << "The TCHEB option can be used only for pressure-dependent reactions with the (+M) option!" << std::endl;
								return false;
							}
							if (chebyshev_temperature_limits_.size() != 0)
							{
								std::cout << "The TCHEB option is used more than once!" << std::endl;
								return false;
							}
							if (iLow_ == true || iHigh_ == true || troe_.size()!=0 || sri_.size()!=0 || iFit1_ == true || iJan_ == true || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The TCHEB option cannot be used together with the following options: LOW | HIGH | TROE | SRI | FIT1 | JAN | LT | PLOG | PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}

							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Temperature limit for Chebyshev Polynomial Rate Expression (TCHEB). ", line, 2, chebyshev_temperature_limits_);
							if (tag == false) return false;
						}
						else if (boost::iequals(keyword, "TDEP"))
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}
						else if (boost::iequals(keyword, "TROE"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The TROE option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}
							if (iTroe_ == true)
							{
								std::cout << "The TROE option is used more than once!" << std::endl;
								return false;
							}
							if (sri_.size()!=0 || iChebyshev_ == true || chebyshev_pressure_limits_.size()!=0 || chebyshev_temperature_limits_.size()!=0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true || iExtLow_ == true)
							{
								std::cout << "The TROE option cannot be used together with the following options: SRI | CHEB | TCHEB | PCHEB | LT | PLOG | PLOGMX | PLOGSP | LOWMX | LOWSP" << std::endl;
								return false;
							}

							std::vector<double> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure-dependent reaction (TROE). ", line, 3, 4, coefficients);
							if (tag == false) return false;
							else
							{
								iTroe_ = true;
								iLindemann_ = false;
								troe_ = coefficients;
							}
						}
						else if (boost::iequals(keyword, "TROEMX"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The TROEMX option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}
							
							if (iLow_ == true || iHigh_ == true || troe_.size() != 0 || sri_.size() != 0 || iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true)
							{
								std::cout << "The TROEMX option cannot be used together with the following options: LOW | HIGH | TROE | SRI | CHEB | TCHEB | PCHEB | LT | PLOG | PLOGMX | PLOGSP" << std::endl;
								return false;
							}

							std::vector<std::string> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure-dependent reaction (TROEMX). ", line, 3, 4, coefficients);
							if (tag == false) return false;
							else
							{
								for (unsigned int i = 0; i < coefficients.size(); i++)
								{
									// Check if they are really numbers
									try
									{
										const double dummy = boost::lexical_cast<double>(coefficients[i]);
									}
									catch (boost::bad_lexical_cast &)
									{
										std::cout << "TROEMX: The coefficients are not written properly (they are not numbers)" << std::endl;
										return false;
									}
								}

								coefficients.insert(coefficients.begin(), "-1");
								coefficients.insert(coefficients.begin(), "Mixture");
								coefficients.insert(coefficients.begin(), "troe");
								extendedfalloff_coefficients_.push_back(coefficients);
								iExtLow_ = true;
							}
						}
						else if (boost::iequals(keyword, "TROESP"))
						{
							if (pressureDependentReaction == "none")
							{
								std::cout << "The TROESP option can be used only for pressure-dependent reactions!" << std::endl;
								return false;
							}

							if (iLow_ == true || iHigh_ == true || troe_.size() != 0 || sri_.size() != 0 || iChebyshev_ == true || chebyshev_pressure_limits_.size() != 0 || chebyshev_temperature_limits_.size() != 0 || iLandauTeller_ == true || iPlog_ == true || iExtPlog_ == true)
							{
								std::cout << "The TROESP option cannot be used together with the following options: LOW | HIGH | TROE | SRI | CHEB | TCHEB | PCHEB | LT | PLOG | PLOGMX | PLOGSP" << std::endl;
								return false;
							}

							std::vector<std::string> coefficients;
							bool tag = OpenSMOKE_Utilities::ReadReactionKeyWordPlusCoefficients("Pressure-dependent reaction (TROESP). ", line, 4, 5, coefficients);
							if (tag == false) return false;
							else
							{
								// Bath species
								{
									std::map<std::string, unsigned int>::const_iterator it = map_of_species.find(coefficients[0]);
									if (it == map_of_species.end())
									{
										std::cout << "TROESP: The following bath species is not available: " << coefficients[0] << std::endl;
										return false;
									}

									coefficients.insert(coefficients.begin() + 1, boost::lexical_cast<std::string>((*it).second));
								}

								for (unsigned int i = 2; i < coefficients.size(); i++)
								{
									// Check if they are really numbers
									try
									{
										const double dummy = boost::lexical_cast<double>(coefficients[i]);
									}
									catch (boost::bad_lexical_cast &)
									{
										std::cout << "TROESP: The coefficients are not written properly (they are not numbers)" << std::endl;
										return false;
									}
								}

								coefficients.insert(coefficients.begin(), "troe");
								extendedfalloff_coefficients_.push_back(coefficients);
								iExtLow_ = true;
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
								else
								{
									std::cout << "Wrong units" << std::endl;
									std::cout << "Available options: CAL || CAL/MOLE || EVOL || EVOLTS || JOUL || JOULES/MOLE || KCAL || KCAL/MOLE || KJOU || KELV || KELVINS || MOLEC || MOLECULES || MOLE || MOLES" << std::endl;
									return false;
								}
							}
						}
						else if (boost::iequals(keyword, "USRPROG"))
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}
						else if (boost::iequals(keyword, "XSMI"))
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}
						else
						{
							std::cout << "Sorry! The " << keyword << " option is not currently available." << std::endl;
							return false;
						}
					}			
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

			// Cross-checks
			if (pressureDependentReaction == "all" && iHigh_ == false && iLow_ == false && iChebyshev_ == false && iExtLow_ == false)
			{
				std::cout << "No pressure dependence option defined." << std::endl;
				return false;
			}

			if (iThirdBody_ == true && pressureDependentReaction != "none") 
			{
				std::cout << "Internal error tag reaction..." << std::endl;
				CheckForFatalError(false);
			}
		
			if (iThirdBody_ == false && pressureDependentReaction == "none") 
				tag_reaction_ = PhysicalConstants::REACTION_SIMPLE;
			else if (iThirdBody_ == true && pressureDependentReaction == "none") 
				tag_reaction_ = PhysicalConstants::REACTION_THIRDBODY;
			else if (pressureDependentReaction != "none") 
			{
				if (iChebyshev_ == true)
					tag_reaction_ = PhysicalConstants::REACTION_CHEBYSHEV;
				else if (iExtLow_ == true)
					tag_reaction_ = PhysicalConstants::REACTION_EXTENDEDFALLOFF;
				else if (iLindemann_ == true && iLow_ == true)
					tag_reaction_ = PhysicalConstants::REACTION_LINDEMANN_FALLOFF;
				else if (iTroe_ == true && iLow_ == true)
					tag_reaction_ = PhysicalConstants::REACTION_TROE_FALLOFF;
				else if (iSRI_ == true && iLow_ == true)
					tag_reaction_ = PhysicalConstants::REACTION_SRI_FALLOFF;
				else if (iLindemann_ == true && iHigh_ == true)
					tag_reaction_ = PhysicalConstants::REACTION_LINDEMANN_CABR;
				else if (iTroe_ == true && iHigh_ == true)
					tag_reaction_ = PhysicalConstants::REACTION_TROE_CABR;
				else if (iSRI_ == true && iHigh_ == true)
					tag_reaction_ = PhysicalConstants::REACTION_SRI_CABR;
			}
			else
			{
				std::cout << "Internal error tag reaction..." << std::endl;
				CheckForFatalError(false);
			}

			// Checking Chebishev coefficients
			if (tag_reaction_ == PhysicalConstants::REACTION_CHEBYSHEV)
			{
				unsigned int NxM = int(chebyshev_coefficients_[0])*int(chebyshev_coefficients_[1]);
				if ( NxM != chebyshev_coefficients_.size()-2 )
				{
					std::cout << "Error in the definition of the Chebyshev coefficients." << std::endl;
					std::cout << "Expected coefficients: " << NxM << " - Given: " << chebyshev_coefficients_.size()-2 << std::endl;
					return false;
				}
				if (chebyshev_pressure_limits_.size() == 0)
				{
					chebyshev_pressure_limits_.resize(2);
					chebyshev_pressure_limits_[0] = 0.001;
					chebyshev_pressure_limits_[1] = 100.;
				}
				if (chebyshev_temperature_limits_.size() == 0)
				{
					chebyshev_temperature_limits_.resize(2);
					chebyshev_temperature_limits_[0] = 300.;
					chebyshev_temperature_limits_[1] = 2500.;
				}
			}

			// Negative frequency factors are allowed only for simple reactions
			{
				if (A_ < 0. && tag_reaction_ != PhysicalConstants::REACTION_SIMPLE)
				{
					std::cout << "Negative frequency factors are allowed only for conventional (SIMPLE) reactions!" << std::endl;
					CheckForFatalError(false);
				}
			}

			ConvertUnits();

			return true;
	}

	void ReactionPolicy_CHEMKIN::ReactionOrders() 
	{
		sumLambdaReactants_ = 0.;
		for(unsigned int i=0;i<reactant_lambda_.size();i++)
			sumLambdaReactants_ += reactant_lambda_[i];
		sumLambdaProducts_ = 0.;
		for(unsigned int i=0;i<product_lambda_.size();i++)
			sumLambdaProducts_ += product_lambda_[i];

		if ( tag_reaction_ == PhysicalConstants::REACTION_THIRDBODY)
		{
			sumLambdaReactants_ += 1.;

			// This includes also the explicitly reversible reactions, which are a particular case of
			// the more general reversible reactions
			if (iReversible_ == true)				sumLambdaProducts_ += 1.;
		}
	}

	void ReactionPolicy_CHEMKIN::ConvertUnits() 
	{
		ConversionFactors();
		ReactionOrders();

		if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || 
			tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF ||
			tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF) 
		{
			E_    *= e_conversion_;																	// J/kmol
			EInf_ *= e_conversion_;																	// J/kmol
			A_    *= std::pow(1.e3, -sumLambdaReactants_)/std::pow(a_conversion_, -sumLambdaReactants_);		// [m, kmol, s]
			AInf_ *= std::pow(1.e3, 1.-sumLambdaReactants_)/std::pow(a_conversion_, 1.-sumLambdaReactants_);	// [m, kmol, s]

			if (iExplicitlyReversible_ == true)
				ErrorMessage("void ReactionPolicy_CHEMKIN::ConvertUnits()", "Explicitly reversible falloff reactions are not allowed!");
		}

		else if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR || 
				 tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
				 tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR) 
		{
			E_    *= e_conversion_;																	// J/kmol
			EInf_ *= e_conversion_;																	// J/kmol
			A_    *= std::pow(1.e3, 1.-sumLambdaReactants_)/std::pow(a_conversion_, 1.-sumLambdaReactants_);		// [m, kmol, s]
			AInf_ *= std::pow(1.e3, 2.-sumLambdaReactants_)/std::pow(a_conversion_, 2.-sumLambdaReactants_);		// [m, kmol, s]

			if (iExplicitlyReversible_ == true)
				ErrorMessage("void ReactionPolicy_CHEMKIN::ConvertUnits()", "Explicitly reversible CABR reactions are not allowed!");
		}

		else if(tag_reaction_ == PhysicalConstants::REACTION_CHEBYSHEV)
		{
			chebyshev_coefficients_.push_back( std::pow(1.e3, 1.-sumLambdaReactants_)/std::pow(a_conversion_, 1.-sumLambdaReactants_) );
		}

		else if (tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
		{
			// Check if the reaction is explicitly reversible
			if (iExplicitlyReversible_ == true)
				ErrorMessage("void ReactionPolicy_CHEMKIN::ConvertUnits()", "Explicitly reversible extended falloff reactions are not allowed!");

			// Third body efficiencies
			{
				std::vector<std::string> third_body_efficiencies(third_body_efficiencies_.size());
				std::vector<std::string> third_body_indices(third_body_indices_.size());
				for (unsigned int i = 0; i < third_body_efficiencies_.size(); i++)
				{
					third_body_efficiencies[i] = boost::lexical_cast<std::string>(third_body_efficiencies_[i]);
					third_body_indices[i] = boost::lexical_cast<std::string>(third_body_indices_[i]);
				}
				extendedfalloff_coefficients_.push_back(third_body_efficiencies);
				extendedfalloff_coefficients_.push_back(third_body_indices);
			}

			// Conversion of low- and high-pressure parameters is done inside the class
			EInf_ = E_;
			AInf_ = A_;	
			betaInf_ = beta_;

			// High pressure kinetic parameters
			{
				std::vector<std::string> high_pressure_kinetic_parameters(3);
				high_pressure_kinetic_parameters[0] = boost::lexical_cast<std::string>(AInf_);
				high_pressure_kinetic_parameters[1] = boost::lexical_cast<std::string>(betaInf_);
				high_pressure_kinetic_parameters[2] = boost::lexical_cast<std::string>(EInf_);
				extendedfalloff_coefficients_.push_back(high_pressure_kinetic_parameters);
			}
			
			// Conversion factors
			{
				std::vector<std::string> conversion_factors(3);
				conversion_factors[0] = boost::lexical_cast<std::string>(std::pow(1.e3, -sumLambdaReactants_) / std::pow(a_conversion_, -sumLambdaReactants_));
				conversion_factors[1] = boost::lexical_cast<std::string>(std::pow(1.e3, 1. - sumLambdaReactants_) / std::pow(a_conversion_, 1. - sumLambdaReactants_));
				conversion_factors[2] = boost::lexical_cast<std::string>(e_conversion_);
				extendedfalloff_coefficients_.push_back(conversion_factors);
			}
		}

		else // simple or third-body
		{
			E_ *= e_conversion_;																	// J/kmol
			A_ *= std::pow(1.e3, 1.-sumLambdaReactants_)/std::pow(a_conversion_, 1.-sumLambdaReactants_);		// [m, kmol, s]

			if (iPlog_ == true)
			{
				plog_coefficients_.push_back( std::pow(1.e3, 1.-sumLambdaReactants_)/std::pow(a_conversion_, 1.-sumLambdaReactants_) );
				plog_coefficients_.push_back( e_conversion_ );
			}

			if (iExtPlog_ == true)
			{
				extendedplog_coefficients_.push_back(boost::lexical_cast<std::string>(std::pow(1.e3, 1. - sumLambdaReactants_) / std::pow(a_conversion_, 1. - sumLambdaReactants_)));
				extendedplog_coefficients_.push_back(boost::lexical_cast<std::string>(e_conversion_));
			}

			if (iExplicitlyReversible_ == true)
			{
				ERev_ *= e_conversion_;																	// J/kmol
				ARev_ *= std::pow(1.e3, 1.-sumLambdaProducts_)/std::pow(a_conversion_, 1.-sumLambdaProducts_);	// [m, kmol, s]
			}
		}
	}

	void ReactionPolicy_CHEMKIN::ConversionFactors() 
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

	void ReactionPolicy_CHEMKIN::GetReactionString(const std::vector<std::string>& list_species, std::string& line_reaction) const
	{
		// Reactants
		for(unsigned int i=0;i<reactant_nu_.size();i++)
		{
			if (reactant_nu_[i] != 1)
			{
                int precision = this->GetDecimalPlaces (reactant_nu_[i]);
				std::stringstream nu;
                nu << std::fixed << std::setprecision(precision) << reactant_nu_[i];
				line_reaction += nu.str();
			}
			line_reaction += list_species[reactant_nu_indices_[i]];
			
			if (tag_reaction_ == PhysicalConstants::REACTION_THIRDBODY)
			{
				if (i == reactant_nu_.size()-1) line_reaction += " +M ";
				else line_reaction += " + ";
			}
			else if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
			{
				if (i == reactant_nu_.size()-1 && pressureDependentSpeciesIndex_ == -1) line_reaction += " (+M) ";
				else if (i == reactant_nu_.size()-1 && pressureDependentSpeciesIndex_ != -1) line_reaction += " (+" + list_species[pressureDependentSpeciesIndex_] + ") ";
				else line_reaction += " + ";
			}
			else
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
                int precision = this->GetDecimalPlaces (product_nu_[i]);
				std::stringstream nu;
                nu << std::fixed << std::setprecision(precision) << product_nu_[i];
				line_reaction += nu.str();
			}
			line_reaction += list_species[product_nu_indices_[i]];

			if (tag_reaction_ == PhysicalConstants::REACTION_THIRDBODY)
			{
				if (i == product_nu_.size()-1) line_reaction += " +M ";
				else line_reaction += " + ";
			}
			else if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR ||
					tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF )
			{
				if (i == product_nu_.size()-1 && pressureDependentSpeciesIndex_ == -1) line_reaction += " (+M) ";
				else if (i == product_nu_.size()-1 && pressureDependentSpeciesIndex_ != -1) line_reaction += " (+" + list_species[pressureDependentSpeciesIndex_] + ") ";
				else line_reaction += " + ";
			}
			else
			{
				if (i == product_nu_.size()-1) line_reaction += " ";
				else line_reaction += " + ";
			}
		}
	}

	void ReactionPolicy_CHEMKIN::GetReactionStringCHEMKIN(	const std::vector<std::string>& list_species, std::stringstream& reaction_data) const
	{
		std::vector<bool> isReducedSpecies(list_species.size());
		for(unsigned int i=0;i<list_species.size();i++)
			isReducedSpecies[i] = true;

		GetReactionStringCHEMKIN( list_species, reaction_data, isReducedSpecies);
	}

	void ReactionPolicy_CHEMKIN::GetReactionStringCHEMKIN(const std::vector<std::string>& list_species,
                std::stringstream& reaction_data, const std::vector<bool>& isReducedSpecies) const
	{
        std::string reaction_string;
        GetReactionString(list_species, reaction_string);
        boost::erase_all(reaction_string, " ");
        reaction_data.precision(4);
        
        if(		Tag() != PhysicalConstants::REACTION_LINDEMANN_FALLOFF && 
				Tag() != PhysicalConstants::REACTION_LINDEMANN_CABR && 
				Tag() != PhysicalConstants::REACTION_TROE_FALLOFF && 
				Tag() != PhysicalConstants::REACTION_TROE_CABR && 
				Tag() != PhysicalConstants::REACTION_SRI_FALLOFF && 
				Tag() != PhysicalConstants::REACTION_SRI_CABR &&
				Tag() != PhysicalConstants::REACTION_EXTENDEDFALLOFF )
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
            
            if(IsPressureLog() == true)
            {
                reaction_data.unsetf(std::ios_base::floatfield);
                reaction_data.precision(6);
                for(unsigned int k = 0; k < plog_coefficients_.size() - 2; k++)
                {
                    if(k % 4 == 0)
                        reaction_data << " PLOG /  ";
                    reaction_data << std::showpoint << std::setw(12) << std::left << plog_coefficients_[k];
                    
                    if((k+1) % 4 == 0 || k == plog_coefficients_.size() - 3)
                        reaction_data << "/" << std::endl;
                }                    
            }

			if (IsExtendedPressureLog() == true)
			{
				// Low-pressure parameters
				OpenSMOKE::ExtendedPressureLogarithmicRateExpression extendedPressureLog;
				extendedPressureLog.Setup(extendedplog_coefficients_);
				extendedPressureLog.WriteCHEMKINOnASCIIFile(reaction_data);
			}
            
            if(IsJanevLanger() == true)
            {
                reaction_data.unsetf(std::ios_base::floatfield);
                reaction_data.precision(6);
                for(unsigned int k = 0; k < janev_langer_coefficients_.size(); k++)
                {
                    if(k % 5 == 0)
                        reaction_data << " JAN /  ";
                    reaction_data << std::showpoint << janev_langer_coefficients_[k]
                            << " ";
                    
                    if((k+1) % 5 == 0 || k == janev_langer_coefficients_.size() - 1)
                        reaction_data << "/" << std::endl;
                }
            }
            
            if(Tag() == PhysicalConstants::REACTION_THIRDBODY)
            {
                for(unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
                {
                    int third_body_index = third_body_indices_[j];
                    if(isReducedSpecies[third_body_index] == true)
                    {
                        reaction_data << list_species[third_body_index] << "/ ";
                        reaction_data.precision(2);
                        reaction_data << std::showpoint << std::fixed << std::left << third_body_efficiencies_[j] << "/ ";
                    }
                }    
                if(third_body_efficiencies().size() != 0)
                    reaction_data << std::endl;
            }
            
            if(Tag() == PhysicalConstants::REACTION_CHEBYSHEV)
            {
                reaction_data.unsetf(std::ios_base::floatfield);          
                
                if( chebyshev_temperature_limits_[0] != 300 || chebyshev_temperature_limits_[1] != 2500 )
                {
                    reaction_data << " TCHEB/ ";
                    reaction_data.precision(1);
                    reaction_data << std::showpoint <<std::fixed << chebyshev_temperature_limits_[0] << " " << chebyshev_temperature_limits_[1];
                    reaction_data << " /" << std::endl;
                }
                
                if(chebyshev_pressure_limits_[0] != 0.001 ||  chebyshev_pressure_limits_[1] != 100)
                {
                    reaction_data << " PCHEB/ ";
                    reaction_data.precision(4);
                    reaction_data << std::showpoint <<std::fixed << std::left << chebyshev_pressure_limits_[0] << " " << chebyshev_pressure_limits_[1];
                    reaction_data << " /" << std::endl;
                }
                
                unsigned int chebyshev_size =	boost::lexical_cast<unsigned int>(chebyshev_coefficients_[0]) * 
												boost::lexical_cast<unsigned int>(chebyshev_coefficients_[1]);
                        
                reaction_data.unsetf(std::ios_base::floatfield);
                reaction_data.precision(6);
                for(unsigned int k=0;k<chebyshev_size+2;k++)
                {                    
                    if(k%6 == 0)
                        reaction_data << " CHEB/ ";
                    
                    if(k < 2)
                        reaction_data << std::noshowpoint << chebyshev_coefficients_[k] << " ";
                        
                    else
                        reaction_data << std::showpoint << chebyshev_coefficients_[k] << " ";
                    
                    if((k+1)%6 == 0)
                        reaction_data << " /" << std::endl;
                    
                    if(k == chebyshev_size + 1 &&
                           (k+1)%6 != 0 )
                        reaction_data << " /" << std::endl;
                }
            }
        }

		else if (Tag() == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
		{
			// High-pressure kinetic parameters
			reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << AInf_;
			reaction_data.precision(3);
			reaction_data.width(9);
			reaction_data << std::fixed << std::right << betaInf_;
			reaction_data.precision(2);
			reaction_data << std::setw(13) << std::right << EInf_ << std::endl;
			reaction_data.unsetf(std::ios_base::floatfield);

			// Low-pressure parameters
			OpenSMOKE::ExtendedFallOff extendedFallOff;
			extendedFallOff.Setup(extendedfalloff_coefficients_);
			extendedFallOff.WriteCHEMKINOnASCIIFile(reaction_data);

			// Add third body efficiencies
			bool iThirdBody_ = false;
			for (unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
			{
				int third_body_index = third_body_indices()[j];
				if (isReducedSpecies[third_body_index] == true)
				{
					reaction_data << list_species[third_body_index] << "/ ";
					reaction_data.precision(2);
					reaction_data << std::fixed << std::showpoint << third_body_efficiencies_[j] << "/  ";
					iThirdBody_ = true;
				}
			}
			if (iThirdBody_ == true)
				reaction_data << std::endl;
		}
        
        else
        {
            reaction_data << std::setw(55) << std::left << reaction_string << " " << std::scientific << AInf_ / A_inf_conversion();
            reaction_data.precision(3);
            reaction_data.width(9);
			reaction_data <<std::fixed << std::right << betaInf_;
            reaction_data.precision(2);
            reaction_data << std::setw(13) << std::right << E_over_R_inf() * PhysicalConstants::R_cal_mol << std::endl;
            
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
            
            
            if(	Tag() == PhysicalConstants::REACTION_TROE_FALLOFF ||
				Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF ||
				Tag() == PhysicalConstants::REACTION_SRI_FALLOFF)
					reaction_data << " LOW/";
            else if(Tag() == PhysicalConstants::REACTION_TROE_CABR ||
					Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR ||
					Tag() == PhysicalConstants::REACTION_SRI_CABR)
					  reaction_data << " HIGH/";
            
            reaction_data.width(11);
            reaction_data.precision(2);
            reaction_data << std::right << std::scientific << A_ / A_conversion();
            reaction_data.precision(3);
            reaction_data.width(11);
			reaction_data <<std::fixed << std::right << beta_;
            reaction_data.precision(1);
            reaction_data << std::setw(13) << std::right << E_over_R() * PhysicalConstants::R_cal_mol << "/" << std::endl;
            
            if(		Tag() == PhysicalConstants::REACTION_TROE_FALLOFF || 
					Tag() == PhysicalConstants::REACTION_TROE_CABR)
            {
                reaction_data << "TROE/";
                reaction_data.width(11);
                reaction_data.unsetf(std::ios_base::floatfield);
                for(unsigned int j = 0; j < troe_.size(); j++)
                {
                    reaction_data.precision(4);
                    reaction_data << std::showpoint << "   " << troe_[j];
                }
                
                reaction_data << "/" << std::endl;                        
            }
            
            if(		Tag() == PhysicalConstants::REACTION_SRI_FALLOFF ||
					Tag() == PhysicalConstants::REACTION_SRI_CABR)
            {
                reaction_data << "SRI/ ";
                for(unsigned int j = 0; j < sri_.size(); j++)
                {
                    reaction_data.precision(4);
                    reaction_data << std::showpoint << "  " << sri_[j];
                }
                
                reaction_data << "/" << std::endl;
            }
            
            bool iThirdBody_ = false;
            for (unsigned int j = 0; j < third_body_efficiencies_.size(); j++)
            {
                int third_body_index = third_body_indices()[j];
                if (isReducedSpecies[third_body_index] == true)
                {
                    reaction_data << list_species[third_body_index] << "/ ";
                    reaction_data.precision(2);
                    reaction_data << std::fixed << std::showpoint << third_body_efficiencies_[j] << "/  ";
                    iThirdBody_ = true;
                }
            }
            
            if(iThirdBody_ == true)
                reaction_data << std::endl;  
        }
        
        if(IsFit1())
        {
            reaction_data.unsetf(std::ios_base::floatfield);
            reaction_data.precision(6);
            for(unsigned int k = 0; k < fit1_coefficients_.size(); k++)
            {
                reaction_data << " FIT1 /  ";
                reaction_data << std::showpoint << fit1_coefficients_[k]
                        << " ";
            }
            reaction_data << "/" << std::endl;
        }
        
        if(IsLandauTeller())
        {
            reaction_data.unsetf(std::ios_base::floatfield);
            reaction_data.precision(3);
            for(unsigned int k = 0; k < landau_teller_coefficients_.size(); k++)
            {
                reaction_data << " LT /  ";
                reaction_data << std::showpoint << landau_teller_coefficients_[k]
                        << " ";
            }
            reaction_data << "/" << std::endl;
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

	double ReactionPolicy_CHEMKIN::A_conversion() const
	{
		double conversion_factor = 0;
		if (	tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || 
				tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF ||
				tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF)
		{
			conversion_factor = std::pow(1.e3, -sumLambdaReactants_)/std::pow(a_conversion_, -sumLambdaReactants_);
		}
            
		else if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR || 
				tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
				tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)
		{
			conversion_factor = std::pow(1.e3, 1.-sumLambdaReactants_)/std::pow(a_conversion_, 1.-sumLambdaReactants_);
		}
            
		else if(tag_reaction_ == PhysicalConstants::REACTION_CHEBYSHEV)
		{
			conversion_factor = 1;
			//ErrorMessage("double A_conversion() const","Conversion of A factor not allowed for Chebyshev reactions");
		}
		else if (tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
		{
			conversion_factor = std::pow(1.e3, -sumLambdaReactants_) / std::pow(a_conversion_, -sumLambdaReactants_);
		}
		else // simple or third-body
		{
			conversion_factor = std::pow(1.e3, 1.-sumLambdaReactants_)/std::pow(a_conversion_, 1.-sumLambdaReactants_);
		}
            
		return conversion_factor;
	}
        
	double ReactionPolicy_CHEMKIN::A_inf_conversion() const
	{
		double conversion_factor = 0;
		if (	tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || 
				tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF ||
				tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF)
		{
			conversion_factor = std::pow(1.e3, 1.-sumLambdaReactants_)/std::pow(a_conversion_, 1.-sumLambdaReactants_);
		}
            
		else if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR || 
				tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
				tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)
		{
			conversion_factor = std::pow(1.e3, 2.-sumLambdaReactants_)/std::pow(a_conversion_, 2.-sumLambdaReactants_);
		}
            
		else if(tag_reaction_ == PhysicalConstants::REACTION_CHEBYSHEV)
		{
			conversion_factor = 1;
			//ErrorMessage("double A_inf_conversion() const","Conversion of A_inf factor not allowed for Chebyshev reactions");
		}

		else if (tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
		{
			conversion_factor = std::pow(1.e3, 1. - sumLambdaReactants_) / std::pow(a_conversion_, 1. - sumLambdaReactants_);
		}

		else // simple or third-body
		{
			ErrorMessage("double A_inf_conversion() const", "Conversion of A_inf factor not allowed for simple"
			" and third body reactions");
		}
            
		return conversion_factor;
	}

	double ReactionPolicy_CHEMKIN::Arev_conversion() const
    {
        double conversion_factor = 0;
        if(	tag_reaction_ == PhysicalConstants::REACTION_SIMPLE || 
            tag_reaction_ == PhysicalConstants::REACTION_THIRDBODY)
        {
            if(iExplicitlyReversible_ == true)
                conversion_factor = std::pow(1.e3, 1.-sumLambdaProducts_)/std::pow(a_conversion_, 1.-sumLambdaProducts_);
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

	void ReactionPolicy_CHEMKIN::WriteAdditionalDataOnASCIIFile(std::ostream& fOut) const
	{
		if (tag_reaction_ == PhysicalConstants::REACTION_SIMPLE)
		{
			if (iPlog_ == true)
			{
				fOut << plog_coefficients_.size() << std::endl;
				for(unsigned int j=0;j<plog_coefficients_.size();j++)
					fOut << plog_coefficients_[j] << " ";
				fOut << std::endl;
			}
			else if (iExtPlog_ == true)
			{
				fOut << extendedplog_coefficients_.size() << std::endl;
				for (unsigned int j = 0; j<extendedplog_coefficients_.size(); j++)
					fOut << extendedplog_coefficients_[j] << " ";
				fOut << std::endl;
			}
			else if (iJan_ == true)
			{
				fOut << janev_langer_coefficients_.size() << std::endl;
				for(unsigned int j=0;j<janev_langer_coefficients_.size();j++)
					fOut << janev_langer_coefficients_[j] << " ";
				fOut << std::endl;
			}
			else if (iFit1_ == true)
			{
				fOut << fit1_coefficients_.size() << std::endl;
				for(unsigned int j=0;j<fit1_coefficients_.size();j++)
					fOut << fit1_coefficients_[j] << " ";
				fOut << std::endl;
			}
			else if (iLandauTeller_ == true)
			{
				fOut << landau_teller_coefficients_.size() << std::endl;
				for(unsigned int j=0;j<landau_teller_coefficients_.size();j++)
					fOut << landau_teller_coefficients_[j] << " ";
				fOut << std::endl;
			}
		}

		else if (tag_reaction_ == PhysicalConstants::REACTION_CHEBYSHEV)
		{
			fOut << chebyshev_coefficients_.size() << std::endl;
			for(unsigned int j=0;j<chebyshev_coefficients_.size();j++)
				fOut << chebyshev_coefficients_[j] << " ";
			fOut << std::endl;
			fOut << chebyshev_pressure_limits_[0] << " " << chebyshev_pressure_limits_[1] << std::endl;
			fOut << chebyshev_temperature_limits_[0] << " " << chebyshev_temperature_limits_[1] << std::endl;
		}

		else if (iExtLow_ == true)
		{
			fOut << extendedfalloff_coefficients_.size() << std::endl;
			for (unsigned int j = 0; j < extendedfalloff_coefficients_.size(); j++)
			{
				fOut << extendedfalloff_coefficients_[j].size() << std::endl;
				for (unsigned int k = 0; k < extendedfalloff_coefficients_[j].size(); k++)
					fOut << extendedfalloff_coefficients_[j][k] << " ";
				fOut << std::endl;
			}
		}
		
	}

	void ReactionPolicy_CHEMKIN::WriteShortSummary(std::ostream& fOut, const std::vector<std::string>& list_species) const
	{
		Eigen::VectorXd reverse_parameters;
		WriteShortSummary(fOut, list_species, reverse_parameters);
	}

	void ReactionPolicy_CHEMKIN::WriteShortSummary(std::ostream& fOut, const std::vector<std::string>& list_species, const Eigen::VectorXd& reverse_parameters) const
	{
		std::string line_reaction;
		GetReactionString(list_species, line_reaction);

		fOut << line_reaction << std::endl;

		// Pressure dependent reactions
		if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF ||
			tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR || 
			tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)	
		{
			fOut << std::setw(9) << " ";
			fOut << TagASCII() << std::endl;
		
			fOut << std::setw(9) << " ";
			fOut << std::setw(9) << std::left << "k0:";
			fOut << std::scientific << std::setprecision(6) << std::right << A_ << "\t";
			fOut << std::setw(8)    << std::setprecision(2) << std::fixed << std::right << beta_;
			fOut << std::setw(14)	  << std::setprecision(2) << std::fixed << std::right << E_/Conversions::J_from_kcal   << std::endl;

			fOut << std::setw(9) << " ";
			fOut << std::setw(9) << std::left << "kInf:";
			fOut << std::scientific << std::setprecision(6) << std::right << AInf_	<< "\t"; 
			fOut << std::setw(8)    << std::setprecision(2) << std::fixed << std::right << betaInf_;
			fOut << std::setw(14)   << std::setprecision(2) << std::fixed << std::right << EInf_/Conversions::J_from_kcal << std::endl;
		}
		else if (tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
		{
			// Do nothing
		}
		else	// Conventional reactions
		{
			fOut << std::setw(9)  << " ";
			fOut << std::setw(9)  << std::left << "k:";
			fOut << std::scientific	<< std::setprecision(6) << std::right << A_ << "\t";
			fOut << std::setw(8) << std::setprecision(2) << std::fixed << std::right << beta_;
			fOut << std::setw(14) << std::setprecision(2) << std::fixed << std::right << E_/Conversions::J_from_kcal << std::endl;

			if (reverse_parameters.size() > 0)
			{
				fOut << std::setw(9) << " ";
				fOut << std::setw(9) << std::left << "kRev:";
				if (reverse_parameters.size() == 3)
				{
					fOut << std::scientific << std::setprecision(6) << std::right << std::exp(reverse_parameters(0)) << "\t";
					fOut << std::setw(8) << std::setprecision(2) << std::fixed << std::right << reverse_parameters(2);
					fOut << std::setw(14) << std::setprecision(2) << std::fixed << std::right << reverse_parameters(1)/Conversions::J_from_kcal << std::endl;
				}
				else if(reverse_parameters.size() == 2)
				{
					fOut << std::scientific << std::setprecision(6) << std::right << std::exp(reverse_parameters(0)) << "\t";
					fOut << std::setw(8) << std::setprecision(2) << std::fixed << std::right << 0.;
					fOut << std::setw(14) << std::setprecision(2) << std::fixed << std::right << reverse_parameters(1)/Conversions::J_from_kcal << std::endl;
				}
			}
		}

		// Troe form
		if (tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR )	
		{
			fOut << std::setw(9) << " "; fOut << "Troe Parameters" << std::endl;
			fOut << std::setw(9) << " "; fOut << std::scientific << "a    " << troe_[0] << std::endl;
			fOut << std::setw(9) << " "; fOut << std::scientific << "T*** " << troe_[1] << std::endl;
			fOut << std::setw(9) << " "; fOut << std::scientific << "T*   " << troe_[2] << std::endl;
			
			if (troe_.size()==3)	
			{	fOut << std::setw(9) << " "; fOut << std::scientific << "T**  " << 0. << std::endl;}

			if (troe_.size()==4)	
			{	fOut << std::setw(9) << " "; fOut << std::scientific << "T**  " << troe_[3] << std::endl;}
		}

		// SRI form
		if (tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR )	
		{
			fOut << std::setw(9) << " "; fOut << "SRI Parameters" << std::endl;
			fOut << std::setw(9) << " "; fOut << std::scientific << "a  " << sri_[0] << std::endl;
			fOut << std::setw(9) << " "; fOut << std::scientific << "b  " << sri_[1] << std::endl;
			fOut << std::setw(9) << " "; fOut << std::scientific << "c  " << sri_[2] << std::endl;
			
			if (sri_.size()==3)	
			{
				fOut << std::setw(9) << " "; fOut << std::scientific << "d  " << 1. << std::endl;
				fOut << std::setw(9) << " "; fOut << std::scientific << "e  " << 0. << std::endl;
			}

			if (sri_.size()==5)	
			{
				fOut << std::setw(9) << " "; fOut << std::scientific << "d  " << sri_[3] << std::endl;
				fOut << std::setw(9) << " "; fOut << std::scientific << "e  " << sri_[4] << std::endl;
			}
		}

		// Landau-Teller
		if (iLandauTeller_ == true)
		{
			fOut << std::setw(9) << " "; fOut << "Landau-Teller Parameters" << std::endl;
			fOut << std::setw(9) << " "; fOut << "B  " << std::scientific << landau_teller_coefficients_[0] << std::endl;
			fOut << std::setw(9) << " "; fOut << "C  " << std::scientific << landau_teller_coefficients_[1] << std::endl;
		}

		// Janev-Langer (TODO units)
		if (iJan_ == true)
		{	
			fOut << std::setw(9) << " "; 
			fOut << "Janev-Langer Parameters" << std::endl;
			for (unsigned int j=0;j<9;j++)
			{	fOut << std::setw(9) << " "; fOut << "b" << j+1 << "  " << std::scientific << janev_langer_coefficients_[j] << std::endl;}
		}

		// Power-Series (FIT1)
		if (iFit1_ == true)
		{	
			fOut << std::setw(9) << " "; 
			fOut << "Power-Series (FIT1) Parameters" << std::endl;
			for (unsigned int j=0;j<4;j++)
			{	fOut << std::setw(9) << " "; fOut << "b" << j+1 << " " << std::scientific << fit1_coefficients_[j] << std::endl;}
		}

		// Chebyshev Polynomials
		if (iChebyshev_ == true)
		{	
			OpenSMOKE::ChebyshevPolynomialRateExpression chebyshev;
			chebyshev.Setup(chebyshev_coefficients_, chebyshev_pressure_limits_, chebyshev_temperature_limits_);
			chebyshev.WriteShortSummaryOnASCIIFile(fOut);	
		}

		// Logarithmic-Pressure Dependence
		if (iPlog_ == true)
		{	
			OpenSMOKE::PressureLogarithmicRateExpression pressureLogarithmic;
			pressureLogarithmic.Setup(plog_coefficients_);
			pressureLogarithmic.WriteShortSummaryOnASCIIFile(fOut);
		}

		// Logarithmic-Pressure Dependence
		if (iExtPlog_ == true)
		{
			OpenSMOKE::ExtendedPressureLogarithmicRateExpression extendedPressureLogarithmic;
			extendedPressureLogarithmic.Setup(extendedplog_coefficients_);
			extendedPressureLogarithmic.WriteShortSummaryOnASCIIFile(fOut);
		}

		// Extended-Falloff Reactions
		if (iExtLow_ == true)
		{
			OpenSMOKE::ExtendedFallOff extendedFallOffReaction;
			extendedFallOffReaction.Setup(extendedfalloff_coefficients_);
			extendedFallOffReaction.WriteShortSummaryOnASCIIFile(fOut);
		}

		// Single species acting as third body
		if (pressureDependentSpeciesIndex_ != -1)
		{
			fOut << std::setw(9) << " ";
			fOut << list_species[pressureDependentSpeciesIndex_] << " is acting as the third body (not the total concentration)";
		}

		// Third Body Efficiencies
		if (third_body_efficiencies_.size() > 0)
		{
			fOut << std::setw(9) << " "; fOut << "Third body efficiencies" << std::endl;

			for(unsigned int i=0;i<third_body_efficiencies_.size();i++)
			{
				fOut << std::setw(9) << " ";
				fOut << std::setw(20) << std::left << list_species[third_body_indices_[i]] << "\tenhanced by\t" << std::fixed 
					 << std::setw(8) << std::right << std::setprecision(3) << third_body_efficiencies_[i] <<std::endl;
			}
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
	
	}

	double ReactionPolicy_CHEMKIN::GetForwardConversionFactor() const
	{
		if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)
						return std::pow(1000.,  1.-(1.-sumLambdaReactants_)-1.);
		else
			 return std::pow(1000., sumLambdaReactants_-1.);
	}

	double ReactionPolicy_CHEMKIN::GetBackwardConversionFactor() const
	{
		return std::pow(1000., sumLambdaProducts_);
	}

	void ReactionPolicy_CHEMKIN::WriteSummary(std::ofstream& fOut, const std::vector<std::string>& list_species, const unsigned int index) const
	{
		std::string line_reaction;
		GetReactionString(list_species, line_reaction);
	
		fOut << "================================================================================================================================" << std::endl;
		fOut << " KINETIC DATA - REACTION  " << index << " " << TagASCII() << std::endl;
		fOut << "  " << line_reaction << std::endl; 
		fOut << "================================================================================================================================" << std::endl;
		fOut << " Change in moles in the reaction = " << sumNuProducts_ - sumNuReactants_ << std::endl;
	
		fOut << " Reaction order (forward)  = " << std::setprecision(3) << sumLambdaReactants_ << std::endl;
		if (iReversible_ == true)
			fOut << " Reaction order (backward) = " << std::setprecision(3) << sumLambdaProducts_ << std::endl;

		if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
					 tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)
		{
			std::string line_kForwardStringSI_High;
			std::string line_kForwardStringCGS_High;
			std::string line_kForwardStringSI_Low;
			std::string line_kForwardStringCGS_Low;

			if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF ||
				tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF )
			{
				OpenSMOKE_Utilities::GetKineticConstantString(A_, beta_, E_, 1.-(-sumLambdaReactants_), line_kForwardStringSI_Low, line_kForwardStringCGS_Low);
				OpenSMOKE_Utilities::GetKineticConstantString(AInf_, betaInf_, EInf_, 1.-(1.-sumLambdaReactants_), line_kForwardStringSI_High, line_kForwardStringCGS_High);
			}
		
			if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
				tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR )
			{
				OpenSMOKE_Utilities::GetKineticConstantString(A_, beta_, E_, 1.-(1.-sumLambdaReactants_), line_kForwardStringSI_Low, line_kForwardStringCGS_Low);
				OpenSMOKE_Utilities::GetKineticConstantString(AInf_, betaInf_, EInf_, 1.-(2.-sumLambdaReactants_), line_kForwardStringSI_High, line_kForwardStringCGS_High);
			}

			fOut << " Low pressure:  " << line_kForwardStringSI_Low   << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsSI(1.-(-sumLambdaReactants_))  << std::endl;
			fOut << " Low pressure:  " << line_kForwardStringCGS_Low  << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsCGS(1.-(-sumLambdaReactants_)) << std::endl;
			fOut << " High pressure: " << line_kForwardStringSI_High  << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsSI(1.-(1.-sumLambdaReactants_))  << std::endl;
			fOut << " High pressure: " << line_kForwardStringCGS_High << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsCGS(1.-(1.-sumLambdaReactants_)) << std::endl;
		}
		else if (tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
		{
			{
				std::string line_kForwardStringSI_High;
				std::string line_kForwardStringCGS_High;
				OpenSMOKE_Utilities::GetKineticConstantString(AInf_*A_inf_conversion(), betaInf_, EInf_*e_conversion_, 1. - (1. - sumLambdaReactants_), line_kForwardStringSI_High, line_kForwardStringCGS_High);
				fOut << " High pressure: " << line_kForwardStringSI_High << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsSI(1. - (1. - sumLambdaReactants_)) << std::endl;
				fOut << " High pressure: " << line_kForwardStringCGS_High << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsCGS(1. - (1. - sumLambdaReactants_)) << std::endl;
			}
			
			OpenSMOKE::ExtendedFallOff extendedFallOff;
			extendedFallOff.Setup(extendedfalloff_coefficients_);
			for (int i=0;i<extendedFallOff.number_of_species();i++)
			{
				std::string line_kForwardStringSI_Low;
				std::string line_kForwardStringCGS_Low;
				OpenSMOKE_Utilities::GetKineticConstantString(extendedFallOff.A0(i), extendedFallOff.Beta0(i), extendedFallOff.E0(i), 1. - (-sumLambdaReactants_), line_kForwardStringSI_Low, line_kForwardStringCGS_Low);
				fOut << " Low pressure (" << extendedFallOff.species(i) << "): " << line_kForwardStringSI_Low << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsSI(1. - (-sumLambdaReactants_)) << std::endl;
				fOut << " Low pressure (" << extendedFallOff.species(i) << "): " << line_kForwardStringCGS_Low << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsCGS(1. - (-sumLambdaReactants_)) << std::endl;
				extendedFallOff.WriteAdditionalParameters(fOut, i);
			}
		}
		else if (tag_reaction_ == PhysicalConstants::REACTION_CHEBYSHEV)
		{
			ChebyshevPolynomialRateExpression chebyshev;
			chebyshev.Setup(chebyshev_coefficients_, chebyshev_pressure_limits_, chebyshev_temperature_limits_);
			chebyshev.WriteStatus(fOut);
		}
		else
		{
			std::string line_kForwardStringSI;
			std::string line_kForwardStringCGS;
			OpenSMOKE_Utilities::GetKineticConstantString(A_, beta_, E_, sumLambdaReactants_, line_kForwardStringSI, line_kForwardStringCGS);
			fOut << " " << line_kForwardStringSI  << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsSI(sumLambdaReactants_)  << " and [J/kmol]" << std::endl;
			fOut << " " << line_kForwardStringCGS << "  " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsCGS(sumLambdaReactants_) << " and [cal/mol]" << std::endl;
		}
	
		if (iReversible_ == true)
			fOut << " Reverse reaction units: " << OpenSMOKE_Utilities::GetUnitsOfKineticConstantsSI(sumLambdaProducts_) << " or " 
												<< OpenSMOKE_Utilities::GetUnitsOfKineticConstantsCGS(sumLambdaProducts_) << std::endl;

		if (tag_reaction_ == PhysicalConstants::REACTION_THIRDBODY ||
			tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR ||
			tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR ||
			tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR ||
			tag_reaction_ == PhysicalConstants::REACTION_CHEBYSHEV ||
			tag_reaction_ == PhysicalConstants::REACTION_EXTENDEDFALLOFF)
		{
			fOut << " Third-body efficiencies:  " << std::endl;
			for(unsigned int i=0;i<third_body_indices_.size();i++)
				fOut << " " << std::setw(16) << std::left << list_species[third_body_indices_[i]] << third_body_efficiencies_[i] << std::endl;
		}

		if (tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR)
		{
			fOut << " Troe parameters:  " << std::endl;
			fOut << std::setw(6) << std::left << " a:" << troe_[0] << std::endl;
			fOut << std::setw(6) << std::left << " T***:" << troe_[1] << std::endl;
			fOut << std::setw(6) << std::left << " T*:"   << troe_[2] << std::endl;
			if (troe_.size() == 4)
				fOut << std::setw(6) << std::left << " T**:"  << troe_[3] << std::endl;
		}

		if (tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)
		{
			fOut << " SRI parameters:  " << std::endl;
			fOut << std::setw(6) << std::left << " a:" << sri_[0] << std::endl;
			fOut << std::setw(6) << std::left << " b:" << sri_[1] << std::endl;
			fOut << std::setw(6) << std::left << " c:" << sri_[2] << std::endl;
			if (sri_.size() == 5)
			{
				fOut << std::setw(6) << std::left << " d:"  << sri_[3] << std::endl;
				fOut << std::setw(6) << std::left << " e:"  << sri_[4] << std::endl;
			}
			else
			{
				fOut << std::setw(6) << std::left << " d:"  << 1. << std::endl;
				fOut << std::setw(6) << std::left << " e:"  << 0. << std::endl;
			}
		}
	}

	void ReactionPolicy_CHEMKIN::WriteKineticsDataOnASCIIFileOldStyle(std::ofstream &fOutput) const
	{
		// 1. Reversible vs non/reversible reactions
		if (iReversible_ == true)	fOutput << 1 << std::endl;
		else                        fOutput << 0 << std::endl;

		// 2. Type of reaction
		if (iThirdBody_	== false && iLandauTeller_ == true)		    fOutput << 100 << std::endl;
		else if (iThirdBody_	== true  && iLandauTeller_ == true)	fOutput << 10  << std::endl;
		else if (iThirdBody_ == false && iJan_ == true)				fOutput << 110 << std::endl;
		else if (iThirdBody_ == true  && iJan_ == true)				fOutput << 11  << std::endl;
		else if (iThirdBody_ == false && iFit1_ == true)			fOutput << 120 << std::endl;
		else if (iThirdBody_ == true  && iFit1_ == true)			fOutput << 12  << std::endl;
		else if (iThirdBody_ == false && iChebyshev_ == true)		fOutput << 130 << std::endl;
	//	else if (iThirdBody_ == true  && iChebyshev_ == true)		ErrorMessage("CHEBISHEV and THREE-BODY option are mutually exclusive...");
		else if (iThirdBody_ == false && iPlog_ == true)			fOutput << 140 << std::endl;
		else if (iExtPlog_ == true)									ErrorMessage("WriteKineticsDataOnASCIIFileOldStyle", "PLOGMX and PLOGSP reactions are not available in the old format...");
		else if (iExtLow_ == true)									ErrorMessage("WriteKineticsDataOnASCIIFileOldStyle", "LOWMX and LOWSP reactions are not available in the old format...");
		//	else if (iThirdBody_ == true  && iPlog_ == true)				ErrorMessage("PLOG and THREE-BODY option are mutually exclusive...");
	//	else if (iThirdBody_ == false && iCollisionEfficiency_ == true)	fOutput << 150;
	//	else if (iThirdBody_ == true  && iCollisionEfficiency_ == true)	ErrorMessage("COLLEFF and THREE-BODY option are mutually exclusive...");

		else if(iLow_ == true)	// Fall-Off Reactions (2,3,4)
		{
			if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF)	fOutput << 2 << std::endl;
			if (tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF)		fOutput << 3 << std::endl;
			if (tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF)		fOutput << 4 << std::endl;
		}
		else if(iHigh_ == true)	// Chemically-Activated Bimolecular Reactions (5,6,7)
		{
			if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR)	fOutput << 5 << std::endl;
			if (tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR)			fOutput << 6 << std::endl;
			if (tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)			fOutput << 7 << std::endl;
		}
		else
		{
			if (iThirdBody_ == true)				fOutput << 1 << std::endl;
			if (iThirdBody_ == false)				fOutput << 0 << std::endl;
		}

		// 3. Third boby efficiencies
		fOutput << third_body_efficiencies_.size() << std::endl;
		if (third_body_efficiencies_.size() != 0)	
		{
			for (unsigned int i=0;i<third_body_efficiencies_.size();i++)
			{
				fOutput << third_body_indices_[i]+1 << std::endl;
				fOutput << third_body_efficiencies_[i] << std::endl;
			}
		}

		// 4. Kinetic parameters
		fOutput << A_ << std::endl;
		fOutput << beta_ << std::endl;
		fOutput << E_/Conversions::J_from_kcal << std::endl;

		if (iLow_ == true || iHigh_ == true)	
		{
			fOutput << AInf_ << std::endl;
			fOutput << betaInf_ << std::endl;
			fOutput << EInf_/Conversions::J_from_kcal << std::endl;
		}

		// 5. Pressure dependent reactions
		if(iLow_ == true || iHigh_==true)
		{
			if (tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || 
				tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR)
			{
				fOutput << troe_[0] << std::endl;
				fOutput << troe_[1] << std::endl;
				fOutput << troe_[2] << std::endl; 
				if (troe_.size()==3)	fOutput << 0. << std::endl << 0. << std::endl;
				if (troe_.size()==4)	fOutput << troe_[3] << std::endl << 0. << std::endl;
			}
		
			else if (tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || 
					 tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR)
			{
				fOutput << sri_[0] << std::endl;
				fOutput << sri_[1] << std::endl;
				fOutput << sri_[2] << std::endl; 
				if (sri_.size()==3)	fOutput << 1. << std::endl << 0. << std::endl;
				if (sri_.size()==5)	fOutput << sri_[3] << std::endl << sri_[4] << std::endl;
			}

			else if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || 
					 tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR)
			{
				fOutput << 0. << std::endl << 0. << std::endl
						<< 0. << std::endl << 0. << std::endl << 0. << std::endl;
			}
		}

		// 6. Additional data
		{
			// Landau-Teller
			if (iLandauTeller_ == true)
				fOutput << landau_teller_coefficients_[0] << std::endl << landau_teller_coefficients_[1] << std::endl;

			// Janev-Langer (TODO units)
			if (iJan_ == true)
				for (unsigned int j=0;j<9;j++)
					fOutput << janev_langer_coefficients_[j] << std::endl;

			// Power-Series
			if (iFit1_ == true)
				for (unsigned int j=0;j<4;j++)
					fOutput << fit1_coefficients_[j] << std::endl;

			// Chebyshev
			if (iChebyshev_ == true)
			{	
				OpenSMOKE::ChebyshevPolynomialRateExpression chebyshev;
				chebyshev.Setup(chebyshev_coefficients_, chebyshev_pressure_limits_, chebyshev_temperature_limits_);
				chebyshev.WriteOnASCIIFileOldStyle(fOutput);
			}

			// Logarithmic-Pressure Dependence
			if (iPlog_ == true)
			{	
				OpenSMOKE::PressureLogarithmicRateExpression pressureLogarithmic;
				pressureLogarithmic.Setup(plog_coefficients_);
				pressureLogarithmic.WriteOnASCIIFileOldStyle(fOutput);
			}

			// Extended Logarithmic-Pressure Dependence
			if (iExtPlog_ == true)
			{
				ErrorMessage("WriteKineticsDataOnASCIIFileOldStyle", "PLOGMX and PLOGSP reactions are not available in the old format...");
			}

			// Efficiency of Collision Frequency
			/*if (iCollisionEfficiency == true)
			{
				if (transport.IsActivated() == true)
				{
					if (indexDirect.Size()>2)	ErrorMessage("COLLEFF option can be applied only to bimolecular reactions!");
					if (indexDirect.Size()==1)
						if (nuDirect[1] != 2.)	ErrorMessage("COLLEFF option can be applied only to bimolecular reactions!");
					if (indexDirect.Size()==2)
						if (nuDirect[1] != 1. || nuDirect[2] != 1.)	ErrorMessage("COLLEFF option can be applied only to bimolecular reactions!");

					if (indexDirect.Size()==1)
					{
						dAB = transport.ReducedDiameter(indexDirect[1], indexDirect[1])*1e-10;
						WAB = thermo.ReducedMolecularWeight(nameDirect[1], nameDirect[1]);
					}
					if (indexDirect.Size()==2)
					{
						dAB = transport.ReducedDiameter(indexDirect[1], indexDirect[2])*1.e-10;
						WAB = thermo.ReducedMolecularWeight(nameDirect[1], nameDirect[2]);
					}

					kStar_collision_frequency = Constants::Nav_kmol*(dAB*dAB)*std::sqrt(8.*Constants::pi*Constants::R_J_kmol/WAB);
					outputFile << kStar_collision_frequency;
				}
				else
					ErrorMessage("COLLEFF option can be used only if transport properties are available...");
			}
			*/
		}
	}

	void ReactionPolicy_CHEMKIN::WriteThirdBodyParametersOnASCIIFile(std::ostream &fOutput) const
	{
			fOutput << third_body_indices_.size() << " ";
			for (unsigned int j=0;j<third_body_indices_.size();j++)
				fOutput << third_body_indices_[j]+1 << " ";
			fOutput << std::endl;
			fOutput << third_body_efficiencies_.size() << " ";
			for (unsigned int j=0;j<third_body_efficiencies_.size();j++)
				fOutput << third_body_efficiencies_[j] << " ";
			fOutput << std::endl;
	}

	void ReactionPolicy_CHEMKIN::WritePressureDependentParametersOnASCIIFile(std::ostream &fOutput) const
	{
		if (tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_LINDEMANN_CABR )
			fOutput << "lindemann" << std::endl;
		else if (tag_reaction_ == PhysicalConstants::REACTION_TROE_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_TROE_CABR )
		{
			fOutput << "troe " << troe_.size() << " ";
			for (unsigned int j=0;j<troe_.size();j++)
				fOutput << troe_[j] << " ";
			fOutput << std::endl;
		}
		else if (tag_reaction_ == PhysicalConstants::REACTION_SRI_FALLOFF || tag_reaction_ == PhysicalConstants::REACTION_SRI_CABR )
		{
			fOutput << "sri " << sri_.size() << " ";
			for (unsigned int j=0;j<sri_.size();j++)
				fOutput << sri_[j] << " ";
			fOutput << std::endl;
		}

		if(pressureDependentSpeciesIndex_ != -1)
		{
			fOutput << "species " << pressureDependentSpeciesIndex_+1 << std::endl;
		}
		else
		{
			fOutput << "thirdbody" << std::endl;
			WriteThirdBodyParametersOnASCIIFile(fOutput);
		}
	}
        
        int ReactionPolicy_CHEMKIN::GetDecimalPlaces(const double nu) const
        {
          int threshold = 15;
          
          double coefficient = nu;
          
          for(int i = 0; i < threshold; i++)
            {
              if(std::fabs(coefficient) / std::fabs(double(int(coefficient))) == 1)
                return i;
              
              else
                coefficient *= 10.;
            }
          
          return threshold;          
          
        }

}



