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
	ReactionPolicy_Solid_CHEMKIN::ReactionPolicy_Solid_CHEMKIN() 
	{
		SetDefaultUnits();
	}

	void ReactionPolicy_Solid_CHEMKIN::SetDefaultUnits()
	{
		iReversible_ = false;
		iExplicitlyReversible_ = false;
		iFord_ = false;
		iRord_ = false;
		iDuplicate_ = false;
		delta_nu_gas_ = 0.;
	}

	ReactionPolicy_Solid_CHEMKIN::ReactionPolicy_Solid_CHEMKIN(const ReactionPolicy_Solid_CHEMKIN& orig) 
	{
	}

	ReactionPolicy_Solid_CHEMKIN::~ReactionPolicy_Solid_CHEMKIN() {
	}

	PhysicalConstants::TAG_REACTION_SURFACE ReactionPolicy_Solid_CHEMKIN::Tag() const
	{
		return tag_reaction_;
	}

	std::string ReactionPolicy_Solid_CHEMKIN::TagASCII() const
	{
		if (tag_reaction_ == PhysicalConstants::REACTION_SURFACE_SIMPLE)
			return "Simple";
		// TODO
		else
			return "Unknown reaction";
	}

	bool ReactionPolicy_Solid_CHEMKIN::SetUnits(const PhysicalConstants::UNITS_REACTION a_units, const PhysicalConstants::UNITS_REACTION e_units, const PhysicalConstants::UNITS_REACTION_COMPOSITION composition_units)
	{
		a_units_ = a_units;
		e_units_ = e_units;
		composition_units_ = composition_units;
		
		return true;
	}


	bool ReactionPolicy_Solid_CHEMKIN::ReadReactionFromStrings(	const std::vector<std::string>& lines, const std::map<std::string, unsigned int>& map_of_species,
																	const unsigned int number_of_gas_species, const unsigned int number_of_solid_species)
	{
		number_of_gas_species_ = number_of_gas_species;
		number_of_solid_species_ = number_of_solid_species;

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
					if (boost::iequals(keyword, "DUP") || boost::iequals(keyword, "DUPLICATE"))
					{
						if (iDuplicate_ == true)
						{
							std::cout << "The DUPLICATE (or DUP) option is used more than once!" << std::endl;
							return false;
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
		}

		if (iReversible_ == true)
		{
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
		
		if (iReversible_ == true && iExplicitlyReversible_ == false)
		{
			double number_of_gas_species_reactant_side	= 0;
			double number_of_gas_species_product_side	= 0;

			for(unsigned int i=0;i<reactant_nu_.size();i++)
			{
				if (reactant_nu_indices_[i] < number_of_gas_species_)	
					number_of_gas_species_reactant_side += reactant_nu_[i];
			}

			for(unsigned int i=0;i<product_nu_.size();i++)
			{
				if (product_nu_indices_[i] < number_of_gas_species_)	
					number_of_gas_species_product_side += product_nu_[i];
			}

			delta_nu_gas_ = number_of_gas_species_product_side - number_of_gas_species_reactant_side;
		}

		tag_reaction_ = PhysicalConstants::REACTION_SURFACE_SIMPLE;

		ConvertUnits();

		return true;
	}

	void ReactionPolicy_Solid_CHEMKIN::ReactionOrders() 
	{

		sumLambdaGasReactants_  = 0.;
		sumLambdaSolidReactants_ = 0.;
		for(unsigned int i=0;i<reactant_lambda_.size();i++)
		{
			if (reactant_lambda_indices_[i] < number_of_gas_species_)								sumLambdaGasReactants_ += reactant_lambda_[i];
			else if (reactant_lambda_indices_[i] < number_of_gas_species_+number_of_solid_species_)	sumLambdaSolidReactants_ += reactant_lambda_[i];
		}

		sumLambdaGasProducts_  = 0.;
		sumLambdaSolidProducts_ = 0.;
		for(unsigned int i=0;i<product_lambda_.size();i++)
		{
			if (product_lambda_indices_[i] < number_of_gas_species_)								sumLambdaGasProducts_  += product_lambda_[i];
			else if (product_lambda_indices_[i] < number_of_gas_species_+number_of_solid_species_)	sumLambdaSolidProducts_ += product_lambda_[i];
		}
	}

	void ReactionPolicy_Solid_CHEMKIN::ConvertUnits() 
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
		{			
			if (composition_units_ == PhysicalConstants::UNITS_STD)
			{
				
				A_ *= std::pow(1.e3, 1.-sumLambdaGasReactants_-sumLambdaSolidReactants_)/
					std::pow(a_conversion_, 1. - sumLambdaGasReactants_ - sumLambdaSolidReactants_);		// [m, kmol, s]
			}
			else 
			{
				double conversion_factor = 1.;

				if (composition_units_ == PhysicalConstants::UNITS_ATM)				conversion_factor = 101325.;
				else if (composition_units_ == PhysicalConstants::UNITS_BAR)		conversion_factor = 100000.;
				else if (composition_units_ == PhysicalConstants::UNITS_PASCALS)	conversion_factor = 1.;
				else if (composition_units_ == PhysicalConstants::UNITS_TORR)		conversion_factor = 101325./760.;
				else if (composition_units_ == PhysicalConstants::UNITS_DYNES)		conversion_factor = 0.1;

				A_ *=	pow(conversion_factor, -sumLambdaGasReactants_)*std::pow(10., 1.-sumLambdaSolidReactants_) /
						pow(a_conversion_, 1.-sumLambdaSolidReactants_);									// [m, kmol, s];
			}

			if (iExplicitlyReversible_ == true)
				ARev_ *=	std::pow(1.e3, 1. - sumLambdaGasProducts_ - sumLambdaSolidProducts_) /
							std::pow(a_conversion_, 1. - sumLambdaGasProducts_ - sumLambdaSolidProducts_);	// [m, kmol, s]
		}
	}

	void ReactionPolicy_Solid_CHEMKIN::ConversionFactors() 
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

	void ReactionPolicy_Solid_CHEMKIN::GetReactionString(const std::vector<std::string>& list_species, std::string& line_reaction) const
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

	void ReactionPolicy_Solid_CHEMKIN::GetReactionStringCHEMKIN(	const std::vector<std::string>& list_species, std::stringstream& reaction_data) const
	{
		std::vector<bool> isReducedSpecies(list_species.size());
		for(unsigned int i=0;i<list_species.size();i++)
			isReducedSpecies[i] = true;

		GetReactionStringCHEMKIN( list_species, reaction_data, isReducedSpecies);
	}

	void ReactionPolicy_Solid_CHEMKIN::GetReactionStringCHEMKIN(const std::vector<std::string>& list_species,
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

	double ReactionPolicy_Solid_CHEMKIN::A_conversion() const
	{
		double conversion_factor = 0;
		// simple
		{
			conversion_factor = std::pow(1.e-3, 1. - sumLambdaGasReactants_ - sumLambdaSolidReactants_) /
								std::pow(1.e-4, 1. - 1.5*sumLambdaGasReactants_ - sumLambdaSolidReactants_) /
								std::pow(a_conversion_, 1. - sumLambdaGasReactants_ - sumLambdaSolidReactants_);
		}
            
		return conversion_factor;
	}
        
	double ReactionPolicy_Solid_CHEMKIN::Arev_conversion() const
    {
        double conversion_factor = 0;
        if(	tag_reaction_ == PhysicalConstants::REACTION_SURFACE_SIMPLE )
        {
            if(iExplicitlyReversible_ == true)
				conversion_factor = std::pow(1.e-3, 1. - sumLambdaGasReactants_ - sumLambdaSolidReactants_) /
									std::pow(1.e-4, 1. - 1.5*sumLambdaGasReactants_ - sumLambdaSolidReactants_) /
									std::pow(a_conversion_, 1. - sumLambdaGasReactants_ - sumLambdaSolidReactants_);
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

	void ReactionPolicy_Solid_CHEMKIN::WriteAdditionalDataOnASCIIFile(std::ostream& fOut) const
	{
		if (tag_reaction_ == PhysicalConstants::REACTION_SURFACE_SIMPLE)
		{
			
		}
	}

	void ReactionPolicy_Solid_CHEMKIN::WriteShortSummary(std::ostream& fOut, const std::vector<std::string>& list_species) const
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
	}

	double ReactionPolicy_Solid_CHEMKIN::GetForwardConversionFactor() const
	{
		return 1. / (	std::pow(1.e-3, 1. - sumLambdaGasReactants_ - sumLambdaSolidReactants_) /
						std::pow(1.e-4, 1. - 1.5*sumLambdaGasReactants_ - sumLambdaSolidReactants_));
	}

	double ReactionPolicy_Solid_CHEMKIN::GetBackwardConversionFactor() const
	{
		return 1. / (	std::pow(1.e-3, 1. - sumLambdaGasReactants_ - sumLambdaSolidReactants_) /
						std::pow(1.e-4, 1. - 1.5*sumLambdaGasReactants_ - sumLambdaSolidReactants_));
	}

	void ReactionPolicy_Solid_CHEMKIN::WriteSummary(std::ofstream& fOut, const std::vector<std::string>& list_species, const unsigned int index) const
	{
		std::string line_reaction;
		GetReactionString(list_species, line_reaction);
	
		fOut << "================================================================================================================================" << std::endl;
		fOut << " KINETIC DATA - REACTION  " << index << " " << TagASCII() << std::endl;
		fOut << "  " << line_reaction << std::endl; 
		fOut << "================================================================================================================================" << std::endl;
		fOut << " Change in moles in the reaction = " << sumNuProducts_ - sumNuReactants_ << std::endl;
	
		fOut << " Reaction order (forward, gas)   = " << std::setprecision(3) << sumLambdaGasReactants_ << std::endl;
		fOut << " Reaction order (forward, solid)  = " << std::setprecision(3) << sumLambdaSolidReactants_ << std::endl;
		if (iReversible_ == true)
		{
			fOut << " Reaction order (backward, gas)  = " << std::setprecision(3) << sumLambdaGasProducts_ << std::endl;
			fOut << " Reaction order (backward, solid) = " << std::setprecision(3) << sumLambdaSolidProducts_ << std::endl;
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

	bool ReactionPolicy_Solid_CHEMKIN::FatalErrorMessage(const std::string message)
	{
		std::cout << message << std::endl;
		return false;
	}
}
