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

#include "math/OpenSMOKEUtilities.h"
#include "ThermodynamicsMap.h"

namespace OpenSMOKE
{
	template<typename map> 
	KineticsMap_Solid_CHEMKIN<map>::KineticsMap_Solid_CHEMKIN(ThermodynamicsMap_Solid_CHEMKIN<map>& thermo, rapidxml::xml_document<>& doc, const unsigned int target, const unsigned int nPoints) :
	thermodynamics_(thermo)
	{
		this->number_of_points_ = nPoints;
		ImportSpeciesFromXMLFile(doc);
		ImportCoefficientsFromXMLFile(doc, target);
		this->T_ = this->P_ = 0.;
	}

	template<typename map> 
	KineticsMap_Solid_CHEMKIN<map>::KineticsMap_Solid_CHEMKIN(ThermodynamicsMap_Solid_CHEMKIN<map>& thermo, rapidxml::xml_document<>& doc, const std::string target, const unsigned int nPoints) :
	thermodynamics_(thermo)
	{
		this->number_of_points_ = nPoints;
		ImportSpeciesFromXMLFile(doc);
		ImportCoefficientsFromXMLFile(doc, target);
		this->T_ = this->P_ = 0.;
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::SetTemperature(const map& T)
	{
		this->T_old_ = this->T_;
		this->T_ = T;
		this->uT_ = 1./this->T_;
		this->logT_ = log(this->T_);
		Patm_over_RT_ = 101325./PhysicalConstants::R_J_kmol/this->T_;
		log_Patm_over_RT_ = log(Patm_over_RT_);

		if (std::fabs(this->T_-this->T_old_)/this->T_>1.e-14)
		{
			arrhenius_kinetic_constants_must_be_recalculated_ = true;
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
			reaction_h_and_s_must_be_recalculated_ = true;
		}
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::SetPressure(const map& P)
	{
		this->P_old_ = this->P_;
		this->P_ = P;
	//	if (std::fabs(this->P_-this->P_old_)/this->P_>1.e-14)
		{
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
		}
	}
/*
	void CheckKeyWord(const std::string read_value, const std::string expected_value)
	{
		if (read_value != expected_value)
		{
			std::cout << "Error in reading the kinetic mechanism: Expected " << expected_value << " - Found: " << read_value << std::endl;
			std::cout << "Please check your kinetic mechanism!" << std::endl;
			std::cout << "Press enter to continue..." << std::endl;
			getchar();
			exit(OPENSMOKE_FATAL_ERROR_EXIT);
		}
	}
*/
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)
	{
		ErrorMessage("void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The solid kinetic map require the user specifies tha material name.");
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc, const std::string target_material_name)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* parent_kinetics_node = opensmoke_node->first_node("Kinetics");

		std::string kinetics_type = parent_kinetics_node->first_attribute("type")->value();
		std::string kinetics_version = parent_kinetics_node->first_attribute("version")->value();
				
		if (kinetics_type != "OpenSMOKE" || kinetics_version != "01-02-2014")
			ErrorMessage("void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current kinetic scheme is not supported.");

		unsigned int target_material_index = 0;
		for (rapidxml::xml_node<> *kinetics_node = parent_kinetics_node->first_node("MaterialKinetics"); kinetics_node; kinetics_node = kinetics_node->next_sibling("MaterialKinetics"))
		{
			if ( target_material_name == kinetics_node->first_attribute("name")->value() )
			{
				target_material_index = boost::lexical_cast<unsigned int>(kinetics_node->first_attribute("index")->value());
				break;
			}
		}

		if (target_material_index == 0)
			ErrorMessage("void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The requested solid material is not available. Please check the name.");
		else
			ImportCoefficientsFromXMLFile(doc, target_material_index);
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc, const unsigned int target_material_index)
	{
		if (target_material_index <=0 || target_material_index > thermodynamics_.number_of_materials())
			ErrorMessage("void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The requested solid material is not available. Please check the name.");

		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* parent_kinetics_node = opensmoke_node->first_node("Kinetics");

		std::string kinetics_type = parent_kinetics_node->first_attribute("type")->value();
		std::string kinetics_version = parent_kinetics_node->first_attribute("version")->value();
				
		if (kinetics_type != "OpenSMOKE" || kinetics_version != "01-02-2014")
			ErrorMessage("void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current kinetic scheme is not supported.");

		for (rapidxml::xml_node<> *kinetics_node = parent_kinetics_node->first_node("MaterialKinetics"); kinetics_node; kinetics_node = kinetics_node->next_sibling("MaterialKinetics"))
		{
			if (target_material_index == boost::lexical_cast<unsigned int>(kinetics_node->first_attribute("index")->value()) )
			{
				if (kinetics_node == 0)
					ErrorMessage("void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "Kinetics tag was not found!");
				
				// Reading the number of species
				{
					std::cout << "Reading number of species..." << std::endl;
					rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
					this->number_of_species_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));
				}

				// Reading the number of reactions
				{
					std::cout << "Reading number of reactions..." << std::endl;
					rapidxml::xml_node<>* number_of_reactions_node = kinetics_node->first_node("NumberOfReactions");
					this->number_of_reactions_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_reactions_node->value())));
				}

				// Irreversible reactions
				{
					std::cout << "Reading irreversible..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("Irreversible");
					std::stringstream fInput;
					fInput << current_node->value();
					indices_of_irreversible_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					number_of_irreversible_reactions_ = indices_of_irreversible_reactions_.Size();
				}
		
				// Reversible reactions
				{
					std::cout << "Reading reversible..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible");
					std::stringstream fInput;
					fInput << current_node->value();
					indices_of_reversible_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					number_of_reversible_reactions_ = indices_of_reversible_reactions_.Size();
				}

				// Thermodynamic Reversible reactions
				{
					std::cout << "Reading reversible thermodynamics..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible-Thermodynamics");
					std::stringstream fInput;
					fInput << current_node->value();
					indices_of_thermodynamic_reversible_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					number_of_thermodynamic_reversible_reactions_ = indices_of_thermodynamic_reversible_reactions_.Size();
				}

				// Explicit Reversible reactions
				{
					std::cout << "Reading reversible explicit..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible-Explicit");
					std::stringstream fInput;
					fInput << current_node->value();
					indices_of_explicitly_reversible_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					number_of_explicitly_reversible_reactions_ = indices_of_explicitly_reversible_reactions_.Size();
				}

				// Reading if the kinetic scheme is conventional
				{
					std::cout << "Reading type of kinetic scheme..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("TypeOfKinetics");
					std::stringstream fInput;
					fInput << current_node->value();
					std::string dummy;
					fInput >> dummy;
					
					if (dummy == "chemkin_conventional")	type_of_solid_kinetics_ = TYPE_OF_SOLID_KINETICS_CHEMKIN_CONVENTIONAL;
					else
					{
						std::cout << "This type of kinetic mechanism is not available: " << dummy << std::endl;
						std::cout << "Press enter to exit..." << std::endl;						
						getchar();
						exit(OPENSMOKE_FATAL_ERROR_EXIT);
					}
				}

				std::cout << " * Reading kinetic parameters of reactions..." << std::endl;	
				{
					rapidxml::xml_node<>* kinetic_parameters_node = kinetics_node->first_node("KineticParameters");

					// Direct side
					{
						rapidxml::xml_node<>* direct_node = kinetic_parameters_node->first_node("Direct");

						// lnA
						{
							rapidxml::xml_node<>* current_node = direct_node->first_node("lnA");
							std::stringstream fInput;
							fInput << current_node->value();
							lnA_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
						}

						// Beta
						{
							rapidxml::xml_node<>* current_node = direct_node->first_node("Beta");
							std::stringstream fInput;
							fInput << current_node->value();
							Beta_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
						}

						// E_over_R
						{
							rapidxml::xml_node<>* current_node = direct_node->first_node("E_over_R");
							std::stringstream fInput;
							fInput << current_node->value();
							E_over_R_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
						}

						// Global kinetic order of forward reaction
						{
							rapidxml::xml_node<>* current_node = direct_node->first_node("ForwardKineticOrder");
							std::stringstream fInput;
							fInput << current_node->value();
							forward_kinetic_order_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
						}
					}

					// Reverse side
					if (number_of_explicitly_reversible_reactions_ != 0)
					{
						rapidxml::xml_node<>* reverse_node = kinetic_parameters_node->first_node("Reverse");

						// lnA
						{
							rapidxml::xml_node<>* current_node = reverse_node->first_node("lnA");
							std::stringstream fInput;
							fInput << current_node->value();
							lnA_reversible_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
						}

						// Beta
						{
							rapidxml::xml_node<>* current_node = reverse_node->first_node("Beta");
							std::stringstream fInput;
							fInput << current_node->value();
							Beta_reversible_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
						}

						// E_over_R
						{
							rapidxml::xml_node<>* current_node = reverse_node->first_node("E_over_R");
							std::stringstream fInput;
							fInput << current_node->value();
							E_over_R_reversible_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
						}
					}		
				}

				// Reactions needing conversion
				{			
					std::string dummy;

					rapidxml::xml_node<>* current_node = kinetics_node->first_node("ReactionsNeedingConversion");
					std::stringstream fInput;
					fInput << current_node->value();

					fInput >> dummy;
					unsigned int number_of_reactions_needing_conversion_ = boost::lexical_cast<unsigned int>(dummy);
					indices_of_reactions_needing_conversion_.resize(number_of_reactions_needing_conversion_);

					for (unsigned int j=0;j<number_of_reactions_needing_conversion_;j++)
					{		
						fInput >> dummy;
						indices_of_reactions_needing_conversion_[j] = boost::lexical_cast<unsigned int>(dummy);
					}
				}

				// Thermodynamic kinetic constants
				{
					std::string dummy;

					rapidxml::xml_node<>* current_node = kinetics_node->first_node("ThermodynamicReversibleReactions");
					std::stringstream fInput;
					fInput << current_node->value();

					fInput >> dummy;
					unsigned int n = boost::lexical_cast<unsigned int>(dummy);

					delta_nu_gas_.resize(n);

					for (unsigned int j=0;j<n;j++)
					{
						fInput >> dummy;
						const unsigned int index_phase = boost::lexical_cast<unsigned int>(dummy);
						fInput >> dummy;
						const double delta_sigma = boost::lexical_cast<double>(dummy);
						fInput >> dummy;
						delta_nu_gas_[j] = boost::lexical_cast<double>(dummy);
					}				
				}

				// Stoichiometry
				{
					stoichiometry_ = new StoichiometricMap(this->number_of_species_, this->number_of_reactions_);
	
					rapidxml::xml_node<>* stoichiometry_node = kinetics_node->first_node("Stoichiometry");

					std::string stoichiometry_type = stoichiometry_node->first_attribute("type")->value();
					std::string stoichiometry_version = stoichiometry_node->first_attribute("version")->value();

					if (stoichiometry_type != "OpenSMOKE" || stoichiometry_version != "01-02-2014")
						ErrorMessage("void KineticsMap_Solid_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current stoichiometric data are not supported.");

					std::stringstream fInput;
					fInput << stoichiometry_node->value();

					stoichiometry_->ReadFromASCIIFile(fInput);
					changeOfMoles_ = stoichiometry_->changeOfMoles();
					{
						OpenSMOKE::OpenSMOKEVectorBool tmp(this->number_of_reactions_);
						tmp = false;
						for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
							tmp[indices_of_thermodynamic_reversible_reactions_[k]] = true;
						stoichiometry_->CompleteChangeOfMoles(tmp);
 					}
			
		

					// Memory allocation			
					ChangeDimensions(this->number_of_species_, &aux_vector_, true);
					ChangeDimensions(this->number_of_species_, &c_, true);

					ChangeDimensions(this->number_of_reactions_, &reaction_s_over_R_, true);
					ChangeDimensions(this->number_of_reactions_, &reaction_h_over_RT_, true);
					ChangeDimensions(this->number_of_reactions_, &kArrhenius_, true);
					ChangeDimensions(this->number_of_reactions_, &kArrheniusModified_, true);
					ChangeDimensions(number_of_thermodynamic_reversible_reactions_, &uKeq_, false);
					ChangeDimensions(number_of_explicitly_reversible_reactions_, &kArrhenius_reversible_, false);

					ChangeDimensions(this->number_of_reactions_, &forwardReactionRates_, true);
					ChangeDimensions(this->number_of_reactions_, &reverseReactionRates_, true);
					ChangeDimensions(this->number_of_reactions_, &netReactionRates_, true);

					ChangeDimensions(this->number_of_reactions_, &isThermodynamicallyReversible_, true);
					ChangeDimensions(this->number_of_reactions_, &isExplicitlyReversible_, true);

					for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
						isThermodynamicallyReversible_[indices_of_thermodynamic_reversible_reactions_[k]] = k;
					for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
						isExplicitlyReversible_[indices_of_explicitly_reversible_reactions_[k]] = k;

			
					// Additional indices for sensitivity analysis
					{
						ChangeDimensions(this->number_of_reactions_, &type_of_reaction_, true);
						ChangeDimensions(this->number_of_reactions_, &local_family_index_, true);

						type_of_reaction_ = PhysicalConstants::REACTION_SIMPLE;
					}

					std::cout << std::endl;
					std::cout << "----------------------------------------------------------------------------" << std::endl;
					std::cout << " Kinetic Mechanism Summary"<< std::endl;
					std::cout << "----------------------------------------------------------------------------" << std::endl;
					std::cout << " Total number of species:          " << this->number_of_species_ << std::endl;
					std::cout << " Total number of reactions:        " << this->number_of_reactions_ << std::endl;
					std::cout << "   Reversible reactions:           " << number_of_reversible_reactions_ << " (" << number_of_reversible_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
					std::cout << "    * by thermodynamics:           " << number_of_thermodynamic_reversible_reactions_ << " (" << number_of_thermodynamic_reversible_reactions_/std::max(1.,double(number_of_reversible_reactions_))*100. << "%)" << std::endl;
					std::cout << "    * by Arrhenius' law:           " << number_of_explicitly_reversible_reactions_ << " (" << number_of_explicitly_reversible_reactions_/std::max(1.,double(number_of_reversible_reactions_))*100. << "%)" << std::endl;
					std::cout << std::endl;

					stoichiometry_->Summary(std::cout);
				}

				break;
			}
		}
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
		try
		{
			this->number_of_species_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));					
		}
		catch(...)
		{
			ErrorMessage("KineticsMap_Solid_CHEMKIN<map>::ImportSpeciesFromXMLFile", "Error in reading the number of species.");
		}
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::ReactionEnthalpiesAndEntropies()
	{
		if (reaction_h_and_s_must_be_recalculated_ == true)
		{
			stoichiometry_->ReactionEnthalpyAndEntropy(	reaction_h_over_RT_, reaction_s_over_R_, 
														thermodynamics_.species_h_over_RT(), thermodynamics_.species_s_over_R() );

			reaction_h_and_s_must_be_recalculated_ = false;
		}
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::KineticConstants()
	{
                ReactionEnthalpiesAndEntropies();
                        
		if (arrhenius_kinetic_constants_must_be_recalculated_ == true)
		{
			// Forward kinetic constants (Arrhenius' Law)
			{
				double *pt_lnA = lnA_.GetHandle();
				double *pt_Beta = Beta_.GetHandle();
				double *pt_E_over_R = E_over_R_.GetHandle();
				double *pt_kArrheniusT = kArrhenius_.GetHandle();
			
				for(unsigned int j=0;j<this->number_of_reactions_;j++)
					*pt_kArrheniusT++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;

				Exp(kArrhenius_, &kArrhenius_);
			}

			// Equilibrium constants (inverse value)
			{
				for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
				{
					unsigned int j = indices_of_thermodynamic_reversible_reactions_[k];
					uKeq_[k] = -reaction_s_over_R_[j] + reaction_h_over_RT_[j] - log_Patm_over_RT_ * delta_nu_gas_[k-1];
				}
				Exp(uKeq_, &uKeq_);
			}

			// Explicit reverse Arrhenius constants
			if (number_of_explicitly_reversible_reactions_ != 0)
			{
				double *pt_lnA = lnA_reversible_.GetHandle();
				double *pt_Beta = Beta_reversible_.GetHandle();
				double *pt_E_over_R = E_over_R_reversible_.GetHandle();
				double *pt_kArrhenius = kArrhenius_reversible_.GetHandle();

				for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
					*pt_kArrhenius++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius_reversible_, &kArrhenius_reversible_);
			}

			arrhenius_kinetic_constants_must_be_recalculated_ = false;
		}

		//if (nonconventional_kinetic_constants_must_be_recalculated_ == true)
		{
			nonconventional_kinetic_constants_must_be_recalculated_ = false;
		}

		// Conversions
		if (indices_of_reactions_needing_conversion_.size() > 0)
		{
			const double R_times_T = PhysicalConstants::R_J_kmol*this->T_;
			for(unsigned int k=0;k<indices_of_reactions_needing_conversion_.size();k++)
			{
				unsigned int j = indices_of_reactions_needing_conversion_[k];
				kArrhenius_[j] *= std::pow(R_times_T, forward_kinetic_order_[j]);
			}
		}

		kArrheniusModified_ = kArrhenius_;
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::ReactionRates(const OpenSMOKEVectorDouble& cGas, const OpenSMOKEVectorDouble& cSolid)
	{
		const double cTot = cGas.SumElements();
		
		unsigned int count = 1;
		for(unsigned int j=0;j<thermodynamics_.number_of_gas_species();j++)
			c_[count++] = cGas[j+1];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_solid_species();j++)
			c_[count++] = cSolid[j+1];
		
		if (type_of_solid_kinetics_ == TYPE_OF_SOLID_KINETICS_CHEMKIN_CONVENTIONAL)
		{
			// 1. Kinetic constants
			KineticConstants();
		}

		// Calculates the product of conenctrations (for forward and reverse reactions)
		// Be careful: the reverseReactionRates_ vector is defined for all the reactions
		// in the kinetic scheme, not only for the reversible reactions. After calling the
		// function reported below the value of reverseReactionRates_ vector for non reversible
		// reactions is put equal to 1.
		stoichiometry_->ProductOfConcentrations(forwardReactionRates_, reverseReactionRates_, c_);

		// Corrects the product of concentrations for reverse reaction by the 
		// thermodynamic equilibrium constant
		for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_thermodynamic_reversible_reactions_[k];
			reverseReactionRates_[j] *= uKeq_[k];
		}

		// Corrects the product of concentrations for reverse reaction by the 
		// explicit Arrhenius kinetic parameters (if provided)
		for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions_[k];
//			reverseReactionRates_[j] *= kArrhenius_reversible_[k]/kArrheniusModified_[j];
			reverseReactionRates_[j] *= kArrhenius_reversible_[k]/kArrhenius_[j];
		}

		// Calculates the net reaction rate
		// Be careful: the netReactionRates_ vector must be multiplied by the effective 
		// forward kinetic constant, to obtain the real reaction rates in [kmol/m3/s]
		netReactionRates_ = forwardReactionRates_;
		for(unsigned int k=1;k<=number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions_[k];
			netReactionRates_[j] -= reverseReactionRates_[j];
		}

		// Multiplies the net reaction rate by the effective kinetic constant (accounting for 
		// third-body effects, fall-off, etc.). At the end of this function the netReactionRates_
		// vector contains the net reaction rates of all the reactions in [kmol/m3/s]
		ElementByElementProduct(netReactionRates_, kArrheniusModified_, &netReactionRates_);
	}

	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::DerivativesOfReactionRatesWithRespectToKineticParameters(const PhysicalConstants::sensitivity_type type, unsigned int jReaction, const OpenSMOKEVectorDouble& c, double& parameter)
	{
		const double cTot = c.SumElements();

		// 1. Kinetic constants
		KineticConstants();
			
		// 2. Calculates the three-body corrections
		ThirdBodyReactions(cTot, c);

		// 3. Correct the effective kinetic constants by three-body coefficients
		for(unsigned int s=1;s<=number_of_thirdbody_reactions_;s++)
		{
			const unsigned int j=indices_of_thirdbody_reactions_[s];
			kArrheniusModified_[j] *= Meff_[s];
		}

		// 4. Correct the effective kinetic constants: Fall-off reactions
		FallOffReactions(cTot, c);
		for(unsigned int s=1;s<=number_of_falloff_reactions_;s++)
		{
			const unsigned int j=indices_of_falloff_reactions_[s];
			kArrheniusModified_[j] *= correction_falloff_[s];
		}

		// 5. Correct the effective kinetic constants: CABR reactions
		ChemicallyActivatedBimolecularReactions(cTot, c);
		for(unsigned int s=1;s<=number_of_cabr_reactions_;s++)
		{
			const unsigned int j=indices_of_cabr_reactions_[s];
			kArrheniusModified_[j] *= correction_cabr_[s];
		}

		if (type == PhysicalConstants::SENSITIVITY_KINETIC_CONSTANT)
		{
			parameter = kArrheniusModified_[jReaction];
			kArrheniusModified_ = 0.;
			kArrheniusModified_[jReaction] = 1.;
		}
		else if (type == PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR)
		{
			if (jReaction <= this->number_of_reactions_)
			{
				const double tmp = kArrheniusModified_[jReaction];
				kArrheniusModified_ = 0.;
				kArrheniusModified_[jReaction] = tmp;
				parameter =  std::exp(lnA_[jReaction]);

				if (type_of_reaction_[jReaction] == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || 
					type_of_reaction_[jReaction] == PhysicalConstants::REACTION_TROE_FALLOFF ||
					type_of_reaction_[jReaction] == PhysicalConstants::REACTION_SRI_FALLOFF)
				{
					unsigned int jFallOff = local_family_index_[jReaction];
					double	kInf = std::exp(lnA_falloff_inf_[jFallOff] + Beta_falloff_inf_[jFallOff]*this->logT_ - E_over_R_falloff_inf_[jFallOff]*this->uT_);
					double F, dF_over_dA0, dF_over_dAInf;this->
					FallOffReactions(jFallOff, cTot, c, F, dF_over_dA0, dF_over_dAInf);

					kArrheniusModified_[jReaction] = kArrheniusModified_[jReaction] / std::exp(lnA_[jReaction]) * (1.-kArrheniusModified_[jReaction]/kInf/F) +
											 kArrheniusModified_[jReaction]/F*dF_over_dA0;			
				}
				else if (type_of_reaction_[jReaction] == PhysicalConstants::REACTION_LINDEMANN_CABR || 
						 type_of_reaction_[jReaction] == PhysicalConstants::REACTION_TROE_CABR ||
						 type_of_reaction_[jReaction] == PhysicalConstants::REACTION_SRI_CABR)
				{
					unsigned int jCABR = local_family_index_[jReaction];
					double	k0   = std::exp(lnA_[jReaction] + Beta_[jReaction]*this->logT_ - E_over_R_[jReaction]*this->uT_);

					double F, dF_over_dA0, dF_over_dAInf;
					ChemicallyActivatedBimolecularReactions(jCABR, cTot, c, F, dF_over_dA0, dF_over_dAInf);
					kArrheniusModified_[jReaction] = kArrheniusModified_[jReaction]*kArrheniusModified_[jReaction] / std::exp(lnA_[jReaction])/k0/F +
											 kArrheniusModified_[jReaction]/F*dF_over_dA0;
				}
				else if (type_of_reaction_[jReaction] == PhysicalConstants::REACTION_CHEBYSHEV)
				{
					parameter = kArrheniusModified_[jReaction];
					kArrheniusModified_[jReaction] = 1.;
				}
				else
				{
					kArrheniusModified_[jReaction] = kArrheniusModified_[jReaction] / std::exp(lnA_[jReaction]);
				}
			}

			else 
			{
				// Fall-off reaction
				if (jReaction <= this->number_of_reactions_ + number_of_falloff_reactions_)
				{
					unsigned int jFallOff = jReaction - this->number_of_reactions_;
					unsigned int index_global = indices_of_falloff_reactions_[jFallOff];
					const double tmp = kArrheniusModified_[index_global];
					kArrheniusModified_ = 0.;
					kArrheniusModified_[index_global] = tmp;
					parameter =  std::exp(lnA_falloff_inf_[jFallOff]);

					double	kInf = std::exp(lnA_falloff_inf_[jFallOff] + Beta_falloff_inf_[jFallOff]*this->logT_ - E_over_R_falloff_inf_[jFallOff]*this->uT_);
					double F, dF_over_dA0, dF_over_dAInf;
					FallOffReactions(jFallOff, cTot, c, F, dF_over_dA0, dF_over_dAInf);

					kArrheniusModified_[index_global] = kArrheniusModified_[index_global] / F * 
											    ( kArrheniusModified_[index_global]/std::exp(lnA_falloff_inf_[jFallOff])/kInf + dF_over_dAInf);		
				}
				// CABR reactions
				else if (jReaction <= this->number_of_reactions_ + number_of_falloff_reactions_ +  number_of_cabr_reactions_)
				{
					unsigned int jCABR = jReaction - this->number_of_reactions_ - number_of_falloff_reactions_;
					unsigned int index_global = indices_of_cabr_reactions_[jCABR];
					const double tmp = kArrheniusModified_[index_global];
					kArrheniusModified_ = 0.;
					kArrheniusModified_[index_global] = tmp;
					parameter =  std::exp(lnA_falloff_inf_[jCABR]);

					double	k0   = std::exp(lnA_[index_global] + Beta_[index_global]*this->logT_ - E_over_R_[index_global]*this->uT_);

					double F, dF_over_dA0, dF_over_dAInf;
					ChemicallyActivatedBimolecularReactions(jCABR, cTot, c, F, dF_over_dA0, dF_over_dAInf);
					kArrheniusModified_[index_global] = kArrheniusModified_[index_global] / F * 
											    ( kArrheniusModified_[index_global]/std::exp(lnA_cabr_inf_[jCABR])/k0 + dF_over_dAInf);	
				}
			}
		}


		// Calculates the product of conenctrations (for forward and reverse reactions)
		// Be careful: the reverseReactionRates_ vector is defined for all the reactions
		// in the kinetic scheme, not only for the reversible reactions. After calling the
		// function reported below the value of reverseReactionRates_ vector for non reversible
		// reactions is put equal to 1.
		stoichiometry_->ProductOfConcentrations(forwardReactionRates_, reverseReactionRates_, c);
		
		// Corrects the product of concentrations for reverse reaction by the 
		// thermodynamic equilibrium constant
		for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_thermodynamic_reversible_reactions_[k];
			reverseReactionRates_[j] *= uKeq_[k];
		}

		// Corrects the product of concentrations for reverse reaction by the 
		// explicit Arrhenius kinetic parameters (if provided)
//		for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
//		{
//			unsigned int j = indices_of_explicitly_reversible_reactions_[k];
		//	reverseReactionRates_[j] *= kArrhenius_reversible_[k]/kArrheniusModified_[j];
//			if (j==jReaction)
//				reverseReactionRates_[j] *= kArrhenius_reversible_[k]/parameter;
//			else
//				reverseReactionRates_[j] = 0.;
//		}

		for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions_[k];
			if (j==jReaction)
				reverseReactionRates_[j] *= kArrhenius_reversible_[k]/kArrhenius_[j];
			else
				reverseReactionRates_[j] = 0.;
		}

		// Calculates the net reaction rate
		// Be careful: the netReactionRates_ vector must be multiplied by the effective 
		// forward kinetic constant, to obtain the real reaction rates in [kmol/m3/s]
		netReactionRates_ = forwardReactionRates_;
		for(unsigned int k=1;k<=number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions_[k];
			netReactionRates_[j] -= reverseReactionRates_[j];

		}

		// Multiplies the net reaction rate by the effective kinetic constant (accounting for 
		// third-body effects, fall-off, etc.). At the end of this function the netReactionRates_
		// vector contains the net reaction rates of all the reactions in [kmol/m3/s]
		ElementByElementProduct(netReactionRates_, kArrheniusModified_, &netReactionRates_);
	}
	*/

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::FormationRates(OpenSMOKEVectorDouble* Rgas, OpenSMOKEVectorDouble* Rsolid)
	{
		OpenSMOKE::OpenSMOKEVectorDouble R(thermodynamics_.NumberOfSpecies());
		stoichiometry_->FormationRatesFromReactionRates(&R, netReactionRates_);
		
		unsigned int count = 1;
		for(unsigned int j=1;j<=thermodynamics_.number_of_gas_species();j++)
			(*Rgas)[j] = R[count++];
		
		for(unsigned int j=1;j<=thermodynamics_.number_of_solid_species();j++)
			(*Rsolid)[j] = R[count++];
	}

	template<typename map> 
	double KineticsMap_Solid_CHEMKIN<map>::HeatRelease(const OpenSMOKEVectorDouble& Rgas, const OpenSMOKEVectorDouble& RSolid)
	{
		unsigned int k = 1;
		for(int j=1;j<=Rgas.Size();j++)
			aux_vector_[k++] = Rgas[j];
		for(int j=1;j<=RSolid.Size();j++)
			aux_vector_[k++] = RSolid[j];

	//	std::cout << "Gas " << std::endl;
	//	for(int j=1;j<=Rgas.Size();j++)
	//	std::cout << j << " " << Rgas[j] << std::endl;
	//	std::cout << "Solid " << std::endl;
	//	for(int j=1;j<=RSolid.Size();j++)
	//	std::cout << j << " " << RSolid[j] << std::endl;
	//	getchar();
		return -Dot(aux_vector_, thermodynamics_.species_h_over_RT()) * PhysicalConstants::R_J_kmol * this->T_;
	}

	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const unsigned int k, const OpenSMOKEVectorDouble& c, OpenSMOKEVectorDouble* Jalfa, double& parameter)
	{
		DerivativesOfReactionRatesWithRespectToKineticParameters(type, k, c, parameter);
		stoichiometry_->FormationRatesFromReactionRates(Jalfa, netReactionRates_);
	}
	*/

	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& parameter)
	{
		DerivativesOfReactionRatesWithRespectToKineticParameters(type, k, c, parameter);
		stoichiometry_->FormationRatesFromReactionRates(Jalfa, netReactionRates_);
		
		JT  = HeatRelease(*Jalfa);
		if (energy_type == CONSTANT_VOLUME_SYSTEM)
		{
			double sumMoleFormationRates = (*Jalfa).SumElements();
			JT  += PhysicalConstants::R_J_kmol*this->T_*sumMoleFormationRates;
		}	
		if (energy_type == PLUGFLOW_SYSTEM)
		{
			// TODO
			//double sumMoleFormationRates = (*Jalfa).SumElements();
		}	
	}
	*/

	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& Jrho, double& parameter)
	{
		DerivativesOfReactionRatesWithRespectToKineticParameters(type, k, c, parameter);
		stoichiometry_->FormationRatesFromReactionRates(Jalfa, netReactionRates_);
		
		double CpMixMolar;
		thermodynamics_.cpMolar_Mixture_From_MoleFractions(CpMixMolar, mole_fractions);
		
		const double dQR_over_parameter = HeatRelease(*Jalfa);
		double sumMoleFormationRates = (*Jalfa).SumElements();
		JT    =  dQR_over_parameter;
		Jrho  = - (dQR_over_parameter/CpMixMolar + this->T_*sumMoleFormationRates);
	}
	*/

	/*
	// This version seems to be wrong
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::ProductionAndDestructionRatesGross(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D)
	{
		OpenSMOKEVectorDouble forward_(this->number_of_reactions_);
		OpenSMOKEVectorDouble backward_(this->number_of_reactions_);

		GetForwardReactionRates(&forward_);
		GetBackwardReactionRates(&backward_);

		stoichiometry_->ProductionAndDestructionRatesFromReactionRatesGross(P, D, forward_, backward_);
	}
	*/

	
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::ProductionAndDestructionRates(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D)
	{
		stoichiometry_->ProductionAndDestructionRatesFromReactionRates(P, D, netReactionRates_);
	}
	

	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::RateOfProductionAnalysis(const bool iNormalize) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, iNormalize);
	}
	*/
	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::RateOfProductionAnalysis(std::ostream& fout) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, false);
		stoichiometry_->WriteRateOfProductionAnalysis(fout);
	}
	*/

	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::RateOfProductionAnalysis(ROPA_Data& ropa) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, false);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}
	*/

	template<typename map> 
	const OpenSMOKEVectorDouble& KineticsMap_Solid_CHEMKIN<map>::GetReactionRates()
	{
		return netReactionRates_;
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::GetReactionRates(OpenSMOKEVectorDouble* r)
	{
		*r = netReactionRates_;
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::GetForwardReactionRates(OpenSMOKEVectorDouble* r)
	{
		ElementByElementProduct(forwardReactionRates_, kArrheniusModified_, r);
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::GetBackwardReactionRates(OpenSMOKEVectorDouble* r)
	{
		*r = 0.;
		for(unsigned int k=1;k<=number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions_[k];
			(*r)[j] = reverseReactionRates_[j]*kArrheniusModified_[j];
		}
		for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions_[k];
			(*r)[j] = reverseReactionRates_[j]*kArrheniusModified_[j];
		}                
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::WriteKineticData(std::ostream& fOut, const unsigned int k)
	{			
		thermodynamics_.SetPressure(101325.);
		thermodynamics_.SetTemperature(298.15);
		SetTemperature(298.15);
		SetPressure(101325.);
		ReactionEnthalpiesAndEntropies();
		
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << (reaction_h_over_RT_[k]-reaction_s_over_R_[k])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_h_over_RT_[k] *PhysicalConstants::R_kcal_mol*this->T_;;						// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_s_over_R_[k] *PhysicalConstants::R_cal_mol;;							// [cal/mol/K]
	}


	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::WriteKineticData(std::ostream& fOut, const unsigned int k, OpenSMOKEVectorDouble& c_bath, const double conversion_forward, const double conversion_backward)
	{			
		const double patm = 101325.;
		SetPressure(patm);
		thermodynamics_.SetPressure(patm);
	
		OpenSMOKEVectorDouble temperatures(5); 
		temperatures[1] = 300.;
		temperatures[2] = 1000.;
		temperatures[3] = 1500.;
		temperatures[4] = 2000.;
		temperatures[5] = 2500.;

		fOut << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOut << "    Temperature   kF            Keq           kR            DG            DH            DS            kF            kR" << std::endl;          
		fOut << "    [K]           [kmol,m3,s]   [-]           [kmol,m3,s]   [kcal/mol]    [kcal/mol]    [cal/mol/K]   [mol,cm3,s]   [mol,cm3,s]" << std::endl;        
		fOut << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;

		for(int i=1;i<=temperatures.Size();i++)
		{
			SetTemperature(temperatures[i]);
			thermodynamics_.SetTemperature(temperatures[i]);

			KineticConstants();
			ReactionRates(c_bath);

			// Temperature
			fOut << "    " << std::setw(14) << std::left << std::setprecision(0) << std::fixed << temperatures[i];
			
			// Forward kinetic constant [kmol, m3, s]
			fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified_[k];

			// Equilibrium and Backward kinetic constant [kmol, m3, s]
			if (isThermodynamicallyReversible_[k] != 0)
			{
				const unsigned int j = isThermodynamicallyReversible_[k];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << 1./uKeq_[j];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified_[k]*uKeq_[j];
			}
			else if (isExplicitlyReversible_[k] != 0)
			{
				const unsigned int j = isExplicitlyReversible_[k];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified_[k]/kArrhenius_reversible_[j];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrhenius_reversible_[j];
			}
			else
			{
				fOut << std::setw(14) << std::left << std::setprecision(0) << std::fixed << 0.;
				fOut << std::setw(14) << std::left << std::setprecision(0) << std::fixed << 0.;
			}

			// Thermodynamic data
			fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << (reaction_h_over_RT_[k]-reaction_s_over_R_[k])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
			fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_h_over_RT_[k] *PhysicalConstants::R_kcal_mol*this->T_;;						// [kcal/mol]
			fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_s_over_R_[k] *PhysicalConstants::R_cal_mol;;							// [cal/mol/K]

			// Forward kinetic constant [mol, cm3, s]
			fOut    << std::setw(14) << std::left << std::setprecision(4) << std::scientific << kArrheniusModified_[k]*conversion_forward;

			// Equilibrium and Backward kinetic constant [mol, cm3, s]
			if (isThermodynamicallyReversible_[k] != 0)
			{
				const unsigned int j = isThermodynamicallyReversible_[k];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified_[k]*uKeq_[j]*conversion_backward;
			}
			else if (isExplicitlyReversible_[k] != 0)
			{
				const unsigned int j = isExplicitlyReversible_[k];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrhenius_reversible_[j]*conversion_backward;
			}
			else
			{
				fOut << std::setw(14) << std::left << std::setprecision(0) << std::fixed << 0.;
			}

			fOut    << std::endl;
		}
		fOut << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOut << std::endl;
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::FittedReverseKineticConstants(OpenSMOKEVectorDouble& x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters)
	{			
		unsigned int npoints = 11;
		OpenSMOKEVectorDouble temperatures(npoints); 
		
		temperatures[1] = 300.;		temperatures[2] = 600.;		temperatures[3] = 900.;		temperatures[4] = 1100.;
		temperatures[5] = 1300.;	temperatures[6] = 1500.;	temperatures[7] = 1700.;	temperatures[8] = 1900.;
		temperatures[9] = 2100.;    temperatures[10] = 2300.;   temperatures[11] = 2500.;
		OpenSMOKEVectorDouble c(x_bath.Size()); 

		Eigen::MatrixXd XTX(nparameters, nparameters);
		Eigen::MatrixXd XT(nparameters, npoints);

		std::cout << "   assembling X and XT matrices..." << std::endl;
		{
			// Assembling X and XT Matrices
			Eigen::MatrixXd X(npoints, nparameters);
		
			for(int i=0;i<temperatures.Size();i++)
			{
				const double T = temperatures[i+1];

				X(i,0) = 1.;
				X(i,1) = -1. / PhysicalConstants::R_J_kmol / T;

				if (nparameters == 3)
					X(i,2) = log(T);
			}

			XT = X.transpose();
			XTX=XT*X;
		}

		const double patm = 101325.;
		SetPressure(patm);
		thermodynamics_.SetPressure(patm);

		Eigen::MatrixXd y(npoints, number_of_thermodynamic_reversible_reactions_);
		Eigen::MatrixXd Y(nparameters, number_of_thermodynamic_reversible_reactions_);

		std::cout << "   evaluating the reaction rates..." << std::endl;
		for(int i=1;i<=temperatures.Size();i++)
		{
			const double ctot = patm/PhysicalConstants::R_J_kmol/temperatures[i];
			Product(ctot, x_bath, &c);

			SetTemperature(temperatures[i]);
			thermodynamics_.SetTemperature(temperatures[i]);

			KineticConstants();
			ReactionRates(c);

			for (unsigned int k=1;k<=this->number_of_reactions_;k++)
			{
				if (isThermodynamicallyReversible_[k] != 0)
				{
					const unsigned int j = isThermodynamicallyReversible_[k];
					y(i-1,j-1) = log(kArrheniusModified_[k]*uKeq_[j]);
				}
			}
		}

		Y=XT*y;

		std::cout << "   solving the linear regressions..." << std::endl;
		fittedKineticParameters = XTX.partialPivLu().solve(Y);
	}

	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters)
	{
		if (isThermodynamicallyReversible_[k] != 0)
		{
			const unsigned int j = isThermodynamicallyReversible_[k];

			fOut << std::setw(18) << std::right << std::scientific << std::setprecision(4) << std::exp(fittedKineticParameters(0, j-1));
			if (fittedKineticParameters.rows() == 2)
				fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << 0.;
			else
				fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << fittedKineticParameters(2, j - 1);
			fOut << std::setw(16) << std::right << std::fixed << std::setprecision(2) << fittedKineticParameters(1,j-1)/Conversions::J_from_kcal;
			fOut << std::setw(5)  << "";
		}
		else if (isExplicitlyReversible_[k] != 0)
		{
			const unsigned int j = isExplicitlyReversible_[k];
				
			fOut << std::setw(18) << std::right << std::scientific << std::setprecision(4) << std::exp(lnA_reversible_[j]);
			fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << Beta_reversible_[j];
			fOut << std::setw(16) << std::right << std::fixed << std::setprecision(2) << E_over_R_reversible_[j];
			fOut << std::setw(5)  << "";
		}
	}

	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* dR_over_domega)
	{
		OpenSMOKE::OpenSMOKEMatrixDouble dc_over_domega(this->number_of_species_, this->number_of_species_);
		OpenSMOKE::OpenSMOKEMatrixDouble dR_over_dc(this->number_of_species_, this->number_of_species_);

		double MW;
		thermodynamics_.MolecularWeight_From_MassFractions(MW, omega);
		double cTot = c.SumElements();

		thermodynamics_.DerivativesOfConcentrationsWithRespectToMassFractions(cTot, MW, omega, &dc_over_domega);
		DerivativesOfFormationRates(c, &dR_over_dc);

		for(unsigned int k=1;k<=this->number_of_species_;k++)
		{
			for(unsigned int i=1;i<=this->number_of_species_;i++)
			{
				double sum = 0.;
				for(unsigned int j=1;j<=this->number_of_species_;j++)
					sum += dR_over_dc[k][j] * dc_over_domega[j][i];
				(*dR_over_domega)[k][i] = sum;
			}
		}
	}
	*/
	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::Derivatives(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* derivatives)
	{
		OpenSMOKE::OpenSMOKEMatrixDouble dc_over_domega(this->number_of_species_, this->number_of_species_);
		OpenSMOKE::OpenSMOKEMatrixDouble derivatives_over_dc(this->number_of_species_+1, this->number_of_species_+1);

		double MW;
		thermodynamics_.MolecularWeight_From_MassFractions(MW, omega);
		double cTot = c.SumElements();

		thermodynamics_.DerivativesOfConcentrationsWithRespectToMassFractions(cTot, MW, omega, &dc_over_domega);
		Derivatives(c, &derivatives_over_dc);

		for(unsigned int k=1;k<=this->number_of_species_;k++)
		{
			for(unsigned int i=1;i<=this->number_of_species_;i++)
			{
				double sum = 0.;
				for(unsigned int j=1;j<=this->number_of_species_;j++)
					sum += derivatives_over_dc[k][j] * dc_over_domega[j][i];
				(*derivatives)[k][i] = sum;
			}
		}

		for(unsigned int i=1;i<=this->number_of_species_;i++)
		{
			double sum = 0.;
			for(unsigned int j=1;j<=this->number_of_species_;j++)
				sum += derivatives_over_dc[this->number_of_species_+1][j] * dc_over_domega[j][i];
			(*derivatives)[this->number_of_species_+1][i] = sum;
		}

		for(unsigned int i=1;i<=this->number_of_species_+1;i++)
			(*derivatives)[i][this->number_of_species_+1] = derivatives_over_dc[i][this->number_of_species_+1];
	}
	*/
	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* dR_over_dC)
	{
		const double ZERO_DER = sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = sqrt(OPENSMOKE_MACH_EPS_DOUBLE);			
		const double TOLR = 100. * OPENSMOKE_MACH_EPS_FLOAT;
		const double TOLA = 1.e-10;
		
		OpenSMOKEVectorDouble c_plus = c;
		OpenSMOKEVectorDouble R_original(this->number_of_species_);
		OpenSMOKEVectorDouble R_plus(this->number_of_species_);

		// Calculates centered values
		ReactionRates(c);
		FormationRates(&R_original);

		// derivata rispetto a ckd
		for(unsigned int kd=1;kd<=this->number_of_species_;kd++)
		{
			if(c[kd] <= 1.e-100)
			{
				double dc = 1.e-10 + 1.e-12 * c_plus[kd];
				double udc = 3.e-8 * Max(c_plus[kd],1./dc);
				udc = std::max(udc,1./dc);
				udc = std::max(udc,1.e-19);
				dc = std::min(udc, 0.001 + 0.001*c_plus[kd]);
				c_plus[kd] += dc;
				udc = 1./dc;

				ReactionRates(c_plus);
				FormationRates(&R_plus);

				for(unsigned int j=1;j<=this->number_of_species_;j++)
					(*dR_over_dC)[j][kd] = (R_plus[j]-R_original[j]) * udc;

				c_plus[kd] = c[kd];
			}
			else
			{
				double hf = 1.e0;
				double error_weight = 1./(TOLA+TOLR*std::fabs(c[kd]));
				double hJ = ETA2 * std::fabs(std::max(c[kd], 1./error_weight));
				double hJf = hf/error_weight;
				hJ = std::max(hJ, hJf);
				hJ = std::max(hJ, ZERO_DER);

				// This is what is done by BzzMath
				double dc = std::min(hJ, 1.e-3 + 1e-3*std::fabs(c[kd]));

				// Thisis what is done in the KPP
				//double dc = TOLR*c[kd];

				double udc = 1. / dc;
				c_plus[kd] += dc;

				ReactionRates(c_plus);
				FormationRates(&R_plus);

				for(unsigned int j=1;j<=this->number_of_species_;j++)
					(*dR_over_dC)[j][kd] = (R_plus[j]-R_original[j]) * udc;

				c_plus[kd] = c[kd];
			}
		}
	}
	*/
	/*
	template<typename map> 
	void KineticsMap_Solid_CHEMKIN<map>::Derivatives(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* derivatives, const bool constant_density)
	{
		const double ZERO_DER = sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = sqrt(OPENSMOKE_MACH_EPS_DOUBLE);			
		const double TOLR = 100. * OPENSMOKE_MACH_EPS_FLOAT;
		const double TOLA = 1.e-10;
	
		OpenSMOKEVectorDouble c_plus = c;
		OpenSMOKEVectorDouble R_original(this->number_of_species_);
		OpenSMOKEVectorDouble R_plus(this->number_of_species_);
		double Q_original;

		SetTemperature(this->T_);
		SetPressure(this->P_);
		thermodynamics_.SetTemperature(this->T_);
		thermodynamics_.SetPressure(this->P_);

		KineticConstants();

		// Calculates centered values
		ReactionRates(c);
		FormationRates(&R_original);
		Q_original = HeatRelease(R_original);

		// derivata rispetto a ckd
		for(unsigned int kd=1;kd<=this->number_of_species_;kd++)
		{
			if(c[kd] <= 1.e-100)
			{
				double dc = 1.e-10 + 1.e-12 * c_plus[kd];
				double udc = 3.e-8 * std::max(c_plus[kd],1./dc);
				udc = std::max(udc,1./dc);
				udc = std::max(udc,1.e-19);
				dc = std::min(udc, 0.001 + 0.001*c_plus[kd]);
				c_plus[kd] += dc;
				udc = 1./dc;

				ReactionRates(c_plus);
				FormationRates(&R_plus);
				double Q_plus = HeatRelease(R_plus);

				for(unsigned int j=1;j<=this->number_of_species_;j++)
					(*derivatives)[j][kd] = (R_plus[j]-R_original[j]) * udc;
				(*derivatives)[this->number_of_species_+1][kd] = (Q_plus-Q_original) * udc;

				c_plus[kd] = c[kd];
			}
			else
			{
				double hf = 1.e0;
				double error_weight = 1./(TOLA+TOLR*std::fabs(c[kd]));
				double hJ = ETA2 * std::fabs(std::max(c[kd], 1./error_weight));
				double hJf = hf/error_weight;
				hJ = std::max(hJ, hJf);
				hJ = std::max(hJ, ZERO_DER);

				// This is what is done by BzzMath
				double dc = std::min(hJ, 1.e-3 + 1e-3*std::fabs(c[kd]));

				// Thisis what is done in the KPP
				//double dc = TOLR*c[kd];
				
				double udc = 1. / dc;
				c_plus[kd] += dc;
				ReactionRates(c_plus);
				FormationRates(&R_plus);
				double Q_plus = HeatRelease(R_plus);

				for(unsigned int j=1;j<=this->number_of_species_;j++)
					(*derivatives)[j][kd] = (R_plus[j]-R_original[j]) * udc;
				(*derivatives)[this->number_of_species_+1][kd] = (Q_plus-Q_original) * udc;

				c_plus[kd] = c[kd];
			}
		}

		// Derivatives with respect to the temperature
		{
			double dT = TOLR*this->T_;
			double T_plus = this->T_+dT;
			double udT = 1. / dT;

			SetTemperature(T_plus);
			SetPressure(this->P_);
			thermodynamics_.SetTemperature(T_plus);
			thermodynamics_.SetPressure(this->P_);

			KineticConstants();

			ReactionRates(c);
			FormationRates(&R_plus);
			double Q_plus = HeatRelease(R_plus);

			for(unsigned int j=1;j<=this->number_of_species_;j++)
				(*derivatives)[j][this->number_of_species_+1] = (R_plus[j]-R_original[j]) * udT;
			(*derivatives)[this->number_of_species_+1][this->number_of_species_+1] = (Q_plus-Q_original) * udT;
		}

		// If density is constant, this correction has to be skipped
		if (constant_density == false)
		{
			SetTemperature(this->T_);
			SetPressure(this->P_);
			thermodynamics_.SetTemperature(this->T_);
			thermodynamics_.SetPressure(this->P_);

			unsigned int index_T = this->number_of_species_+1;
			for(unsigned int k=1;k<=this->number_of_species_+1;k++)
			{
				double sum = 0.;
				for(unsigned int j=1;j<=this->number_of_species_;j++)
					sum += c[j]/this->T_ * (*derivatives)[k][j];
				(*derivatives)[k][index_T] -= sum;
			}
		}
	}
	*/
}

