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
	KineticsMap_Solid_CHEMKIN::KineticsMap_Solid_CHEMKIN(ThermodynamicsMap_Solid_CHEMKIN& thermo, rapidxml::xml_document<>& doc, const unsigned int target) :
	thermodynamics_(thermo)
	{
		ImportSpeciesFromXMLFile(doc);
		ImportCoefficientsFromXMLFile(doc, target);
		this->T_ = this->P_ = 0.;
	}

	KineticsMap_Solid_CHEMKIN::KineticsMap_Solid_CHEMKIN(ThermodynamicsMap_Solid_CHEMKIN& thermo, rapidxml::xml_document<>& doc, const std::string target) :
	thermodynamics_(thermo)
	{
		ImportSpeciesFromXMLFile(doc);
		ImportCoefficientsFromXMLFile(doc, target);
		this->T_ = this->P_ = 0.;
	}

	void KineticsMap_Solid_CHEMKIN::SetTemperature(const double& T)
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

	void KineticsMap_Solid_CHEMKIN::SetPressure(const double& P)
	{
		this->P_old_ = this->P_;
		this->P_ = P;
	//	if (std::fabs(this->P_-this->P_old_)/this->P_>1.e-14)
		{
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
		}
	}

	void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)
	{
		ErrorMessage("void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The solid kinetic map require the user specifies tha material name.");
	}

	void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc, const std::string target_material_name)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* parent_kinetics_node = opensmoke_node->first_node("Kinetics");

		std::string kinetics_type = parent_kinetics_node->first_attribute("type")->value();
		std::string kinetics_version = parent_kinetics_node->first_attribute("version")->value();
				
		if (kinetics_type != "OpenSMOKE" || kinetics_version != "01-02-2014")
			ErrorMessage("void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current kinetic scheme is not supported.");

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
			ErrorMessage("void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The requested solid material is not available. Please check the name.");
		else
			ImportCoefficientsFromXMLFile(doc, target_material_index);
	}
 
	void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc, const unsigned int target_material_index)
	{
		if (target_material_index <=0 || target_material_index > thermodynamics_.number_of_materials())
			ErrorMessage("void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The requested solid material is not available. Please check the name.");

		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* parent_kinetics_node = opensmoke_node->first_node("Kinetics");

		std::string kinetics_type = parent_kinetics_node->first_attribute("type")->value();
		std::string kinetics_version = parent_kinetics_node->first_attribute("version")->value();
				
		if (kinetics_type != "OpenSMOKE" || kinetics_version != "01-02-2014")
			ErrorMessage("void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current kinetic scheme is not supported.");

		for (rapidxml::xml_node<> *kinetics_node = parent_kinetics_node->first_node("MaterialKinetics"); kinetics_node; kinetics_node = kinetics_node->next_sibling("MaterialKinetics"))
		{
			if (target_material_index == boost::lexical_cast<unsigned int>(kinetics_node->first_attribute("index")->value()) )
			{
				if (kinetics_node == 0)
					ErrorMessage("void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "Kinetics tag was not found!");
				
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
						ErrorMessage("void KineticsMap_Solid_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current stoichiometric data are not supported.");

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

	void KineticsMap_Solid_CHEMKIN::ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
		try
		{
			this->number_of_species_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));					
		}
		catch(...)
		{
			ErrorMessage("KineticsMap_Solid_CHEMKIN::ImportSpeciesFromXMLFile", "Error in reading the number of species.");
		}
	}

	void KineticsMap_Solid_CHEMKIN::ReactionEnthalpiesAndEntropies()
	{
		if (reaction_h_and_s_must_be_recalculated_ == true)
		{
			stoichiometry_->ReactionEnthalpyAndEntropy(	reaction_h_over_RT_, reaction_s_over_R_, 
														thermodynamics_.species_h_over_RT(), thermodynamics_.species_s_over_R() );

			reaction_h_and_s_must_be_recalculated_ = false;
		}
	}

	void KineticsMap_Solid_CHEMKIN::KineticConstants()
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

	void KineticsMap_Solid_CHEMKIN::ReactionRates(const OpenSMOKEVectorDouble& cGas, const OpenSMOKEVectorDouble& cSolid)
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
 
	void KineticsMap_Solid_CHEMKIN::FormationRates(OpenSMOKEVectorDouble* Rgas, OpenSMOKEVectorDouble* Rsolid)
	{
		OpenSMOKE::OpenSMOKEVectorDouble R(thermodynamics_.NumberOfSpecies());
		stoichiometry_->FormationRatesFromReactionRates(&R, netReactionRates_);
		
		unsigned int count = 1;
		for(unsigned int j=1;j<=thermodynamics_.number_of_gas_species();j++)
			(*Rgas)[j] = R[count++];
		
		for(unsigned int j=1;j<=thermodynamics_.number_of_solid_species();j++)
			(*Rsolid)[j] = R[count++];
	}
 
	double KineticsMap_Solid_CHEMKIN::HeatRelease(const OpenSMOKEVectorDouble& Rgas, const OpenSMOKEVectorDouble& RSolid)
	{
		unsigned int k = 1;
		for(int j=1;j<=Rgas.Size();j++)
			aux_vector_[k++] = Rgas[j];
		for(int j=1;j<=RSolid.Size();j++)
			aux_vector_[k++] = RSolid[j];

		return -Dot(aux_vector_, thermodynamics_.species_h_over_RT()) * PhysicalConstants::R_J_kmol * this->T_;
	}

	void KineticsMap_Solid_CHEMKIN::ProductionAndDestructionRates(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D)
	{
		stoichiometry_->ProductionAndDestructionRatesFromReactionRates(P, D, netReactionRates_);
	}
	
	const OpenSMOKEVectorDouble& KineticsMap_Solid_CHEMKIN::GetReactionRates()
	{
		return netReactionRates_;
	}

	void KineticsMap_Solid_CHEMKIN::GetReactionRates(OpenSMOKEVectorDouble* r)
	{
		*r = netReactionRates_;
	}

	void KineticsMap_Solid_CHEMKIN::GetForwardReactionRates(OpenSMOKEVectorDouble* r)
	{
		ElementByElementProduct(forwardReactionRates_, kArrheniusModified_, r);
	}
 
	void KineticsMap_Solid_CHEMKIN::GetBackwardReactionRates(OpenSMOKEVectorDouble* r)
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

	void KineticsMap_Solid_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k)
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
 
	void KineticsMap_Solid_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k, OpenSMOKEVectorDouble& c_bath, const double conversion_forward, const double conversion_backward)
	{			
		// TODO
	}

	void KineticsMap_Solid_CHEMKIN::FittedReverseKineticConstants(OpenSMOKEVectorDouble& x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters)
	{			
		// TODO
	}

	void KineticsMap_Solid_CHEMKIN::FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters)
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
}

