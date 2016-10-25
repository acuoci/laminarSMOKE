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
	KineticsMap_Surface_CHEMKIN<map>::KineticsMap_Surface_CHEMKIN(ThermodynamicsMap_Surface_CHEMKIN<map>& thermo, const unsigned int nSpecies, const unsigned int nPoints) :
	thermodynamics_(thermo)
	{
		this->number_of_species_ = nSpecies;
		this->number_of_points_ = nPoints;
		this->T_ = this->P_ = 0.;
	}

	template<typename map> 
	KineticsMap_Surface_CHEMKIN<map>::KineticsMap_Surface_CHEMKIN(ThermodynamicsMap_Surface_CHEMKIN<map>& thermo, rapidxml::xml_document<>& doc, const unsigned int nPoints) :
	thermodynamics_(thermo)
	{
		this->number_of_points_ = nPoints;
		ImportSpeciesFromXMLFile(doc);
		ImportCoefficientsFromXMLFile(doc);
		this->T_ = this->P_ = 0.;
	}

	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::SetTemperature(const map& T)
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
	void KineticsMap_Surface_CHEMKIN<map>::SetPressure(const map& P)
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
	void KineticsMap_Surface_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* parent_kinetics_node = opensmoke_node->first_node("Kinetics");

		std::string kinetics_type = parent_kinetics_node->first_attribute("type")->value();
		std::string kinetics_version = parent_kinetics_node->first_attribute("version")->value();
				
		if (kinetics_type != "OpenSMOKE" || kinetics_version != "01-02-2014")
			ErrorMessage("void KineticsMap_Surface_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current kinetic scheme is not supported.");

		
		unsigned int target_material = 1;
		for (rapidxml::xml_node<> *kinetics_node = parent_kinetics_node->first_node("MaterialKinetics"); kinetics_node; kinetics_node = kinetics_node->next_sibling("MaterialKinetics"))
		{
			if (target_material == boost::lexical_cast<unsigned int>(kinetics_node->first_attribute("index")->value()) )
			{
				
				if (kinetics_node == 0)
					ErrorMessage("void KineticsMap_Surface_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "Kinetics tag was not found!");
				
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

				// Stick reactions
				{
					std::cout << "Reading stick reactions..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("Stick");
					std::stringstream fInput;
					fInput << current_node->value();
					indices_of_stick_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					number_of_stick_reactions_ = indices_of_stick_reactions_.Size();
				}

				// Cov dependent reactions
				{
					std::cout << "Reading coverage dependent reactions..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("CoverageDependent");
					std::stringstream fInput;
					fInput << current_node->value();
					indices_of_coverage_dependent_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					number_of_coverage_dependent_reactions_ = indices_of_coverage_dependent_reactions_.Size();
				}

				// Langmuir-Hinshelwood reactions
				{
					std::cout << "Reading Langmuir-Hinshelwood reactions..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("Langmuir");
					std::stringstream fInput;
					fInput << current_node->value();
					indices_of_langmuir_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					number_of_langmuir_reactions_ = indices_of_langmuir_reactions_.Size();
				}

				// Lumped reactions
				{
					std::cout << "Reading Lumped reactions..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("Lumped");
					std::stringstream fInput;
					fInput << current_node->value();
					indices_of_lumped_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					number_of_lumped_reactions_ = indices_of_lumped_reactions_.Size();
				}

				// Reading if the kinetic scheme is conventional or UBI-QEP
				{
					std::cout << "Reading type of kinetic scheme..." << std::endl;
					rapidxml::xml_node<>* current_node = kinetics_node->first_node("TypeOfKinetics");
					std::stringstream fInput;
					fInput << current_node->value();
					std::string dummy;
					fInput >> dummy;
					
					if (dummy == "UBIQEP")						type_of_kinetics_ = TYPE_OF_KINETICS_UBI_QEP;
					else if (dummy == "chemkin_conventional")	type_of_kinetics_ = TYPE_OF_KINETICS_CHEMKIN_CONVENTIONAL;
					else
					{
						std::cout << "This type of kinetic mechanism is not available: " << dummy << std::endl;
						std::cout << "Press enter to exit..." << std::endl;						
						getchar();
						exit(OPENSMOKE_FATAL_ERROR_EXIT);
					}
				}

				// Reading UBI
				if (type_of_kinetics_ == TYPE_OF_KINETICS_UBI_QEP)
				{
					ubiqep_submechanism_ = new UBIQEP_SubMechanism();
					ubiqep_submechanism_->ReadFromXMLFile(kinetics_node, thermodynamics_.MW());
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

					std::cout << " * Reading kinetic parameters of stick reactions..." << std::endl;
					if (number_of_stick_reactions_ != 0)
					{
						stick_power_.resize(number_of_stick_reactions_);
						stick_motz_wise_.resize(number_of_stick_reactions_);
						stick_constant_coefficient_.resize(number_of_stick_reactions_);

						rapidxml::xml_node<>* stick_node =  kinetic_parameters_node->first_node("StickParameters");
						
						std::stringstream fparameters;
						fparameters << stick_node->value();
						
						for(unsigned int i=0;i<number_of_stick_reactions_;i++)
						{
							std::string dummy;
							fparameters >> dummy;
							stick_motz_wise_[i] = boost::lexical_cast<bool>(dummy);
							fparameters >> dummy;
							stick_power_[i] = boost::lexical_cast<double>(dummy);
							fparameters >> dummy;
							const unsigned int stick_gas_species = boost::lexical_cast<unsigned int>(dummy);

							fparameters >> dummy;
							const unsigned int n = boost::lexical_cast<unsigned int>(dummy);
							
							double sigma_power = 1.;
							for (unsigned int ii=0;ii<n;ii++)
							{
								fparameters >> dummy;
								const unsigned int index = boost::lexical_cast<unsigned int>(dummy);
								fparameters >> dummy;
								const double exponent = boost::lexical_cast<double>(dummy);

								
								const std::string name_site_species = thermodynamics_.vector_names_site_species()[index-1];
								for(unsigned int j=0;j<thermodynamics_.matrix_occupancies_site_species()[0].size();j++)
									for(unsigned int jj=0;jj<thermodynamics_.matrix_occupancies_site_species()[0][j].size();jj++)
										if (thermodynamics_.matrix_names_site_species()[0][j][jj] == name_site_species)
										{	
											const double sigma = thermodynamics_.matrix_occupancies_site_species()[0][j][jj];
											sigma_power *= std::pow(sigma, exponent);
											break;
										}
							}

							stick_constant_coefficient_[i] = sigma_power * 
															sqrt(PhysicalConstants::R_J_kmol/2./PhysicalConstants::pi/thermodynamics_.MW()[stick_gas_species]);
						}
					}

					std::cout << " * Reading kinetic parameters of coverage dependent reactions..." << std::endl;
					if (number_of_coverage_dependent_reactions_ != 0)
					{
						coverage_dependent_species_site_type_.resize(number_of_coverage_dependent_reactions_);
						coverage_dependent_species_index_.resize(number_of_coverage_dependent_reactions_);
						coverage_dependent_eta_.resize(number_of_coverage_dependent_reactions_);
						coverage_dependent_mu_.resize(number_of_coverage_dependent_reactions_);
						coverage_dependent_epsilon_.resize(number_of_coverage_dependent_reactions_);

						rapidxml::xml_node<>* coverage_dependent_node =  kinetic_parameters_node->first_node("CoverageDependentParameters");
						
						std::stringstream fparameters;
						fparameters << coverage_dependent_node->value();

						for(unsigned int i=0;i<number_of_coverage_dependent_reactions_;i++)
						{
							std::string dummy;
							fparameters >> dummy;
							const unsigned int n = boost::lexical_cast<unsigned int>(dummy);

							coverage_dependent_species_site_type_[i].resize(n);
							coverage_dependent_species_index_[i].resize(n);
							coverage_dependent_eta_[i].resize(n);
							coverage_dependent_mu_[i].resize(n);
							coverage_dependent_epsilon_[i].resize(n);

							for(unsigned int j=0;j<n;j++)
							{
								fparameters >> dummy;
								coverage_dependent_species_site_type_[i][j] = boost::lexical_cast<bool>(dummy);
								fparameters >> dummy;
								coverage_dependent_species_index_[i][j] = boost::lexical_cast<unsigned int>(dummy);
								fparameters >> dummy;
								coverage_dependent_eta_[i][j] = boost::lexical_cast<double>(dummy);
								fparameters >> dummy;
								coverage_dependent_mu_[i][j] = boost::lexical_cast<double>(dummy);
								fparameters >> dummy;
								coverage_dependent_epsilon_[i][j] = boost::lexical_cast<double>(dummy);
							}
						}
					}

					std::cout << " * Reading kinetic parameters of Langmuir-Hinshelwood reactions..." << std::endl;
					if (number_of_langmuir_reactions_ != 0)
					{
						langmuir_denominator_order_.resize(number_of_langmuir_reactions_);
						langmuir_units_.resize(number_of_langmuir_reactions_);
						langmuir_species_index_.resize(number_of_langmuir_reactions_);
						langmuir_numerator_species_.resize(number_of_langmuir_reactions_);
						langmuir_lnA_.resize(number_of_langmuir_reactions_);
						langmuir_Beta_.resize(number_of_langmuir_reactions_);
						langmuir_H_over_R_.resize(number_of_langmuir_reactions_);
						langmuir_order_.resize(number_of_langmuir_reactions_);

						rapidxml::xml_node<>* langmuir_node =  kinetic_parameters_node->first_node("LangmuirParameters");
						
						std::stringstream fparameters;
						fparameters << langmuir_node->value();

						for(unsigned int i=0;i<number_of_langmuir_reactions_;i++)
						{
							std::string dummy;
							fparameters >> dummy;
							const unsigned int n = boost::lexical_cast<unsigned int>(dummy);

							fparameters >> dummy;
							langmuir_denominator_order_[i] = boost::lexical_cast<double>(dummy);


							fparameters >> dummy;
							langmuir_units_[i] = 
								static_cast<PhysicalConstants::UNITS_REACTION_COMPOSITION>(boost::lexical_cast<unsigned int>(dummy));

							langmuir_species_index_[i].resize(n);
							langmuir_numerator_species_[i].resize(n);
							langmuir_lnA_[i].resize(n);
							langmuir_Beta_[i].resize(n);
							langmuir_H_over_R_[i].resize(n);
							langmuir_order_[i].resize(n);

							for(unsigned int j=0;j<n;j++)
							{
								fparameters >> dummy;
								langmuir_species_index_[i][j] = boost::lexical_cast<unsigned int>(dummy);
								
								fparameters >> dummy;
								langmuir_numerator_species_[i][j] = boost::lexical_cast<bool>(dummy);
								
								fparameters >> dummy;
								langmuir_lnA_[i][j] = boost::lexical_cast<double>(dummy);
								if (langmuir_lnA_[i][j] > 0.)
									langmuir_lnA_[i][j] = log(langmuir_lnA_[i][j]);

								fparameters >> dummy;
								langmuir_Beta_[i][j] = boost::lexical_cast<double>(dummy);
								
								fparameters >> dummy;
								langmuir_H_over_R_[i][j] = boost::lexical_cast<double>(dummy);
								
								fparameters >> dummy;
								langmuir_order_[i][j] = boost::lexical_cast<double>(dummy);
							}
						}
					}

					if (number_of_lumped_reactions_ != 0)
					{
						names_of_lumped_functions_.resize(number_of_lumped_reactions_);

						rapidxml::xml_node<>* lumped_node = kinetic_parameters_node->first_node("LumpedParameters");

						std::stringstream fparameters;
						fparameters << lumped_node->value();

						for (unsigned int i = 0; i<number_of_lumped_reactions_; i++)
							fparameters >> names_of_lumped_functions_[i];
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

				// Reactions with non conservation of sites
				{			
					std::string dummy;

					rapidxml::xml_node<>* current_node = kinetics_node->first_node("ConservationOfSites");
					std::stringstream fInput;
					fInput << current_node->value();

					fInput >> dummy;
					unsigned int number_of_reactions_with_non_conservation_of_sites = boost::lexical_cast<unsigned int>(dummy);
					non_conservation_of_sites_indices_of_reactions_.resize(number_of_reactions_with_non_conservation_of_sites);
					non_conservation_of_sites_phase_of_reactions_.resize(number_of_reactions_with_non_conservation_of_sites);
					non_conservation_of_sites_delta_sigma_.resize(number_of_reactions_with_non_conservation_of_sites);

					for (unsigned int j=0;j<number_of_reactions_with_non_conservation_of_sites;j++)
					{		
						fInput >> dummy;
						non_conservation_of_sites_indices_of_reactions_[j] = boost::lexical_cast<unsigned int>(dummy);
						fInput >> dummy;
						non_conservation_of_sites_phase_of_reactions_[j] = boost::lexical_cast<unsigned int>(dummy);
						fInput >> dummy;
						non_conservation_of_sites_delta_sigma_[j] = boost::lexical_cast<double>(dummy);
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

					delta_sigma_times_log_Gamma0_.resize(n);
					delta_nu_gas_.resize(n);

					for (unsigned int j=0;j<n;j++)
					{
						fInput >> dummy;
						const unsigned int index_phase = boost::lexical_cast<unsigned int>(dummy);
						fInput >> dummy;
						const double delta_sigma = boost::lexical_cast<double>(dummy);
						fInput >> dummy;
						delta_nu_gas_[j] = boost::lexical_cast<double>(dummy);

						delta_sigma_times_log_Gamma0_[j] = delta_sigma * log(thermodynamics_.matrix_densities_site_phases()[0][index_phase-1]);
					}				
				}

				// Stoichiometry
				{
					stoichiometry_ = new StoichiometricMap(this->number_of_species_, this->number_of_reactions_);
	
					rapidxml::xml_node<>* stoichiometry_node = kinetics_node->first_node("Stoichiometry");

					std::string stoichiometry_type = stoichiometry_node->first_attribute("type")->value();
					std::string stoichiometry_version = stoichiometry_node->first_attribute("version")->value();

					if (stoichiometry_type != "OpenSMOKE" || stoichiometry_version != "01-02-2014")
						ErrorMessage("void KineticsMap_Surface_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current stoichiometric data are not supported.");

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
					ChangeDimensions(thermodynamics_.number_of_site_species(), &cSites_, true);
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

					//	for(unsigned int k=1;k<=number_of_thirdbody_reactions_;k++)
				//		{
					//		type_of_reaction_[indices_of_thirdbody_reactions_[k]] = PhysicalConstants::REACTION_THIRDBODY;
					//		local_family_index_[indices_of_thirdbody_reactions_[k]] = k;
					//	}
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
					std::cout << "   Stick reactions:                " << number_of_stick_reactions_ << " (" << number_of_stick_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
					std::cout << "   Coverage dependent reactions:   " << number_of_coverage_dependent_reactions_ << " (" << number_of_coverage_dependent_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
					std::cout << "   Langmuir-Hinshelwood reactions: " << number_of_langmuir_reactions_ << " (" << number_of_langmuir_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
					std::cout << "   Lumped reactions:               " << number_of_lumped_reactions_ << " (" << number_of_lumped_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
					std::cout << std::endl;

					stoichiometry_->Summary(std::cout);
				}
			}
		}
	}

	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
		try
		{
			this->number_of_species_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));					
		}
		catch(...)
		{
			ErrorMessage("KineticsMap_Surface_CHEMKIN<map>::ImportSpeciesFromXMLFile", "Error in reading the number of species.");
		}
	}

	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::ReactionEnthalpiesAndEntropies()
	{
		if (reaction_h_and_s_must_be_recalculated_ == true)
		{
			stoichiometry_->ReactionEnthalpyAndEntropy(	reaction_h_over_RT_, reaction_s_over_R_, 
														thermodynamics_.species_h_over_RT(), thermodynamics_.species_s_over_R() );

			reaction_h_and_s_must_be_recalculated_ = false;
		}
	}

	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::ArrheniusKineticConstants()
	{
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
					uKeq_[k] = -reaction_s_over_R_[j] + reaction_h_over_RT_[j] - log_Patm_over_RT_ * delta_nu_gas_[k-1] - delta_sigma_times_log_Gamma0_[k-1];
				}
				Exp(uKeq_, &uKeq_);

				/* 
				// Alternative method as suggested by Lu & Law for calculating the equilibrium constants
				// This method can result in better efficiency if the number of species is much smaller than
				// the number of reversible reactions
				OpenSMOKEVectorDouble Kp_(this->number_of_reactions_);
				OpenSMOKEVectorDouble species_exp_g_over_RT_(this->number_of_species_);

				const double t3=OpenSMOKE::OpenSMOKEGetCpuTime();
				Exp(thermodynamics_.species_g_over_RT_, &species_exp_g_over_RT_);
				stoichiometry_->EquilibriumConstants(Kp_, species_exp_g_over_RT_, Patm_over_RT_);
				for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
				{
					unsigned int j = indices_of_thermodynamic_reversible_reactions_[k];
					uKeq_[k] = Kp_[j];
				}
				*/
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
	void KineticsMap_Surface_CHEMKIN<map>::ReactionRates(const OpenSMOKEVectorDouble& cGas, const OpenSMOKEVectorDouble& z, const OpenSMOKEVectorDouble& a, const OpenSMOKEVectorDouble& Gamma)
	{
		const double cTot = cGas.SumElements();

		for(int j=1;j<=z.Size();j++)
			cSites_[j] = z[j]*Gamma[thermodynamics_.vector_site_phases_belonging()[j-1]+1] / 
								thermodynamics_.vector_occupancies_site_species()[j-1];
		
		unsigned int count = 1;
		for(unsigned int j=0;j<thermodynamics_.number_of_gas_species();j++)
			c_[count++] = cGas[j+1];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_site_species();j++)
			c_[count++] = cSites_[j+1];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_bulk_species();j++)
			c_[count++] = a[j+1];
		
		if (type_of_kinetics_ == TYPE_OF_KINETICS_CHEMKIN_CONVENTIONAL)
		{
			// 1. Kinetic constants
			ArrheniusKineticConstants();
		
			// 2. Correct the effective kinetic constants by stick reactions
			if (number_of_stick_reactions_>0)
			{
				double total_site_density = Gamma.SumElements();

				for(unsigned int s=1;s<=number_of_stick_reactions_;s++)
				{
					const unsigned int j=indices_of_stick_reactions_[s];

					const double gamma = std::min(1., kArrhenius_[j]);
					if (stick_motz_wise_[s-1] == false)	kArrheniusModified_[j] = gamma;
					else                                kArrheniusModified_[j] = gamma/(1.-gamma/2.);

					kArrheniusModified_[j] *= stick_constant_coefficient_[s-1]*std::sqrt(this->T_)/std::pow(total_site_density, stick_power_[s-1]);
				}
			}

			// 3. Correct the effective kinetic constants by coverage dependent reactions
			if (number_of_coverage_dependent_reactions_>0)
			{
				for(unsigned int s=1;s<=number_of_coverage_dependent_reactions_;s++)
				{
					double correction = 0.;
					for(unsigned int k=0;k<coverage_dependent_species_site_type_[s-1].size();k++)
					{
						double value;
						if (coverage_dependent_species_site_type_[s-1][k] == true)
							value = z[coverage_dependent_species_index_[s-1][k]];
						else
							value = a[coverage_dependent_species_index_[s-1][k]];

						const double eps_value = 1.e-20;
						correction +=	PhysicalConstants::ln_10*coverage_dependent_eta_[s-1][k]*value + 
										coverage_dependent_mu_[s-1][k]*log(value+eps_value) -
										coverage_dependent_epsilon_[s-1][k]*value/(PhysicalConstants::R_J_kmol*this->T_);
						
					}
					correction = std::exp(correction);
				
					const unsigned int j=indices_of_coverage_dependent_reactions_[s];
					kArrheniusModified_[j] *= correction;
				}
			}

			// 4. Correct the effective kinetic constants by Langmuir-Hinshelwood reactions
			if (number_of_langmuir_reactions_>0)
			{
				for(unsigned int s=1;s<=number_of_langmuir_reactions_;s++)
				{
					double correction_denominator = 1.;
					double correction_numerator = 1.;
					for(unsigned int k=0;k<langmuir_species_index_[s-1].size();k++)
					{
						const double  K = std::exp( langmuir_lnA_[s-1][k] + langmuir_Beta_[s-1][k]*log(this->T_) - langmuir_H_over_R_[s-1][k]/this->T_ );

						double value = cGas[langmuir_species_index_[s-1][k]];
					
						if (langmuir_units_[s-1] != PhysicalConstants::UNITS_STD)
							value *= PhysicalConstants::R_J_kmol*this->T_;

						correction_denominator += K*std::pow(value, langmuir_order_[s-1][k]);

						if (langmuir_numerator_species_[s-1][k] == true)
							correction_numerator *= K;
					}

					const unsigned int j=indices_of_langmuir_reactions_[s];
					kArrheniusModified_[j] *= correction_numerator/std::pow(correction_denominator,langmuir_denominator_order_[s-1]);
				}
			}
		}
		
		else if (type_of_kinetics_ == TYPE_OF_KINETICS_UBI_QEP)
		{
			ubiqep_submechanism_->CalculateDissociationEnergies(thermodynamics_.species_h_over_RT());
			ubiqep_submechanism_->CalculateChemisorptionHeats(this->T_, z);
			ubiqep_submechanism_->CalculateSurfaceEnthalpies(this->T_);
			ubiqep_submechanism_->CalculateActivationEnergies();

			double total_site_density = Gamma.SumElements();
			ubiqep_submechanism_->ForwardKineticConstants(this->T_, total_site_density, kArrheniusModified_);
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

		// User defined reaction rates
		if (number_of_lumped_reactions_ != 0)
			UserDefinedReactionRates(cGas, z, a, Gamma);
	}

	template<typename map>
	void KineticsMap_Surface_CHEMKIN<map>::UserDefinedReactionRates(const OpenSMOKEVectorDouble& cGas, const OpenSMOKEVectorDouble& z, const OpenSMOKEVectorDouble& a, const OpenSMOKEVectorDouble& Gamma)
	{
		FatalErrorMessage("KineticsMap_Surface_CHEMKIN<map>::UserDefinedReactionRates: No user defined reaction rates are provided by the user!");
	}

	/*
	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::DerivativesOfReactionRatesWithRespectToKineticParameters(const PhysicalConstants::sensitivity_type type, unsigned int jReaction, const OpenSMOKEVectorDouble& c, double& parameter)
	{
		const double cTot = c.SumElements();

		// 1. Kinetic constants
		ArrheniusKineticConstants();
			
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
	void KineticsMap_Surface_CHEMKIN<map>::FormationRates(OpenSMOKEVectorDouble* Rgas, OpenSMOKEVectorDouble* Rsite, OpenSMOKEVectorDouble* Rbulk, OpenSMOKEVectorDouble* RsitePhases)
	{
		OpenSMOKE::OpenSMOKEVectorDouble R(thermodynamics_.NumberOfSpecies());
		stoichiometry_->FormationRatesFromReactionRates(&R, netReactionRates_);
		
		unsigned int count = 1;
		for(unsigned int j=1;j<=thermodynamics_.number_of_gas_species();j++)
			(*Rgas)[j] = R[count++];
		
		for(unsigned int j=1;j<=thermodynamics_.number_of_site_species();j++)
			(*Rsite)[j] = R[count++];
		
		for(unsigned int j=1;j<=thermodynamics_.number_of_bulk_species();j++)
			(*Rbulk)[j] = R[count++];
		
		*RsitePhases = 0.;
		for(unsigned int j=0;j<non_conservation_of_sites_indices_of_reactions_.size();j++)
			(*RsitePhases)[non_conservation_of_sites_phase_of_reactions_[j]] += 
				netReactionRates_[non_conservation_of_sites_indices_of_reactions_[j]] * non_conservation_of_sites_delta_sigma_[j];
	}

	template<typename map> 
	double KineticsMap_Surface_CHEMKIN<map>::HeatRelease(const OpenSMOKEVectorDouble& Rgas, const OpenSMOKEVectorDouble& Rsurface, const OpenSMOKEVectorDouble& Rbulk)
	{
		unsigned int k = 1;
		for(int j=1;j<=Rgas.Size();j++)
			aux_vector_[k++] = Rgas[j];
		for(int j=1;j<=Rsurface.Size();j++)
			aux_vector_[k++] = Rsurface[j];
		for(int j=1;j<=Rbulk.Size();j++)
			aux_vector_[k++] = Rbulk[j];

		return -Dot(aux_vector_, thermodynamics_.species_h_over_RT()) * PhysicalConstants::R_J_kmol * this->T_;
	}

	/*
	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const unsigned int k, const OpenSMOKEVectorDouble& c, OpenSMOKEVectorDouble* Jalfa, double& parameter)
	{
		DerivativesOfReactionRatesWithRespectToKineticParameters(type, k, c, parameter);
		stoichiometry_->FormationRatesFromReactionRates(Jalfa, netReactionRates_);
	}
	*/

	/*
	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& parameter)
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
	void KineticsMap_Surface_CHEMKIN<map>::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& Jrho, double& parameter)
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
	void KineticsMap_Surface_CHEMKIN<map>::ProductionAndDestructionRatesGross(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D)
	{
		OpenSMOKEVectorDouble forward_(this->number_of_reactions_);
		OpenSMOKEVectorDouble backward_(this->number_of_reactions_);

		GetForwardReactionRates(&forward_);
		GetBackwardReactionRates(&backward_);

		stoichiometry_->ProductionAndDestructionRatesFromReactionRatesGross(P, D, forward_, backward_);
	}
	*/

	/*
	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::ProductionAndDestructionRates(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D)
	{
		stoichiometry_->ProductionAndDestructionRatesFromReactionRates(P, D, netReactionRates_);
	}
	*/

	/*
	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::RateOfProductionAnalysis(const bool iNormalize) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, iNormalize);
	}
	*/
	/*
	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::RateOfProductionAnalysis(std::ostream& fout) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, false);
		stoichiometry_->WriteRateOfProductionAnalysis(fout);
	}
	*/

	/*
	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::RateOfProductionAnalysis(ROPA_Data& ropa) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, false);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}
	*/

	template<typename map> 
	const OpenSMOKEVectorDouble& KineticsMap_Surface_CHEMKIN<map>::GetReactionRates()
	{
		return netReactionRates_;
	}

	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::GetReactionRates(OpenSMOKEVectorDouble* r)
	{
		*r = netReactionRates_;
	}

	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::GetForwardReactionRates(OpenSMOKEVectorDouble* r)
	{
		ElementByElementProduct(forwardReactionRates_, kArrheniusModified_, r);
	}

	template<typename map> 
	void KineticsMap_Surface_CHEMKIN<map>::GetBackwardReactionRates(OpenSMOKEVectorDouble* r)
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
	void KineticsMap_Surface_CHEMKIN<map>::WriteKineticData(std::ostream& fOut, const unsigned int k)
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
	void KineticsMap_Surface_CHEMKIN<map>::WriteKineticData(std::ostream& fOut, const unsigned int k, OpenSMOKEVectorDouble& c_bath, const double conversion_forward, const double conversion_backward)
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

			ReactionEnthalpiesAndEntropies();
			ArrheniusKineticConstants();
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
	void KineticsMap_Surface_CHEMKIN<map>::FittedReverseKineticConstants(OpenSMOKEVectorDouble& x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters)
	{			
		const unsigned int npoints = 11;
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

				X(i, 0) = 1.;
				X(i, 1) = -1. / PhysicalConstants::R_J_kmol / T;

				if (nparameters == 3)
					X(i, 2) = log(T);
			}

			XT = X.transpose();
			XTX=XT*X;
		}

		const double patm = 101325.;
		SetPressure(patm);
		thermodynamics_.SetPressure(patm);

		Eigen::MatrixXd y(11, number_of_thermodynamic_reversible_reactions_);
		Eigen::MatrixXd Y(3, number_of_thermodynamic_reversible_reactions_);

		std::cout << "   evaluating the reaction rates..." << std::endl;
		for(int i=1;i<=temperatures.Size();i++)
		{
			const double ctot = patm/PhysicalConstants::R_J_kmol/temperatures[i];
			Product(ctot, x_bath, &c);

			SetTemperature(temperatures[i]);
			thermodynamics_.SetTemperature(temperatures[i]);
			ReactionEnthalpiesAndEntropies();
			ArrheniusKineticConstants();
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
	void KineticsMap_Surface_CHEMKIN<map>::FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters)
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
	void KineticsMap_Surface_CHEMKIN<map>::DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* dR_over_domega)
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
	void KineticsMap_Surface_CHEMKIN<map>::Derivatives(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* derivatives)
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
	void KineticsMap_Surface_CHEMKIN<map>::DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* dR_over_dC)
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
	void KineticsMap_Surface_CHEMKIN<map>::Derivatives(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* derivatives, const bool constant_density)
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
		ReactionEnthalpiesAndEntropies();
		ArrheniusKineticConstants();

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
			ReactionEnthalpiesAndEntropies();
			ArrheniusKineticConstants();

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

	void UBIQEP_SubMechanism::ReadFromXMLFile(rapidxml::xml_node<> *xml_node, const OpenSMOKE::OpenSMOKEVectorDouble& MW)
	{
		std::cout << "Reading Heats of chemisorption..." << std::endl;
		{
			std::string dummy;
			rapidxml::xml_node<>* current_node = xml_node->first_node("HeatsOfChemisorption");
			std::stringstream fInput;
			fInput << current_node->value();

			fInput >> dummy;
			number_of_gas_species_ = boost::lexical_cast<unsigned int>(dummy);

			fInput >> dummy;
			number_of_site_species_ = boost::lexical_cast<unsigned int>(dummy);

			fInput >> dummy;
			unsigned int number_of_gas_species_chemisorption_heats = boost::lexical_cast<unsigned int>(dummy);

			fInput >> dummy;
			reference_temperature_ = boost::lexical_cast<double>(dummy);

			unsigned int size = number_of_site_species_+number_of_gas_species_chemisorption_heats;
			ChangeDimensions(size, &chemisorption_heats_constant_coefficient_, true);
			ChangeDimensions(size, &chemisorption_heats_temperature_coefficient_, true);
			ChangeDimensions(size, &QStar_, true);
			ChangeDimensions(number_of_gas_species_chemisorption_heats, &chemisorption_heats_gas_indices_, true);

			chemisorption_heats_coefficients_ = new OpenSMOKE::OpenSMOKEVectorDouble[size+1];
			chemisorption_heats_indices_ = new OpenSMOKE::OpenSMOKEVectorUnsignedInt[size+1];
			for (unsigned int j=1;j<=number_of_site_species_+number_of_gas_species_chemisorption_heats;j++)
			{		
				if (j > number_of_site_species_)
				{
					fInput >> dummy;
					chemisorption_heats_gas_indices_[j-number_of_site_species_] = boost::lexical_cast<unsigned int>(dummy);
				}

				fInput >> dummy;
				chemisorption_heats_temperature_coefficient_[j] = boost::lexical_cast<double>(dummy);

				fInput >> dummy;
				chemisorption_heats_constant_coefficient_[j] = boost::lexical_cast<double>(dummy);

				fInput >> dummy;
				unsigned int n = boost::lexical_cast<unsigned int>(dummy);

				ChangeDimensions(n, &chemisorption_heats_coefficients_[j], true);
				ChangeDimensions(n, &chemisorption_heats_indices_[j], true);

				for (unsigned int k=1;k<=n;k++)
				{
					fInput >> dummy;
					chemisorption_heats_indices_[j][k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					chemisorption_heats_coefficients_[j][k] = boost::lexical_cast<double>(dummy);
				}
			}
		}

		std::cout << "Reading reaction parameters (direct reactions)" << std::endl;
		{
			std::string dummy;
			rapidxml::xml_node<>* parent_node = xml_node->first_node("UBIParameters");
			rapidxml::xml_node<>* current_node = parent_node->first_node("UBIDirect");

			std::stringstream fInput;
			fInput << current_node->value();

			fInput >> dummy;
			half_number_of_ubiqep_reactions_ = boost::lexical_cast<unsigned int>(dummy);
			
			ChangeDimensions(half_number_of_ubiqep_reactions_, &dissociation_energies_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &surface_enthalpies_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &E_forward_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &E_backward_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &adsorption_coefficient_, true);			

			ChangeDimensions(half_number_of_ubiqep_reactions_, &sigma_direct_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &Beta_direct_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &lnA_direct_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &lambda_, true);
			
			ChangeDimensions(half_number_of_ubiqep_reactions_, &ubiqep_class_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &ubiqep_type_direct_, true);

			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_A_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_B_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_C_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_D_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_Star_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_A2_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_A2_Chemisorption_Heats_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_AB_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_AB_Chemisorption_Heats_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_ABStar_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_AStar_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_BStar_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_CStar_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &index_DStar_, true);
			
			for(unsigned int k=1;k<=half_number_of_ubiqep_reactions_;k++)
			{
				fInput >> dummy;
				lnA_direct_[k] = boost::lexical_cast<double>(dummy);
				lnA_direct_[k] = lnA_direct_[k] > 0. ? log(lnA_direct_[k]) : lnA_direct_[k];
				fInput >> dummy;
				Beta_direct_[k] = boost::lexical_cast<double>(dummy);
				fInput >> dummy;
				sigma_direct_[k] = boost::lexical_cast<double>(dummy);

				fInput >> dummy;
				ubiqep_class_[k] = boost::lexical_cast<unsigned int>(dummy);
				
				fInput >> dummy;
				ubiqep_type_direct_[k] = static_cast<PhysicalConstants::UBIQEP_TYPE>(boost::lexical_cast<unsigned int>(dummy));

				if (ubiqep_type_direct_[k] == PhysicalConstants::UBIQEP_TYPE_ADSORPTION)
				{
					fInput >> dummy;
					const unsigned int gas_index = boost::lexical_cast<unsigned int>(dummy);
					adsorption_coefficient_[k] = sqrt(PhysicalConstants::R_J_kmol/2/PhysicalConstants::pi/MW[gas_index]);
				}

				if (ubiqep_class_[k] == 0)
				{	
				}
				else if (ubiqep_class_[k] == 1)
				{
					fInput >> dummy;
					index_A_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_Star_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_AStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					lambda_[k] = 1.;
				}
				else if (ubiqep_class_[k] == 2 || ubiqep_class_[k] == 3)
				{
					fInput >> dummy;
					index_A2_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_Star_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_AStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_A_[k] = boost::lexical_cast<unsigned int>(dummy);

					lambda_[k] = 2.;

					// Gas phase
					for (int j=1;j<=chemisorption_heats_gas_indices_.Size();j++)
						if (index_A2_[k] == chemisorption_heats_gas_indices_[j])
						{
							index_A2_Chemisorption_Heats_[k] = number_of_site_species_+j;
							break;
						}
				}
				else if (ubiqep_class_[k] == 4)
				{
					fInput >> dummy;
					index_AB_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_Star_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_AStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_BStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_A_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_B_[k] = boost::lexical_cast<unsigned int>(dummy);


					lambda_[k] = 2.;

					// Gas phase
					for (int j=1;j<=chemisorption_heats_gas_indices_.Size();j++)
						if (index_AB_[k] == chemisorption_heats_gas_indices_[j])
						{
							index_AB_Chemisorption_Heats_[k] = number_of_site_species_+j;
							break;
						}
				}
				else if (ubiqep_class_[k] == 5)
				{
					fInput >> dummy;
					index_ABStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_Star_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_AStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_BStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_A_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_B_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_AB_[k] = boost::lexical_cast<unsigned int>(dummy);

					lambda_[k] = 2.;
				}
				else if (ubiqep_class_[k] == 6)
				{
					fInput >> dummy;
					index_CStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_DStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_AStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_BStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_A_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_B_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_C_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_D_[k] = boost::lexical_cast<unsigned int>(dummy);

					lambda_[k] = 2.;
				}
				else if (ubiqep_class_[k] == 7)
				{
					fInput >> dummy;
					index_CStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_AStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_BStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_A_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_B_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_C_[k] = boost::lexical_cast<unsigned int>(dummy);

					lambda_[k] = 2.;
				}
				else if (ubiqep_class_[k] == 8)
				{
					fInput >> dummy;
					index_CStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_DStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_AStar_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_A_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_C_[k] = boost::lexical_cast<unsigned int>(dummy);

					fInput >> dummy;
					index_D_[k] = boost::lexical_cast<unsigned int>(dummy);

					lambda_[k] = 2.;
				}
			}
		}

		std::cout << "Reading reaction parameters (reverse reactions)" << std::endl;
		{
			std::string dummy;
			rapidxml::xml_node<>* parent_node = xml_node->first_node("UBIParameters");
			rapidxml::xml_node<>* current_node = parent_node->first_node("UBIReverse");

			std::stringstream fInput;
			fInput << current_node->value();

			fInput >> dummy;
			half_number_of_ubiqep_reactions_ = boost::lexical_cast<unsigned int>(dummy);
			
			ChangeDimensions(half_number_of_ubiqep_reactions_, &sigma_reverse_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &Beta_reverse_, true);
			ChangeDimensions(half_number_of_ubiqep_reactions_, &lnA_reverse_, true);		
			ChangeDimensions(half_number_of_ubiqep_reactions_, &ubiqep_type_reverse_, true);

			for(unsigned int k=1;k<=half_number_of_ubiqep_reactions_;k++)
			{
				fInput >> dummy;
				lnA_reverse_[k] = boost::lexical_cast<double>(dummy);
				lnA_reverse_[k] = lnA_reverse_[k] > 0. ? log(lnA_reverse_[k]) : lnA_reverse_[k];
				fInput >> dummy;
				Beta_reverse_[k] = boost::lexical_cast<double>(dummy);
				fInput >> dummy;
				sigma_reverse_[k] = boost::lexical_cast<double>(dummy);

				fInput >> dummy;	// class of reaction
				fInput >> dummy;    // 
				ubiqep_type_reverse_[k] = static_cast<PhysicalConstants::UBIQEP_TYPE>(boost::lexical_cast<unsigned int>(dummy));
			}
		}
	}

	void UBIQEP_SubMechanism::CalculateChemisorptionHeats(const double T, const OpenSMOKE::OpenSMOKEVectorDouble& Z)
	{
		const double deltaT_times_R = (T-reference_temperature_)*PhysicalConstants::R_kcal_mol;
		
		for (int j=1;j<=QStar_.Size();j++)
		{
			QStar_[j] = chemisorption_heats_constant_coefficient_[j] - 
						chemisorption_heats_temperature_coefficient_[j] * deltaT_times_R;

			for (int k=1;k<=chemisorption_heats_indices_[j].Size();k++)
				QStar_[j] += chemisorption_heats_coefficients_[j][k] * Z[chemisorption_heats_indices_[j][k]];
		}
	}

	void UBIQEP_SubMechanism::CalculateDissociationEnergies(const OpenSMOKE::OpenSMOKEVectorDouble& H)
	{
		for (unsigned int k=1;k<=half_number_of_ubiqep_reactions_;k++)
		{
			switch(ubiqep_class_[k])
			{
			case 1:
				dissociation_energies_[k] = 0;
				
				break;
			case 2:
				dissociation_energies_[k] = 2.*H[index_A_[k]] - H[index_A2_[k]];
				break;
			case 3:
				dissociation_energies_[k] = 2.*H[index_A_[k]] - H[index_A2_[k]];
				break;
			case 4:
				dissociation_energies_[k] = H[index_A_[k]] + H[index_B_[k]] - H[index_AB_[k]];
				break;
			case 5:
				dissociation_energies_[k] = H[index_A_[k]] + H[index_B_[k]] - H[index_AB_[k]];
				break;
			case 6:
				dissociation_energies_[k] = H[index_A_[k]] + H[index_B_[k]] - H[index_C_[k]]- H[index_D_[k]];
				break;
			case 7:
				dissociation_energies_[k] = H[index_A_[k]] + H[index_B_[k]] - 2.*H[index_C_[k]];
				break;
			case 8:
				dissociation_energies_[k] = 2.*H[index_A_[k]] - H[index_C_[k]] - H[index_D_[k]];
				break;
			}

			dissociation_energies_[k] = std::fabs(dissociation_energies_[k]);
		}
	}

	void UBIQEP_SubMechanism::CalculateSurfaceEnthalpies(const double T)
	{
		const double R_times_T = PhysicalConstants::R_kcal_mol * T;

		for (unsigned int k=1;k<=half_number_of_ubiqep_reactions_;k++)
		{
			switch(ubiqep_class_[k])
			{
			case 1:
				surface_enthalpies_[k] = ( QStar_[index_Star_[k]] - QStar_[index_AStar_[k]] );
				break;
			case 2:
				surface_enthalpies_[k] = 2. * ( QStar_[index_Star_[k]] - QStar_[index_AStar_[k]] );
				break;
			case 3:
				surface_enthalpies_[k] = 2. * ( QStar_[index_Star_[k]] - QStar_[index_AStar_[k]] );
				break;
			case 4:
				surface_enthalpies_[k] = 2.*QStar_[index_Star_[k]] - (QStar_[index_AStar_[k]] + QStar_[index_BStar_[k]]);
				break;
			case 5:
				surface_enthalpies_[k] = QStar_[index_ABStar_[k]] + QStar_[index_Star_[k]] -
					                     (QStar_[index_AStar_[k]] + QStar_[index_BStar_[k]] );
				break;
			case 6:
				surface_enthalpies_[k] = QStar_[index_CStar_[k]] + QStar_[index_DStar_[k]] -
					                     (QStar_[index_AStar_[k]] + QStar_[index_BStar_[k]] );
				break;
			case 7:
				surface_enthalpies_[k] = 2.*QStar_[index_CStar_[k]] -
					                     (QStar_[index_AStar_[k]] + QStar_[index_BStar_[k]] );
				break;
			case 8:
				surface_enthalpies_[k] = QStar_[index_CStar_[k]] + QStar_[index_DStar_[k]] - 
										 2.*QStar_[index_AStar_[k]];
				break;
			}

			surface_enthalpies_[k] += dissociation_energies_[k] * R_times_T;
		}
	}

	void UBIQEP_SubMechanism::CalculateActivationEnergies()
	{
		for (unsigned int k = 1; k <= half_number_of_ubiqep_reactions_; k++)
		{
			switch (ubiqep_class_[k])
			{
			case 1:
				E_forward_[k] = 0.;
				break;
			case 2:
				E_forward_[k] = 0.;
				break;
			case 3:
				E_forward_[k] = sigma_direct_[k] * (surface_enthalpies_[k] +
					0.5*QStar_[index_AStar_[k]] - QStar_[index_A2_Chemisorption_Heats_[k]]);
				break;
			case 4:
				E_forward_[k] = sigma_direct_[k] * (surface_enthalpies_[k] +
					QStar_[index_AStar_[k]] * QStar_[index_BStar_[k]] / (QStar_[index_AStar_[k]] + QStar_[index_BStar_[k]]) -
					QStar_[index_AB_Chemisorption_Heats_[k]]);
				break;
			case 5:
				E_forward_[k] = sigma_direct_[k] * (surface_enthalpies_[k] +
					QStar_[index_AStar_[k]] * QStar_[index_BStar_[k]] / (QStar_[index_AStar_[k]] + QStar_[index_BStar_[k]]));
				break;
			case 6:
				E_forward_[k] = sigma_direct_[k] * (surface_enthalpies_[k] +
					QStar_[index_AStar_[k]] * QStar_[index_BStar_[k]] / (QStar_[index_AStar_[k]] + QStar_[index_BStar_[k]]));
				break;
			case 7:
				E_forward_[k] = sigma_direct_[k] * (surface_enthalpies_[k] +
					QStar_[index_AStar_[k]] * QStar_[index_BStar_[k]] / (QStar_[index_AStar_[k]] + QStar_[index_BStar_[k]]));
				break;
			case 8:
				E_forward_[k] = sigma_direct_[k] * (surface_enthalpies_[k] + 0.5*QStar_[index_AStar_[k]]);
				break;
			}

			E_backward_[k] = E_forward_[k] - surface_enthalpies_[k];

			// Activation energy check
			if (E_forward_[k] < 0.)
			{
				E_forward_[k] = 0.;
				E_backward_[k] = std::fabs(surface_enthalpies_[k]);
			}
			else if (E_backward_[k] < 0.)
			{
				E_backward_[k] = 0.;
				E_forward_[k] = std::fabs(surface_enthalpies_[k]);
			}
		}
	}

	void UBIQEP_SubMechanism::ForwardKineticConstants(const double T, const double total_site_density, OpenSMOKE::OpenSMOKEVectorDouble& kForward)
	{
		const double R_times_T = PhysicalConstants::R_kcal_mol * T;
		const double ln_T_over_T0 = log(T/reference_temperature_);
		const double ln_density = log(total_site_density);
		const double sqrt_T = sqrt(T);

		unsigned int j=1;
		for (unsigned int k=1;k<=half_number_of_ubiqep_reactions_;k++)
		{
			if (ubiqep_type_direct_[k] == PhysicalConstants::UBIQEP_TYPE_ADSORPTION)
				kForward[j]   = std::exp( lnA_direct_[k]  -lambda_[k]*ln_density + Beta_direct_[k]*ln_T_over_T0 - E_forward_[k]/R_times_T ) * adsorption_coefficient_[k] * sqrt_T;
			else
				kForward[j]   = std::exp( lnA_direct_[k]  + (1.-lambda_[k])*ln_density + Beta_direct_[k]*ln_T_over_T0 - E_forward_[k]/R_times_T );

			if (ubiqep_type_reverse_[k] == PhysicalConstants::UBIQEP_TYPE_ADSORPTION)
				kForward[j+1] = std::exp( lnA_reverse_[k] -lambda_[k]*ln_density + Beta_reverse_[k]*ln_T_over_T0 - E_backward_[k]/R_times_T ) * adsorption_coefficient_[k] * sqrt_T;
			else
				kForward[j+1] = std::exp( lnA_reverse_[k] + (1.-lambda_[k])*ln_density + Beta_reverse_[k]*ln_T_over_T0 - E_backward_[k]/R_times_T );

			j+=2;
		}
	}
}

