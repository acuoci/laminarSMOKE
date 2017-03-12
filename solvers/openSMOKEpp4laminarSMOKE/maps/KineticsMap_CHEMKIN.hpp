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
	KineticsMap_CHEMKIN::KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN& thermo, const unsigned int nSpecies, const unsigned int nPoints) :
	thermodynamics_(thermo)
	{
		this->number_of_species_ = nSpecies;
		this->number_of_points_ = nPoints;
        this->verbose_output_ = true;
        this->isJacobianSparsityMapAvailable_ = false;
                
		this->T_ = this->P_ = 0.;
	}
 
	KineticsMap_CHEMKIN::KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN& thermo, rapidxml::xml_document<>& doc, const unsigned int nPoints) :
	thermodynamics_(thermo)
	{
		this->number_of_points_ = nPoints;
                
        this->verbose_output_ = true;
        this->isJacobianSparsityMapAvailable_ = false;
                
		ImportSpeciesFromXMLFile(doc);
		ImportCoefficientsFromXMLFile(doc);
		this->T_ = this->P_ = 0.;
	}
        
	KineticsMap_CHEMKIN::KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN& thermo, rapidxml::xml_document<>& doc, bool verbose, const unsigned int nPoints) :
	thermodynamics_(thermo)
	{
		this->number_of_points_ = nPoints;
                
        this->verbose_output_ = verbose;
        this->isJacobianSparsityMapAvailable_ = false;
                
		ImportSpeciesFromXMLFile(doc);
		ImportCoefficientsFromXMLFile(doc);
		this->T_ = this->P_ = 0.;
	}
         
	KineticsMap_CHEMKIN::KineticsMap_CHEMKIN( const KineticsMap_CHEMKIN& rhs) :
    thermodynamics_(rhs.thermodynamics_)
    {
        CopyFromMap(rhs);
    }
        
	KineticsMap_CHEMKIN::KineticsMap_CHEMKIN( const KineticsMap_CHEMKIN& rhs, ThermodynamicsMap_CHEMKIN& thermo) :
    thermodynamics_(thermo)
    {
        CopyFromMap(rhs);
    }
        
	void KineticsMap_CHEMKIN::CopyFromMap( const KineticsMap_CHEMKIN& rhs)
    {
        this->number_of_species_ = rhs.number_of_species_;
        this->number_of_reactions_ = rhs.number_of_reactions_;
		this->number_of_points_ = rhs.number_of_points_;
            
        this->reaction_entropy_over_R_ = rhs.reaction_entropy_over_R_ ;
        this->reaction_enthalpy_over_RT_ = rhs.reaction_enthalpy_over_RT_ ;

        this->indices_of_irreversible_reactions_ = rhs.indices_of_irreversible_reactions_ ;
        this->indices_of_reversible_reactions_ = rhs.indices_of_reversible_reactions_ ;
        this->indices_of_thermodynamic_reversible_reactions_ = rhs.indices_of_thermodynamic_reversible_reactions_ ;
        this->indices_of_explicitly_reversible_reactions_ = rhs.indices_of_explicitly_reversible_reactions_ ;
        this->indices_of_thirdbody_reactions_ = rhs.indices_of_thirdbody_reactions_ ;
        this->indices_of_falloff_reactions_ = rhs.indices_of_falloff_reactions_ ;
		this->indices_of_extendedfalloff_reactions_ = rhs.indices_of_extendedfalloff_reactions_;
        this->indices_of_cabr_reactions_ = rhs.indices_of_cabr_reactions_ ;
        this->indices_of_chebyshev_reactions_ = rhs.indices_of_chebyshev_reactions_ ;
        this->indices_of_pressurelog_reactions_ = rhs.indices_of_pressurelog_reactions_ ;
		this->indices_of_extendedpressurelog_reactions_ = rhs.indices_of_extendedpressurelog_reactions_;
        this->indices_of_fit1_reactions_ = rhs.indices_of_fit1_reactions_ ;
        this->indices_of_janevlanger_reactions_ = rhs.indices_of_janevlanger_reactions_ ;
        this->indices_of_landauteller_reactions_ = rhs.indices_of_landauteller_reactions_ ;

        this->indices_of_falloff_lindemann_reactions_ = rhs.indices_of_falloff_lindemann_reactions_ ;
        this->indices_of_cabr_lindemann_reactions_ = rhs.indices_of_cabr_lindemann_reactions_ ;
        this->indices_of_falloff_troe_reactions_ = rhs.indices_of_falloff_troe_reactions_ ;
        this->indices_of_cabr_troe_reactions_ = rhs.indices_of_cabr_troe_reactions_ ;
        this->indices_of_falloff_sri_reactions_ = rhs.indices_of_falloff_sri_reactions_ ;
        this->indices_of_cabr_sri_reactions_ = rhs.indices_of_cabr_sri_reactions_ ;

        this->number_of_irreversible_reactions_ = rhs.number_of_irreversible_reactions_ ;
        this->number_of_reversible_reactions_ = rhs.number_of_reversible_reactions_ ;
        this->number_of_thermodynamic_reversible_reactions_ = rhs.number_of_thermodynamic_reversible_reactions_ ;
        this->number_of_explicitly_reversible_reactions_ = rhs.number_of_explicitly_reversible_reactions_ ;
        this->number_of_thirdbody_reactions_ = rhs.number_of_thirdbody_reactions_ ;
		this->number_of_falloff_reactions_ = rhs.number_of_falloff_reactions_;
		this->number_of_extendedfalloff_reactions_ = rhs.number_of_extendedfalloff_reactions_;
		this->number_of_cabr_reactions_ = rhs.number_of_cabr_reactions_ ;
        this->number_of_chebyshev_reactions_ = rhs.number_of_chebyshev_reactions_ ;
        this->number_of_pressurelog_reactions_ = rhs.number_of_pressurelog_reactions_ ;
		this->number_of_extendedpressurelog_reactions_ = rhs.number_of_extendedpressurelog_reactions_;
        this->number_of_fit1_reactions_ = rhs.number_of_fit1_reactions_ ;
        this->number_of_janevlanger_reactions_ = rhs.number_of_janevlanger_reactions_ ;
        this->number_of_landauteller_reactions_ = rhs.number_of_landauteller_reactions_ ;

        this->number_of_falloff_lindemann_reactions_ = rhs.number_of_falloff_lindemann_reactions_ ;
        this->number_of_cabr_lindemann_reactions_ = rhs.number_of_cabr_lindemann_reactions_ ;
        this->number_of_falloff_troe_reactions_ = rhs.number_of_falloff_troe_reactions_ ;
        this->number_of_cabr_troe_reactions_ = rhs.number_of_cabr_troe_reactions_ ;
        this->number_of_falloff_sri_reactions_ = rhs.number_of_falloff_sri_reactions_ ;
        this->number_of_cabr_sri_reactions_ = rhs.number_of_cabr_sri_reactions_ ;

        this->verbose_output_ = rhs.verbose_output_ ;

        this->lnA_ = rhs.lnA_ ;
        this->Beta_ = rhs.Beta_ ;
        this->E_over_R_ = rhs.E_over_R_ ;
		this->negative_lnA_ = rhs.negative_lnA_;
		this->sign_lnA_ = rhs.sign_lnA_;

        this->lnA_reversible_ = rhs.lnA_reversible_ ;
        this->Beta_reversible_ = rhs.Beta_reversible_ ;
        this->E_over_R_reversible_ = rhs.E_over_R_reversible_ ;

        this->Meff_ = rhs.Meff_ ;
        this->indices_of_thirdbody_species_ = rhs.indices_of_thirdbody_species_ ;
        this->indices_of_thirdbody_efficiencies_ = rhs.indices_of_thirdbody_efficiencies_ ;

        this->lnA_falloff_inf_ = rhs.lnA_falloff_inf_ ;
        this->Beta_falloff_inf_ = rhs.Beta_falloff_inf_ ;
        this->E_over_R_falloff_inf_ = rhs.E_over_R_falloff_inf_ ;

        this->falloff_indices_of_thirdbody_species_ = rhs.falloff_indices_of_thirdbody_species_ ;
        this->falloff_indices_of_thirdbody_efficiencies_ = rhs.falloff_indices_of_thirdbody_efficiencies_ ;
        this->falloff_index_of_single_thirdbody_species_ = rhs.falloff_index_of_single_thirdbody_species_ ;

        this->falloff_reaction_type_ = rhs.falloff_reaction_type_ ;
        this->a_falloff_ = rhs.a_falloff_ ;
        this->b_falloff_ = rhs.b_falloff_ ;
        this->c_falloff_ = rhs.c_falloff_ ;
        this->d_falloff_ = rhs.d_falloff_ ;
        this->e_falloff_ = rhs.e_falloff_ ;
        this->logFcent_falloff_ = rhs.logFcent_falloff_ ;

        this->lnA_cabr_inf_ = rhs.lnA_cabr_inf_ ;
        this->Beta_cabr_inf_ = rhs.Beta_cabr_inf_ ;
        this->E_over_R_cabr_inf_ = rhs.E_over_R_cabr_inf_ ;

        this->cabr_indices_of_thirdbody_species_ = rhs.cabr_indices_of_thirdbody_species_ ;
        this->cabr_indices_of_thirdbody_efficiencies_ = rhs.cabr_indices_of_thirdbody_efficiencies_ ;
        this->cabr_index_of_single_thirdbody_species_ = rhs.cabr_index_of_single_thirdbody_species_ ;

        this->cabr_reaction_type_ = rhs.cabr_reaction_type_ ;
        this->a_cabr_ = rhs.a_cabr_ ;
        this->b_cabr_ = rhs.b_cabr_ ;
        this->c_cabr_ = rhs.c_cabr_ ;
        this->d_cabr_ = rhs.d_cabr_ ;
        this->e_cabr_ = rhs.e_cabr_ ;
        this->logFcent_cabr_ = rhs.logFcent_cabr_ ;

        this->changeOfMoles_ = rhs.changeOfMoles_ ;
            
        stoichiometry_ = new StoichiometricMap(*rhs.stoichiometry_);

        this->arrhenius_kinetic_constants_must_be_recalculated_ = true;
        this->nonconventional_kinetic_constants_must_be_recalculated_ = true;
        this->reaction_h_and_s_must_be_recalculated_ = true;
        this->isJacobianSparsityMapAvailable_ = false;

        ChangeDimensions(rhs.reaction_s_over_R_.Size(), &reaction_s_over_R_, true);
        ChangeDimensions(rhs.reaction_h_over_RT_.Size(), &reaction_h_over_RT_, true);
        ChangeDimensions(rhs.kArrheniusModified_.Size(), &kArrheniusModified_, true);
        ChangeDimensions(rhs.kArrhenius_.Size(), &kArrhenius_, true);
        ChangeDimensions(rhs.uKeq_.Size(), &uKeq_, true);
        ChangeDimensions(rhs.kArrhenius_reversible_.Size(), &kArrhenius_reversible_, true);
        ChangeDimensions(rhs.kArrhenius_falloff_inf_.Size(), &kArrhenius_falloff_inf_, true);
        ChangeDimensions(rhs.kArrhenius_cabr_inf_.Size(), &kArrhenius_cabr_inf_, true);

        ChangeDimensions(rhs.forwardReactionRates_.Size(), &forwardReactionRates_, true);
        ChangeDimensions(rhs.reverseReactionRates_.Size(), &reverseReactionRates_, true);
        ChangeDimensions(rhs.netReactionRates_.Size(), &netReactionRates_, true);

        ChangeDimensions(rhs.correction_falloff_.Size(), &correction_falloff_, true);
        ChangeDimensions(rhs.correction_cabr_.Size(), &correction_cabr_, true);
            
        this->chebyshev_reactions_ = new ChebyshevPolynomialRateExpression[this->number_of_chebyshev_reactions_];
		for(unsigned int j=0;j<this->number_of_chebyshev_reactions_;j++)
            this->chebyshev_reactions_[j] = rhs.chebyshev_reactions_[j];
            
        this->pressurelog_reactions_ = new PressureLogarithmicRateExpression[this->number_of_pressurelog_reactions_];
		for(unsigned int j=0;j<this->number_of_pressurelog_reactions_;j++)
            this->pressurelog_reactions_[j] = rhs.pressurelog_reactions_[j];

		this->extendedpressurelog_reactions_ = new ExtendedPressureLogarithmicRateExpression[this->number_of_extendedpressurelog_reactions_];
		for (unsigned int j = 0; j<this->number_of_extendedpressurelog_reactions_; j++)
			this->extendedpressurelog_reactions_[j] = rhs.extendedpressurelog_reactions_[j];

		this->extendedfalloff_reactions_ = new ExtendedFallOff[this->number_of_extendedfalloff_reactions_];
		for (unsigned int j = 0; j<this->number_of_extendedfalloff_reactions_; j++)
			this->extendedfalloff_reactions_[j] = rhs.extendedfalloff_reactions_[j];

        //this->fit1_reactions_ = new Fit1RateExpression[this->number_of_fit1_reactions_];
		//for(unsigned int j=0;j<this->number_of_fit1_reactions_;j++)
        //    this->fit1_reactions_[j] = rhs.fit1_reactions_[j];
            
        //this->janevlanger_reactions_ = new JanevLangerRateExpression[this->number_of_janevlanger_reactions_];
		//for(unsigned int j=0;j<this->number_of_janevlanger_reactions_;j++)
        //    this->janevlanger_reactions_[j] = rhs.janevlanger_reactions_[j];         
            
        //this->landauteller_reactions_ = new JanevLangerRateExpression[this->number_of_landauteller_reactions_];
		//for(unsigned int j=0;j<this->number_of_landauteller_reactions_;j++)
        //    this->landauteller_reactions_[j] = rhs.landauteller_reactions_[j];    
                        
        this->isThermodynamicallyReversible_ = rhs.isThermodynamicallyReversible_ ;
        this->isExplicitlyReversible_ = rhs.isExplicitlyReversible_ ;

        this->type_of_reaction_ = rhs.type_of_reaction_ ;
        this->local_family_index_ = rhs.local_family_index_ ;
    }

	KineticsMap_CHEMKIN::~KineticsMap_CHEMKIN()
	{
		indices_of_thirdbody_species_.clear();			
		indices_of_thirdbody_efficiencies_.clear();		
		falloff_indices_of_thirdbody_species_.clear();
		falloff_indices_of_thirdbody_efficiencies_.clear();
		cabr_indices_of_thirdbody_species_.clear();
		cabr_indices_of_thirdbody_efficiencies_.clear();
                
        if (number_of_chebyshev_reactions_ > 0)
            delete [] chebyshev_reactions_;
        if (number_of_pressurelog_reactions_ > 0)
            delete [] pressurelog_reactions_;
		if (number_of_extendedpressurelog_reactions_ > 0)
			delete[] extendedpressurelog_reactions_;
		if (number_of_extendedfalloff_reactions_ > 0)
			delete[] extendedfalloff_reactions_;
		delete stoichiometry_;
        if(isJacobianSparsityMapAvailable_)
        	delete jacobian_sparsity_pattern_map_;
	}

	void KineticsMap_CHEMKIN::SetTemperature(const double& T)
	{
		this->T_old_ = this->T_;
		this->T_ = T;
		this->uT_ = 1./this->T_;
		this->logT_ = std::log(this->T_);
		Patm_over_RT_ = 101325./PhysicalConstants::R_J_kmol/this->T_;
		log_Patm_over_RT_ = std::log(Patm_over_RT_);

		const double epsilon = 0.;
		if (std::fabs(this->T_-this->T_old_)/this->T_>epsilon)
		{
			arrhenius_kinetic_constants_must_be_recalculated_ = true;
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
			reaction_h_and_s_must_be_recalculated_ = true;
		}
	}

	void KineticsMap_CHEMKIN::SetPressure(const double& P)
	{
		this->P_old_ = this->P_;
		this->P_ = P;
	//	if (std::fabs(this->P_-this->P_old_)/this->P_>1.e-14)
		{
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
		}
	}

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
 
	void KineticsMap_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* kinetics_node = opensmoke_node->first_node("Kinetics");
		
		
		if (kinetics_node == 0)
			ErrorMessage("void KineticsMap_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "Kinetics tag was not found!");

		std::string kinetics_type = kinetics_node->first_attribute("type")->value();
		std::string kinetics_version = kinetics_node->first_attribute("version")->value();

		if (kinetics_type != "OpenSMOKE" || kinetics_version != "04-22-2013")
			ErrorMessage("void KineticsMap_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current kinetic scheme is not supported.");

		// Reading the number of species
		{
			rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
			this->number_of_species_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));
		}

		// Reading the number of reactions
		{
			rapidxml::xml_node<>* number_of_reactions_node = kinetics_node->first_node("NumberOfReactions");
			this->number_of_reactions_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_reactions_node->value())));
		}

		// Irreversible reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("Irreversible");
			std::stringstream fInput;
			fInput << current_node->value();
			indices_of_irreversible_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_irreversible_reactions_ = indices_of_irreversible_reactions_.Size();
		}
		
		// Reversible reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible");
			std::stringstream fInput;
			fInput << current_node->value();
			indices_of_reversible_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_reversible_reactions_ = indices_of_reversible_reactions_.Size();
		}

		// Thermodynamic Reversible reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible-Thermodynamics");
			std::stringstream fInput;
			fInput << current_node->value();
			indices_of_thermodynamic_reversible_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_thermodynamic_reversible_reactions_ = indices_of_thermodynamic_reversible_reactions_.Size();
		}

		// Explicit Reversible reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible-Explicit");
			std::stringstream fInput;
			fInput << current_node->value();
			indices_of_explicitly_reversible_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_explicitly_reversible_reactions_ = indices_of_explicitly_reversible_reactions_.Size();
		}

		// Three-body reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("ThreeBody");
			std::stringstream fInput;
			fInput << current_node->value();
			indices_of_thirdbody_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_thirdbody_reactions_  = indices_of_thirdbody_reactions_.Size();
		}

		// FallOff
		{
			rapidxml::xml_node<>* falloff_node = kinetics_node->first_node("FallOff");

			{
				std::stringstream fInput;
				fInput << falloff_node->value();
				indices_of_falloff_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_falloff_reactions_  = indices_of_falloff_reactions_.Size();
			}

			// FallOff Lindemann
			{
				rapidxml::xml_node<>* current_node = falloff_node->first_node("Lindemann");
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_falloff_lindemann_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_falloff_lindemann_reactions_  = indices_of_falloff_lindemann_reactions_.Size();
			}

			// FallOff Troe
			{
				rapidxml::xml_node<>* current_node = falloff_node->first_node("Troe");
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_falloff_troe_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_falloff_troe_reactions_  = indices_of_falloff_troe_reactions_.Size();
			}

			// FallOff SRI
			{
				rapidxml::xml_node<>* current_node = falloff_node->first_node("SRI");
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_falloff_sri_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_falloff_sri_reactions_  = indices_of_falloff_sri_reactions_.Size();
			}
		}	

		// CABR
		{
			rapidxml::xml_node<>* cabr_node = kinetics_node->first_node("CABR");

			{
				std::stringstream fInput;
				fInput << cabr_node->value();
				indices_of_cabr_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_cabr_reactions_  = indices_of_cabr_reactions_.Size();
			}

			// CABR Lindemann
			{
				rapidxml::xml_node<>* current_node = cabr_node->first_node("Lindemann");
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_cabr_lindemann_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_cabr_lindemann_reactions_  = indices_of_cabr_lindemann_reactions_.Size();
			}

			// CABR Troe
			{
				rapidxml::xml_node<>* current_node = cabr_node->first_node("Troe");
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_cabr_troe_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_cabr_troe_reactions_  = indices_of_cabr_troe_reactions_.Size();
			}

			// CABR SRI
			{
				rapidxml::xml_node<>* current_node = cabr_node->first_node("SRI");
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_cabr_sri_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_cabr_sri_reactions_  = indices_of_cabr_sri_reactions_.Size();
			}
		}	

		// Chebyshev reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("Chebyshev");
			number_of_chebyshev_reactions_ = 0;

			if (current_node != 0)
			{
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_chebyshev_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_chebyshev_reactions_ = indices_of_chebyshev_reactions_.Size();
			}
		}

		// PressureLog reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("PressureLog");
			number_of_pressurelog_reactions_ = 0;
			
			if (current_node != 0)
			{
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_pressurelog_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_pressurelog_reactions_ = indices_of_pressurelog_reactions_.Size();
			}
		}

		// ExtendedPressureLog reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("ExtendedPressureLog");
			number_of_extendedpressurelog_reactions_ = 0;
			
			if (current_node != 0)
			{
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_extendedpressurelog_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_extendedpressurelog_reactions_ = indices_of_extendedpressurelog_reactions_.Size();
			}
		}

		// ExtendedFalloff reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("ExtendedFallOff");
			number_of_extendedfalloff_reactions_ = 0;
			
			if (current_node != 0)
			{
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_extendedfalloff_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_extendedfalloff_reactions_ = indices_of_extendedfalloff_reactions_.Size();
			}
		}

		// FIT1 reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("FIT1");
			number_of_fit1_reactions_ = 0;
			
			if (current_node != 0)
			{
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_fit1_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_fit1_reactions_ = indices_of_fit1_reactions_.Size();
			}
		}

		// JanevLanger reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("JanevLanger");
			number_of_janevlanger_reactions_ = 0;

			if (current_node != 0)
			{
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_janevlanger_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_janevlanger_reactions_ = indices_of_janevlanger_reactions_.Size();
			}
		}

		// LandauTeller reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("LandauTeller");
			number_of_landauteller_reactions_ = 0;

			if (current_node != 0)
			{
				std::stringstream fInput;
				fInput << current_node->value();
				indices_of_landauteller_reactions_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_landauteller_reactions_ = indices_of_landauteller_reactions_.Size();
			}
		}

        if(verbose_output_ == true)
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

				// negative-lnA
				{
					// List of reactions with negative frequency factors
					rapidxml::xml_node<>* current_node = direct_node->first_node("negative-lnA");
					if (current_node != 0)
					{
						std::stringstream fInput;
						fInput << current_node->value();
						negative_lnA_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
					}

					// Sign of frequency factors
					ChangeDimensions(lnA_.Size(), &sign_lnA_, true);
					sign_lnA_ = 1;
					for (int i = 1; i <= negative_lnA_.Size(); i++)
						sign_lnA_[negative_lnA_[i]] = -1;
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
                        
                        if(verbose_output_ == true)
			    std::cout << " * Reading kinetic parameters of third body reactions..." << std::endl;
			if (number_of_thirdbody_reactions_ != 0)
			{
				indices_of_thirdbody_species_.resize(number_of_thirdbody_reactions_);
				indices_of_thirdbody_efficiencies_.resize(number_of_thirdbody_reactions_);

				rapidxml::xml_node<>* threebody_node = kinetic_parameters_node->first_node("ThreeBody");
			
				unsigned int i = 0;
				for (rapidxml::xml_node<> *current_node = threebody_node->first_node("reaction"); current_node; current_node = current_node->next_sibling())
				{
					std::stringstream fInput;
					fInput << current_node->value();

					indices_of_thirdbody_species_[i].Load(fInput, OPENSMOKE_FORMATTED_FILE);
					indices_of_thirdbody_efficiencies_[i].Load(fInput, OPENSMOKE_FORMATTED_FILE);
					indices_of_thirdbody_efficiencies_[i] -= 1.;

					i++;
				}
			}

                        if(verbose_output_ == true)
			    std::cout << " * Reading kinetic parameters of pressure-dependent reactions..." << std::endl;
			if (number_of_falloff_reactions_ != 0)
			{
				rapidxml::xml_node<>* falloff_node = kinetic_parameters_node->first_node("FallOff");

				// lnA
				{
					rapidxml::xml_node<>* current_node = falloff_node->first_node("lnA");
					std::stringstream fInput;
					fInput << current_node->value();
					lnA_falloff_inf_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Beta
				{
					rapidxml::xml_node<>* current_node = falloff_node->first_node("Beta");
					std::stringstream fInput;
					fInput << current_node->value();
					Beta_falloff_inf_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// E_over_R
				{
					rapidxml::xml_node<>* current_node = falloff_node->first_node("E_over_R");
					std::stringstream fInput;
					fInput << current_node->value();
					E_over_R_falloff_inf_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Parameters
				{
					ChangeDimensions(number_of_falloff_reactions_, &falloff_reaction_type_, true);
					ChangeDimensions(number_of_falloff_reactions_, &falloff_index_of_single_thirdbody_species_, true);

					ChangeDimensions(number_of_falloff_reactions_, &a_falloff_, true);
					ChangeDimensions(number_of_falloff_reactions_, &b_falloff_, true);
					ChangeDimensions(number_of_falloff_reactions_, &c_falloff_, true);
					ChangeDimensions(number_of_falloff_reactions_, &d_falloff_, true);
					ChangeDimensions(number_of_falloff_reactions_, &e_falloff_, true);

					falloff_indices_of_thirdbody_species_.resize(number_of_falloff_reactions_);
					falloff_indices_of_thirdbody_efficiencies_.resize(number_of_falloff_reactions_);

					rapidxml::xml_node<>* parameters_node = falloff_node->first_node("Parameters");

					unsigned int i = 1;
					for (rapidxml::xml_node<> *current_node = parameters_node->first_node("reaction"); current_node; current_node = current_node->next_sibling())
					{
						std::stringstream fInput;
						fInput << current_node->value();

						std::string dummy;
						fInput >> dummy;
						
						if (dummy == "lindemann")
						{
							falloff_reaction_type_[i] = PhysicalConstants::REACTION_LINDEMANN_FALLOFF;
						}
						if (dummy == "troe")
						{
							falloff_reaction_type_[i] = PhysicalConstants::REACTION_TROE_FALLOFF;

							OpenSMOKEVectorDouble coefficients;
							coefficients.Load(fInput, OPENSMOKE_FORMATTED_FILE);
							a_falloff_[i] = coefficients[1];
							b_falloff_[i] = coefficients[2];
							c_falloff_[i] = coefficients[3];
							if (coefficients.Size() == 4)
								d_falloff_[i] = coefficients[4];
						}
						else if (dummy == "sri")
						{
							falloff_reaction_type_[i] = PhysicalConstants::REACTION_SRI_FALLOFF;
					
							OpenSMOKEVectorDouble coefficients;
							coefficients.Load(fInput, OPENSMOKE_FORMATTED_FILE);
							a_falloff_[i] = coefficients[1];
							b_falloff_[i] = coefficients[2];
							c_falloff_[i] = coefficients[3];
							if (coefficients.Size() == 5)
							{
								d_falloff_[i] = coefficients[4];
								e_falloff_[i] = coefficients[5];
							}
							else
							{
								d_falloff_[i] = 1.;
								e_falloff_[i] = 0.;
							}
						}

						fInput >> dummy;
						if (dummy == "species")
						{
							unsigned int species;
							fInput >> species;
							falloff_index_of_single_thirdbody_species_[i] = species;
						}
						else
						{
							falloff_indices_of_thirdbody_species_[i-1].Load(fInput, OPENSMOKE_FORMATTED_FILE);
							falloff_indices_of_thirdbody_efficiencies_[i-1].Load(fInput, OPENSMOKE_FORMATTED_FILE);
							falloff_indices_of_thirdbody_efficiencies_[i-1] -= 1.;
						}

						i++;
					}
				}
			}

			if (number_of_cabr_reactions_ != 0)
			{
				rapidxml::xml_node<>* cabr_node = kinetic_parameters_node->first_node("CABR");

				// lnA
				{
					rapidxml::xml_node<>* current_node = cabr_node->first_node("lnA");
					std::stringstream fInput;
					fInput << current_node->value();
					lnA_cabr_inf_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Beta
				{
					rapidxml::xml_node<>* current_node = cabr_node->first_node("Beta");
					std::stringstream fInput;
					fInput << current_node->value();
					Beta_cabr_inf_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// E_over_R
				{
					rapidxml::xml_node<>* current_node = cabr_node->first_node("E_over_R");
					std::stringstream fInput;
					fInput << current_node->value();
					E_over_R_cabr_inf_.Load(fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Parameters
				{
					ChangeDimensions(number_of_cabr_reactions_, &cabr_reaction_type_, true);
					ChangeDimensions(number_of_cabr_reactions_, &cabr_index_of_single_thirdbody_species_, true);

					ChangeDimensions(number_of_cabr_reactions_, &a_cabr_, true);
					ChangeDimensions(number_of_cabr_reactions_, &b_cabr_, true);
					ChangeDimensions(number_of_cabr_reactions_, &c_cabr_, true);
					ChangeDimensions(number_of_cabr_reactions_, &d_cabr_, true);
					ChangeDimensions(number_of_cabr_reactions_, &e_cabr_, true);

					cabr_indices_of_thirdbody_species_.resize(number_of_cabr_reactions_);
					cabr_indices_of_thirdbody_efficiencies_.resize(number_of_cabr_reactions_);

					rapidxml::xml_node<>* parameters_node = cabr_node->first_node("Parameters");

					unsigned int i = 1;
					for (rapidxml::xml_node<> *current_node = parameters_node->first_node("reaction"); current_node; current_node = current_node->next_sibling())
					{
						std::stringstream fInput;
						fInput << current_node->value();

						std::string dummy;
						fInput >> dummy;
						if (dummy == "lindemann")
						{
							cabr_reaction_type_[i] = PhysicalConstants::REACTION_LINDEMANN_CABR;
						}
						if (dummy == "troe")
						{
							cabr_reaction_type_[i] = PhysicalConstants::REACTION_TROE_CABR;

							OpenSMOKEVectorDouble coefficients;
							coefficients.Load(fInput, OPENSMOKE_FORMATTED_FILE);
							a_cabr_[i] = coefficients[1];
							b_cabr_[i] = coefficients[2];
							c_cabr_[i] = coefficients[3];
							if (coefficients.Size() == 4)
								d_cabr_[i] = coefficients[4];
						}
						else if (dummy == "sri")
						{
							cabr_reaction_type_[i] = PhysicalConstants::REACTION_SRI_CABR;
					
							OpenSMOKEVectorDouble coefficients;
							coefficients.Load(fInput, OPENSMOKE_FORMATTED_FILE);
							a_cabr_[i] = coefficients[1];
							b_cabr_[i] = coefficients[2];
							c_cabr_[i] = coefficients[3];
							if (coefficients.Size() == 5)
							{
								d_cabr_[i] = coefficients[4];
								e_cabr_[i] = coefficients[5];
							}
							else
							{
								d_cabr_[i] = 1.;
								e_cabr_[i] = 0.;
							}
						}

						fInput >> dummy;
						if (dummy == "species")
						{
							unsigned int species;
							fInput >> species;
							cabr_index_of_single_thirdbody_species_[i] = species;
						}
						else
						{
							cabr_indices_of_thirdbody_species_[i-1].Load(fInput, OPENSMOKE_FORMATTED_FILE);
							cabr_indices_of_thirdbody_efficiencies_[i-1].Load(fInput, OPENSMOKE_FORMATTED_FILE);
							cabr_indices_of_thirdbody_efficiencies_[i-1] -= 1.;
						}

						i++;
					}
				}
			}
                        
                        if(verbose_output_ == true)
			    std::cout << " * Reading kinetic parameters of additional reactions..." << std::endl;
			{
				// Chebyshev reaction parameters
				if (number_of_chebyshev_reactions_ != 0)
				{
					rapidxml::xml_node<>* current_node = kinetic_parameters_node->first_node("Chebyshev");
					std::stringstream fInput;
					fInput << current_node->value();

					chebyshev_reactions_ = new ChebyshevPolynomialRateExpression[number_of_chebyshev_reactions_];
					for(unsigned int j=0;j<number_of_chebyshev_reactions_;j++)
						chebyshev_reactions_[j].ReadFromASCIIFile(fInput);
				}

				// PressureLog reaction parameters
				if (number_of_pressurelog_reactions_ != 0)
				{
					rapidxml::xml_node<>* current_node = kinetic_parameters_node->first_node("PressureLog");
					std::stringstream fInput;
					fInput << current_node->value();

					pressurelog_reactions_ = new PressureLogarithmicRateExpression[number_of_pressurelog_reactions_];
					for(unsigned int j=0;j<number_of_pressurelog_reactions_;j++)
						pressurelog_reactions_[j].ReadFromASCIIFile(fInput);
				}

				// Extended PressureLog reaction parameters
				if (number_of_extendedpressurelog_reactions_ != 0)
				{
					rapidxml::xml_node<>* current_node = kinetic_parameters_node->first_node("ExtendedPressureLog");
					
					std::stringstream fInput;
					fInput << current_node->value();

					extendedpressurelog_reactions_ = new ExtendedPressureLogarithmicRateExpression[number_of_extendedpressurelog_reactions_];
					for (unsigned int j = 0; j < number_of_extendedpressurelog_reactions_; j++)
						extendedpressurelog_reactions_[j].ReadFromASCIIFile(fInput);
				}

				// Extended PressureLog reaction parameters
				if (number_of_extendedfalloff_reactions_ != 0)
				{
					rapidxml::xml_node<>* current_node = kinetic_parameters_node->first_node("ExtendedFallOff");

					std::stringstream fInput;
					fInput << current_node->value();

					extendedfalloff_reactions_ = new ExtendedFallOff[number_of_extendedfalloff_reactions_];
					for (unsigned int j = 0; j < number_of_extendedfalloff_reactions_; j++)
						extendedfalloff_reactions_[j].ReadFromASCIIFile(fInput);
				}

				// FIT1 reaction parameters
				if (number_of_fit1_reactions_ != 0)
				{
					rapidxml::xml_node<>* current_node = kinetic_parameters_node->first_node("FIT1");
					std::stringstream fInput;
					fInput << current_node->value();

					//	fit1_reactions_ = new Fit1RateExpression[number_of_fit1_reactions_];
					//	for(unsigned int j=0;j<number_of_fit1_reactions_;j++)
					//		fit1_reactions_[j].ReadFromASCIIFile(fInput);
				}

				// JanevLanger reaction parameters
				if (number_of_janevlanger_reactions_ != 0)
				{
					rapidxml::xml_node<>* current_node = kinetic_parameters_node->first_node("JanevLanger");
					std::stringstream fInput;
					fInput << current_node->value();

					//	janevlanger_reactions_ = new JanevLangerRateExpression[number_of_janevlanger_reactions_];
					//	for(unsigned int j=0;j<number_of_janevlanger_reactions_;j++)
					//		janevlanger_reactions_[j].ReadFromASCIIFile(fInput);
				}

				// LandauTeller reaction parameters
				if (number_of_landauteller_reactions_ != 0)
				{
					rapidxml::xml_node<>* current_node = kinetic_parameters_node->first_node("LandauTeller");
					std::stringstream fInput;
					fInput << current_node->value();

					//	landauteller_reactions_ = new LandauTellerRateExpression[number_of_landauteller_reactions_];
					//	for(unsigned int j=0;j<number_of_landauteller_reactions_;j++)
					//		landauteller_reactions_[j].ReadFromASCIIFile(fInput);
				}
			}
		} // kinetic parameters

		// Stoichiometry
		{
			stoichiometry_ = new StoichiometricMap(this->number_of_species_, this->number_of_reactions_, verbose_output_);
	
			rapidxml::xml_node<>* stoichiometry_node = kinetics_node->first_node("Stoichiometry");

			std::string stoichiometry_type = stoichiometry_node->first_attribute("type")->value();
			std::string stoichiometry_version = stoichiometry_node->first_attribute("version")->value();

			if (stoichiometry_type != "OpenSMOKE" || stoichiometry_version != "04-22-2013")
				ErrorMessage("void KineticsMap_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current stoichiometric data are not supported.");

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
			ChangeDimensions(this->number_of_reactions_, &reaction_s_over_R_, true);
			ChangeDimensions(this->number_of_reactions_, &reaction_h_over_RT_, true);
			ChangeDimensions(this->number_of_reactions_, &kArrhenius_, true);
			ChangeDimensions(this->number_of_reactions_, &kArrheniusModified_, true);
			ChangeDimensions(number_of_thermodynamic_reversible_reactions_, &uKeq_, false);
			ChangeDimensions(number_of_explicitly_reversible_reactions_, &kArrhenius_reversible_, false);
		
			ChangeDimensions(number_of_falloff_reactions_, &kArrhenius_falloff_inf_, false);
			ChangeDimensions(number_of_falloff_reactions_, &logFcent_falloff_, false);
			ChangeDimensions(number_of_cabr_reactions_, &kArrhenius_cabr_inf_, false);
			ChangeDimensions(number_of_cabr_reactions_, &logFcent_cabr_, false);

			ChangeDimensions(number_of_thirdbody_reactions_, &Meff_, false);
			ChangeDimensions(number_of_falloff_reactions_, &correction_falloff_, false);
			ChangeDimensions(number_of_cabr_reactions_, &correction_cabr_, false);

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

				for(unsigned int k=1;k<=number_of_thirdbody_reactions_;k++)
				{
					type_of_reaction_[indices_of_thirdbody_reactions_[k]] = PhysicalConstants::REACTION_THIRDBODY;
					local_family_index_[indices_of_thirdbody_reactions_[k]] = k;
				}
				for(unsigned int k=1;k<=number_of_chebyshev_reactions_;k++)
				{
					type_of_reaction_[indices_of_chebyshev_reactions_[k]] = PhysicalConstants::REACTION_CHEBYSHEV;
					local_family_index_[indices_of_chebyshev_reactions_[k]] = k;
				}

				for(unsigned int k=1;k<=number_of_falloff_lindemann_reactions_;k++)
					type_of_reaction_[indices_of_falloff_lindemann_reactions_[k]] = PhysicalConstants::REACTION_LINDEMANN_FALLOFF;	
				for(unsigned int k=1;k<=number_of_falloff_troe_reactions_;k++)
					type_of_reaction_[indices_of_falloff_troe_reactions_[k]] = PhysicalConstants::REACTION_TROE_FALLOFF;
				for(unsigned int k=1;k<=number_of_falloff_sri_reactions_;k++)
					type_of_reaction_[indices_of_falloff_sri_reactions_[k]] = PhysicalConstants::REACTION_SRI_FALLOFF;
				for(unsigned int k=1;k<=number_of_falloff_reactions_;k++)
					local_family_index_[indices_of_falloff_reactions_[k]] = k;

				for(unsigned int k=1;k<=number_of_cabr_lindemann_reactions_;k++)
					type_of_reaction_[indices_of_cabr_lindemann_reactions_[k]] = PhysicalConstants::REACTION_LINDEMANN_CABR;
				for(unsigned int k=1;k<=number_of_cabr_troe_reactions_;k++)
					type_of_reaction_[indices_of_cabr_troe_reactions_[k]] = PhysicalConstants::REACTION_TROE_CABR;
				for(unsigned int k=1;k<=number_of_cabr_sri_reactions_;k++)
					type_of_reaction_[indices_of_cabr_sri_reactions_[k]] = PhysicalConstants::REACTION_SRI_CABR;
				for(unsigned int k=1;k<=number_of_cabr_reactions_;k++)
					local_family_index_[indices_of_cabr_reactions_[k]] = k;
			}

			// Sparsity Analysis
			const bool sparsity_features_enabled_ = true;
			if (sparsity_features_enabled_ == true)
			{
				jacobian_sparsity_pattern_map_ = new JacobianSparsityPatternMap<KineticsMap_CHEMKIN>(*this);
                isJacobianSparsityMapAvailable_ = true;
			}
                       
			// Output
			if (verbose_output_ == true)
			{
				std::cout << std::endl;
				std::cout << "----------------------------------------------------------------------------" << std::endl;
				std::cout << " Kinetic Mechanism Summary" << std::endl;
				std::cout << "----------------------------------------------------------------------------" << std::endl;
				std::cout << " Total number of species:        " << this->number_of_species_ << std::endl;
				std::cout << " Total number of reactions:      " << this->number_of_reactions_ << std::endl;
				std::cout << "   Reversible reactions:         " << number_of_reversible_reactions_ << " (" << number_of_reversible_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << "    * by thermodynamics:         " << number_of_thermodynamic_reversible_reactions_ << " (" << number_of_thermodynamic_reversible_reactions_ / std::max(1., double(number_of_reversible_reactions_))*100. << "%)" << std::endl;
				std::cout << "    * by Arrhenius' law:         " << number_of_explicitly_reversible_reactions_ << " (" << number_of_explicitly_reversible_reactions_ / std::max(1., double(number_of_reversible_reactions_))*100. << "%)" << std::endl;
				std::cout << "   Pressure dependent reactions: " << number_of_falloff_reactions_ + number_of_cabr_reactions_ << " (" << (number_of_falloff_reactions_ + number_of_cabr_reactions_) / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << "    * fall-off reactions:        " << number_of_falloff_reactions_ << " (" << number_of_falloff_reactions_ / std::max(1., double(number_of_falloff_reactions_ + number_of_cabr_reactions_))*100. << "%)" << std::endl;
				std::cout << "      ** lindemann form:         " << number_of_falloff_lindemann_reactions_ << " (" << number_of_falloff_lindemann_reactions_ / std::max(1., double(number_of_falloff_reactions_))*100. << "%)" << std::endl;
				std::cout << "      ** troe form:              " << number_of_falloff_troe_reactions_ << " (" << number_of_falloff_troe_reactions_ / std::max(1., double(number_of_falloff_reactions_))*100. << "%)" << std::endl;
				std::cout << "      ** sri form:               " << number_of_falloff_sri_reactions_ << " (" << number_of_falloff_sri_reactions_ / std::max(1., double(number_of_falloff_reactions_))*100. << "%)" << std::endl;
				std::cout << "    * cabr reactions:            " << number_of_cabr_reactions_ << " (" << number_of_cabr_reactions_ / std::max(1., double(number_of_falloff_reactions_ + number_of_cabr_reactions_))*100. << "%)" << std::endl;
				std::cout << "   Chebyshev reactions:          " << number_of_chebyshev_reactions_ << " (" << number_of_chebyshev_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << "   Pressure-Log reactions:       " << number_of_pressurelog_reactions_ << " (" << number_of_pressurelog_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << "   Ext-Pressure-Log reactions:   " << number_of_extendedpressurelog_reactions_ << " (" << number_of_extendedpressurelog_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << "   Ext-Falloff reactions:        " << number_of_extendedfalloff_reactions_ << " (" << number_of_extendedfalloff_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << "   Fit1 reactions:               " << number_of_fit1_reactions_ << " (" << number_of_fit1_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << "   Janev-Langer reactions:       " << number_of_janevlanger_reactions_ << " (" << number_of_janevlanger_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << "   Landau-Teller reactions:      " << number_of_landauteller_reactions_ << " (" << number_of_landauteller_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << " Negative frequency factors:     " << negative_lnA_.Size() << " (" << negative_lnA_.Size() / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
				std::cout << std::endl;

				stoichiometry_->Summary(std::cout);
			}
		}
	}

	void KineticsMap_CHEMKIN::ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
		try
		{
			this->number_of_species_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));					
		}
		catch(...)
		{
			ErrorMessage("KineticsMap_CHEMKIN::ImportSpeciesFromXMLFile", "Error in reading the number of species.");
		}
	}
 
	void KineticsMap_CHEMKIN::ReactionEnthalpiesAndEntropies()
	{
		if (reaction_h_and_s_must_be_recalculated_ == true)
		{
			stoichiometry_->ReactionEnthalpyAndEntropy(	reaction_h_over_RT_, reaction_s_over_R_, 
									thermodynamics_.species_h_over_RT(), thermodynamics_.species_s_over_R() );

			reaction_h_and_s_must_be_recalculated_ = false;
		}
	}

	void KineticsMap_CHEMKIN::KineticConstants()
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

				// Negative frequency factors: reactions must be reversed
				for (int j = 1; j <= negative_lnA_.Size(); j++)
					kArrhenius_[negative_lnA_[j]] *= -1.;
			}

			// Equilibrium constants (inverse value)
			{
				for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
				{
					unsigned int j = indices_of_thermodynamic_reversible_reactions_[k];
					uKeq_[k] = -reaction_s_over_R_[j] + reaction_h_over_RT_[j] - log_Patm_over_RT_ * changeOfMoles_[j];
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

			// Fall-off high temperature region kinetic constants
			if (number_of_falloff_reactions_ != 0)
			{
				double *pt_lnA = lnA_falloff_inf_.GetHandle();
				double *pt_Beta = Beta_falloff_inf_.GetHandle();
				double *pt_E_over_R = E_over_R_falloff_inf_.GetHandle();
				double *pt_kArrhenius = kArrhenius_falloff_inf_.GetHandle();

				for(unsigned int k=0;k<number_of_falloff_reactions_;k++)
					*pt_kArrhenius++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius_falloff_inf_, &kArrhenius_falloff_inf_);

				for(unsigned int k=1;k<=number_of_falloff_reactions_;k++)
				{
					switch(falloff_reaction_type_[k])
					{
						case PhysicalConstants::REACTION_TROE_FALLOFF:
							
							logFcent_falloff_[k] = (d_falloff_[k] != 0.) ? (1.-a_falloff_[k])*std::exp(-this->T_/b_falloff_[k]) + a_falloff_[k]*std::exp(-this->T_/c_falloff_[k]) + std::exp(-d_falloff_[k]/this->T_) :
																		 (1.-a_falloff_[k])*std::exp(-this->T_/b_falloff_[k]) + a_falloff_[k]*std::exp(-this->T_/c_falloff_[k]);
							
							if(logFcent_falloff_[k] < 1.e-300)	logFcent_falloff_[k] = -300.;
							else								logFcent_falloff_[k] = std::log10(logFcent_falloff_[k]);

							break;
						
						case PhysicalConstants::REACTION_SRI_FALLOFF:

							logFcent_falloff_[k] = ( a_falloff_[k]*std::exp(-b_falloff_[k]/this->T_) +  std::exp(-this->T_/c_falloff_[k]));
							
							break;
					}
				}
			}

			// Cabr high temperature region kinetic constants
			if (number_of_cabr_reactions_ != 0)
			{
				double *pt_lnA = lnA_cabr_inf_.GetHandle();
				double *pt_Beta = Beta_cabr_inf_.GetHandle();
				double *pt_E_over_R = E_over_R_cabr_inf_.GetHandle();
				double *pt_kArrhenius = kArrhenius_cabr_inf_.GetHandle();

				for(unsigned int k=0;k<number_of_cabr_reactions_;k++)
					*pt_kArrhenius++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius_cabr_inf_, &kArrhenius_cabr_inf_);

				for(unsigned int k=1;k<=number_of_cabr_reactions_;k++)
				{
					switch(cabr_reaction_type_[k])
					{
						case PhysicalConstants::REACTION_TROE_CABR:
							
							logFcent_cabr_[k] = (d_cabr_[k] != 0.) ? (1.-a_cabr_[k])*std::exp(-this->T_/b_cabr_[k]) + a_cabr_[k]*std::exp(-this->T_/c_cabr_[k]) + std::exp(-d_cabr_[k]/this->T_) :
																		 (1.-a_cabr_[k])*std::exp(-this->T_/b_cabr_[k]) + a_cabr_[k]*std::exp(-this->T_/c_cabr_[k]);
							
							if(logFcent_cabr_[k] < 1.e-300)	logFcent_cabr_[k] = -300.;
							else							logFcent_cabr_[k] = std::log10(logFcent_cabr_[k]);
							
							break;
						
						case PhysicalConstants::REACTION_SRI_CABR:

							logFcent_cabr_[k] = ( a_cabr_[k]*std::exp(-b_cabr_[k]/this->T_) +  std::exp(-this->T_/c_cabr_[k]));
							
							break;
					}
				}
			}

			arrhenius_kinetic_constants_must_be_recalculated_ = false;
		}

		//if (nonconventional_kinetic_constants_must_be_recalculated_ == true)
		{
			// Chebishev-Polynomials reactions
			{
				for(unsigned int k=1;k<=number_of_chebyshev_reactions_;k++)
				{
					const unsigned int j = indices_of_chebyshev_reactions_[k];
					kArrhenius_[j] = chebyshev_reactions_[k-1].KineticConstant(this->T_, this->P_);
				}
			}

			// Pressure logarithmic interpolated reactions
			{
				for(unsigned int k=1;k<=number_of_pressurelog_reactions_;k++)
				{
					const unsigned int j = indices_of_pressurelog_reactions_[k];
					kArrhenius_[j] = pressurelog_reactions_[k-1].KineticConstant(this->T_,this->P_);
				}
			}

			// Extended pressure logarithmic interpolated reactions
			// Extended falloff reactions
			// They are calculated separately

			// Fit1 reactions
			{
				for(unsigned int k=1;k<=number_of_fit1_reactions_;k++)
				{
					const unsigned int j = indices_of_fit1_reactions_[k];
				//	kArrhenius_[j] = fit1_reactions_[k-1].KineticConstant(this->T_,this->P_);
				}
			}

			// Janev-Langer reactions
			{
				for(unsigned int k=1;k<=number_of_janevlanger_reactions_;k++)
				{
					const unsigned int j = indices_of_janevlanger_reactions_[k];
				//	kArrhenius_[j] = janevlanger_reactions_[k-1].KineticConstant(this->T_,this->P_);
				}
			}

			// Landau-Teller reactions
			{
				for(unsigned int k=1;k<=number_of_landauteller_reactions_;k++)
				{
					const unsigned int j = indices_of_landauteller_reactions_[k];
				//	kArrhenius_[j] = landauteller_reactions_[k-1].KineticConstant(this->T_,this->P_);
				}
			}

			nonconventional_kinetic_constants_must_be_recalculated_ = false;
		}

		kArrheniusModified_ = kArrhenius_;
	}

	void KineticsMap_CHEMKIN::ThirdBodyReactions(const double cTot, const OpenSMOKEVectorDouble& c)
	{	
		for(unsigned int s=1;s<=number_of_thirdbody_reactions_;s++)
		{
			const unsigned int j=indices_of_thirdbody_reactions_[s];

			Meff_[s] = cTot;
			for(int k=1;k<=indices_of_thirdbody_species_[s-1].Size();k++)
				Meff_[s] += c[ indices_of_thirdbody_species_[s-1][k] ] * indices_of_thirdbody_efficiencies_[s-1][k];
		}
	}

	void KineticsMap_CHEMKIN::WeakThirdBodyConcentrationEfficiencies(std::vector<unsigned int>& reactions, std::vector<unsigned int>& species)
	{
		reactions.resize(0);
		species.resize(0);

		// Third-body third-body
		for (unsigned int s = 1; s <= number_of_thirdbody_reactions_; s++)
		{
			const unsigned int j = indices_of_thirdbody_reactions_[s];

			for (int k = 1; k <= indices_of_thirdbody_species_[s - 1].Size(); k++)
			{
				reactions.push_back(j);
				species.push_back(indices_of_thirdbody_species_[s - 1][k]);
			}
		}
	}

	void KineticsMap_CHEMKIN::WeakFallOffConcentrationEfficiencies(std::vector<unsigned int>& reactions, std::vector<unsigned int>& species)
	{
		reactions.resize(0);
		species.resize(0);

		// Fall-off reactions
		for (unsigned int k = 1; k <= number_of_falloff_reactions_; k++)
		{
			const unsigned int j = indices_of_falloff_reactions_[k];

			if (falloff_index_of_single_thirdbody_species_[k] == 0)
			{
				for (int s = 1; s <= falloff_indices_of_thirdbody_species_[k - 1].Size(); s++)
				{
					reactions.push_back(j);
					species.push_back(falloff_indices_of_thirdbody_species_[k - 1][s]);
				}
			}
			else
			{
				reactions.push_back(j);
				species.push_back(falloff_index_of_single_thirdbody_species_[k]);
			}
		}
	}

	void KineticsMap_CHEMKIN::WeakCABRConcentrationEfficiencies(std::vector<unsigned int>& reactions, std::vector<unsigned int>& species)
	{
		reactions.resize(0);
		species.resize(0);

		// CABR
		for (unsigned int k = 1; k <= number_of_cabr_reactions_; k++)
		{
			const unsigned int j = indices_of_cabr_reactions_[k];

			if (cabr_index_of_single_thirdbody_species_[k] == 0)
			{
				for (int s = 1; s <= cabr_indices_of_thirdbody_species_[k - 1].Size(); s++)
				{
					reactions.push_back(j);
					species.push_back(cabr_indices_of_thirdbody_species_[k - 1][s]);
				}
			}
			else
			{
				reactions.push_back(j);
				species.push_back(cabr_index_of_single_thirdbody_species_[k]);
			}
		}
	}

	void KineticsMap_CHEMKIN::StrongConcentrationEfficiencies(std::vector<unsigned int>& reactions)
	{
		reactions.resize(0);

		// Third-body third-body
		for (unsigned int s = 1; s <= number_of_thirdbody_reactions_; s++)
		{
			const unsigned int j = indices_of_thirdbody_reactions_[s];
			reactions.push_back(j);
		}

		// Fall-off reactions
		for (unsigned int k = 1; k <= number_of_falloff_reactions_; k++)
		{
			const unsigned int j = indices_of_falloff_reactions_[k];
			reactions.push_back(j);
		}

		// CABR
		for (unsigned int k = 1; k <= number_of_cabr_reactions_; k++)
		{
			const unsigned int j = indices_of_cabr_reactions_[k];
			reactions.push_back(j);
		}
	}
	
	void KineticsMap_CHEMKIN::FallOffReactions(const double cTot, const OpenSMOKEVectorDouble& c)
	{
		//if (arrhenius_kinetic_constants_must_be_corrected_ == true)
		{
			for(unsigned int k=1;k<=number_of_falloff_reactions_;k++)
			{
				double M;
				if (falloff_index_of_single_thirdbody_species_[k] == 0)
				{
					M = cTot;
					for(int s=1;s<=falloff_indices_of_thirdbody_species_[k-1].Size();s++)
						M += c[ falloff_indices_of_thirdbody_species_[k-1][s] ] * falloff_indices_of_thirdbody_efficiencies_[k-1][s];
				}
				else
				{
					const double epsilon = 1.e-16;
					M = c[falloff_index_of_single_thirdbody_species_[k]] + epsilon;
				}
			
				const unsigned int j=indices_of_falloff_reactions_[k];
				const double Pr = kArrhenius_[j] * M / kArrhenius_falloff_inf_[k];
		
				double wF = 1.;
				switch(falloff_reaction_type_[k])
				{
					case PhysicalConstants::REACTION_TROE_FALLOFF:

                        if (Pr > 1.e-32)
                        {
                            const double nTroe = 0.75-1.27*logFcent_falloff_[k];
                            const double cTroe = -0.4-0.67*logFcent_falloff_[k];
                            const double sTroe = std::log10(Pr) + cTroe;
                            wF = std::pow(10., logFcent_falloff_[k]/(1. + boost::math::pow<2>(sTroe/(nTroe-0.14*sTroe)))); 
                        }
                        else
                        {
                            // Asymptotic value for wF when sTroe --> -Inf
                            wF = std::pow(10., logFcent_falloff_[k]/(1. + boost::math::pow<2>(1./0.14)));
                        }
                                                
						break;

					case PhysicalConstants::REACTION_SRI_FALLOFF:

						const double xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
						wF = std::pow(logFcent_falloff_[k], xSRI) * d_falloff_[k];
						if(e_falloff_[k] != 0.)
							wF *= std::pow(this->T_, e_falloff_[k]);
						break;
				}
	
				correction_falloff_[k] = kArrhenius_falloff_inf_[k] * (Pr/(1.+Pr)) * wF / kArrhenius_[j];
			}
		}
	}

	double KineticsMap_CHEMKIN::FallOffReactionsCorrection(const unsigned int local_k, const double cTot, const OpenSMOKEVectorDouble& c)
	{
		double M;
		if (falloff_index_of_single_thirdbody_species_[local_k] == 0)
		{
			M = cTot;
			for (int s = 1; s <= falloff_indices_of_thirdbody_species_[local_k - 1].Size(); s++)
				M += c[falloff_indices_of_thirdbody_species_[local_k - 1][s]] * falloff_indices_of_thirdbody_efficiencies_[local_k - 1][s];
		}
		else
		{
			const double epsilon = 1.e-16;
			M = c[falloff_index_of_single_thirdbody_species_[local_k]] + epsilon;
		}

		const unsigned int j = indices_of_falloff_reactions_[local_k];
		const double Pr = kArrhenius_[j] * M / kArrhenius_falloff_inf_[local_k];

		double wF = 1.;
		switch (falloff_reaction_type_[local_k])
		{
			case PhysicalConstants::REACTION_TROE_FALLOFF:

				if (Pr > 1.e-32)
				{
					const double nTroe = 0.75 - 1.27*logFcent_falloff_[local_k];
					const double cTroe = -0.4 - 0.67*logFcent_falloff_[local_k];
					const double sTroe = std::log10(Pr) + cTroe;
					wF = std::pow(10., logFcent_falloff_[local_k] / (1. + boost::math::pow<2>(sTroe / (nTroe - 0.14*sTroe))));
				}
				else
				{
					// Asymptotic value for wF when sTroe --> -Inf
					wF = std::pow(10., logFcent_falloff_[local_k] / (1. + boost::math::pow<2>(1. / 0.14)));
				}

				break;

			case PhysicalConstants::REACTION_SRI_FALLOFF:

				if (Pr > 1.e-32)
				{
					const double xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
					wF = std::pow(logFcent_falloff_[local_k], xSRI) * d_falloff_[local_k];
					if (e_falloff_[local_k] != 0.)
						wF *= std::pow(this->T_, e_falloff_[local_k]);
				}
				else
				{
					const double xSRI = 0.;
					wF = std::pow(logFcent_falloff_[local_k], xSRI) * d_falloff_[local_k];
					if (e_falloff_[local_k] != 0.)
						wF *= std::pow(this->T_, e_falloff_[local_k]);
				}

				break;
		}

		return ( kArrhenius_falloff_inf_[local_k] * (Pr / (1. + Pr)) * wF / kArrhenius_[j] );
	}
 
	void KineticsMap_CHEMKIN::ChemicallyActivatedBimolecularReactions(const double cTot, const OpenSMOKEVectorDouble& c)
	{
		//if (arrhenius_kinetic_constants_must_be_corrected_ == true)
		{
			for(unsigned int k=1;k<=number_of_cabr_reactions_;k++)
			{
				double M;
				if (cabr_index_of_single_thirdbody_species_[k] == 0)
				{
					M = cTot;
					for(int s=1;s<=cabr_indices_of_thirdbody_species_[k-1].Size();s++)
						M += c[ cabr_indices_of_thirdbody_species_[k-1][s] ] * cabr_indices_of_thirdbody_efficiencies_[k-1][s];
				}
				else
				{
					const double epsilon = 1.e-16;
					M = c[cabr_index_of_single_thirdbody_species_[k]] + epsilon;
				}
			
				const unsigned int j=indices_of_cabr_reactions_[k];
				const double Pr = kArrhenius_[j] * M / kArrhenius_cabr_inf_[k];
		
				double wF = 1.;
				double nTroe, cTroe, sTroe, xSRI;
				switch(cabr_reaction_type_[k])
				{
					case PhysicalConstants::REACTION_TROE_CABR:

						nTroe = 0.75-1.27*logFcent_cabr_[k];
						cTroe = -0.4-0.67*logFcent_cabr_[k];
						sTroe = std::log10(Pr) + cTroe;
						wF = std::pow(10., logFcent_cabr_[k]/(1. + boost::math::pow<2>(sTroe/(nTroe-0.14*sTroe))));
						
						break;
					
					case PhysicalConstants::REACTION_SRI_CABR:
						xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
						wF = std::pow(logFcent_cabr_[k], xSRI) * d_cabr_[k];
						if(e_cabr_[k] != 0.)
							wF *= std::pow(this->T_, e_cabr_[k]);
						break;
				}
	
				correction_cabr_[k] = (1./(1.+Pr)) * wF;
			}
		}
	}

	// Extended Pressure logarithmic interpolated reactions
	void KineticsMap_CHEMKIN::ExtendedPressureLogReactions(const double cTot, const OpenSMOKEVectorDouble& c)
	{
		for (unsigned int k = 1; k <= number_of_extendedpressurelog_reactions_; k++)
		{
			const unsigned int j = indices_of_extendedpressurelog_reactions_[k];
			kArrhenius_[j] = extendedpressurelog_reactions_[k-1].KineticConstant(this->T_, this->P_, cTot, c.GetHandle());
			kArrheniusModified_[j] = kArrhenius_[j];
		}
	}

	// Extended Falloff reactions
	void KineticsMap_CHEMKIN::ExtendedFallOffReactions(const double cTot, const OpenSMOKEVectorDouble& c)
	{
		for (unsigned int k = 1; k <= number_of_extendedfalloff_reactions_; k++)
		{
			const unsigned int j = indices_of_extendedfalloff_reactions_[k];
			kArrhenius_[j] = extendedfalloff_reactions_[k - 1].KineticConstant(this->T_, this->P_, cTot, c.GetHandle());
			kArrheniusModified_[j] = kArrhenius_[j];
		}
	}
	
	void KineticsMap_CHEMKIN::ReactionRates(const OpenSMOKEVectorDouble& c)
	{
		const double cTot = c.SumElements();
		ReactionRates(c, cTot);
	}

	void KineticsMap_CHEMKIN::ReactionRates(const OpenSMOKEVectorDouble& c, const double cTot)
	{
		// 1. Kinetic constants
		KineticConstants();

		// 1.bis Extended pressure log reactions
		ExtendedPressureLogReactions(cTot, c);

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

		// 4.bis Extended fallExtendedFallOffReactionsoff reactions
		ExtendedFallOffReactions(cTot, c);

		// 5. Correct the effective kinetic constants: CABR reactions
		ChemicallyActivatedBimolecularReactions(cTot, c);
		for(unsigned int s=1;s<=number_of_cabr_reactions_;s++)
		{
			const unsigned int j=indices_of_cabr_reactions_[s];
			kArrheniusModified_[j] *= correction_cabr_[s];
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

	void KineticsMap_CHEMKIN::DerivativesOfReactionRatesWithRespectToKineticParameters(const PhysicalConstants::sensitivity_type type, unsigned int jReaction, const OpenSMOKEVectorDouble& c, double& parameter)
	{
		const double cTot = c.SumElements();

		// 1. Kinetic constants
		KineticConstants();

		// 1.bis Extended pressure log reactions
		ExtendedPressureLogReactions(cTot, c);
			
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

		// 4.bis Extended fallExtendedFallOffReactionsoff reactions
		ExtendedFallOffReactions(cTot, c);

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
				parameter =  sign_lnA_[jReaction]*std::exp(lnA_[jReaction]);

				if (type_of_reaction_[jReaction] == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || 
					type_of_reaction_[jReaction] == PhysicalConstants::REACTION_TROE_FALLOFF ||
					type_of_reaction_[jReaction] == PhysicalConstants::REACTION_SRI_FALLOFF)
				{
					unsigned int jFallOff = local_family_index_[jReaction];
					double	kInf = std::exp(lnA_falloff_inf_[jFallOff] + Beta_falloff_inf_[jFallOff]*this->logT_ - E_over_R_falloff_inf_[jFallOff]*this->uT_);
					double F, dF_over_dA0, dF_over_dAInf;
					this->FallOffReactions(jFallOff, cTot, c, F, dF_over_dA0, dF_over_dAInf);

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
					kArrheniusModified_[jReaction] = kArrheniusModified_[jReaction] / (sign_lnA_[jReaction]*std::exp(lnA_[jReaction]));
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
					
					double F = 0.;
					double dF_over_dA0 = 0.; 
					double dF_over_dAInf = 0.;
					FallOffReactions(jFallOff, cTot, c, F, dF_over_dA0, dF_over_dAInf);

					kArrheniusModified_[index_global] = kArrheniusModified_[index_global] / F * 
														( kArrheniusModified_[index_global]/std::exp(lnA_falloff_inf_[jFallOff])/kInf + dF_over_dAInf);	

					if (parameter == 0.) std::cout << "parameter " << jReaction - this->number_of_reactions_ << std::endl;
					if (F == 0.) std::cout << "Falloff inv F " << jReaction - this->number_of_reactions_ << std::endl;
					if (kInf == 0.) std::cout << "Falloff inv kInf " << jReaction - this->number_of_reactions_ << std::endl;
					if (std::exp(lnA_falloff_inf_[jFallOff]) == 0.) std::cout << "Falloff inv exp " << jReaction - this->number_of_reactions_ << std::endl;

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

					double F = 0.;
					double dF_over_dA0 = 0.;
					double dF_over_dAInf = 0.;
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

	void KineticsMap_CHEMKIN::FormationRates(OpenSMOKEVectorDouble* R)
	{
		stoichiometry_->FormationRatesFromReactionRates(R, netReactionRates_);
	}

	double KineticsMap_CHEMKIN::HeatRelease(const OpenSMOKEVectorDouble& R)
	{
		return -Dot(R, thermodynamics_.species_h_over_RT()) * PhysicalConstants::R_J_kmol * this->T_;
	}

	void KineticsMap_CHEMKIN::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const unsigned int k, const OpenSMOKEVectorDouble& c, OpenSMOKEVectorDouble* Jalfa, double& parameter)
	{
		DerivativesOfReactionRatesWithRespectToKineticParameters(type, k, c, parameter);
		stoichiometry_->FormationRatesFromReactionRates(Jalfa, netReactionRates_);
	}

	void KineticsMap_CHEMKIN::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& parameter)
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
 
	void KineticsMap_CHEMKIN::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& Jrho, double& parameter)
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


	// This version seems to be wrong
	void KineticsMap_CHEMKIN::ProductionAndDestructionRatesGross(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D)
	{
		OpenSMOKEVectorDouble forward_(this->number_of_reactions_);
		OpenSMOKEVectorDouble backward_(this->number_of_reactions_);

		GetForwardReactionRates(&forward_);
		GetBackwardReactionRates(&backward_);

		stoichiometry_->ProductionAndDestructionRatesFromReactionRatesGross(P, D, forward_, backward_);
	}
	
	void KineticsMap_CHEMKIN::ProductionAndDestructionRates(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D)
	{
		stoichiometry_->ProductionAndDestructionRatesFromReactionRates(P, D, netReactionRates_);
	}
 
	void KineticsMap_CHEMKIN::RateOfProductionAnalysis(const bool iNormalize) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, iNormalize);
	}

	void KineticsMap_CHEMKIN::RateOfProductionAnalysis(std::ostream& fout) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, false);
		stoichiometry_->WriteRateOfProductionAnalysis(fout);
	}

	void KineticsMap_CHEMKIN::RateOfProductionAnalysis(ROPA_Data& ropa) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates_, false);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}

	void KineticsMap_CHEMKIN::RateOfProductionAnalysis(ROPA_Data& ropa, const OpenSMOKE::OpenSMOKEVectorDouble& rf, const OpenSMOKE::OpenSMOKEVectorDouble& rb) const
	{
		stoichiometry_->RateOfProductionAnalysis(rf, rb);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}

	const OpenSMOKEVectorDouble& KineticsMap_CHEMKIN::GetReactionRates()
	{
		return netReactionRates_;
	}

	void KineticsMap_CHEMKIN::GetReactionRates(OpenSMOKEVectorDouble* r)
	{
		*r = netReactionRates_;
	}

	void KineticsMap_CHEMKIN::GetForwardReactionRates(OpenSMOKEVectorDouble* r)
	{
		ElementByElementProduct(forwardReactionRates_, kArrheniusModified_, r);
	}

	void KineticsMap_CHEMKIN::GetBackwardReactionRates(OpenSMOKEVectorDouble* r)
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

	void KineticsMap_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k)
	{			
		thermodynamics_.SetPressure(101325.);
		thermodynamics_.SetTemperature(298.15);
		SetTemperature(298.15);
		SetPressure(101325.);
		ReactionEnthalpiesAndEntropies();
		
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << (reaction_h_over_RT_[k]-reaction_s_over_R_[k])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_h_over_RT_[k] *PhysicalConstants::R_kcal_mol*this->T_;;							// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_s_over_R_[k] *PhysicalConstants::R_cal_mol;;									// [cal/mol/K]
	}

	void KineticsMap_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k, OpenSMOKEVectorDouble& c_bath, const std::vector<double> list_of_temperatures, const double conversion_forward, const double conversion_backward)
	{			
		const double patm = 101325.;
		SetPressure(patm);
		thermodynamics_.SetPressure(patm);
	
		OpenSMOKEVectorDouble temperatures;
		
		if (list_of_temperatures.size() == 0)
		{
			OpenSMOKE::ChangeDimensions(5, &temperatures, true);
			temperatures[1] = 300.;
			temperatures[2] = 1000.;
			temperatures[3] = 1500.;
			temperatures[4] = 2000.;
			temperatures[5] = 2500.;
		}
		else
		{
			OpenSMOKE::ChangeDimensions(list_of_temperatures.size(), &temperatures, true);
			for (int i = 1; i <= temperatures.Size(); i++)
				temperatures[i] = list_of_temperatures[i - 1];
		}

		fOut << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOut << "    Temperature   kF            Keq           kR            DG            DH            DS            kF            kR" << std::endl;          
		fOut << "    [K]           [kmol,m3,s]   [-]           [kmol,m3,s]   [kcal/mol]    [kcal/mol]    [cal/mol/K]   [mol,cm3,s]   [mol,cm3,s]" << std::endl;        
		fOut << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;

		for(int i=1;i<=temperatures.Size();i++)
		{
			SetTemperature(temperatures[i]);
			thermodynamics_.SetTemperature(temperatures[i]);

			ReactionEnthalpiesAndEntropies();
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

	void KineticsMap_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k, const double T, const double P_Pa, OpenSMOKEVectorDouble& c)
	{
		SetPressure(P_Pa);
		thermodynamics_.SetPressure(P_Pa);
		SetTemperature(T);
		thermodynamics_.SetTemperature(T);

		ReactionEnthalpiesAndEntropies();
		KineticConstants();
		ReactionRates(c);

		// Label
		if (k == 1)
		{
			fOut << std::setw(16) << std::left << "Index";
			fOut << std::setw(16) << std::left << "kF[kmol,m3,s]";
			fOut << std::setw(16) << std::left << "Keq[-]";
			fOut << std::setw(16) << std::left << "kR[kmol,m3,s]";
			fOut << std::setw(16) << std::left << "DG[kcal/mol]";
			fOut << std::setw(16) << std::left << "DH[kcal/mol]";
			fOut << std::setw(16) << std::left << "DS[kcal/mol/K]";
			fOut << std::setw(16) << std::left << "kInf[kmol,m3,s]";
			fOut << std::setw(16) << std::left << "F[-]";
			fOut << std::endl;
		}

		// Index
		fOut << std::setw(16) << std::left << k;

		// Forward kinetic constant [kmol, m3, s]
		fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kArrheniusModified_[k];

		// Equilibrium and Backward kinetic constant [kmol, m3, s]
		if (isThermodynamicallyReversible_[k] != 0)
		{
			const unsigned int j = isThermodynamicallyReversible_[k];
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << 1. / uKeq_[j];
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kArrheniusModified_[k] * uKeq_[j];
		}
		else if (isExplicitlyReversible_[k] != 0)
		{
			const unsigned int j = isExplicitlyReversible_[k];
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kArrheniusModified_[k] / kArrhenius_reversible_[j];
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kArrhenius_reversible_[j];
		}
		else
		{
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::fixed << 0.;
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::fixed << 0.;
		}

		// Thermodynamic data
		fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << (reaction_h_over_RT_[k] - reaction_s_over_R_[k])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
		fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << reaction_h_over_RT_[k] * PhysicalConstants::R_kcal_mol*this->T_;;						// [kcal/mol]
		fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << reaction_s_over_R_[k] * PhysicalConstants::R_kcal_mol;;							// [kcal/mol/K]


		if (type_of_reaction_[k] == PhysicalConstants::REACTION_LINDEMANN_FALLOFF ||
			type_of_reaction_[k] == PhysicalConstants::REACTION_TROE_FALLOFF ||
			type_of_reaction_[k] == PhysicalConstants::REACTION_SRI_FALLOFF)
		{
			unsigned int jFallOff = local_family_index_[k];
			double	kInf = std::exp(lnA_falloff_inf_[jFallOff] + Beta_falloff_inf_[jFallOff] * this->logT_ - E_over_R_falloff_inf_[jFallOff] * this->uT_);
			
			double F = 0.;
			double dF_over_dA0 = 0.;
			double dF_over_dAInf = 0.;
			this->FallOffReactions(jFallOff, c.SumElements(), c, F, dF_over_dA0, dF_over_dAInf);

			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kInf;
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << F;
		}
		else
		{
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << 0.;
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << 0.;
		}

		fOut << std::endl;
	}

	void KineticsMap_CHEMKIN::FittedReverseKineticConstants(const OpenSMOKEVectorDouble& x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters, const bool only_reversible)
	{			
		unsigned int npoints = 11;
		Eigen::VectorXd T(npoints);
		T(0) = 300.;	T(1) = 600.;	T(2) = 900.;	T(3) = 1100.;
		T(4) = 1300.;	T(5) = 1500.;	T(6) = 1700.;	T(7) = 1900.;
		T(8) = 2100.;   T(9) = 2300.;   T(10) = 2500.;

		Eigen::MatrixXd XT(nparameters, npoints);
		Eigen::PartialPivLU<Eigen::MatrixXd> LU;
		
        if(verbose_output_ == true)
		    std::cout << "   assembling X and XT matrices..." << std::endl;
		{
			// Assembling X and XT Matrices
			Eigen::MatrixXd X(npoints, nparameters);
			Eigen::MatrixXd XTX(nparameters, nparameters);
			
			for(unsigned int i=0;i<npoints;i++)
			{
				X(i, 0) = 1.;
				X(i, 1) = -1. / PhysicalConstants::R_J_kmol / T(i);

				// In case of 3 parameters we add also the exponent of temperature
				if (nparameters == 3)
					X(i, 2) = std::log(T(i));
			}

			XT = X.transpose();
			XTX=XT*X;

			// Factorize
			LU.compute(XTX);
		}

		if (verbose_output_ == true)
			std::cout << "   evaluating the reaction rates..." << std::endl;

		const double patm = 101325.;
		SetPressure(patm);
		thermodynamics_.SetPressure(patm);

		unsigned int nr = number_of_thermodynamic_reversible_reactions_;
		if (only_reversible == false)
			nr = this->number_of_reactions_;
		
		Eigen::MatrixXd Y(nparameters, nr);
		Y.setConstant(0.);

		// Building reaction rates
		{
			Eigen::MatrixXd y(npoints, nr);
			y.setConstant(0.);

			OpenSMOKEVectorDouble c(x_bath.Size());
			for (unsigned int i = 1; i <= npoints; i++)
			{
				const double ctot = patm / PhysicalConstants::R_J_kmol / T(i - 1);
				Product(ctot, x_bath, &c);

				SetTemperature(T(i - 1));
				thermodynamics_.SetTemperature(T(i - 1));
				ReactionEnthalpiesAndEntropies();
				KineticConstants();
				ReactionRates(c);

				for (unsigned int k = 1; k <= this->number_of_reactions_; k++)
				{
					if (only_reversible == false)
					{
						double one_over_Keq = -reaction_s_over_R_[k] + reaction_h_over_RT_[k] - log_Patm_over_RT_ * changeOfMoles_[k];	
						one_over_Keq = std::exp(one_over_Keq);

						y(i - 1, k - 1) = std::log(kArrheniusModified_[k]* one_over_Keq);
					}
					else
					{
						if (isThermodynamicallyReversible_[k] != 0)
						{
							const unsigned int j = isThermodynamicallyReversible_[k];
							y(i - 1, j - 1) = std::log(kArrheniusModified_[k] * uKeq_[j]);
						}
					}
				}
			}

			Y = XT*y;
		}

		// Solving linear system
		{
			if (verbose_output_ == true)
				std::cout << "   solving the linear regressions (full pivoting LU)..." << std::endl;
			fittedKineticParameters = LU.solve(Y);
		}
	}

	void KineticsMap_CHEMKIN::FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters)
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

	void KineticsMap_CHEMKIN::DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* dR_over_domega)
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

	void KineticsMap_CHEMKIN::Derivatives(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* derivatives)
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

	void KineticsMap_CHEMKIN::DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* dR_over_dC)
	{
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OPENSMOKE_MACH_EPS_DOUBLE);			
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

	void KineticsMap_CHEMKIN::Derivatives(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* derivatives, const bool constant_density)
	{
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OPENSMOKE_MACH_EPS_DOUBLE);			
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

	void KineticsMap_CHEMKIN::FallOffReactions(const unsigned int k, const double cTot, const OpenSMOKEVectorDouble& c,
													double &F, double &dF_over_dA0, double &dF_over_dAInf)
	{
		double M;
		if (falloff_index_of_single_thirdbody_species_[k] == 0)
		{
			M = cTot;
			for(int s=1;s<=falloff_indices_of_thirdbody_species_[k-1].Size();s++)
				M += c[ falloff_indices_of_thirdbody_species_[k-1][s] ] * falloff_indices_of_thirdbody_efficiencies_[k-1][s];
		}
		else
		{
			const double epsilon = 1.e-16;
			M = c[falloff_index_of_single_thirdbody_species_[k]] + epsilon;
		}
			
		const unsigned int j=indices_of_falloff_reactions_[k];
		const double Pr = std::exp(lnA_[j] + Beta_[j]*this->logT_-E_over_R_[j]*this->uT_) * M / 
						  std::exp(lnA_falloff_inf_[k] + Beta_falloff_inf_[k]*this->logT_-E_over_R_falloff_inf_[k]*this->uT_);
		
		F = 1.;
		double dF_over_dPr = 0.;
		double nTroe, cTroe, sTroe, xSRI;
		switch(falloff_reaction_type_[k])
		{
			case PhysicalConstants::REACTION_TROE_FALLOFF:

				nTroe = 0.75-1.27*logFcent_falloff_[k];
				cTroe = -0.4-0.67*logFcent_falloff_[k];
				sTroe = std::log10(Pr) + cTroe;
				F = std::pow(10., logFcent_falloff_[k]/(1. + boost::math::pow<2>(sTroe/(nTroe-0.14*sTroe))));

				// Calculates the derivative 
				if (Pr > 1.e-32)
				{
					const double a = logFcent_falloff_[k];
					const double LOG10 = std::log(10.);
					const double LOGPR = std::log(Pr);

					dF_over_dPr =
					(-0.868589*std::pow(10., a / (1. + (9.62305*std::pow(cTroe*LOG10 + LOGPR, 2.)) / std::pow(cTroe - 7.14286*nTroe + 0.434294*LOGPR, 2.)))*a*nTroe*
					(-0.14*cTroe + nTroe - 0.0608012*LOGPR)*(cTroe*LOG10 + LOGPR)) / (Pr*std::pow(1.0196*std::pow(cTroe, 2.) - 
					0.28*cTroe*nTroe + std::pow(nTroe, 2.) +
					0.885613*cTroe*LOGPR - 0.121602*nTroe*LOGPR + 0.192308*std::pow(LOGPR, 2.), 2.));
				}

				break;

			case PhysicalConstants::REACTION_SRI_FALLOFF:

				// Calculates the derivative 
				if (Pr > 1.e-32)
				{
					xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
					F = std::pow(logFcent_falloff_[k], xSRI) * d_falloff_[k];
					if (e_falloff_[k] != 0.)
						F *= std::pow(this->T_, e_falloff_[k]);

					const double eps = 0.01;
					const double PrPlus = (1.+eps)*Pr;
					const double xSRIPlus = 1. / (1. + boost::math::pow<2>(std::log10(PrPlus)));
					double FPlus = std::pow(logFcent_falloff_[k], xSRIPlus) * d_falloff_[k];
					if(e_falloff_[k] != 0.)
						FPlus *= std::pow(this->T_, e_falloff_[k]);
					dF_over_dPr = (FPlus-F)/(eps*Pr);
				}
				else
				{
					xSRI = 0.;
					F = std::pow(logFcent_falloff_[k], xSRI) * d_falloff_[k];
					if (e_falloff_[k] != 0.)
						F *= std::pow(this->T_, e_falloff_[k]);
					dF_over_dPr = 0.;
				}

				break;
		}

		dF_over_dA0   =  Pr/std::exp(lnA_[j]) * dF_over_dPr;
		dF_over_dAInf = -Pr/std::exp(lnA_falloff_inf_[k]) * dF_over_dPr;
	}

	void KineticsMap_CHEMKIN::ChemicallyActivatedBimolecularReactions(const unsigned int k, const double cTot, const OpenSMOKEVectorDouble& c,
																		   double &F, double &dF_over_dA0, double &dF_over_dAInf)
	{
		double M;
		if (cabr_index_of_single_thirdbody_species_[k] == 0)
		{
			M = cTot;
			for(int s=1;s<=cabr_indices_of_thirdbody_species_[k-1].Size();s++)
				M += c[ cabr_indices_of_thirdbody_species_[k-1][s] ] * cabr_indices_of_thirdbody_efficiencies_[k-1][s];
		}
		else
		{
			const double epsilon = 1.e-16;
			M = c[cabr_index_of_single_thirdbody_species_[k]] + epsilon;
		}
			
		const unsigned int j=indices_of_cabr_reactions_[k];
		const double Pr = std::exp(lnA_[j] + Beta_[j]*this->logT_-E_over_R_[j]*this->uT_) * M / 
						  std::exp(lnA_cabr_inf_[k] + Beta_cabr_inf_[k]*this->logT_-E_over_R_cabr_inf_[k]*this->uT_);
		
		F = 1.;
		double dF_over_dPr = 0.;
		double nTroe, cTroe, sTroe, xSRI;
		switch(cabr_reaction_type_[k])
		{
			case PhysicalConstants::REACTION_TROE_CABR:

				nTroe = 0.75-1.27*logFcent_cabr_[k];
				cTroe = -0.4-0.67*logFcent_cabr_[k];
				sTroe = std::log10(Pr) + cTroe;
				F = std::pow(10., logFcent_cabr_[k]/(1. + boost::math::pow<2>(sTroe/(nTroe-0.14*sTroe))));

				// Calculates the derivative 
				if (Pr > 1.e-32)
				{
					const double eps = 0.01;
					const double PrPlus = (1.+eps)*Pr;
					const double nTroePlus = 0.75-1.27*logFcent_cabr_[k];
					const double cTroePlus = -0.4-0.67*logFcent_cabr_[k];
					const double sTroePlus = std::log10(PrPlus) + cTroePlus;
					const double FPlus = std::pow(10., logFcent_cabr_[k]/(1. + boost::math::pow<2>(sTroePlus/(nTroePlus-0.14*sTroePlus))));
					dF_over_dPr = (FPlus-F)/(eps*Pr);
				}

				break;

			case PhysicalConstants::REACTION_SRI_CABR:

				xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
				F = std::pow(logFcent_cabr_[k], xSRI) * d_cabr_[k];
				if(e_cabr_[k] != 0.)
					F *= std::pow(this->T_, e_cabr_[k]);

				// Calculates the derivative
				if (Pr > 1.e-32)
				{
					const double eps = 0.01;
					const double PrPlus = (1.+eps)*Pr;
					const double xSRIPlus = 1. / (1. + boost::math::pow<2>(std::log10(PrPlus)));
					double FPlus = std::pow(logFcent_cabr_[k], xSRIPlus) * d_cabr_[k];
					if(e_cabr_[k] != 0.)
						FPlus *= std::pow(this->T_, e_cabr_[k]);
					dF_over_dPr = (FPlus-F)/(eps*Pr);
				}

				break;
		}

		dF_over_dA0 = Pr/std::exp(lnA_[j]) * dF_over_dPr;
		dF_over_dAInf = -Pr/std::exp(lnA_falloff_inf_[k]) * dF_over_dPr;
	}
}

