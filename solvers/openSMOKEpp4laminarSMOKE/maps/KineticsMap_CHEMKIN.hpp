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

#include "math/OpenSMOKEUtilities.h"

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

        this->indices_of_irreversible_reactions__ = rhs.indices_of_irreversible_reactions__ ;
        this->indices_of_reversible_reactions__ = rhs.indices_of_reversible_reactions__ ;
        this->indices_of_thermodynamic_reversible_reactions__ = rhs.indices_of_thermodynamic_reversible_reactions__ ;
        this->indices_of_explicitly_reversible_reactions__ = rhs.indices_of_explicitly_reversible_reactions__ ;
        this->indices_of_thirdbody_reactions__ = rhs.indices_of_thirdbody_reactions__ ;
        this->indices_of_falloff_reactions__ = rhs.indices_of_falloff_reactions__ ;
		this->indices_of_extendedfalloff_reactions__ = rhs.indices_of_extendedfalloff_reactions__;
        this->indices_of_cabr_reactions__ = rhs.indices_of_cabr_reactions__ ;
        this->indices_of_chebyshev_reactions__ = rhs.indices_of_chebyshev_reactions__ ;
        this->indices_of_pressurelog_reactions__ = rhs.indices_of_pressurelog_reactions__ ;
		this->indices_of_extendedpressurelog_reactions__ = rhs.indices_of_extendedpressurelog_reactions__;
        this->indices_of_fit1_reactions__ = rhs.indices_of_fit1_reactions__ ;
        this->indices_of_janevlanger_reactions__ = rhs.indices_of_janevlanger_reactions__ ;
        this->indices_of_landauteller_reactions__ = rhs.indices_of_landauteller_reactions__ ;

        this->indices_of_falloff_lindemann_reactions__ = rhs.indices_of_falloff_lindemann_reactions__ ;
        this->indices_of_cabr_lindemann_reactions__ = rhs.indices_of_cabr_lindemann_reactions__ ;
        this->indices_of_falloff_troe_reactions__ = rhs.indices_of_falloff_troe_reactions__ ;
        this->indices_of_cabr_troe_reactions__ = rhs.indices_of_cabr_troe_reactions__ ;
        this->indices_of_falloff_sri_reactions__ = rhs.indices_of_falloff_sri_reactions__ ;
        this->indices_of_cabr_sri_reactions__ = rhs.indices_of_cabr_sri_reactions__ ;

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

        this->lnA__ = rhs.lnA__;
        this->Beta__ = rhs.Beta__;
        this->E_over_R__ = rhs.E_over_R__;
		this->negative_lnA__ = rhs.negative_lnA__;
		this->sign_lnA__ = rhs.sign_lnA__;

        this->lnA_reversible__ = rhs.lnA_reversible__;
        this->Beta_reversible__ = rhs.Beta_reversible__;
        this->E_over_R_reversible__ = rhs.E_over_R_reversible__;

        this->Meff__ = rhs.Meff__ ;
        this->indices_of_thirdbody_species__ = rhs.indices_of_thirdbody_species__ ;
        this->indices_of_thirdbody_efficiencies__ = rhs.indices_of_thirdbody_efficiencies__ ;

        this->lnA_falloff_inf__ = rhs.lnA_falloff_inf__;
        this->Beta_falloff_inf__ = rhs.Beta_falloff_inf__;
        this->E_over_R_falloff_inf__ = rhs.E_over_R_falloff_inf__;

        this->falloff_indices_of_thirdbody_species__ = rhs.falloff_indices_of_thirdbody_species__ ;
        this->falloff_indices_of_thirdbody_efficiencies__ = rhs.falloff_indices_of_thirdbody_efficiencies__ ;
        this->falloff_index_of_single_thirdbody_species__ = rhs.falloff_index_of_single_thirdbody_species__ ;

        this->falloff_reaction_type__ = rhs.falloff_reaction_type__ ;
        this->a_falloff__ = rhs.a_falloff__ ;
        this->b_falloff__ = rhs.b_falloff__ ;
        this->c_falloff__ = rhs.c_falloff__ ;
        this->d_falloff__ = rhs.d_falloff__ ;
        this->e_falloff__ = rhs.e_falloff__ ;
        this->logFcent_falloff__ = rhs.logFcent_falloff__ ;

        this->lnA_cabr_inf__ = rhs.lnA_cabr_inf__ ;
        this->Beta_cabr_inf__ = rhs.Beta_cabr_inf__ ;
        this->E_over_R_cabr_inf__ = rhs.E_over_R_cabr_inf__ ;

        this->cabr_indices_of_thirdbody_species__ = rhs.cabr_indices_of_thirdbody_species__ ;
        this->cabr_indices_of_thirdbody_efficiencies__ = rhs.cabr_indices_of_thirdbody_efficiencies__ ;
        this->cabr_index_of_single_thirdbody_species__ = rhs.cabr_index_of_single_thirdbody_species__ ;

        this->cabr_reaction_type__ = rhs.cabr_reaction_type__ ;
        this->a_cabr__ = rhs.a_cabr__ ;
        this->b_cabr__ = rhs.b_cabr__ ;
        this->c_cabr__ = rhs.c_cabr__ ;
        this->d_cabr__ = rhs.d_cabr__ ;
        this->e_cabr__ = rhs.e_cabr__ ;
        this->logFcent_cabr__ = rhs.logFcent_cabr__ ;

        this->changeOfMoles__ = rhs.changeOfMoles__ ;
            
        stoichiometry_ = new StoichiometricMap(*rhs.stoichiometry_);

        this->arrhenius_kinetic_constants_must_be_recalculated_ = true;
        this->nonconventional_kinetic_constants_must_be_recalculated_ = true;
        this->reaction_h_and_s_must_be_recalculated_ = true;
        this->isJacobianSparsityMapAvailable_ = false;

		reaction_s_over_R__.resize(rhs.reaction_s_over_R__.size());
		reaction_h_over_RT__.resize(rhs.reaction_h_over_RT__.size());
		kArrheniusModified__.resize(rhs.kArrheniusModified__.size());
		kArrhenius__.resize(rhs.kArrhenius__.size());
		uKeq__.resize(rhs.uKeq__.size());
		kArrhenius_reversible__.resize(rhs.kArrhenius_reversible__.size());
		kArrhenius_falloff_inf__.resize(rhs.kArrhenius_falloff_inf__.size());
		kArrhenius_cabr_inf__.resize(rhs.kArrhenius_cabr_inf__.size());

		forwardReactionRates__.resize(rhs.forwardReactionRates__.size());
		reverseReactionRates__.resize(rhs.reverseReactionRates__.size());
		netReactionRates__.resize(rhs.netReactionRates__.size());

		correction_falloff__.resize(rhs.correction_falloff__.size());
		correction_cabr__.resize(rhs.correction_cabr__.size());
            
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
                        
        this->isThermodynamicallyReversible__ = rhs.isThermodynamicallyReversible__ ;
        this->isExplicitlyReversible__ = rhs.isExplicitlyReversible__ ;

        this->type_of_reaction__ = rhs.type_of_reaction__ ;
        this->local_family_index__ = rhs.local_family_index__ ;
    }

	KineticsMap_CHEMKIN::~KineticsMap_CHEMKIN()
	{
		indices_of_thirdbody_species__.clear();			
		indices_of_thirdbody_efficiencies__.clear();		
		falloff_indices_of_thirdbody_species__.clear();
		falloff_indices_of_thirdbody_efficiencies__.clear();
		cabr_indices_of_thirdbody_species__.clear();
		cabr_indices_of_thirdbody_efficiencies__.clear();
                
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
		
		// Violation of Chebyshev Polynomials
		bool is_violation_allowed_in_chebyshev_polynomials = false;
		{
			rapidxml::xml_node<>* violation_chebyshev_node = opensmoke_node->first_node("ViolationChebyshev");
			if (violation_chebyshev_node != 0)
			{
				std::string flag = violation_chebyshev_node->first_attribute("flag")->value();
				if (flag == "true")
					is_violation_allowed_in_chebyshev_polynomials = true;
			}
		}
		
		// Kinetics
		rapidxml::xml_node<>* kinetics_node = opensmoke_node->first_node("Kinetics");
		if (kinetics_node == 0)
			ErrorMessage("void KineticsMap_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "Kinetics tag was not found!");

		std::string kinetics_type = kinetics_node->first_attribute("type")->value();
		std::string kinetics_version = kinetics_node->first_attribute("version")->value();

		if (kinetics_type != "OpenSMOKE" || (kinetics_version != "04-22-2013" && kinetics_version != "01-02-2014"))
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
			Load(indices_of_irreversible_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_irreversible_reactions_ = indices_of_irreversible_reactions__.size();
		}
		
		// Reversible reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible");
			std::stringstream fInput;
			fInput << current_node->value();
			Load(indices_of_reversible_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_reversible_reactions_ = indices_of_reversible_reactions__.size();
		}

		// Thermodynamic Reversible reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible-Thermodynamics");
			std::stringstream fInput;
			fInput << current_node->value();
			Load(indices_of_thermodynamic_reversible_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_thermodynamic_reversible_reactions_ = indices_of_thermodynamic_reversible_reactions__.size();
		}

		// Explicit Reversible reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("Reversible-Explicit");
			std::stringstream fInput;
			fInput << current_node->value();
			Load(indices_of_explicitly_reversible_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_explicitly_reversible_reactions_ = indices_of_explicitly_reversible_reactions__.size();
		}

		// Three-body reactions
		{
			rapidxml::xml_node<>* current_node = kinetics_node->first_node("ThreeBody");
			std::stringstream fInput;
			fInput << current_node->value();
			Load(indices_of_thirdbody_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
			number_of_thirdbody_reactions_  = indices_of_thirdbody_reactions__.size();
		}

		// FallOff
		{
			rapidxml::xml_node<>* falloff_node = kinetics_node->first_node("FallOff");

			{
				std::stringstream fInput;
				fInput << falloff_node->value();
				Load(indices_of_falloff_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_falloff_reactions_  = indices_of_falloff_reactions__.size();
			}

			// FallOff Lindemann
			{
				rapidxml::xml_node<>* current_node = falloff_node->first_node("Lindemann");
				std::stringstream fInput;
				fInput << current_node->value();
				Load(indices_of_falloff_lindemann_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_falloff_lindemann_reactions_  = indices_of_falloff_lindemann_reactions__.size();
			}

			// FallOff Troe
			{
				rapidxml::xml_node<>* current_node = falloff_node->first_node("Troe");
				std::stringstream fInput;
				fInput << current_node->value();
				Load(indices_of_falloff_troe_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_falloff_troe_reactions_  = indices_of_falloff_troe_reactions__.size();
			}

			// FallOff SRI
			{
				rapidxml::xml_node<>* current_node = falloff_node->first_node("SRI");
				std::stringstream fInput;
				fInput << current_node->value();
				Load(indices_of_falloff_sri_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_falloff_sri_reactions_  = indices_of_falloff_sri_reactions__.size();
			}
		}	

		// CABR
		{
			rapidxml::xml_node<>* cabr_node = kinetics_node->first_node("CABR");

			{
				std::stringstream fInput;
				fInput << cabr_node->value();
				Load(indices_of_cabr_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_cabr_reactions_  = indices_of_cabr_reactions__.size();
			}

			// CABR Lindemann
			{
				rapidxml::xml_node<>* current_node = cabr_node->first_node("Lindemann");
				std::stringstream fInput;
				fInput << current_node->value();
				Load(indices_of_cabr_lindemann_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_cabr_lindemann_reactions_  = indices_of_cabr_lindemann_reactions__.size();
			}

			// CABR Troe
			{
				rapidxml::xml_node<>* current_node = cabr_node->first_node("Troe");
				std::stringstream fInput;
				fInput << current_node->value();
				Load(indices_of_cabr_troe_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_cabr_troe_reactions_  = indices_of_cabr_troe_reactions__.size();
			}

			// CABR SRI
			{
				rapidxml::xml_node<>* current_node = cabr_node->first_node("SRI");
				std::stringstream fInput;
				fInput << current_node->value();
				Load(indices_of_cabr_sri_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_cabr_sri_reactions_  = indices_of_cabr_sri_reactions__.size();
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
				Load(indices_of_chebyshev_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_chebyshev_reactions_ = indices_of_chebyshev_reactions__.size();
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
				Load(indices_of_pressurelog_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_pressurelog_reactions_ = indices_of_pressurelog_reactions__.size();
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
				Load(indices_of_extendedpressurelog_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_extendedpressurelog_reactions_ = indices_of_extendedpressurelog_reactions__.size();
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
				Load(indices_of_extendedfalloff_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_extendedfalloff_reactions_ = indices_of_extendedfalloff_reactions__.size();
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
				Load(indices_of_fit1_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_fit1_reactions_ = indices_of_fit1_reactions__.size();
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
				Load(indices_of_janevlanger_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_janevlanger_reactions_ = indices_of_janevlanger_reactions__.size();
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
				Load(indices_of_landauteller_reactions__, fInput, OPENSMOKE_FORMATTED_FILE);
				number_of_landauteller_reactions_ = indices_of_landauteller_reactions__.size();
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
					Load(lnA__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// negative-lnA
				{
					// List of reactions with negative frequency factors
					rapidxml::xml_node<>* current_node = direct_node->first_node("negative-lnA");
					if (current_node != 0)
					{
						std::stringstream fInput;
						fInput << current_node->value();
						Load(negative_lnA__, fInput, OPENSMOKE_FORMATTED_FILE);
					}

					// Sign of frequency factors
					sign_lnA__.resize(lnA__.size());
					std::fill(sign_lnA__.begin(), sign_lnA__.end(), 1);
					for (unsigned int i = 0; i < negative_lnA__.size(); i++)
						sign_lnA__[negative_lnA__[i]-1] = -1;
				}

				// Beta
				{
					rapidxml::xml_node<>* current_node = direct_node->first_node("Beta");
					std::stringstream fInput;
					fInput << current_node->value();
					Load(Beta__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// E_over_R
				{
					rapidxml::xml_node<>* current_node = direct_node->first_node("E_over_R");
					std::stringstream fInput;
					fInput << current_node->value();
					Load(E_over_R__, fInput, OPENSMOKE_FORMATTED_FILE);
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
					Load(lnA_reversible__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Beta
				{
					rapidxml::xml_node<>* current_node = reverse_node->first_node("Beta");
					std::stringstream fInput;
					fInput << current_node->value();
					Load(Beta_reversible__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// E_over_R
				{
					rapidxml::xml_node<>* current_node = reverse_node->first_node("E_over_R");
					std::stringstream fInput;
					fInput << current_node->value();
					Load(E_over_R_reversible__, fInput, OPENSMOKE_FORMATTED_FILE);
				}
			}
                        
                        if(verbose_output_ == true)
			    std::cout << " * Reading kinetic parameters of third body reactions..." << std::endl;
			if (number_of_thirdbody_reactions_ != 0)
			{
				indices_of_thirdbody_species__.resize(number_of_thirdbody_reactions_);
				indices_of_thirdbody_efficiencies__.resize(number_of_thirdbody_reactions_);

				rapidxml::xml_node<>* threebody_node = kinetic_parameters_node->first_node("ThreeBody");
			
				unsigned int i = 0;
				for (rapidxml::xml_node<> *current_node = threebody_node->first_node("reaction"); current_node; current_node = current_node->next_sibling())
				{
					std::stringstream fInput;
					fInput << current_node->value();

					Load(indices_of_thirdbody_species__[i], fInput, OPENSMOKE_FORMATTED_FILE);
					Load(indices_of_thirdbody_efficiencies__[i], fInput, OPENSMOKE_FORMATTED_FILE);

					for (std::vector<double>::iterator it = indices_of_thirdbody_efficiencies__[i].begin(); it != indices_of_thirdbody_efficiencies__[i].end(); ++it)
						*it -= 1.;

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
					Load(lnA_falloff_inf__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Beta
				{
					rapidxml::xml_node<>* current_node = falloff_node->first_node("Beta");
					std::stringstream fInput;
					fInput << current_node->value();
					Load(Beta_falloff_inf__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// E_over_R
				{
					rapidxml::xml_node<>* current_node = falloff_node->first_node("E_over_R");
					std::stringstream fInput;
					fInput << current_node->value();
					Load(E_over_R_falloff_inf__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Parameters
				{
					falloff_reaction_type__.resize(number_of_falloff_reactions_);

					falloff_index_of_single_thirdbody_species__.resize(number_of_falloff_reactions_);
					std::fill(falloff_index_of_single_thirdbody_species__.begin(), falloff_index_of_single_thirdbody_species__.end(), 0);

					a_falloff__.resize(number_of_falloff_reactions_);
					std::fill(a_falloff__.begin(), a_falloff__.end(), 0.);

					b_falloff__.resize(number_of_falloff_reactions_);
					std::fill(b_falloff__.begin(), b_falloff__.end(), 0.);

					c_falloff__.resize(number_of_falloff_reactions_);
					std::fill(c_falloff__.begin(), c_falloff__.end(), 0.);

					d_falloff__.resize(number_of_falloff_reactions_);
					std::fill(d_falloff__.begin(), d_falloff__.end(), 0.);

					e_falloff__.resize(number_of_falloff_reactions_);
					std::fill(e_falloff__.begin(), e_falloff__.end(), 0.);

					falloff_indices_of_thirdbody_species__.resize(number_of_falloff_reactions_);
					falloff_indices_of_thirdbody_efficiencies__.resize(number_of_falloff_reactions_);

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
							falloff_reaction_type__[i-1] = PhysicalConstants::REACTION_LINDEMANN_FALLOFF;
						}
						if (dummy == "troe")
						{
							falloff_reaction_type__[i-1] = PhysicalConstants::REACTION_TROE_FALLOFF;

							OpenSMOKEVectorDouble coefficients;
							coefficients.Load(fInput, OPENSMOKE_FORMATTED_FILE);
							a_falloff__[i-1] = coefficients[1];
							b_falloff__[i-1] = coefficients[2];
							c_falloff__[i-1] = coefficients[3];
							if (coefficients.Size() == 4)
								d_falloff__[i-1] = coefficients[4];
						}
						else if (dummy == "sri")
						{
							falloff_reaction_type__[i-1] = PhysicalConstants::REACTION_SRI_FALLOFF;
					
							OpenSMOKEVectorDouble coefficients;
							coefficients.Load(fInput, OPENSMOKE_FORMATTED_FILE);
							a_falloff__[i-1] = coefficients[1];
							b_falloff__[i-1] = coefficients[2];
							c_falloff__[i-1] = coefficients[3];
							if (coefficients.Size() == 5)
							{
								d_falloff__[i-1] = coefficients[4];
								e_falloff__[i-1] = coefficients[5];
							}
							else
							{
								d_falloff__[i-1] = 1.;
								e_falloff__[i-1] = 0.;
							}
						}

						fInput >> dummy;
						if (dummy == "species")
						{
							unsigned int species;
							fInput >> species;
							falloff_index_of_single_thirdbody_species__[i-1] = species;
						}
						else
						{
							Load(falloff_indices_of_thirdbody_species__[i - 1], fInput, OPENSMOKE_FORMATTED_FILE);
							Load(falloff_indices_of_thirdbody_efficiencies__[i - 1], fInput, OPENSMOKE_FORMATTED_FILE);
							
							for (std::vector<double>::iterator it = falloff_indices_of_thirdbody_efficiencies__[i - 1].begin(); it != falloff_indices_of_thirdbody_efficiencies__[i - 1].end(); ++it)
								*it -= 1.;
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
					Load(lnA_cabr_inf__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Beta
				{
					rapidxml::xml_node<>* current_node = cabr_node->first_node("Beta");
					std::stringstream fInput;
					fInput << current_node->value();
					Load(Beta_cabr_inf__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// E_over_R
				{
					rapidxml::xml_node<>* current_node = cabr_node->first_node("E_over_R");
					std::stringstream fInput;
					fInput << current_node->value();
					Load(E_over_R_cabr_inf__, fInput, OPENSMOKE_FORMATTED_FILE);
				}

				// Parameters
				{
					cabr_reaction_type__.resize(number_of_cabr_reactions_);

					cabr_index_of_single_thirdbody_species__.resize(number_of_cabr_reactions_);
					std::fill(cabr_index_of_single_thirdbody_species__.begin(), cabr_index_of_single_thirdbody_species__.end(), 0);

					a_cabr__.resize(number_of_cabr_reactions_);
					std::fill(a_cabr__.begin(), a_cabr__.end(), 0.);

					b_cabr__.resize(number_of_cabr_reactions_);
					std::fill(b_cabr__.begin(), b_cabr__.end(), 0.);

					c_cabr__.resize(number_of_cabr_reactions_);
					std::fill(c_cabr__.begin(), c_cabr__.end(), 0.);

					d_cabr__.resize(number_of_cabr_reactions_);
					std::fill(d_cabr__.begin(), d_cabr__.end(), 0.);

					e_cabr__.resize(number_of_cabr_reactions_);
					std::fill(e_cabr__.begin(), e_cabr__.end(), 0.);

					cabr_indices_of_thirdbody_species__.resize(number_of_cabr_reactions_);
					cabr_indices_of_thirdbody_efficiencies__.resize(number_of_cabr_reactions_);

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
							cabr_reaction_type__[i-1] = PhysicalConstants::REACTION_LINDEMANN_CABR;
						}
						if (dummy == "troe")
						{
							cabr_reaction_type__[i-1] = PhysicalConstants::REACTION_TROE_CABR;

							OpenSMOKEVectorDouble coefficients;
							coefficients.Load(fInput, OPENSMOKE_FORMATTED_FILE);
							a_cabr__[i-1] = coefficients[1];
							b_cabr__[i-1] = coefficients[2];
							c_cabr__[i-1] = coefficients[3];
							if (coefficients.Size() == 4)
								d_cabr__[i-1] = coefficients[4];
						}
						else if (dummy == "sri")
						{
							cabr_reaction_type__[i-1] = PhysicalConstants::REACTION_SRI_CABR;
					
							OpenSMOKEVectorDouble coefficients;
							coefficients.Load(fInput, OPENSMOKE_FORMATTED_FILE);
							a_cabr__[i-1] = coefficients[1];
							b_cabr__[i-1] = coefficients[2];
							c_cabr__[i-1] = coefficients[3];
							if (coefficients.Size() == 5)
							{
								d_cabr__[i-1] = coefficients[4];
								e_cabr__[i-1] = coefficients[5];
							}
							else
							{
								d_cabr__[i-1] = 1.;
								e_cabr__[i-1] = 0.;
							}
						}

						fInput >> dummy;
						if (dummy == "species")
						{
							unsigned int species;
							fInput >> species;
							cabr_index_of_single_thirdbody_species__[i-1] = species;
						}
						else
						{
							Load(cabr_indices_of_thirdbody_species__[i - 1], fInput, OPENSMOKE_FORMATTED_FILE);
							Load(cabr_indices_of_thirdbody_efficiencies__[i - 1], fInput, OPENSMOKE_FORMATTED_FILE);
							
							for (std::vector<double>::iterator it = cabr_indices_of_thirdbody_efficiencies__[i - 1].begin(); it != cabr_indices_of_thirdbody_efficiencies__[i - 1].end(); ++it)
								*it -= 1.;
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

					if (is_violation_allowed_in_chebyshev_polynomials == true)
					{
						std::cout << " * WARNING: violation of T and/or P limits in Chebyshev Polynomials is turned on" << std::endl;
						for (unsigned int j = 0; j < number_of_chebyshev_reactions_; j++)
							chebyshev_reactions_[j].SetViolationAllowed(true);
					}
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

			if (stoichiometry_type != "OpenSMOKE" || (stoichiometry_version != "04-22-2013" && stoichiometry_version != "01-02-2014"))
				ErrorMessage("void KineticsMap_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)", "The current stoichiometric data are not supported.");

			std::stringstream fInput;
			fInput << stoichiometry_node->value();

			stoichiometry_->ReadFromASCIIFile(fInput);
			changeOfMoles__ = stoichiometry_->ChangeOfMoles();
			{
				OpenSMOKE::OpenSMOKEVectorBool tmp(this->number_of_reactions_);
				tmp = false;
				for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
					tmp[indices_of_thermodynamic_reversible_reactions__[k-1]] = true;
				stoichiometry_->CompleteChangeOfMoles(tmp.GetHandle());
 			}
			
			// Memory allocation
			reaction_s_over_R__.resize(this->number_of_reactions_);
			std::fill(reaction_s_over_R__.begin(), reaction_s_over_R__.end(), 0.);
			
			reaction_h_over_RT__.resize(this->number_of_reactions_);
			std::fill(reaction_h_over_RT__.begin(), reaction_h_over_RT__.end(), 0.);
			
			kArrhenius__.resize(this->number_of_reactions_);
			std::fill(kArrhenius__.begin(), kArrhenius__.end(), 0.);
			
			kArrheniusModified__.resize(this->number_of_reactions_);
			std::fill(kArrheniusModified__.begin(), kArrheniusModified__.end(), 0.);
			
			uKeq__.resize(this->number_of_thermodynamic_reversible_reactions_);
			std::fill(uKeq__.begin(), uKeq__.end(), 0.);
		
			kArrhenius_reversible__.resize(this->number_of_explicitly_reversible_reactions_);
			std::fill(kArrhenius_reversible__.begin(), kArrhenius_reversible__.end(), 0.);
		
			kArrhenius_falloff_inf__.resize(this->number_of_falloff_reactions_);
			std::fill(kArrhenius_falloff_inf__.begin(), kArrhenius_falloff_inf__.end(), 0.);

			logFcent_falloff__.resize(number_of_falloff_reactions_);
			std::fill(logFcent_falloff__.begin(), logFcent_falloff__.end(), 0.);

			kArrhenius_cabr_inf__.resize(number_of_cabr_reactions_);
			std::fill(kArrhenius_cabr_inf__.begin(), kArrhenius_cabr_inf__.end(), 0.);
			
			logFcent_cabr__.resize(number_of_cabr_reactions_);
			std::fill(logFcent_cabr__.begin(), logFcent_cabr__.end(), 0.);

			Meff__.resize(number_of_thirdbody_reactions_);
			std::fill(Meff__.begin(), Meff__.end(), 0.);

			correction_falloff__.resize(number_of_falloff_reactions_);
			std::fill(correction_falloff__.begin(), correction_falloff__.end(), 0.);

			correction_cabr__.resize(number_of_cabr_reactions_);
			std::fill(correction_cabr__.begin(), correction_cabr__.end(), 0.);

			forwardReactionRates__.resize(this->number_of_reactions_);
			std::fill(forwardReactionRates__.begin(), forwardReactionRates__.end(), 0.);

			reverseReactionRates__.resize(this->number_of_reactions_);
			std::fill(reverseReactionRates__.begin(), reverseReactionRates__.end(), 0.);

			netReactionRates__.resize(this->number_of_reactions_);
			std::fill(netReactionRates__.begin(), netReactionRates__.end(), 0.);

			isThermodynamicallyReversible__.resize(this->number_of_reactions_);
			std::fill(isThermodynamicallyReversible__.begin(), isThermodynamicallyReversible__.end(), 0);
			
			isExplicitlyReversible__.resize(this->number_of_reactions_);
			std::fill(isExplicitlyReversible__.begin(), isExplicitlyReversible__.end(), 0);

			for(unsigned int k=1;k<=number_of_thermodynamic_reversible_reactions_;k++)
				isThermodynamicallyReversible__[indices_of_thermodynamic_reversible_reactions__[k-1]-1] = k;
			for(unsigned int k=1;k<=number_of_explicitly_reversible_reactions_;k++)
				isExplicitlyReversible__[indices_of_explicitly_reversible_reactions__[k-1]-1] = k;

			// Additional indices for sensitivity analysis
			{
				type_of_reaction__.resize(this->number_of_reactions_);
				std::fill(type_of_reaction__.begin(), type_of_reaction__.end(), PhysicalConstants::REACTION_SIMPLE);

				local_family_index__.resize(this->number_of_reactions_);
				std::fill(local_family_index__.begin(), local_family_index__.end(), 0);

				for(unsigned int k=1;k<=number_of_thirdbody_reactions_;k++)
				{
					type_of_reaction__[indices_of_thirdbody_reactions__[k-1]-1] = PhysicalConstants::REACTION_THIRDBODY;
					local_family_index__[indices_of_thirdbody_reactions__[k-1]-1] = k;
				}
				for(unsigned int k=1;k<=number_of_chebyshev_reactions_;k++)
				{
					type_of_reaction__[indices_of_chebyshev_reactions__[k-1]-1] = PhysicalConstants::REACTION_CHEBYSHEV;
					local_family_index__[indices_of_chebyshev_reactions__[k-1]-1] = k;
				}

				for(unsigned int k=1;k<=number_of_falloff_lindemann_reactions_;k++)
					type_of_reaction__[indices_of_falloff_lindemann_reactions__[k-1]-1] = PhysicalConstants::REACTION_LINDEMANN_FALLOFF;	
				for(unsigned int k=1;k<=number_of_falloff_troe_reactions_;k++)
					type_of_reaction__[indices_of_falloff_troe_reactions__[k-1]-1] = PhysicalConstants::REACTION_TROE_FALLOFF;
				for(unsigned int k=1;k<=number_of_falloff_sri_reactions_;k++)
					type_of_reaction__[indices_of_falloff_sri_reactions__[k-1]-1] = PhysicalConstants::REACTION_SRI_FALLOFF;
				for(unsigned int k=1;k<=number_of_falloff_reactions_;k++)
					local_family_index__[indices_of_falloff_reactions__[k-1]-1] = k;

				for(unsigned int k=1;k<=number_of_cabr_lindemann_reactions_;k++)
					type_of_reaction__[indices_of_cabr_lindemann_reactions__[k-1]-1] = PhysicalConstants::REACTION_LINDEMANN_CABR;
				for(unsigned int k=1;k<=number_of_cabr_troe_reactions_;k++)
					type_of_reaction__[indices_of_cabr_troe_reactions__[k-1]-1] = PhysicalConstants::REACTION_TROE_CABR;
				for(unsigned int k=1;k<=number_of_cabr_sri_reactions_;k++)
					type_of_reaction__[indices_of_cabr_sri_reactions__[k-1]-1] = PhysicalConstants::REACTION_SRI_CABR;
				for(unsigned int k=1;k<=number_of_cabr_reactions_;k++)
					local_family_index__[indices_of_cabr_reactions__[k-1]-1] = k;
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
				std::cout << " Negative frequency factors:     " << negative_lnA__.size() << " (" << negative_lnA__.size() / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
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
			stoichiometry_->ReactionEnthalpyAndEntropy(	reaction_h_over_RT__, reaction_s_over_R__, 
														thermodynamics_.Species_H_over_RT(), thermodynamics_.Species_S_over_R() );

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
				double *pt_lnA = lnA__.data();
				double *pt_Beta = Beta__.data();
				double *pt_E_over_R = E_over_R__.data();
				double *pt_kArrheniusT = kArrhenius__.data();
			
				for(unsigned int j=0;j<this->number_of_reactions_;j++)
					*pt_kArrheniusT++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius__, &kArrhenius__);

				// Negative frequency factors: reactions must be reversed
				for (unsigned int j = 0; j < negative_lnA__.size(); j++)
					kArrhenius__[negative_lnA__[j]-1] *= -1.;
			}

			// Equilibrium constants (inverse value)
			{
				for(unsigned int k=0;k<number_of_thermodynamic_reversible_reactions_;k++)
				{
					unsigned int j = indices_of_thermodynamic_reversible_reactions__[k]-1;
					uKeq__[k] = -reaction_s_over_R__[j] + reaction_h_over_RT__[j] - log_Patm_over_RT_ * changeOfMoles__[j];
				}
				Exp(uKeq__, &uKeq__);

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
				double *pt_lnA = lnA_reversible__.data();
				double *pt_Beta = Beta_reversible__.data();
				double *pt_E_over_R = E_over_R_reversible__.data();
				double *pt_kArrhenius = kArrhenius_reversible__.data();

				for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
					*pt_kArrhenius++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius_reversible__, &kArrhenius_reversible__);
			}

			// Fall-off high temperature region kinetic constants
			if (number_of_falloff_reactions_ != 0)
			{
				double *pt_lnA = lnA_falloff_inf__.data();
				double *pt_Beta = Beta_falloff_inf__.data();
				double *pt_E_over_R = E_over_R_falloff_inf__.data();
				double *pt_kArrhenius = kArrhenius_falloff_inf__.data();

				for(unsigned int k=0;k<number_of_falloff_reactions_;k++)
					*pt_kArrhenius++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius_falloff_inf__, &kArrhenius_falloff_inf__);

				for(unsigned int k=0;k<number_of_falloff_reactions_;k++)
				{
					switch(falloff_reaction_type__[k])
					{
						case PhysicalConstants::REACTION_TROE_FALLOFF:
							
							logFcent_falloff__[k] = (d_falloff__[k] != 0.) ? (1.-a_falloff__[k])*std::exp(-this->T_/b_falloff__[k]) + a_falloff__[k]*std::exp(-this->T_/c_falloff__[k]) + std::exp(-d_falloff__[k]/this->T_) :
																		 (1.-a_falloff__[k])*std::exp(-this->T_/b_falloff__[k]) + a_falloff__[k]*std::exp(-this->T_/c_falloff__[k]);
							
							if(logFcent_falloff__[k] < 1.e-300)	logFcent_falloff__[k] = -300.;
							else								logFcent_falloff__[k] = std::log10(logFcent_falloff__[k]);

							break;
						
						case PhysicalConstants::REACTION_SRI_FALLOFF:

							logFcent_falloff__[k] = ( a_falloff__[k]*std::exp(-b_falloff__[k]/this->T_) +  std::exp(-this->T_/c_falloff__[k]));
							
							break;
					}
				}
			}

			// Cabr high temperature region kinetic constants
			if (number_of_cabr_reactions_ != 0)
			{
				double *pt_lnA = lnA_cabr_inf__.data();
				double *pt_Beta = Beta_cabr_inf__.data();
				double *pt_E_over_R = E_over_R_cabr_inf__.data();
				double *pt_kArrhenius = kArrhenius_cabr_inf__.data();

				for(unsigned int k=0;k<number_of_cabr_reactions_;k++)
					*pt_kArrhenius++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius_cabr_inf__, &kArrhenius_cabr_inf__);

				for(unsigned int k=0;k<number_of_cabr_reactions_;k++)
				{
					switch(cabr_reaction_type__[k])
					{
						case PhysicalConstants::REACTION_TROE_CABR:
							
							logFcent_cabr__[k] = (d_cabr__[k] != 0.) ? (1.-a_cabr__[k])*std::exp(-this->T_/b_cabr__[k]) + a_cabr__[k]*std::exp(-this->T_/c_cabr__[k]) + std::exp(-d_cabr__[k]/this->T_) :
																		 (1.-a_cabr__[k])*std::exp(-this->T_/b_cabr__[k]) + a_cabr__[k]*std::exp(-this->T_/c_cabr__[k]);
							
							if(logFcent_cabr__[k] < 1.e-300)	logFcent_cabr__[k] = -300.;
							else							logFcent_cabr__[k] = std::log10(logFcent_cabr__[k]);
							
							break;
						
						case PhysicalConstants::REACTION_SRI_CABR:

							logFcent_cabr__[k] = ( a_cabr__[k]*std::exp(-b_cabr__[k]/this->T_) +  std::exp(-this->T_/c_cabr__[k]));
							
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
				for(unsigned int k=0;k<number_of_chebyshev_reactions_;k++)
				{
					const unsigned int j = indices_of_chebyshev_reactions__[k];
					kArrhenius__[j-1] = chebyshev_reactions_[k].KineticConstant(this->T_, this->P_);
				}
			}

			// Pressure logarithmic interpolated reactions
			{
				for(unsigned int k=0;k<number_of_pressurelog_reactions_;k++)
				{
					const unsigned int j = indices_of_pressurelog_reactions__[k];
					kArrhenius__[j-1] = pressurelog_reactions_[k].KineticConstant(this->T_,this->P_);
				}
			}

			// Extended pressure logarithmic interpolated reactions
			// Extended falloff reactions
			// They are calculated separately

			// Fit1 reactions
			{
				for(unsigned int k=0;k<number_of_fit1_reactions_;k++)
				{
					const unsigned int j = indices_of_fit1_reactions__[k];
				//	kArrhenius_[j] = fit1_reactions_[k].KineticConstant(this->T_,this->P_);
				}
			}

			// Janev-Langer reactions
			{
				for(unsigned int k=0;k<number_of_janevlanger_reactions_;k++)
				{
					const unsigned int j = indices_of_janevlanger_reactions__[k];
				//	kArrhenius_[j] = janevlanger_reactions_[k].KineticConstant(this->T_,this->P_);
				}
			}

			// Landau-Teller reactions
			{
				for(unsigned int k=0;k<number_of_landauteller_reactions_;k++)
				{
					const unsigned int j = indices_of_landauteller_reactions__[k];
				//	kArrhenius_[j] = landauteller_reactions_[k].KineticConstant(this->T_,this->P_);
				}
			}

			nonconventional_kinetic_constants_must_be_recalculated_ = false;
		}

		kArrheniusModified__ = kArrhenius__;
	}

	void KineticsMap_CHEMKIN::ThirdBodyReactions(const double cTot, const double* c)
	{	
		for(unsigned int s=0;s<number_of_thirdbody_reactions_;s++)
		{
			const unsigned int j=indices_of_thirdbody_reactions__[s];

			Meff__[s] = cTot;
			for(unsigned int k=0;k<indices_of_thirdbody_species__[s].size();k++)
				Meff__[s] += c[indices_of_thirdbody_species__[s][k] - 1] * indices_of_thirdbody_efficiencies__[s][k];
		}
	}

	void KineticsMap_CHEMKIN::WeakThirdBodyConcentrationEfficiencies(std::vector<unsigned int>& reactions, std::vector<unsigned int>& species)
	{
		reactions.resize(0);
		species.resize(0);

		// Third-body third-body
		for (unsigned int s = 0; s < number_of_thirdbody_reactions_; s++)
		{
			const unsigned int j = indices_of_thirdbody_reactions__[s];

			for (unsigned int k = 0; k < indices_of_thirdbody_species__[s].size(); k++)
			{
				reactions.push_back(j);
				species.push_back(indices_of_thirdbody_species__[s][k]);
			}
		}
	}

	void KineticsMap_CHEMKIN::WeakFallOffConcentrationEfficiencies(std::vector<unsigned int>& reactions, std::vector<unsigned int>& species)
	{
		reactions.resize(0);
		species.resize(0);

		// Fall-off reactions
		for (unsigned int k = 0; k < number_of_falloff_reactions_; k++)
		{
			const unsigned int j = indices_of_falloff_reactions__[k];

			if (falloff_index_of_single_thirdbody_species__[k] == 0)
			{
				for (unsigned int s = 0; s < falloff_indices_of_thirdbody_species__[k].size(); s++)
				{
					reactions.push_back(j);
					species.push_back(falloff_indices_of_thirdbody_species__[k][s]);
				}
			}
			else
			{
				reactions.push_back(j);
				species.push_back(falloff_index_of_single_thirdbody_species__[k]);
			}
		}
	}

	void KineticsMap_CHEMKIN::WeakCABRConcentrationEfficiencies(std::vector<unsigned int>& reactions, std::vector<unsigned int>& species)
	{
		reactions.resize(0);
		species.resize(0);

		// CABR
		for (unsigned int k = 0; k < number_of_cabr_reactions_; k++)
		{
			const unsigned int j = indices_of_cabr_reactions__[k];

			if (cabr_index_of_single_thirdbody_species__[k] == 0)
			{
				for (unsigned int s = 0; s < cabr_indices_of_thirdbody_species__[k].size(); s++)
				{
					reactions.push_back(j);
					species.push_back(cabr_indices_of_thirdbody_species__[k][s]);
				}
			}
			else
			{
				reactions.push_back(j);
				species.push_back(cabr_index_of_single_thirdbody_species__[k]);
			}
		}
	}

	void KineticsMap_CHEMKIN::StrongConcentrationEfficiencies(std::vector<unsigned int>& reactions)
	{
		reactions.resize(0);

		// Third-body third-body
		for (unsigned int s = 0; s < number_of_thirdbody_reactions_; s++)
		{
			const unsigned int j = indices_of_thirdbody_reactions__[s];
			reactions.push_back(j);
		}

		// Fall-off reactions
		for (unsigned int k = 0; k < number_of_falloff_reactions_; k++)
		{
			const unsigned int j = indices_of_falloff_reactions__[k];
			reactions.push_back(j);
		}

		// CABR
		for (unsigned int k = 0; k < number_of_cabr_reactions_; k++)
		{
			const unsigned int j = indices_of_cabr_reactions__[k];
			reactions.push_back(j);
		}
	}
	
	void KineticsMap_CHEMKIN::FallOffReactions(const double cTot, const double* c)
	{
		//if (arrhenius_kinetic_constants_must_be_corrected_ == true)
		{
			for(unsigned int k=0;k<number_of_falloff_reactions_;k++)
			{
				double M;
				if (falloff_index_of_single_thirdbody_species__[k] == 0)
				{
					M = cTot;
					for(unsigned int s=0;s<falloff_indices_of_thirdbody_species__[k].size();s++)
						M += c[falloff_indices_of_thirdbody_species__[k][s]-1] * falloff_indices_of_thirdbody_efficiencies__[k][s];
				}
				else
				{
					const double epsilon = 1.e-16;
					M = c[falloff_index_of_single_thirdbody_species__[k]-1] + epsilon;
				}
			
				const unsigned int j=indices_of_falloff_reactions__[k];
				const double Pr = kArrhenius__[j-1] * M / kArrhenius_falloff_inf__[k];
		
				double wF = 1.;
				switch(falloff_reaction_type__[k])
				{
					case PhysicalConstants::REACTION_TROE_FALLOFF:

                        if (Pr > 1.e-32)
                        {
                            const double nTroe = 0.75-1.27*logFcent_falloff__[k];
                            const double cTroe = -0.4-0.67*logFcent_falloff__[k];
                            const double sTroe = std::log10(Pr) + cTroe;
                            wF = std::pow(10., logFcent_falloff__[k]/(1. + boost::math::pow<2>(sTroe/(nTroe-0.14*sTroe)))); 
                        }
                        else
                        {
                            // Asymptotic value for wF when sTroe --> -Inf
                            wF = std::pow(10., logFcent_falloff__[k]/(1. + boost::math::pow<2>(1./0.14)));
                        }
                                                
						break;

					case PhysicalConstants::REACTION_SRI_FALLOFF:

						const double xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
						wF = std::pow(logFcent_falloff__[k], xSRI) * d_falloff__[k];
						if(e_falloff__[k] != 0.)
							wF *= std::pow(this->T_, e_falloff__[k]);
						break;
				}
	
				correction_falloff__[k] = kArrhenius_falloff_inf__[k] * (Pr/(1.+Pr)) * wF / kArrhenius__[j-1];
			}
		}
	}

	double KineticsMap_CHEMKIN::FallOffReactionsCorrection(const unsigned int local_k, const double cTot, const double* c)
	{
		double M;
		if (falloff_index_of_single_thirdbody_species__[local_k-1] == 0)
		{
			M = cTot;
			for (unsigned int s = 0; s < falloff_indices_of_thirdbody_species__[local_k - 1].size(); s++)
				M += c[falloff_indices_of_thirdbody_species__[local_k - 1][s]-1] * falloff_indices_of_thirdbody_efficiencies__[local_k - 1][s];
		}
		else
		{
			const double epsilon = 1.e-16;
			M = c[falloff_index_of_single_thirdbody_species__[local_k-1]-1] + epsilon;
		}

		const unsigned int j = indices_of_falloff_reactions__[local_k-1];
		const double Pr = kArrhenius__[j-1] * M / kArrhenius_falloff_inf__[local_k-1];

		double wF = 1.;
		switch (falloff_reaction_type__[local_k-1])
		{
			case PhysicalConstants::REACTION_TROE_FALLOFF:

				if (Pr > 1.e-32)
				{
					const double nTroe = 0.75 - 1.27*logFcent_falloff__[local_k-1];
					const double cTroe = -0.4 - 0.67*logFcent_falloff__[local_k-1];
					const double sTroe = std::log10(Pr) + cTroe;
					wF = std::pow(10., logFcent_falloff__[local_k-1] / (1. + boost::math::pow<2>(sTroe / (nTroe - 0.14*sTroe))));
				}
				else
				{
					// Asymptotic value for wF when sTroe --> -Inf
					wF = std::pow(10., logFcent_falloff__[local_k-1] / (1. + boost::math::pow<2>(1. / 0.14)));
				}

				break;

			case PhysicalConstants::REACTION_SRI_FALLOFF:

				if (Pr > 1.e-32)
				{
					const double xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
					wF = std::pow(logFcent_falloff__[local_k-1], xSRI) * d_falloff__[local_k-1];
					if (e_falloff__[local_k-1] != 0.)
						wF *= std::pow(this->T_, e_falloff__[local_k-1]);
				}
				else
				{
					const double xSRI = 0.;
					wF = std::pow(logFcent_falloff__[local_k-1], xSRI) * d_falloff__[local_k-1];
					if (e_falloff__[local_k-1] != 0.)
						wF *= std::pow(this->T_, e_falloff__[local_k-1]);
				}

				break;
		}

		return ( kArrhenius_falloff_inf__[local_k-1] * (Pr / (1. + Pr)) * wF / kArrhenius__[j-1] );
	}
 
	void KineticsMap_CHEMKIN::ChemicallyActivatedBimolecularReactions(const double cTot, const double* c)
	{
		//if (arrhenius_kinetic_constants_must_be_corrected_ == true)
		{
			for(unsigned int k=0;k<number_of_cabr_reactions_;k++)
			{
				double M;
				if (cabr_index_of_single_thirdbody_species__[k] == 0)
				{
					M = cTot;
					for(unsigned int s=0;s<cabr_indices_of_thirdbody_species__[k].size();s++)
						M += c[cabr_indices_of_thirdbody_species__[k][s]-1] * cabr_indices_of_thirdbody_efficiencies__[k][s];
				}
				else
				{
					const double epsilon = 1.e-16;
					M = c[cabr_index_of_single_thirdbody_species__[k]-1] + epsilon;
				}
			
				const unsigned int j=indices_of_cabr_reactions__[k];
				const double Pr = kArrhenius__[j-1] * M / kArrhenius_cabr_inf__[k];
		
				double wF = 1.;
				double nTroe, cTroe, sTroe, xSRI;
				switch(cabr_reaction_type__[k])
				{
					case PhysicalConstants::REACTION_TROE_CABR:

						nTroe = 0.75-1.27*logFcent_cabr__[k];
						cTroe = -0.4-0.67*logFcent_cabr__[k];
						sTroe = std::log10(Pr) + cTroe;
						wF = std::pow(10., logFcent_cabr__[k]/(1. + boost::math::pow<2>(sTroe/(nTroe-0.14*sTroe))));
						
						break;
					
					case PhysicalConstants::REACTION_SRI_CABR:
						xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
						wF = std::pow(logFcent_cabr__[k], xSRI) * d_cabr__[k];
						if(e_cabr__[k] != 0.)
							wF *= std::pow(this->T_, e_cabr__[k]);
						break;
				}
	
				correction_cabr__[k] = (1./(1.+Pr)) * wF;
			}
		}
	}

	// Extended Pressure logarithmic interpolated reactions
	void KineticsMap_CHEMKIN::ExtendedPressureLogReactions(const double cTot, const double* c)
	{
		for (unsigned int k = 0; k < number_of_extendedpressurelog_reactions_; k++)
		{
			const unsigned int j = indices_of_extendedpressurelog_reactions__[k]-1;
			kArrhenius__[j] = extendedpressurelog_reactions_[k].KineticConstant(this->T_, this->P_, cTot, c);
			kArrheniusModified__[j] = kArrhenius__[j];
		}
	}

	// Extended Falloff reactions
	void KineticsMap_CHEMKIN::ExtendedFallOffReactions(const double cTot, const double* c)
	{
		for (unsigned int k = 0; k < number_of_extendedfalloff_reactions_; k++)
		{
			const unsigned int j = indices_of_extendedfalloff_reactions__[k]-1;
			kArrhenius__[j] = extendedfalloff_reactions_[k].KineticConstant(this->T_, this->P_, cTot, c);
			kArrheniusModified__[j] = kArrhenius__[j];
		}
	}
	
	void KineticsMap_CHEMKIN::ReactionRates(const double* c)
	{
		double cTot = 0.;
		for(unsigned int i=0;i<this->number_of_species_;i++)
			cTot += c[i];

		ReactionRates(c, cTot);
	}

	void KineticsMap_CHEMKIN::ReactionRates(const double* c, const double cTot)
	{
		// 1. Kinetic constants
		KineticConstants();

		// 1.bis Extended pressure log reactions
		ExtendedPressureLogReactions(cTot, c);

		// 2. Calculates the three-body corrections
		ThirdBodyReactions(cTot, c);

		// 3. Correct the effective kinetic constants by three-body coefficients
		for(unsigned int s=0;s<number_of_thirdbody_reactions_;s++)
		{
			const unsigned int j=indices_of_thirdbody_reactions__[s];
			kArrheniusModified__[j-1] *= Meff__[s];
		}

		// 4. Correct the effective kinetic constants: Fall-off reactions
		FallOffReactions(cTot, c);
		for(unsigned int s=0;s<number_of_falloff_reactions_;s++)
		{
			const unsigned int j=indices_of_falloff_reactions__[s];
			kArrheniusModified__[j-1] *= correction_falloff__[s];
		}

		// 4.bis Extended fallExtendedFallOffReactionsoff reactions
		ExtendedFallOffReactions(cTot, c);

		// 5. Correct the effective kinetic constants: CABR reactions
		ChemicallyActivatedBimolecularReactions(cTot, c);
		for(unsigned int s=0;s<number_of_cabr_reactions_;s++)
		{
			const unsigned int j=indices_of_cabr_reactions__[s];
			kArrheniusModified__[j-1] *= correction_cabr__[s];
		}

		// Calculates the product of conenctrations (for forward and reverse reactions)
		// Be careful: the reverseReactionRates_ vector is defined for all the reactions
		// in the kinetic scheme, not only for the reversible reactions. After calling the
		// function reported below the value of reverseReactionRates_ vector for non reversible
		// reactions is put equal to 1.
		stoichiometry_->ProductOfConcentrations(forwardReactionRates__, reverseReactionRates__, c);

		// Corrects the product of concentrations for reverse reaction by the 
		// thermodynamic equilibrium constant
		for(unsigned int k=0;k<number_of_thermodynamic_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_thermodynamic_reversible_reactions__[k];
			reverseReactionRates__[j-1] *= uKeq__[k];
		}

		// Corrects the product of concentrations for reverse reaction by the 
		// explicit Arrhenius kinetic parameters (if provided)
		for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions__[k];
			reverseReactionRates__[j-1] *= kArrhenius_reversible__[k]/kArrhenius__[j-1];
		}

		// Calculates the net reaction rate
		// Be careful: the netReactionRates_ vector must be multiplied by the effective 
		// forward kinetic constant, to obtain the real reaction rates in [kmol/m3/s]
		netReactionRates__ = forwardReactionRates__;
		for(unsigned int k=0;k<number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions__[k]-1;
			netReactionRates__[j] -= reverseReactionRates__[j];
		}

		// Multiplies the net reaction rate by the effective kinetic constant (accounting for 
		// third-body effects, fall-off, etc.). At the end of this function the netReactionRates_
		// vector contains the net reaction rates of all the reactions in [kmol/m3/s]
		// TODO ElementByElementProduct(netReactionRates_, kArrheniusModified__, &netReactionRates_); 
		for (unsigned int i = 0; i < netReactionRates__.size(); i++)
			netReactionRates__[i] *= kArrheniusModified__[i];
	}

	void KineticsMap_CHEMKIN::DerivativesOfReactionRatesWithRespectToKineticParameters(const PhysicalConstants::sensitivity_type type, unsigned int jReaction, const double* c, double& parameter)
	{
		// Total concentration
		double cTot = 0.;
		for (unsigned int i = 0; i<this->number_of_species_; i++)
			cTot += c[i];

		// 1. Kinetic constants
		KineticConstants();

		// 1.bis Extended pressure log reactions
		ExtendedPressureLogReactions(cTot, c);
			
		// 2. Calculates the three-body corrections
		ThirdBodyReactions(cTot, c);

		// 3. Correct the effective kinetic constants by three-body coefficients
		for(unsigned int s=0;s<number_of_thirdbody_reactions_;s++)
		{
			const unsigned int j=indices_of_thirdbody_reactions__[s];
			kArrheniusModified__[j-1] *= Meff__[s];
		}

		// 4. Correct the effective kinetic constants: Fall-off reactions
		FallOffReactions(cTot, c);
		for(unsigned int s=0;s<number_of_falloff_reactions_;s++)
		{
			const unsigned int j=indices_of_falloff_reactions__[s];
			kArrheniusModified__[j-1] *= correction_falloff__[s];
		}

		// 4.bis Extended fallExtendedFallOffReactionsoff reactions
		ExtendedFallOffReactions(cTot, c);

		// 5. Correct the effective kinetic constants: CABR reactions
		ChemicallyActivatedBimolecularReactions(cTot, c);
		for(unsigned int s=0;s<number_of_cabr_reactions_;s++)
		{
			const unsigned int j=indices_of_cabr_reactions__[s];
			kArrheniusModified__[j-1] *= correction_cabr__[s];
		}

		if (type == PhysicalConstants::SENSITIVITY_KINETIC_CONSTANT)
		{
			parameter = kArrheniusModified__[jReaction-1];
			std::fill(kArrheniusModified__.begin(), kArrheniusModified__.end(), 0.);
			kArrheniusModified__[jReaction-1] = 1.;
		}
		else if (type == PhysicalConstants::SENSITIVITY_FREQUENCY_FACTOR)
		{
			if (jReaction <= this->number_of_reactions_)
			{
				const double tmp = kArrheniusModified__[jReaction-1];
				std::fill(kArrheniusModified__.begin(), kArrheniusModified__.end(), 0.);
				kArrheniusModified__[jReaction-1] = tmp;
				parameter =  sign_lnA__[jReaction-1]*std::exp(lnA__[jReaction-1]);

				if (type_of_reaction__[jReaction-1] == PhysicalConstants::REACTION_LINDEMANN_FALLOFF || 
					type_of_reaction__[jReaction-1] == PhysicalConstants::REACTION_TROE_FALLOFF ||
					type_of_reaction__[jReaction-1] == PhysicalConstants::REACTION_SRI_FALLOFF)
				{
					unsigned int jFallOff = local_family_index__[jReaction-1];
					double	kInf = std::exp(lnA_falloff_inf__[jFallOff-1] + Beta_falloff_inf__[jFallOff-1]*this->logT_ - E_over_R_falloff_inf__[jFallOff-1]*this->uT_);
					double F, dF_over_dA0, dF_over_dAInf;
					this->FallOffReactions(jFallOff, cTot, c, F, dF_over_dA0, dF_over_dAInf);

					kArrheniusModified__[jReaction-1] = kArrheniusModified__[jReaction-1] / std::exp(lnA__[jReaction-1]) * (1.-kArrheniusModified__[jReaction-1]/kInf/F) +
														kArrheniusModified__[jReaction-1]/F*dF_over_dA0;			
				}
				else if (type_of_reaction__[jReaction-1] == PhysicalConstants::REACTION_LINDEMANN_CABR || 
						 type_of_reaction__[jReaction-1] == PhysicalConstants::REACTION_TROE_CABR ||
						 type_of_reaction__[jReaction-1] == PhysicalConstants::REACTION_SRI_CABR)
				{
					unsigned int jCABR = local_family_index__[jReaction-1];
					double	k0   = std::exp(lnA__[jReaction-1] + Beta__[jReaction-1]*this->logT_ - E_over_R__[jReaction-1]*this->uT_);

					double F, dF_over_dA0, dF_over_dAInf;
					ChemicallyActivatedBimolecularReactions(jCABR, cTot, c, F, dF_over_dA0, dF_over_dAInf);
					kArrheniusModified__[jReaction-1] = kArrheniusModified__[jReaction-1]*kArrheniusModified__[jReaction-1] / std::exp(lnA__[jReaction-1])/k0/F +
														kArrheniusModified__[jReaction-1]/F*dF_over_dA0;
				}
				else if (type_of_reaction__[jReaction-1] == PhysicalConstants::REACTION_CHEBYSHEV)
				{
					parameter = kArrheniusModified__[jReaction-1];
					kArrheniusModified__[jReaction-1] = 1.;
				}
				else
				{
					kArrheniusModified__[jReaction-1] = kArrheniusModified__[jReaction-1] / (sign_lnA__[jReaction-1]*std::exp(lnA__[jReaction-1]));
				}
			}

			else 
			{
				// Fall-off reaction
				if (jReaction <= this->number_of_reactions_ + number_of_falloff_reactions_)
				{
					unsigned int jFallOff = jReaction - this->number_of_reactions_;
					unsigned int index_global = indices_of_falloff_reactions__[jFallOff-1];
					const double tmp = kArrheniusModified__[index_global-1];
					std::fill(kArrheniusModified__.begin(), kArrheniusModified__.end(), 0.);
					kArrheniusModified__[index_global-1] = tmp;
					parameter =  std::exp(lnA_falloff_inf__[jFallOff-1]);

					double	kInf = std::exp(lnA_falloff_inf__[jFallOff-1] + Beta_falloff_inf__[jFallOff-1]*this->logT_ - E_over_R_falloff_inf__[jFallOff-1]*this->uT_);
					
					double F = 0.;
					double dF_over_dA0 = 0.; 
					double dF_over_dAInf = 0.;
					FallOffReactions(jFallOff, cTot, c, F, dF_over_dA0, dF_over_dAInf);

					kArrheniusModified__[index_global-1] = kArrheniusModified__[index_global-1] / F * 
														 ( kArrheniusModified__[index_global-1]/std::exp(lnA_falloff_inf__[jFallOff-1])/kInf + dF_over_dAInf);	

					if (parameter == 0.) std::cout << "parameter " << jReaction - this->number_of_reactions_ << std::endl;
					if (F == 0.) std::cout << "Falloff inv F " << jReaction - this->number_of_reactions_ << std::endl;
					if (kInf == 0.) std::cout << "Falloff inv kInf " << jReaction - this->number_of_reactions_ << std::endl;
					if (std::exp(lnA_falloff_inf__[jFallOff-1]) == 0.) std::cout << "Falloff inv exp " << jReaction - this->number_of_reactions_ << std::endl;

				}
				// CABR reactions
				else if (jReaction <= this->number_of_reactions_ + number_of_falloff_reactions_ +  number_of_cabr_reactions_)
				{
					unsigned int jCABR = jReaction - this->number_of_reactions_ - number_of_falloff_reactions_;
					unsigned int index_global = indices_of_cabr_reactions__[jCABR-1];
					const double tmp = kArrheniusModified__[index_global-1];
					std::fill(kArrheniusModified__.begin(), kArrheniusModified__.end(), 0.);
					kArrheniusModified__[index_global-1] = tmp;
					parameter =  std::exp(lnA_falloff_inf__[jCABR-1]);

					double	k0   = std::exp(lnA__[index_global-1] + Beta__[index_global-1]*this->logT_ - E_over_R__[index_global-1]*this->uT_);

					double F = 0.;
					double dF_over_dA0 = 0.;
					double dF_over_dAInf = 0.;
					ChemicallyActivatedBimolecularReactions(jCABR, cTot, c, F, dF_over_dA0, dF_over_dAInf);
					
					kArrheniusModified__[index_global-1] = kArrheniusModified__[index_global-1] / F * 
														 ( kArrheniusModified__[index_global-1]/std::exp(lnA_cabr_inf__[jCABR-1])/k0 + dF_over_dAInf);	
				}
			}
		}


		// Calculates the product of conenctrations (for forward and reverse reactions)
		// Be careful: the reverseReactionRates_ vector is defined for all the reactions
		// in the kinetic scheme, not only for the reversible reactions. After calling the
		// function reported below the value of reverseReactionRates_ vector for non reversible
		// reactions is put equal to 1.
		stoichiometry_->ProductOfConcentrations(forwardReactionRates__, reverseReactionRates__, c);
		
		// Corrects the product of concentrations for reverse reaction by the 
		// thermodynamic equilibrium constant
		for(unsigned int k=0;k<number_of_thermodynamic_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_thermodynamic_reversible_reactions__[k];
			reverseReactionRates__[j-1] *= uKeq__[k];
		}

		for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions__[k];
			if (j==jReaction)
				reverseReactionRates__[j-1] *= kArrhenius_reversible__[k]/kArrhenius__[j-1];
			else
				reverseReactionRates__[j-1] = 0.;
		}

		// Calculates the net reaction rate
		// Be careful: the netReactionRates_ vector must be multiplied by the effective 
		// forward kinetic constant, to obtain the real reaction rates in [kmol/m3/s]
		netReactionRates__ = forwardReactionRates__;
		for(unsigned int k=0;k<number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions__[k]-1;
			netReactionRates__[j] -= reverseReactionRates__[j];
		}

		// Multiplies the net reaction rate by the effective kinetic constant (accounting for 
		// third-body effects, fall-off, etc.). At the end of this function the netReactionRates_
		// vector contains the net reaction rates of all the reactions in [kmol/m3/s]
		// TODO ElementByElementProduct(netReactionRates_, kArrheniusModified__, &netReactionRates_);
		for (unsigned int i = 0; i < netReactionRates__.size(); i++)
			netReactionRates__[i] *= kArrheniusModified__[i];
	}

	void KineticsMap_CHEMKIN::FormationRates(double* R)
	{
		stoichiometry_->FormationRatesFromReactionRates(R, netReactionRates__.data());
	}

	double KineticsMap_CHEMKIN::HeatRelease(const double* R)
	{
		return -Dot(this->number_of_species_, R, thermodynamics_.Species_H_over_RT().data()) * PhysicalConstants::R_J_kmol * this->T_;
	}

	void KineticsMap_CHEMKIN::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const unsigned int k, const double* c, double* Jalfa, double& parameter)
	{
		DerivativesOfReactionRatesWithRespectToKineticParameters(type, k, c, parameter);
		stoichiometry_->FormationRatesFromReactionRates(Jalfa, netReactionRates__.data());
	}

	void KineticsMap_CHEMKIN::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const double* c, const double* mole_fractions, double* Jalfa, double& JT, double& parameter)
	{
		DerivativesOfReactionRatesWithRespectToKineticParameters(type, k, c, parameter);
		stoichiometry_->FormationRatesFromReactionRates(Jalfa, netReactionRates__.data());
		
		JT  = HeatRelease(Jalfa);
		if (energy_type == CONSTANT_VOLUME_SYSTEM)
		{
			double sumMoleFormationRates = 0.;
			for(unsigned int i=0;i<this->number_of_species_;i++)
				sumMoleFormationRates += Jalfa[i];

			JT  += PhysicalConstants::R_J_kmol*this->T_*sumMoleFormationRates;
		}	
		if (energy_type == PLUGFLOW_SYSTEM)
		{
			// TODO
			//double sumMoleFormationRates = (*Jalfa).SumElements();
		}	
	}
 
	void KineticsMap_CHEMKIN::SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const double* c, const double* mole_fractions, double* Jalfa, double& JT, double& Jrho, double& parameter)
	{
		DerivativesOfReactionRatesWithRespectToKineticParameters(type, k, c, parameter);
		stoichiometry_->FormationRatesFromReactionRates(Jalfa, netReactionRates__.data());
		
		const double CpMixMolar = thermodynamics_.cpMolar_Mixture_From_MoleFractions(mole_fractions);
		
		const double dQR_over_parameter = HeatRelease(Jalfa);
		
		double sumMoleFormationRates = 0.;
		for (unsigned int i = 0; i<this->number_of_species_; i++)
			sumMoleFormationRates += Jalfa[i];

		JT    =  dQR_over_parameter;
		Jrho  = - (dQR_over_parameter/CpMixMolar + this->T_*sumMoleFormationRates);
	}


	// This version seems to be wrong
	void KineticsMap_CHEMKIN::ProductionAndDestructionRatesGross(double* P, double* D)
	{
		std::vector<double> forward_(this->number_of_reactions_);
		std::vector<double> backward_(this->number_of_reactions_);

		GetForwardReactionRates(forward_.data());
		GetBackwardReactionRates(backward_.data());

		stoichiometry_->ProductionAndDestructionRatesFromReactionRatesGross(P, D, forward_.data(), backward_.data());
	}
	
	void KineticsMap_CHEMKIN::ProductionAndDestructionRates(double* P, double* D)
	{
		stoichiometry_->ProductionAndDestructionRatesFromReactionRates(P, D, netReactionRates__.data());
	}
 
	void KineticsMap_CHEMKIN::RateOfProductionAnalysis(const bool iNormalize) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates__.data(), iNormalize);
	}

	void KineticsMap_CHEMKIN::RateOfProductionAnalysis(std::ostream& fout) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates__.data(), false);
		stoichiometry_->WriteRateOfProductionAnalysis(fout);
	}

	void KineticsMap_CHEMKIN::RateOfProductionAnalysis(ROPA_Data& ropa) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates__.data(), false);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}

	void KineticsMap_CHEMKIN::RateOfProductionAnalysis(ROPA_Data& ropa, const double* rf, const double* rb) const
	{
		stoichiometry_->RateOfProductionAnalysis(rf, rb);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}

	const std::vector<double>& KineticsMap_CHEMKIN::GiveMeReactionRates()
	{
		return netReactionRates__;
	}

	void KineticsMap_CHEMKIN::GiveMeReactionRates(double* r)
	{
		for(unsigned int i=0;i<netReactionRates__.size();i++)
			r[i] = netReactionRates__[i];
	}

	void KineticsMap_CHEMKIN::GetForwardReactionRates(double* r)
	{
		// TODO ElementByElementProduct(forwardReactionRates_, kArrheniusModified__, r);

		for (unsigned int i = 0; i < kArrheniusModified__.size(); i++)
			r[i] = kArrheniusModified__[i]*forwardReactionRates__[i];
	}

	void KineticsMap_CHEMKIN::GetBackwardReactionRates(double* r)
	{
		*r = 0.;
		for(unsigned int k=0;k<number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions__[k]-1;
			r[j] = reverseReactionRates__[j]*kArrheniusModified__[j];
		}
		for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions__[k]-1;
			r[j] = reverseReactionRates__[j]*kArrheniusModified__[j];
		}                
	}

	void KineticsMap_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k)
	{			
		thermodynamics_.SetPressure(101325.);
		thermodynamics_.SetTemperature(298.15);
		SetTemperature(298.15);
		SetPressure(101325.);
		ReactionEnthalpiesAndEntropies();
		
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << (reaction_h_over_RT__[k-1]-reaction_s_over_R__[k-1])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_h_over_RT__[k-1] *PhysicalConstants::R_kcal_mol*this->T_;;							// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_s_over_R__[k-1] *PhysicalConstants::R_cal_mol;;									// [cal/mol/K]
	}

	void KineticsMap_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k, const double* c_bath, const std::vector<double> list_of_temperatures, const double conversion_forward, const double conversion_backward)
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
			fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified__[k-1];

			// Equilibrium and Backward kinetic constant [kmol, m3, s]
			if (isThermodynamicallyReversible__[k-1] != 0)
			{
				const unsigned int j = isThermodynamicallyReversible__[k-1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << 1./uKeq__[j-1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified__[k-1]*uKeq__[j-1];
			}
			else if (isExplicitlyReversible__[k-1] != 0)
			{
				const unsigned int j = isExplicitlyReversible__[k-1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified__[k-1]/kArrhenius_reversible__[j-1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrhenius_reversible__[j-1];
			}
			else
			{
				fOut << std::setw(14) << std::left << std::setprecision(0) << std::fixed << 0.;
				fOut << std::setw(14) << std::left << std::setprecision(0) << std::fixed << 0.;
			}

			// Thermodynamic data
			fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << (reaction_h_over_RT__[k-1]-reaction_s_over_R__[k-1])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
			fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_h_over_RT__[k-1] *PhysicalConstants::R_kcal_mol*this->T_;;						// [kcal/mol]
			fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_s_over_R__[k-1] *PhysicalConstants::R_cal_mol;;							// [cal/mol/K]

			// Forward kinetic constant [mol, cm3, s]
			fOut    << std::setw(14) << std::left << std::setprecision(4) << std::scientific << kArrheniusModified__[k-1]*conversion_forward;

			// Equilibrium and Backward kinetic constant [mol, cm3, s]
			if (isThermodynamicallyReversible__[k-1] != 0)
			{
				const unsigned int j = isThermodynamicallyReversible__[k-1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrheniusModified__[k-1]*uKeq__[j-1]*conversion_backward;
			}
			else if (isExplicitlyReversible__[k-1] != 0)
			{
				const unsigned int j = isExplicitlyReversible__[k-1];
				fOut << std::setw(14) << std::left << std::setprecision(3) << std::scientific << kArrhenius_reversible__[j-1]*conversion_backward;
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

	void KineticsMap_CHEMKIN::WriteCollisionRateConstantForBimolecularReactions(
		TransportPropertiesMap_CHEMKIN& transport, 
		std::ostream& fOut, 
		const std::vector<std::string>& reaction_names,
		const std::vector<double>& list_of_temperatures)
	{
		// Find bimolecular direct reactions
		Eigen::VectorXd sum_forward;
		Eigen::VectorXd sum_backward;
		stoichiometry_->GetSumOfStoichiometricCoefficientsOfReactants(sum_forward);
		stoichiometry_->GetSumOfStoichiometricCoefficientsOfProducts(sum_backward);

		// Select only forward bimolecular reactions
		std::vector<unsigned int> indices_forward;
		for (int i = 0; i < sum_forward.size(); i++)
			if (sum_forward[i] == 2.)
			{
				std::vector<unsigned int>::iterator location;
				location = std::find(indices_of_thirdbody_reactions__.begin(), indices_of_thirdbody_reactions__.end(), i+1);
				if (location == indices_of_thirdbody_reactions__.end())
					indices_forward.push_back(i);
			}

		// Select only backward bimolecular reactions if reversible
		std::vector<unsigned int> indices_backward;
		for (int i = 0; i < sum_backward.size(); i++)
			if (sum_backward[i] == 2.)
			{
				if (isExplicitlyReversible__[i] != 0)
				{
					std::vector<unsigned int>::iterator location;
					location = std::find(indices_of_thirdbody_reactions__.begin(), indices_of_thirdbody_reactions__.end(), i + 1);
					if (location == indices_of_thirdbody_reactions__.end())
						indices_backward.push_back(i);
				}
				//if (isThermodynamicallyReversible__[i] != 0 && stoichiometry_->ChangeOfMoles()[i]==0.)	// original version (Oct 2017)
				if (isThermodynamicallyReversible__[i] != 0)	// Hai Wang version (private communication, Nov 2017)
				{
					std::vector<unsigned int>::iterator location;
					location = std::find(indices_of_thirdbody_reactions__.begin(), indices_of_thirdbody_reactions__.end(), i + 1);
					if (location == indices_of_thirdbody_reactions__.end())
						indices_backward.push_back(i);
				}
			}

		std::vector<std::vector<unsigned int>> indices_species_forward(indices_forward.size());
		for (int k = 0; k<stoichiometry_->stoichiometric_matrix_reactants().outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometry_->stoichiometric_matrix_reactants(), k); it; ++it)
			{
				std::vector<unsigned int>::iterator location;
				location = std::find(indices_forward.begin(), indices_forward.end(), it.row());

				if (location != indices_forward.end())
					indices_species_forward[location - indices_forward.begin()].push_back(it.col());
			}
		}

		std::vector<std::vector<unsigned int>> indices_species_backward(indices_backward.size());
		for (int k = 0; k<stoichiometry_->stoichiometric_matrix_products().outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometry_->stoichiometric_matrix_products(), k); it; ++it)
			{
				std::vector<unsigned int>::iterator location;
				location = std::find(indices_backward.begin(), indices_backward.end(), it.row());

				if (location != indices_backward.end())
					indices_species_backward[location - indices_backward.begin()].push_back(it.col());
			}
		}

		// Check and adjust (forward)
		for (unsigned int i = 0; i < indices_forward.size(); i++)
		{
			if (indices_species_forward[i].size() != 2)
			{
				if (indices_species_forward[i].size() == 1)
					indices_species_forward[i].push_back(indices_species_forward[i][0]);
				else
				{
					std::cout << "Rate of collision analysis for bimolecular reactions" << std::endl;
					std::cout << "Analyzing forward reaction #" << indices_forward[i] + 1 << std::endl;
					OpenSMOKE::FatalErrorMessage("The reaction was expected to be bimolecular");
				}
			}
		}

		// Check and adjust (backward)
		for (unsigned int i = 0; i < indices_backward.size(); i++)
		{
			if (indices_species_backward[i].size() != 2)
			{
				if (indices_species_backward[i].size() == 1)
					indices_species_backward[i].push_back(indices_species_backward[i][0]);
				else
				{
					std::cout << "Rate of collision analysis for bimolecular reactions" << std::endl;
					std::cout << "Analyzing backward reaction #" << indices_backward[i] + 1 << std::endl;
					OpenSMOKE::FatalErrorMessage("The reaction was expected to be bimolecular");
				}
			}
		}

		const double P = 1e12;	// [Pa] super-high pressure for fall-off reactions
		Eigen::VectorXd c(thermodynamics_.NumberOfSpecies());
		this->SetPressure(P);
		thermodynamics_.SetPressure(P);

		unsigned int ni = 23;
		std::vector<double> temperatures(ni);
		temperatures[0] = 300.;
		for (unsigned int k = 1; k < ni; k++)
			temperatures[k] = temperatures[k-1] + 100.;

		if (list_of_temperatures.size() != 0)
		{
			ni = list_of_temperatures.size();
			temperatures = list_of_temperatures;
		}

		Eigen::MatrixXd kForward(indices_forward.size(), ni);
		Eigen::MatrixXd kCollisionForward(indices_forward.size(), ni);
		Eigen::MatrixXd kStarForward(indices_forward.size(),ni);
		Eigen::MatrixXd kBackward(indices_backward.size(), ni);
		Eigen::MatrixXd kCollisionBackward(indices_backward.size(),ni);
		Eigen::MatrixXd kStarBackward(indices_backward.size(),ni);
		for (unsigned int k = 0; k < ni; k++)
		{
			const double T = temperatures[k];
			const double cTot = P / PhysicalConstants::R_J_kmol / T;
			c.setConstant(cTot / double(thermodynamics_.NumberOfSpecies()));

			this->SetTemperature(T);
			thermodynamics_.SetTemperature(T);

			this->ReactionEnthalpiesAndEntropies();
			this->KineticConstants();
			this->ReactionRates(c.data());

			// Forward kinetic constants
			for (unsigned int i = 0; i < indices_forward.size(); i++)
			{
				kForward(i, k) = kArrheniusModified__[indices_forward[i]];
				kCollisionForward(i, k) = transport.kCollision(indices_species_forward[i][0], indices_species_forward[i][1], T);
				kStarForward(i, k) = kForward(i, k) / kCollisionForward(i, k);
			}

			// Backward kinetic constants
			for (unsigned int i = 0; i < indices_backward.size(); i++)
			{
				kCollisionBackward(i, k) = transport.kCollision(indices_species_backward[i][0], indices_species_backward[i][1], T);

				if (isThermodynamicallyReversible__[indices_backward[i]] != 0)
				{
					const unsigned int j = isThermodynamicallyReversible__[indices_backward[i]] - 1;
					kBackward(i, k) = (kArrheniusModified__[indices_backward[i]] * uKeq__[j]);
					kStarBackward(i, k) = kBackward(i, k) / kCollisionBackward(i, k);
				}
				else if (isExplicitlyReversible__[indices_backward[i]] != 0)
				{
					const unsigned int j = isExplicitlyReversible__[indices_backward[i]] - 1;
					kBackward(i, k) = kArrhenius_reversible__[j];
					kStarBackward(i, k) = kBackward(i, k) / kCollisionBackward(i, k);
				}
			}
		}

		// General data
		{
			fOut << "----------------------------------------------------------------" << std::endl;
			fOut << " General data                                                   " << std::endl;
			fOut << "----------------------------------------------------------------" << std::endl;
			fOut << "Bimolecular reactions (total):    " << indices_forward.size() + indices_backward.size() << std::endl;
			fOut << "Bimolecular reactions (forward):  " << indices_forward.size() << std::endl;
			fOut << "Bimolecular reactions (backward): " << indices_backward.size() << std::endl;
			fOut << "----------------------------------------------------------------" << std::endl;
			fOut << std::endl;
		}

		// List of unfeasible reactions
		Eigen::VectorXi unfeasible_forward_reactions(indices_forward.size());
		Eigen::VectorXi unfeasible_backward_reactions(indices_backward.size());
		{
			unfeasible_forward_reactions.setZero();
			for (unsigned int i = 0; i < indices_forward.size(); i++)
				for (unsigned int k = 0; k < ni; k++)
					if (kStarForward(i, k) >= 1.)	unfeasible_forward_reactions(i)++;

			unfeasible_backward_reactions.setZero();
			for (unsigned int i = 0; i < indices_backward.size(); i++)
				for (unsigned int k = 0; k < ni; k++)
					if (kStarBackward(i, k) >= 1.)	unfeasible_backward_reactions(i)++;

			fOut << "-----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(11) << std::left << "Forward";
			fOut << std::setw(11) << std::left << "Violations";
			fOut << std::setw(11) << std::left << "Max";
			fOut << std::setw(80) << std::left << "Reaction";
			fOut << std::setw(11) << std::left << "T>=300K";
			fOut << std::setw(11) << std::left << "T>=500K";
			fOut << std::setw(11) << std::left << "T>=700K";
			fOut << std::setw(11) << std::left << "T>=1000K";
			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			for (unsigned int i = 0; i < indices_forward.size(); i++)
				if (unfeasible_forward_reactions(i) > 0)
				{
					fOut << std::setw(11) << std::left << std::fixed << indices_forward[i] + 1;
					fOut << std::setw(11) << std::left << std::fixed << unfeasible_forward_reactions(i);
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarForward.row(i).maxCoeff();
					fOut << std::setw(80) << std::left << reaction_names[indices_forward[i]];

					double kStar300  = 0.;
					double kStar500 = 0.;
					double kStar700 = 0.;
					double kStar1000 = 0.;
					for (unsigned int k = 0; k < ni; k++)
					{
						if (temperatures[k] >= 300.)
							if (kStarForward(i, k) > kStar300)	kStar300 = kStarForward(i, k);
						if (temperatures[k] >= 500.)
							if (kStarForward(i, k) > kStar500)	kStar500 = kStarForward(i, k);
						if (temperatures[k] >= 700.)
							if (kStarForward(i, k) > kStar700)	kStar700 = kStarForward(i, k);
						if (temperatures[k] >= 1000.)
							if (kStarForward(i, k) > kStar1000)	kStar1000 = kStarForward(i, k);
					}
					
					if (kStar300>1.)	fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStar300;
					else                fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << "*";
					if (kStar500>1.)	fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStar500;
					else                fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << "*";
					if (kStar700>1.)	fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStar700;
					else                fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << "*";
					if (kStar1000>1.)	fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStar1000;
					else                fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << "*";

					fOut << std::endl;
				}
			fOut << "-----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			fOut << std::endl;

			fOut << "-----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			fOut << std::setw(11) << std::left << "Backward";
			fOut << std::setw(11) << std::left << "Violations";
			fOut << std::setw(11) << std::left << "Max";
			fOut << std::setw(80) << std::left << "Reaction";
			fOut << std::setw(11) << std::left << "T>=300K";
			fOut << std::setw(11) << std::left << "T>=500K";
			fOut << std::setw(11) << std::left << "T>=700K";
			fOut << std::setw(11) << std::left << "T>=1000K";
			fOut << std::endl;
			fOut << "-----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			for (unsigned int i = 0; i < indices_backward.size(); i++)
				if (unfeasible_backward_reactions(i) > 0)
				{
					fOut << std::setw(11) << std::left << std::fixed << indices_backward[i] + 1;
					fOut << std::setw(11) << std::left << std::fixed << unfeasible_backward_reactions(i);
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarBackward.row(i).maxCoeff();
					fOut << std::setw(80) << std::left << reaction_names[indices_backward[i]];

					double kStar300 = 0.;
					double kStar500 = 0.;
					double kStar700 = 0.;
					double kStar1000 = 0.;
					for (unsigned int k = 0; k < ni; k++)
					{
						if (temperatures[k] >= 300.)
							if (kStarForward(i, k) > kStar300)	kStar300 = kStarBackward(i, k);
						if (temperatures[k] >= 500.)
							if (kStarForward(i, k) > kStar500)	kStar500 = kStarBackward(i, k);
						if (temperatures[k] >= 700.)
							if (kStarForward(i, k) > kStar700)	kStar700 = kStarBackward(i, k);
						if (temperatures[k] >= 1000.)
							if (kStarForward(i, k) > kStar1000)	kStar1000 = kStarBackward(i, k);
					}

					if (kStar300>1.)	fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStar300;
					else                fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << "*";
					if (kStar500>1.)	fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStar500;
					else                fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << "*";
					if (kStar700>1.)	fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStar700;
					else                fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << "*";
					if (kStar1000>1.)	fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStar1000;
					else                fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << "*";

					fOut << std::endl;
				}
			fOut << "-----------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			fOut << std::endl;
		}

		// Analysis of violations
		{
			Eigen::VectorXi violations_forward(ni);
			violations_forward.setZero();
			for (unsigned int i = 0; i < indices_forward.size(); i++)
				for (unsigned int k = 0; k < ni; k++)
					if (kStarForward(i, k) >= 1.)	violations_forward(k)++;

			Eigen::VectorXi violations_backward(ni);
			violations_backward.setZero();
			for (unsigned int i = 0; i < indices_backward.size(); i++)
				for (unsigned int k = 0; k < ni; k++)
					if (kStarBackward(i, k) >= 1.)	violations_backward(k)++;

			fOut << std::string(275, '-') << std::endl;
			fOut << std::setw(11) << std::left << "Violations";
			fOut << std::setw(11) << std::left << "Total";
			for (unsigned int k = 0; k < ni; k++)
				fOut << std::setw(11) << std::left << std::setprecision(0) << std::fixed << temperatures[k];
			fOut << std::endl;
			fOut << std::string(275, '-') << std::endl;

			fOut << std::setw(11) << std::left << "Forward";
			fOut << std::setw(11) << std::left << std::fixed << violations_forward.sum();
			for (unsigned int k = 0; k < ni; k++)
				fOut << std::setw(11) << std::left << std::fixed << violations_forward(k);
			fOut << std::endl;

			fOut << std::setw(11) << std::left << "Backward";
			fOut << std::setw(11) << std::left << std::fixed << violations_backward.sum();
			for (unsigned int k = 0; k < ni; k++)
				fOut << std::setw(11) << std::left << std::fixed << violations_backward(k);
			fOut << std::endl;

			fOut << std::setw(11) << std::left << "Total";
			fOut << std::setw(11) << std::left << std::fixed << violations_forward.sum() + violations_backward.sum();
			for (unsigned int k = 0; k < ni; k++)
				fOut << std::setw(11) << std::left << std::fixed << violations_forward(k) + violations_backward(k);
			fOut << std::endl;

			fOut << std::endl;
		}

		// Details about unfeasible reactions
		{
			for (unsigned int i = 0; i < indices_forward.size(); i++)
				if (unfeasible_forward_reactions(i) > 0)
				{
					fOut << std::string(100, '-') << std::endl;
					fOut << "Forward reaction: " << indices_forward[i] + 1 << std::endl;
					fOut << std::string(100, '-') << std::endl;
					fOut << "Stoichiometry:    " << reaction_names[indices_forward[i]] << std::endl;
					fOut << std::endl;

					const unsigned int i1 = indices_species_forward[i][0];
					const unsigned int i2 = indices_species_forward[i][1];

					fOut << std::setw(20) << std::left << "Species";
					fOut << std::setw(11) << std::left << "mass[kg]";
					fOut << std::setw(11) << std::left << "sigma[A]";
					fOut << std::setw(11) << std::left << "eps/kb[K]";
					fOut << std::endl;
					fOut << std::string(55, '-') << std::endl;

					fOut << std::setw(20) << std::left << std::fixed << thermodynamics_.NamesOfSpecies()[i1];
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.mu()[i1];
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.sigma()[i1]*1.e10;
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.epsilon_over_kb()[i1];
					fOut << std::endl;

					fOut << std::setw(20) << std::left << std::fixed << thermodynamics_.NamesOfSpecies()[i2];
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.mu()[i2];
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.sigma()[i2]*1.e10;
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.epsilon_over_kb()[i2];
					fOut << std::endl;

					fOut << std::setw(20) << std::left << std::fixed << "Mean";
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << (transport.mu()[i1] * transport.mu()[i2])/(transport.mu()[i1] + transport.mu()[i2]);
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << 0.5*(transport.sigma()[i1]+transport.sigma()[i2])*1.e10;
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << std::sqrt(transport.epsilon_over_kb()[i1]*transport.epsilon_over_kb()[i2]);
					fOut << std::endl;
					fOut << std::string(55, '-') << std::endl;
					fOut << std::endl;

					fOut << std::setw(10) << std::left << std::fixed << "T[K]";
					fOut << std::setw(11) << std::left << std::fixed << "TStar[-]";
					fOut << std::setw(13) << std::left << std::fixed << "Omega11*[-]";
					fOut << std::setw(15) << std::left << std::fixed << "k[m3/kmol/s]";
					fOut << std::setw(16) << std::left << std::fixed << "kCol[m3/kmol/s]";
					fOut << std::setw(11) << std::left << std::fixed << "k*[-]";
					fOut << std::endl;
					fOut << std::string(75, '-') << std::endl;

					const double eps_over_kb = std::sqrt(transport.epsilon_over_kb()[i1] * transport.epsilon_over_kb()[i2]);
					for (unsigned int j = 0; j < ni; j++)
					{
						const double T = temperatures[j];
						const double TStar = T / eps_over_kb;
						const double Omega11Star = transport.Omega11(TStar);

						fOut << std::setw(10) << std::left << std::setprecision(0) << std::fixed << T;
						fOut << std::setw(11) << std::left << std::setprecision(3) << std::scientific << TStar;
						fOut << std::setw(13) << std::left << std::setprecision(3) << std::scientific << Omega11Star;
						fOut << std::setw(15) << std::left << std::setprecision(3) << std::scientific << kForward(i, j);
						fOut << std::setw(16) << std::left << std::setprecision(3) << std::scientific << kCollisionForward(i, j);
						fOut << std::setw(11) << std::left << std::setprecision(3) << std::scientific << kForward(i, j)/kCollisionForward(i, j);
						fOut << std::endl;
					}
					fOut << std::string(75, '-') << std::endl;
					fOut << std::endl;
				}

			for (unsigned int i = 0; i < indices_backward.size(); i++)
				if (unfeasible_backward_reactions(i) > 0)
				{
					fOut << std::string(100, '-') << std::endl;
					fOut << "Backward reaction: " << indices_backward[i] + 1 << std::endl;
					fOut << std::string(100, '-') << std::endl;
					fOut << "Stoichiometry:    " << reaction_names[indices_backward[i]] << std::endl;
					fOut << std::endl;

					const unsigned int i1 = indices_species_backward[i][0];
					const unsigned int i2 = indices_species_backward[i][1];

					fOut << std::setw(20) << std::left << "Species";
					fOut << std::setw(11) << std::left << "mass[kg]";
					fOut << std::setw(11) << std::left << "sigma[A]";
					fOut << std::setw(11) << std::left << "eps/kb[K]";
					fOut << std::endl;
					fOut << std::string(55, '-') << std::endl;

					fOut << std::setw(20) << std::left << std::fixed << thermodynamics_.NamesOfSpecies()[i1];
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.mu()[i1];
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.sigma()[i1] * 1.e10;
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.epsilon_over_kb()[i1];
					fOut << std::endl;

					fOut << std::setw(20) << std::left << std::fixed << thermodynamics_.NamesOfSpecies()[i2];
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.mu()[i2];
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.sigma()[i2] * 1.e10;
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << transport.epsilon_over_kb()[i2];
					fOut << std::endl;

					fOut << std::setw(20) << std::left << std::fixed << "Mean";
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << (transport.mu()[i1] * transport.mu()[i2]) / (transport.mu()[i1] + transport.mu()[i2]);
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << 0.5*(transport.sigma()[i1] + transport.sigma()[i2])*1.e10;
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << std::sqrt(transport.epsilon_over_kb()[i1] * transport.epsilon_over_kb()[i2]);
					fOut << std::endl;
					fOut << std::string(55, '-') << std::endl;
					fOut << std::endl;

					fOut << std::setw(10) << std::left << std::fixed << "T[K]";
					fOut << std::setw(11) << std::left << std::fixed << "TStar[-]";
					fOut << std::setw(13) << std::left << std::fixed << "Omega11*[-]";
					fOut << std::setw(15) << std::left << std::fixed << "k[m3/kmol/s]";
					fOut << std::setw(16) << std::left << std::fixed << "kCol[m3/kmol/s]";
					fOut << std::setw(11) << std::left << std::fixed << "k*[-]";
					fOut << std::endl;
					fOut << std::string(75, '-') << std::endl;

					const double eps_over_kb = std::sqrt(transport.epsilon_over_kb()[i1] * transport.epsilon_over_kb()[i2]);
					for (unsigned int j = 0; j < ni; j++)
					{
						const double T = temperatures[j];
						const double TStar = T / eps_over_kb;
						const double Omega11Star = transport.Omega11(TStar);

						fOut << std::setw(10) << std::left << std::setprecision(0) << std::fixed << T;
						fOut << std::setw(11) << std::left << std::setprecision(3) << std::scientific << TStar;
						fOut << std::setw(13) << std::left << std::setprecision(3) << std::scientific << Omega11Star;
						fOut << std::setw(15) << std::left << std::setprecision(3) << std::scientific << kBackward(i, j);
						fOut << std::setw(16) << std::left << std::setprecision(3) << std::scientific << kCollisionBackward(i, j);
						fOut << std::setw(11) << std::left << std::setprecision(3) << std::scientific << kBackward(i, j) / kCollisionBackward(i, j);
						fOut << std::endl;
					}
					fOut << std::string(75, '-') << std::endl;
					fOut << std::endl;
				}
		}

		// Complete matrices
		{
			fOut << std::string(275, '-') << std::endl;
			fOut << std::setw(11) << std::left << "Reaction";
			fOut << std::setw(11) << std::left << "Max";
			for (unsigned int k = 0; k < ni; k++)
				fOut << std::setw(11) << std::left << std::setprecision(0) << std::fixed << temperatures[k];
			fOut << std::endl;
			fOut << std::string(275, '-') << std::endl;

			fOut << std::setw(11) << std::left << std::fixed << "Max";
			fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarForward.maxCoeff();
			for (unsigned int k = 0; k < ni; k++)
				fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarForward.col(k).maxCoeff();
			fOut << std::endl;
			for (unsigned int i = 0; i < indices_forward.size(); i++)
			{
				fOut << std::setw(11) << std::left << std::fixed << indices_forward[i] + 1;
				fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarForward.row(i).maxCoeff();
				for (unsigned int k = 0; k < ni; k++)
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarForward(i, k);
				fOut << std::endl;
			}
			fOut << std::string(275, '-') << std::endl;
			fOut << std::endl;

			fOut << std::string(275, '-') << std::endl;
			fOut << std::setw(11) << std::left << "Reaction";
			fOut << std::setw(11) << std::left << "Max";
			for (unsigned int k = 0; k < ni; k++)
				fOut << std::setw(11) << std::left << std::setprecision(0) << std::fixed << temperatures[k];
			fOut << std::endl;
			fOut << std::string(275, '-') << std::endl;

			fOut << std::setw(11) << std::left << std::fixed << "Max";
			fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarBackward.maxCoeff();
			for (unsigned int k = 0; k < ni; k++)
				fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarBackward.col(k).maxCoeff();
			fOut << std::endl;
			for (unsigned int i = 0; i < indices_backward.size(); i++)
			{
				fOut << std::setw(11) << std::left << std::fixed << indices_backward[i] + 1;
				fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarBackward.row(i).maxCoeff();
				for (unsigned int k = 0; k < ni; k++)
					fOut << std::setw(11) << std::left << std::setprecision(2) << std::scientific << kStarBackward(i, k);
				fOut << std::endl;
			}
			fOut << std::string(275, '-') << std::endl;
			fOut << std::endl;
		}
	}

	void KineticsMap_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k, const double T, const double P_Pa, const double* c)
	{
		double cTot = 0.;
		for (unsigned int i = 0; i < this->number_of_species_; i++)
			cTot += c[i];

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
		fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kArrheniusModified__[k-1];

		// Equilibrium and Backward kinetic constant [kmol, m3, s]
		if (isThermodynamicallyReversible__[k-1] != 0)
		{
			const unsigned int j = isThermodynamicallyReversible__[k-1];
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << 1. / uKeq__[j-1];
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kArrheniusModified__[k-1] * uKeq__[j-1];
		}
		else if (isExplicitlyReversible__[k-1] != 0)
		{
			const unsigned int j = isExplicitlyReversible__[k-1];
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kArrheniusModified__[k-1] / kArrhenius_reversible__[j-1];
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << kArrhenius_reversible__[j-1];
		}
		else
		{
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::fixed << 0.;
			fOut << std::setw(16) << std::left << std::setprecision(6) << std::fixed << 0.;
		}

		// Thermodynamic data
		fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << (reaction_h_over_RT__[k-1] - reaction_s_over_R__[k-1])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
		fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << reaction_h_over_RT__[k-1] * PhysicalConstants::R_kcal_mol*this->T_;;						// [kcal/mol]
		fOut << std::setw(16) << std::left << std::setprecision(6) << std::scientific << reaction_s_over_R__[k-1] * PhysicalConstants::R_kcal_mol;;							// [kcal/mol/K]


		if (type_of_reaction__[k-1] == PhysicalConstants::REACTION_LINDEMANN_FALLOFF ||
			type_of_reaction__[k-1] == PhysicalConstants::REACTION_TROE_FALLOFF ||
			type_of_reaction__[k-1] == PhysicalConstants::REACTION_SRI_FALLOFF)
		{
			unsigned int jFallOff = local_family_index__[k-1];
			double	kInf = std::exp(lnA_falloff_inf__[jFallOff-1] + Beta_falloff_inf__[jFallOff-1] * this->logT_ - E_over_R_falloff_inf__[jFallOff-1] * this->uT_);
			
			double F = 0.;
			double dF_over_dA0 = 0.;
			double dF_over_dAInf = 0.;
			this->FallOffReactions(jFallOff, cTot, c, F, dF_over_dA0, dF_over_dAInf);

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

	void KineticsMap_CHEMKIN::FittedReverseKineticConstants(const double* x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters, const bool only_reversible)
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

			std::vector<double> c(this->number_of_species_);
			for (unsigned int i = 1; i <= npoints; i++)
			{
				const double ctot = patm / PhysicalConstants::R_J_kmol / T(i - 1);
				Prod(this->number_of_species_, ctot, x_bath, c.data());

				SetTemperature(T(i - 1));
				thermodynamics_.SetTemperature(T(i - 1));
				ReactionEnthalpiesAndEntropies();
				KineticConstants();
				ReactionRates(c.data());

				for (unsigned int k = 1; k <= this->number_of_reactions_; k++)
				{
					if (only_reversible == false)
					{
						double one_over_Keq = -reaction_s_over_R__[k-1] + reaction_h_over_RT__[k-1] - log_Patm_over_RT_ * changeOfMoles__[k-1];	
						one_over_Keq = std::exp(one_over_Keq);

						y(i - 1, k - 1) = std::log(kArrheniusModified__[k-1]* one_over_Keq);
					}
					else
					{
						if (isThermodynamicallyReversible__[k-1] != 0)
						{
							const unsigned int j = isThermodynamicallyReversible__[k-1];
							y(i - 1, j - 1) = std::log(kArrheniusModified__[k-1] * uKeq__[j-1]);
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
		if (isThermodynamicallyReversible__[k-1] != 0)
		{
			const unsigned int j = isThermodynamicallyReversible__[k-1];

			fOut << std::setw(18) << std::right << std::scientific << std::setprecision(4) << std::exp(fittedKineticParameters(0, j-1));
			if (fittedKineticParameters.rows() == 2)
				fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << 0.;
			else
				fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << fittedKineticParameters(2, j - 1);
			fOut << std::setw(16) << std::right << std::fixed << std::setprecision(2) << fittedKineticParameters(1,j-1)/Conversions::J_from_kcal;
			fOut << std::setw(5)  << "";
		}
		else if (isExplicitlyReversible__[k-1] != 0)
		{
			const unsigned int j = isExplicitlyReversible__[k-1];
				
			fOut << std::setw(18) << std::right << std::scientific << std::setprecision(4) << std::exp(lnA_reversible__[j-1]);
			fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << Beta_reversible__[j-1];
			fOut << std::setw(16) << std::right << std::fixed << std::setprecision(2) << E_over_R_reversible__[j-1];
			fOut << std::setw(5)  << "";
		}
	}

	void KineticsMap_CHEMKIN::DerivativesOfFormationRates(const double* c, const double* omega, Eigen::MatrixXd* dR_over_domega)
	{
		Eigen::MatrixXd dc_over_domega(this->number_of_species_, this->number_of_species_);
		Eigen::MatrixXd dR_over_dc(this->number_of_species_, this->number_of_species_);

		// Molecular weight [kg/kmol]
		const double MW = thermodynamics_.MolecularWeight_From_MassFractions(omega);
		
		// Total concentration [kmol/m3]
		double cTot = 0.;
		for (unsigned int i = 0; i < this->number_of_species_; i++)
			cTot += c[i];

		thermodynamics_.DerivativesOfConcentrationsWithRespectToMassFractions(cTot, MW, omega, &dc_over_domega);
		DerivativesOfFormationRates(c, &dR_over_dc);

		for(unsigned int k=1;k<=this->number_of_species_;k++)
		{
			for(unsigned int i=1;i<=this->number_of_species_;i++)
			{
				double sum = 0.;
				for(unsigned int j=1;j<=this->number_of_species_;j++)
					sum += dR_over_dc(k-1,j-1) * dc_over_domega(j-1,i-1);
				(*dR_over_domega)(k-1,i-1) = sum;
			}
		}
	}

	void KineticsMap_CHEMKIN::Derivatives(const double* c, const double* omega, Eigen::MatrixXd* derivatives)
	{
		Eigen::MatrixXd dc_over_domega(this->number_of_species_, this->number_of_species_);
		Eigen::MatrixXd derivatives_over_dc(this->number_of_species_+1, this->number_of_species_+1);

		// Molecular weight [kg/kmol]
		const double MW = thermodynamics_.MolecularWeight_From_MassFractions(omega);
		
		// Total concentration [kmol/m3]
		double cTot = 0.;
		for (unsigned int i = 0; i < this->number_of_species_; i++)
			cTot += c[i];

		thermodynamics_.DerivativesOfConcentrationsWithRespectToMassFractions(cTot, MW, omega, &dc_over_domega);
		Derivatives(c, &derivatives_over_dc);

		for(unsigned int k=1;k<=this->number_of_species_;k++)
		{
			for(unsigned int i=1;i<=this->number_of_species_;i++)
			{
				double sum = 0.;
				for(unsigned int j=1;j<=this->number_of_species_;j++)
					sum += derivatives_over_dc(k-1,j-1) * dc_over_domega(j-1,i-1);
				(*derivatives)(k-1,i-1) = sum;
			}
		}

		for(unsigned int i=1;i<=this->number_of_species_;i++)
		{
			double sum = 0.;
			for(unsigned int j=1;j<=this->number_of_species_;j++)
				sum += derivatives_over_dc(this->number_of_species_,j-1) * dc_over_domega(j-1,i-1);
			(*derivatives)(this->number_of_species_, i-1) = sum;
		}

		for(unsigned int i=1;i<=this->number_of_species_+1;i++)
			(*derivatives)(i-1, this->number_of_species_) = derivatives_over_dc(i-1,this->number_of_species_);
	}

	void KineticsMap_CHEMKIN::DerivativesOfFormationRates(const double* c, Eigen::MatrixXd* dR_over_dC)
	{
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OPENSMOKE_MACH_EPS_DOUBLE);			
		const double TOLR = 100. * OPENSMOKE_MACH_EPS_FLOAT;
		const double TOLA = 1.e-10;
		
		OpenSMOKEVectorDouble c_plus(this->number_of_species_);
		OpenSMOKEVectorDouble R_original(this->number_of_species_);
		OpenSMOKEVectorDouble R_plus(this->number_of_species_);

		for (unsigned int i = 0; i < this->number_of_species_; i++)
			c_plus[i + 1] = c[i];

		// Calculates centered values
		ReactionRates(c);
		FormationRates(R_original.GetHandle());

		// derivata rispetto a ckd
		for(unsigned int kd=1;kd<=this->number_of_species_;kd++)
		{
			if(c[kd-1] <= 1.e-100)
			{
				double dc = 1.e-10 + 1.e-12 * c_plus[kd];
				double udc = 3.e-8 * Max(c_plus[kd],1./dc);
				udc = std::max(udc,1./dc);
				udc = std::max(udc,1.e-19);
				dc = std::min(udc, 0.001 + 0.001*c_plus[kd]);
				c_plus[kd] += dc;
				udc = 1./dc;

				ReactionRates(c_plus.GetHandle());
				FormationRates(R_plus.GetHandle());

				for(unsigned int j=1;j<=this->number_of_species_;j++)
					(*dR_over_dC)(j-1,kd-1) = (R_plus[j]-R_original[j]) * udc;

				c_plus[kd] = c[kd-1];
			}
			else
			{
				double hf = 1.e0;
				double error_weight = 1./(TOLA+TOLR*std::fabs(c[kd-1]));
				double hJ = ETA2 * std::fabs(std::max(c[kd-1], 1./error_weight));
				double hJf = hf/error_weight;
				hJ = std::max(hJ, hJf);
				hJ = std::max(hJ, ZERO_DER);

				// This is what is done by BzzMath
				double dc = std::min(hJ, 1.e-3 + 1e-3*std::fabs(c[kd-1]));

				// Thisis what is done in the KPP
				//double dc = TOLR*c[kd-1];

				double udc = 1. / dc;
				c_plus[kd] += dc;

				ReactionRates(c_plus.GetHandle());
				FormationRates(R_plus.GetHandle());

				for(unsigned int j=1;j<=this->number_of_species_;j++)
					(*dR_over_dC)(j-1,kd-1) = (R_plus[j]-R_original[j]) * udc;

				c_plus[kd] = c[kd-1];
			}
		}
	}

	void KineticsMap_CHEMKIN::Derivatives(const double* c, Eigen::MatrixXd* derivatives, const bool constant_density)
	{
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OPENSMOKE_MACH_EPS_DOUBLE);			
		const double TOLR = 100. * OPENSMOKE_MACH_EPS_FLOAT;
		const double TOLA = 1.e-10;
	
		double Q_original;
		OpenSMOKEVectorDouble c_plus(this->number_of_species_);
		OpenSMOKEVectorDouble R_original(this->number_of_species_);
		OpenSMOKEVectorDouble R_plus(this->number_of_species_);

		for (unsigned int i = 0; i < this->number_of_species_; i++)
			c_plus[i + 1] = c[i];
		
		SetTemperature(this->T_);
		SetPressure(this->P_);
		thermodynamics_.SetTemperature(this->T_);
		thermodynamics_.SetPressure(this->P_);
		
		KineticConstants();

		// Calculates centered values
		ReactionRates(c);
		FormationRates(R_original.GetHandle());
		Q_original = HeatRelease(R_original.GetHandle());

		// derivata rispetto a ckd
		for(unsigned int kd=1;kd<=this->number_of_species_;kd++)
		{
			if(c[kd-1] <= 1.e-100)
			{
				double dc = 1.e-10 + 1.e-12 * c_plus[kd];
				double udc = 3.e-8 * std::max(c_plus[kd],1./dc);
				udc = std::max(udc,1./dc);
				udc = std::max(udc,1.e-19);
				dc = std::min(udc, 0.001 + 0.001*c_plus[kd]);
				c_plus[kd] += dc;
				udc = 1./dc;

				ReactionRates(c_plus.GetHandle());
				FormationRates(R_plus.GetHandle());
				double Q_plus = HeatRelease(R_plus.GetHandle());

				for(unsigned int j=1;j<=this->number_of_species_;j++)
					(*derivatives)(j-1,kd-1) = (R_plus[j]-R_original[j]) * udc;
				(*derivatives)(this->number_of_species_,kd-1) = (Q_plus-Q_original) * udc;

				c_plus[kd] = c[kd-1];
			}
			else
			{
				double hf = 1.e0;
				double error_weight = 1./(TOLA+TOLR*std::fabs(c[kd-1]));
				double hJ = ETA2 * std::fabs(std::max(c[kd-1], 1./error_weight));
				double hJf = hf/error_weight;
				hJ = std::max(hJ, hJf);
				hJ = std::max(hJ, ZERO_DER);

				// This is what is done by BzzMath
				double dc = std::min(hJ, 1.e-3 + 1e-3*std::fabs(c[kd-1]));

				// Thisis what is done in the KPP
				//double dc = TOLR*c[kd-1];
				
				double udc = 1. / dc;
				c_plus[kd] += dc;
				ReactionRates(c_plus.GetHandle());
				FormationRates(R_plus.GetHandle());
				double Q_plus = HeatRelease(R_plus.GetHandle());

				for(unsigned int j=1;j<=this->number_of_species_;j++)
					(*derivatives)(j-1,kd-1) = (R_plus[j]-R_original[j]) * udc;
				(*derivatives)(this->number_of_species_, kd-1) = (Q_plus-Q_original) * udc;

				c_plus[kd] = c[kd-1];
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
			FormationRates(R_plus.GetHandle());
			double Q_plus = HeatRelease(R_plus.GetHandle());

			for(unsigned int j=1;j<=this->number_of_species_;j++)
				(*derivatives)(j-1, this->number_of_species_) = (R_plus[j]-R_original[j]) * udT;
			(*derivatives)(this->number_of_species_, this->number_of_species_) = (Q_plus-Q_original) * udT;
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
					sum += c[j-1]/this->T_ * (*derivatives)(k-1,j-1);
				(*derivatives)(k-1,index_T-1) -= sum;
			}
		}
	}

	void KineticsMap_CHEMKIN::FallOffReactions(	const unsigned int k, const double cTot, const double* c,
												double &F, double &dF_over_dA0, double &dF_over_dAInf)
	{
		double M;
		if (falloff_index_of_single_thirdbody_species__[k-1] == 0)
		{
			M = cTot;
			for(unsigned int s=0;s<falloff_indices_of_thirdbody_species__[k-1].size();s++)
				M += c[falloff_indices_of_thirdbody_species__[k-1][s]-1] * falloff_indices_of_thirdbody_efficiencies__[k-1][s];
		}
		else
		{
			const double epsilon = 1.e-16;
			M = c[falloff_index_of_single_thirdbody_species__[k-1]-1] + epsilon;
		}
			
		const unsigned int j=indices_of_falloff_reactions__[k-1];
		const double Pr = std::exp(lnA__[j-1] + Beta__[j-1]*this->logT_-E_over_R__[j-1]*this->uT_) * M / 
						  std::exp(lnA_falloff_inf__[k-1] + Beta_falloff_inf__[k-1]*this->logT_-E_over_R_falloff_inf__[k-1]*this->uT_);
		
		F = 1.;
		double dF_over_dPr = 0.;
		double nTroe, cTroe, sTroe, xSRI;
		switch(falloff_reaction_type__[k-1])
		{
			case PhysicalConstants::REACTION_TROE_FALLOFF:

				nTroe = 0.75-1.27*logFcent_falloff__[k-1];
				cTroe = -0.4-0.67*logFcent_falloff__[k-1];
				sTroe = std::log10(Pr) + cTroe;
				F = std::pow(10., logFcent_falloff__[k-1]/(1. + boost::math::pow<2>(sTroe/(nTroe-0.14*sTroe))));

				// Calculates the derivative 
				if (Pr > 1.e-32)
				{
					const double a = logFcent_falloff__[k-1];
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
					F = std::pow(logFcent_falloff__[k-1], xSRI) * d_falloff__[k-1];
					if (e_falloff__[k-1] != 0.)
						F *= std::pow(this->T_, e_falloff__[k-1]);

					const double eps = 0.01;
					const double PrPlus = (1.+eps)*Pr;
					const double xSRIPlus = 1. / (1. + boost::math::pow<2>(std::log10(PrPlus)));
					double FPlus = std::pow(logFcent_falloff__[k-1], xSRIPlus) * d_falloff__[k-1];
					if(e_falloff__[k-1] != 0.)
						FPlus *= std::pow(this->T_, e_falloff__[k-1]);
					dF_over_dPr = (FPlus-F)/(eps*Pr);
				}
				else
				{
					xSRI = 0.;
					F = std::pow(logFcent_falloff__[k-1], xSRI) * d_falloff__[k-1];
					if (e_falloff__[k-1] != 0.)
						F *= std::pow(this->T_, e_falloff__[k-1]);
					dF_over_dPr = 0.;
				}

				break;
		}

		dF_over_dA0   =  Pr/std::exp(lnA__[j-1]) * dF_over_dPr;
		dF_over_dAInf = -Pr/std::exp(lnA_falloff_inf__[k-1]) * dF_over_dPr;
	}

	void KineticsMap_CHEMKIN::ChemicallyActivatedBimolecularReactions(	const unsigned int k, const double cTot, const double* c,
																		double &F, double &dF_over_dA0, double &dF_over_dAInf)
	{
		double M;
		if (cabr_index_of_single_thirdbody_species__[k-1] == 0)
		{
			M = cTot;
			for(unsigned int s=0;s<cabr_indices_of_thirdbody_species__[k-1].size();s++)
				M += c[cabr_indices_of_thirdbody_species__[k-1][s]-1] * cabr_indices_of_thirdbody_efficiencies__[k-1][s];
		}
		else
		{
			const double epsilon = 1.e-16;
			M = c[cabr_index_of_single_thirdbody_species__[k-1]-1] + epsilon;
		}
			
		const unsigned int j=indices_of_cabr_reactions__[k-1];
		const double Pr = std::exp(lnA__[j-1] + Beta__[j-1]*this->logT_-E_over_R__[j-1]*this->uT_) * M / 
						  std::exp(lnA_cabr_inf__[k-1] + Beta_cabr_inf__[k-1]*this->logT_-E_over_R_cabr_inf__[k-1]*this->uT_);
		
		F = 1.;
		double dF_over_dPr = 0.;
		double nTroe, cTroe, sTroe, xSRI;
		switch(cabr_reaction_type__[k-1])
		{
			case PhysicalConstants::REACTION_TROE_CABR:

				nTroe = 0.75-1.27*logFcent_cabr__[k-1];
				cTroe = -0.4-0.67*logFcent_cabr__[k-1];
				sTroe = std::log10(Pr) + cTroe;
				F = std::pow(10., logFcent_cabr__[k-1]/(1. + boost::math::pow<2>(sTroe/(nTroe-0.14*sTroe))));

				// Calculates the derivative 
				if (Pr > 1.e-32)
				{
					const double eps = 0.01;
					const double PrPlus = (1.+eps)*Pr;
					const double nTroePlus = 0.75-1.27*logFcent_cabr__[k-1];
					const double cTroePlus = -0.4-0.67*logFcent_cabr__[k-1];
					const double sTroePlus = std::log10(PrPlus) + cTroePlus;
					const double FPlus = std::pow(10., logFcent_cabr__[k-1]/(1. + boost::math::pow<2>(sTroePlus/(nTroePlus-0.14*sTroePlus))));
					dF_over_dPr = (FPlus-F)/(eps*Pr);
				}

				break;

			case PhysicalConstants::REACTION_SRI_CABR:

				xSRI = 1. / (1. + boost::math::pow<2>(std::log10(Pr)));
				F = std::pow(logFcent_cabr__[k-1], xSRI) * d_cabr__[k-1];
				if(e_cabr__[k-1] != 0.)
					F *= std::pow(this->T_, e_cabr__[k-1]);

				// Calculates the derivative
				if (Pr > 1.e-32)
				{
					const double eps = 0.01;
					const double PrPlus = (1.+eps)*Pr;
					const double xSRIPlus = 1. / (1. + boost::math::pow<2>(std::log10(PrPlus)));
					double FPlus = std::pow(logFcent_cabr__[k-1], xSRIPlus) * d_cabr__[k-1];
					if(e_cabr__[k-1] != 0.)
						FPlus *= std::pow(this->T_, e_cabr__[k-1]);
					dF_over_dPr = (FPlus-F)/(eps*Pr);
				}

				break;
		}

		dF_over_dA0 = Pr/std::exp(lnA__[j-1]) * dF_over_dPr;
		dF_over_dAInf = -Pr/std::exp(lnA_falloff_inf__[k-1]) * dF_over_dPr;
	}

	double KineticsMap_CHEMKIN::ThirdBody(const unsigned int j, const unsigned int k)
	{
		for (unsigned int s = 0; s < number_of_thirdbody_reactions_; s++)
		{
			if (indices_of_thirdbody_reactions__[s] == j + 1)
			{
				for (unsigned int i = 0; i < indices_of_thirdbody_species__[s].size(); i++)
				{
					if (indices_of_thirdbody_species__[s][i] == k + 1)
					{
						return indices_of_thirdbody_efficiencies__[s][i]+1.;
						break;
					}
				}

				break;
			}

			if (indices_of_thirdbody_reactions__[s] > j + 1)
				break;
		}

		for (unsigned int s = 0; s < number_of_falloff_reactions_; s++)
		{
			if (indices_of_falloff_reactions__[s] == j + 1)
			{
				for (unsigned int i = 0; i < falloff_indices_of_thirdbody_species__[s].size(); i++)
				{
					if (falloff_indices_of_thirdbody_species__[s][i] == k + 1)
					{
						return falloff_indices_of_thirdbody_efficiencies__[s][i] + 1.;
						break;
					}
				}

				break;
			}

			if (indices_of_falloff_reactions__[s] > j + 1)
				break;
		}

		return 1.;
	}

	void KineticsMap_CHEMKIN::Set_ThirdBody(const unsigned int j, const unsigned int k, const double efficiency)
	{
		for (unsigned int s = 0; s < number_of_thirdbody_reactions_; s++)
		{
			if (indices_of_thirdbody_reactions__[s] == j+1)
			{
				for (unsigned int i = 0; i < indices_of_thirdbody_species__[s].size(); i++)
				{
					if (indices_of_thirdbody_species__[s][i] == k + 1)
					{
						indices_of_thirdbody_efficiencies__[s][i] = efficiency-1.;
						break;
					}
				}

				break;
			}

			if (indices_of_thirdbody_reactions__[s] > j + 1)
				break;
		}

		for (unsigned int s = 0; s < number_of_falloff_reactions_; s++)
		{
			if (indices_of_falloff_reactions__[s] == j + 1)
			{
				for (unsigned int i = 0; i < falloff_indices_of_thirdbody_species__[s].size(); i++)
				{
					if (falloff_indices_of_thirdbody_species__[s][i] == k + 1)
					{
						falloff_indices_of_thirdbody_efficiencies__[s][i] = efficiency - 1.;
						break;
					}
				}

				break;
			}

			if (indices_of_falloff_reactions__[s] > j + 1)
				break;
		}
	}
}

