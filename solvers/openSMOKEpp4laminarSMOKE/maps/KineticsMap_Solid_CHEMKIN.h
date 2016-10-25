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

#ifndef OpenSMOKE_KineticsMap_Solid_CHEMKIN_CHEMKIN_H
#define OpenSMOKE_KineticsMap_Solid_CHEMKIN_CHEMKIN_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include "KineticsMap.h"
#include "StoichiometricMap.h"
#include "rapidxml.hpp"

namespace OpenSMOKE
{
	typedef OpenSMOKEVector<PhysicalConstants::TAG_REACTION, OpenSMOKE::OneIndexPolicy > OpenSMOKEVectorReactionTag;

	enum TYPE_OF_SOLID_KINETICS { TYPE_OF_SOLID_KINETICS_CHEMKIN_CONVENTIONAL };

	//!  A class to efficiently evaluate the reaction and formation rates, to be used in production codes
	/*!
		 This class provides the tools to calculate in a very efficient way the reaction rates and the
		 formation rates. In order to ensure a good efficiency a map is created to store all the data
		 depending on the temperature. Inthis way they are recalculated only if strictly needed, i.e. only
		 if the temperature changes
	*/

	template<typename map> 
	class KineticsMap_Solid_CHEMKIN : public KineticsMap<map>
	{

	public:

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file  
		*@param target solid material index (from 1 to the number of materials)
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_Solid_CHEMKIN(ThermodynamicsMap_Solid_CHEMKIN<map>& thermo, rapidxml::xml_document<>& doc, const unsigned int target, const unsigned int nPoints = 1);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file  
		*@param target solid material name
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_Solid_CHEMKIN(ThermodynamicsMap_Solid_CHEMKIN<map>& thermo, rapidxml::xml_document<>& doc, const std::string target, const unsigned int nPoints = 1);

		/**
		*@brief Set the temperature at which the properties have to be evaluated
		*@param T the temperature value in K
		*/
		virtual void SetTemperature(const map& T);

		/**
		*@brief Set the pressure at which the properties have to be evaluated
		*@param P the pressure value in Pa
		*/
		virtual void SetPressure(const map& P);		

		/**
		*@brief Returns the names of the species
		        Please remember that the returned vector is 0-index based
		*/
		const std::vector<std::string>& NamesOfSpecies() const { return thermodynamics_.names(); }

		/**
		*@brief Imports the kinetic schemes from a file in XML format
		*/
		virtual void ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc);

		/**
		*@brief Imports the kinetic schemes from a file in XML format
		*/
		void ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc, const std::string target_material_name);

		/**
		*@brief Imports the kinetic schemes from a file in XML format
		*/
		void ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc, const unsigned int target_material_index);

		/**
		*@brief Imports the list of species from a file in XML format
		*/
		virtual void ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc);

		/**
		*@brief Calculates the kinetic constants of the reverse reactions
		*/
		void FittedReverseKineticConstants(OpenSMOKEVectorDouble& x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters);

		/**
		*@brief Calculates the kinetic constants of the reverse reactions
		*/
		void FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters);

		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k, OpenSMOKEVectorDouble& c_bath, const double conversion_forward=1., const double conversion_backward=1.);
		
		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k);

		/**
		*@brief Calculates the formation rates for all the species in the kinetic mechanism
		*@param Rgas the formation rates of all the gaseous species in [kmol/m3/s]
		*@param Rsolid the formation rates of all the solid-phase spacies in [kmol/m3/s]
		*/
		void FormationRates(OpenSMOKEVectorDouble* Rgas, OpenSMOKEVectorDouble* Rsolid);

		/**
		*@brief Calculates the heat release
		*@param Rgas the formation rates of all the gaseous species in [kmol/m3/s] (as returned by the FormationRates function)
		*@param Rsolid the formation rates of all the solid-phase spacies in [kmol/m3/s] (as returned by the FormationRates function)
		*@return the heat release due to the gas/solid reactions in [J/m3/s]
		*/
		double HeatRelease(const OpenSMOKEVectorDouble& Rgas, const OpenSMOKEVectorDouble& RSolid);

		/**
		*@brief Returns the forward reaction rates for all the reactions in the kinetic scheme
		*@param r the reaction rates of all the gaseous species in [kmol/m3/s]
		*/
		void GetForwardReactionRates(OpenSMOKEVectorDouble* r);

		/**
		*@brief Returns the backward reaction rates for all the reactions in the kinetic scheme
		        If a reaction is irreversible, it returns zero
		*@param r the reaction rates of all the gaseous species in [kmol/m3/s]
		*/
		void GetBackwardReactionRates(OpenSMOKEVectorDouble* r);

		/**
		*@brief Returns the net reaction rates in [kmol/m3/s]
		*/
		const OpenSMOKEVectorDouble& GetReactionRates();

		/**
		*@brief Returns the net reaction rates in [kmol/m3/s]
		*@param r the net reaction rates of all the reactions in [kmol/m3/s]
		*/
		void GetReactionRates(OpenSMOKEVectorDouble* r);

		/**
		*@brief Calculates the production and the destruction rates for all the species in the kinetic mechanism
		*/
		void ProductionAndDestructionRates(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D);
//      void ProductionAndDestructionRatesGross(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D);

		/**
		*@brief Calculates the reaction rates for all the reactions in the kinetic scheme
		*@param cGas the concentrations of all the gas-phase species in [kmol/m3/s]
		*@param cSolid the concentrations of all the solid-phase species in [kmol/m3/s]
		*/
		void ReactionRates(const OpenSMOKEVectorDouble& cGas, const OpenSMOKEVectorDouble& cSolid);

		/**
		*@brief Returns the indices of the reversible reactions
		*/
		const OpenSMOKEVectorUnsignedInt& indices_of_reversible_reactions() const { return indices_of_reversible_reactions_; }

		/**
		*@brief Returns the stoichiometric map
		*@brief the stoichiometric map
		*/
		StoichiometricMap& stoichiometry() { return *stoichiometry_; }
                
                /**
		*@brief Calculates the reaction enthalpies and entropies (to be used for the kinetic constants)
		*/
		void ReactionEnthalpiesAndEntropies();

		/**
		*@brief Calculates the kinetic constants
		*/
		void KineticConstants();


	private:

		/**
		*@brief Calculates the modified Arrhenius constants
		*/
		const OpenSMOKEVectorDouble& kArrheniusModified() const { return kArrheniusModified_; }

		/**
		*@brief Calculates the Arrhenius constants
		*/
		const OpenSMOKEVectorDouble& kArrhenius() const { return kArrhenius_; }

	private:

		ThermodynamicsMap_Solid_CHEMKIN<map>& thermodynamics_;		//!< reference to the thermodynamics
		OpenSMOKEVectorDouble c_;
		OpenSMOKEVectorDouble aux_vector_;

		OpenSMOKEVectorDouble reaction_entropy_over_R_;			//!< list of reaction entropies
		OpenSMOKEVectorDouble reaction_enthalpy_over_RT_;		//!< list of reaction enthalpies

		OpenSMOKEVectorUnsignedInt indices_of_irreversible_reactions_;				//!< indices of irreversible reactions
		OpenSMOKEVectorUnsignedInt indices_of_reversible_reactions_;				//!< indices of reversible reactions
		OpenSMOKEVectorUnsignedInt indices_of_thermodynamic_reversible_reactions_;	//!< indices of reversible (thermodynamic) reactions
		OpenSMOKEVectorUnsignedInt indices_of_explicitly_reversible_reactions_;		//!< indices of reversible (explicit) reactions


		unsigned int number_of_irreversible_reactions_;
		unsigned int number_of_reversible_reactions_;
		unsigned int number_of_thermodynamic_reversible_reactions_;
		unsigned int number_of_explicitly_reversible_reactions_;

		OpenSMOKEVectorDouble lnA_;								//!< frequency factors (log)
		OpenSMOKEVectorDouble Beta_;							//!< temperature exponents
		OpenSMOKEVectorDouble E_over_R_;						//!< activation temperatures
		OpenSMOKEVectorDouble forward_kinetic_order_;			//!< global kinetic order for forward reactions

		std::vector<unsigned int> indices_of_reactions_needing_conversion_;

		OpenSMOKEVectorDouble lnA_reversible_;					//!< frequency factors (log) for explicitly reversible reactions
		OpenSMOKEVectorDouble Beta_reversible_;					//!< temperature exponents for explicitly reversible reactions
		OpenSMOKEVectorDouble E_over_R_reversible_;				//!< activation temperatures for explicitly reversible reactions

		OpenSMOKEVectorDouble changeOfMoles_;		//!< list of change of moles

		StoichiometricMap* stoichiometry_;			//!< pointer to the stoichiometry

		bool arrhenius_kinetic_constants_must_be_recalculated_;
		bool nonconventional_kinetic_constants_must_be_recalculated_;
		bool reaction_h_and_s_must_be_recalculated_;

		OpenSMOKEVectorDouble reaction_s_over_R_;
		OpenSMOKEVectorDouble reaction_h_over_RT_;
		OpenSMOKEVectorDouble kArrheniusModified_;
		OpenSMOKEVectorDouble kArrhenius_;
		OpenSMOKEVectorDouble kArrhenius_reversible_;
		OpenSMOKEVectorDouble uKeq_;

		OpenSMOKEVectorDouble forwardReactionRates_;
		OpenSMOKEVectorDouble reverseReactionRates_;
		OpenSMOKEVectorDouble netReactionRates_;

		double Patm_over_RT_;
		double log_Patm_over_RT_;

		OpenSMOKEVectorUnsignedInt isThermodynamicallyReversible_;		//!< vector containing the local index of thermodynamically reversible reactions
		OpenSMOKEVectorUnsignedInt isExplicitlyReversible_;				//!< vector containing the local index of explicitly reversible reactions

		OpenSMOKEVector<PhysicalConstants::TAG_REACTION> type_of_reaction_;
		OpenSMOKEVectorUnsignedInt local_family_index_;

		// Thermodynamic reversible reactions
		std::vector<double> delta_nu_gas_;

		// Type of kinetics
		TYPE_OF_SOLID_KINETICS type_of_solid_kinetics_;
	};
}

#include "KineticsMap_Solid_CHEMKIN.hpp"

#endif /* OpenSMOKE_KineticsMap_Solid_CHEMKIN_H */
