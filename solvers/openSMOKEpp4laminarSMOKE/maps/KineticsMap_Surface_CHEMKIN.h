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

#ifndef OpenSMOKE_KineticsMap_Surface_CHEMKIN_CHEMKIN_H
#define OpenSMOKE_KineticsMap_Surface_CHEMKIN_CHEMKIN_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include "KineticsMap.h"
#include "StoichiometricMap.h"
#include "rapidxml.hpp"

namespace OpenSMOKE
{
	class UBIQEP_SubMechanism
	{
	public:

		void ReadFromXMLFile(rapidxml::xml_node<> *xml_node, const OpenSMOKE::OpenSMOKEVectorDouble& MW);
		void CalculateChemisorptionHeats(const double T, const OpenSMOKE::OpenSMOKEVectorDouble& Z);
		void CalculateDissociationEnergies(const OpenSMOKE::OpenSMOKEVectorDouble& H);
		void CalculateSurfaceEnthalpies(const double T);
		void CalculateActivationEnergies();
		void ForwardKineticConstants(const double T, const double total_site_density, OpenSMOKE::OpenSMOKEVectorDouble& kForward);

	private:
		unsigned int number_of_gas_species_;
		unsigned int number_of_site_species_;
		unsigned int half_number_of_ubiqep_reactions_;
		double reference_temperature_;

		OpenSMOKE::OpenSMOKEVectorUnsignedInt chemisorption_heats_gas_indices_;
		OpenSMOKE::OpenSMOKEVectorDouble chemisorption_heats_constant_coefficient_;
		OpenSMOKE::OpenSMOKEVectorDouble chemisorption_heats_temperature_coefficient_;
		OpenSMOKE::OpenSMOKEVectorDouble* chemisorption_heats_coefficients_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt* chemisorption_heats_indices_;
		OpenSMOKE::OpenSMOKEVectorDouble QStar_;
		OpenSMOKE::OpenSMOKEVectorDouble dissociation_energies_;
		OpenSMOKE::OpenSMOKEVectorDouble surface_enthalpies_;
		OpenSMOKE::OpenSMOKEVectorDouble E_forward_;
		OpenSMOKE::OpenSMOKEVectorDouble E_backward_;
		OpenSMOKE::OpenSMOKEVectorDouble adsorption_coefficient_;


		OpenSMOKE::OpenSMOKEVectorDouble sigma_direct_;
		OpenSMOKE::OpenSMOKEVectorDouble Beta_direct_;
		OpenSMOKE::OpenSMOKEVectorDouble lnA_direct_;

		OpenSMOKE::OpenSMOKEVectorDouble sigma_reverse_;
		OpenSMOKE::OpenSMOKEVectorDouble Beta_reverse_;
		OpenSMOKE::OpenSMOKEVectorDouble lnA_reverse_;

		OpenSMOKE::OpenSMOKEVectorDouble lambda_;
			
		OpenSMOKE::OpenSMOKEVectorUnsignedInt ubiqep_class_;
		OpenSMOKE::OpenSMOKEVector<PhysicalConstants::UBIQEP_TYPE> ubiqep_type_direct_;
		OpenSMOKE::OpenSMOKEVector<PhysicalConstants::UBIQEP_TYPE> ubiqep_type_reverse_;

		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_A_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_B_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_C_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_D_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_Star_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_A2_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_AB_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_ABStar_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_AStar_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_BStar_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_CStar_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_DStar_;

		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_A2_Chemisorption_Heats_;
		OpenSMOKE::OpenSMOKEVectorUnsignedInt index_AB_Chemisorption_Heats_;
	};


	typedef OpenSMOKEVector<PhysicalConstants::TAG_REACTION, OpenSMOKE::OneIndexPolicy > OpenSMOKEVectorReactionTag;

	enum TYPE_OF_KINETICS { TYPE_OF_KINETICS_CHEMKIN_CONVENTIONAL, TYPE_OF_KINETICS_UBI_QEP };

	//!  A class to efficiently evaluate the reaction and formation rates, to be used in production codes
	/*!
		 This class provides the tools to calculate in a very efficient way the reaction rates and the
		 formation rates. In order to ensure a good efficiency a map is created to store all the data
		 depending on the temperature. Inthis way they are recalculated only if strictly needed, i.e. only
		 if the temperature changes
	*/

	template<typename map> 
	class KineticsMap_Surface_CHEMKIN : public KineticsMap<map>
	{

	public:

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties (obsolete, TOREMOVE)
		*@param thermo the thermodynamic map
		*@param nSpecies number of species 
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_Surface_CHEMKIN(ThermodynamicsMap_Surface_CHEMKIN<map>& thermo, const unsigned int nSpecies, const unsigned int nPoints = 1);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file  
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_Surface_CHEMKIN(ThermodynamicsMap_Surface_CHEMKIN<map>& thermo, rapidxml::xml_document<>& doc, const unsigned int nPoints = 1);

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
		*/
		const std::vector<std::string>& NamesOfSpecies() const { return thermodynamics_.names(); }

		/**
		*@brief Imports the kinetic schemes from a file in XML format
		*/
		virtual void ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc);

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
		*/
		void FormationRates(OpenSMOKEVectorDouble* Rgas, OpenSMOKEVectorDouble* Rsite, OpenSMOKEVectorDouble* Rbulk, OpenSMOKEVectorDouble* RsitePhases);

		/**
		*@brief Calculates the heat release
		*/
		double HeatRelease(const OpenSMOKEVectorDouble& Rgas, const OpenSMOKEVectorDouble& Rsurface, const OpenSMOKEVectorDouble& Rbulk);

		/**
		*@brief Calculates the production and the destruction rates for all the species in the kinetic mechanism
		*/
//		void ProductionAndDestructionRates(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D);
//      void ProductionAndDestructionRatesGross(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D);


		/**
		*@brief Returns the forward reaction rates for all the reactions in the kinetic scheme
		*/
		void GetForwardReactionRates(OpenSMOKEVectorDouble* r);

		/**
		*@brief Returns the backward reaction rates for all the reactions in the kinetic scheme
		        If a reaction is irreversible, it returns zero
		*/
		void GetBackwardReactionRates(OpenSMOKEVectorDouble* r);

		/**
		*@brief Calculates the reaction rates for all the reactions in the kinetic scheme
		*/
		void ReactionRates(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& z, const OpenSMOKEVectorDouble& a, const OpenSMOKEVectorDouble& gamma);

		/**
		*@brief Calculates the reaction rates for all the lumped reactions in the kinetic scheme
		*/
		virtual void UserDefinedReactionRates(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& z, const OpenSMOKEVectorDouble& a, const OpenSMOKEVectorDouble& gamma);


//		void DerivativesOfReactionRatesWithRespectToKineticParameters(const PhysicalConstants::sensitivity_type type, unsigned int jReaction, const OpenSMOKEVectorDouble& c, double& parameter);

		/**
		*@brief Returns the indices of the reversible reactions
		*/
		const OpenSMOKEVectorUnsignedInt& indices_of_reversible_reactions() const { return indices_of_reversible_reactions_; }

		/**
		*@brief Calculates the reaction enthalpies and entropies (to be used for the kinetic constants)
		*/
		void ReactionEnthalpiesAndEntropies();

		/**
		*@brief Calculates the kinetic constants
		*/
		void ArrheniusKineticConstants();

		/**
		*@brief Return the net reaction rates in [kmol/m2/s]
		*/
		const OpenSMOKEVectorDouble& GetReactionRates();

		/**
		*@brief Return the net reaction rates in [kmol/m2/s]
		*/
		void GetReactionRates(OpenSMOKEVectorDouble* r);

//		void RateOfProductionAnalysis(const bool iNormalize) const;

//		void RateOfProductionAnalysis(std::ostream& fout) const;
//		void RateOfProductionAnalysis(ROPA_Data& ropa) const;

		const OpenSMOKEVectorDouble& kArrheniusModified() const { return kArrheniusModified_; }
		const OpenSMOKEVectorDouble& kArrhenius() const { return kArrhenius_; }

//		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const unsigned int k, const OpenSMOKEVectorDouble& c, OpenSMOKEVectorDouble* Jalfa, double& parameter);
//		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& parameter);
//		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& Jrho, double& parameter);


//		void DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* dR_over_dC);
//		void DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* dR_over_domega);

//		void Derivatives(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* derivatives, const bool constant_density = false);
//		void Derivatives(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* derivatives);


		StoichiometricMap& stoichiometry() { return *stoichiometry_; }

	protected:

		ThermodynamicsMap_Surface_CHEMKIN<map>& thermodynamics_;		//!< reference to the thermodynamics

		OpenSMOKEVectorDouble cSites_;
		OpenSMOKEVectorDouble c_;
		OpenSMOKEVectorDouble aux_vector_;

		OpenSMOKEVectorDouble reaction_entropy_over_R_;			//!< list of reaction entropies
		OpenSMOKEVectorDouble reaction_enthalpy_over_RT_;		//!< list of reaction enthalpies

		OpenSMOKEVectorUnsignedInt indices_of_irreversible_reactions_;				//!< indices of irreversible reactions
		OpenSMOKEVectorUnsignedInt indices_of_reversible_reactions_;				//!< indices of reversible reactions
		OpenSMOKEVectorUnsignedInt indices_of_thermodynamic_reversible_reactions_;	//!< indices of reversible (thermodynamic) reactions
		OpenSMOKEVectorUnsignedInt indices_of_explicitly_reversible_reactions_;		//!< indices of reversible (explicit) reactions
		OpenSMOKEVectorUnsignedInt indices_of_stick_reactions_;						//!< indices of stick reactions
		OpenSMOKEVectorUnsignedInt indices_of_coverage_dependent_reactions_;		//!< indices of coverage dependent reactions
		OpenSMOKEVectorUnsignedInt indices_of_langmuir_reactions_;					//!< indices of coverage dependent reactions
		OpenSMOKEVectorUnsignedInt indices_of_lumped_reactions_;					//!< indices of lumped reactions

		unsigned int number_of_irreversible_reactions_;
		unsigned int number_of_reversible_reactions_;
		unsigned int number_of_thermodynamic_reversible_reactions_;
		unsigned int number_of_explicitly_reversible_reactions_;
		unsigned int number_of_stick_reactions_;
		unsigned int number_of_coverage_dependent_reactions_;
		unsigned int number_of_langmuir_reactions_;
		unsigned int number_of_lumped_reactions_;

		OpenSMOKEVectorDouble lnA_;								//!< frequency factors (log)
		OpenSMOKEVectorDouble Beta_;							//!< temperature exponents
		OpenSMOKEVectorDouble E_over_R_;						//!< activation temperatures
		OpenSMOKEVectorDouble forward_kinetic_order_;			//!< global kinetic order for forward reactions

	//	std::vector<PhysicalConstants::UNITS_REACTION_COMPOSITION> units_of_reactions_needing_conversion_;
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

		// Stick reactions
		std::vector<double> stick_constant_coefficient_;
		std::vector<double> stick_power_;
		std::vector<bool> stick_motz_wise_;

		// Coverage dependent reactions
		std::vector< std::vector<bool> >			coverage_dependent_species_site_type_;
		std::vector< std::vector<unsigned int> >	coverage_dependent_species_index_;
		std::vector< std::vector<double> >			coverage_dependent_eta_;
		std::vector< std::vector<double> >			coverage_dependent_mu_;
		std::vector< std::vector<double> >			coverage_dependent_epsilon_;

		// Langmuir-Hinshelwood reactions
		std::vector< double >											langmuir_denominator_order_;
		std::vector< PhysicalConstants::UNITS_REACTION_COMPOSITION >	langmuir_units_;
		std::vector< std::vector<unsigned int> >						langmuir_species_index_;
		std::vector< std::vector<bool> >								langmuir_numerator_species_;
		std::vector< std::vector<double> >								langmuir_lnA_;
		std::vector< std::vector<double> >								langmuir_Beta_;
		std::vector< std::vector<double> >								langmuir_H_over_R_;
		std::vector< std::vector<double> >								langmuir_order_;

		// Lumped reactions
		std::vector< std::string >										names_of_lumped_functions_;

		// Non conservation of sites
		std::vector<unsigned int>	non_conservation_of_sites_indices_of_reactions_;
		std::vector<unsigned int>	non_conservation_of_sites_phase_of_reactions_;
		std::vector<double>			non_conservation_of_sites_delta_sigma_;

		// Thermodynamic reversible reactions
		std::vector<double> delta_sigma_times_log_Gamma0_;
		std::vector<double> delta_nu_gas_;

		// UBI
		TYPE_OF_KINETICS type_of_kinetics_;
		UBIQEP_SubMechanism* ubiqep_submechanism_;


		//void FallOffReactions(const unsigned int k, const double cTot, const OpenSMOKEVectorDouble& c, double &F, double &dF_over_dA0, double &dF_over_dAInf);
		//void ChemicallyActivatedBimolecularReactions(const unsigned int k, const double cTot, const OpenSMOKEVectorDouble& c, double &F, double &dF_over_dA0, double &dF_over_dAInf);
	};


}

#include "KineticsMap_Surface_CHEMKIN.hpp"

#endif /* OpenSMOKE_KineticsMap_Surface_CHEMKIN_H */
