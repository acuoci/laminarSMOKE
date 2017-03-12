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

#ifndef OpenSMOKE_ReactionPolicy_Surface_CHEMKIN_H
#define	OpenSMOKE_ReactionPolicy_Surface_CHEMKIN_H

namespace OpenSMOKE
{
	//!  A class to manage a single reaction according to the CHEMKIN standard
	/*!
			This class provides the tools to manage a single reaction according to the CHEMKIN
			standard. In particular it allows to read and interpret a reaction from a CHEMKIN input,
			write the reaction in a readable format and evaluates the reaction rate. The 
			functions are coded only for analysis purposes and not for intensive calculations
			(i.e. they are not efficient)
	*/

	class ReactionPolicy_Surface_CHEMKIN 
	{
	public:

		/**
		* Default constructor
		*/
		ReactionPolicy_Surface_CHEMKIN();

		/**
		* Copy constructor
		*/
		ReactionPolicy_Surface_CHEMKIN(const ReactionPolicy_Surface_CHEMKIN& orig);
    
		/**
		* Default destructor
		*/
		virtual ~ReactionPolicy_Surface_CHEMKIN();

		/**
		*@brief Returns the type of reaction
		*/
		PhysicalConstants::TAG_REACTION_SURFACE Tag() const; 

		/**
		*@brief Returns the type of reaction in a readable format
		*/
		std::string TagASCII() const;

		/**
		*@brief Reads and interprets the reaction from a vector of lines extracted from a kinetic file in CHEMKIN format
		*/
		bool ReadReactionFromStrings(	const std::vector<std::string>& lines, const std::map<std::string, unsigned int>& map_of_species,
										const unsigned int number_of_gas_species, const unsigned int number_of_site_species, const unsigned int number_of_bulk_species,
										const std::vector<double>& occupancy_site_species,
										const std::vector<unsigned int>& site_phase_membership, const std::vector<unsigned int>& bulk_phase_membership);

		/**
		*@brief Sets the units for the pre-exponential factor and the activation energy, according to what reported in 
		        the kinetic input file in CHEMKIN file
		*/
		bool SetUnits(const PhysicalConstants::UNITS_REACTION a_units, const PhysicalConstants::UNITS_REACTION e_units, const PhysicalConstants::UNITS_REACTION_COMPOSITION composition_units);

		/**
		*@brief Sets the Motz-Wise correction for sticking-coefficient reactions
		*/
		bool SetMotzWiseCorrection(const bool iMotzWiseCorrection);

		/**
		*@brief Sets the global non conservation of sites
		*/
		bool SetGlobalNonConservationOfSites(const bool global_non_conservation_of_sites);
	
		/**
		*@brief Returns the stoichiometric coefficients of reactants
		*/
		const std::vector<double>& reactant_nu() const { return reactant_nu_;}

		/**
		*@brief Returns the stoichiometric coefficients of products
		*/
		const std::vector<double>& product_nu() const { return product_nu_;}

		/**
		*@brief Returns the indices of reactant species (in the stoichiometric sense)
		*/
		const std::vector<unsigned int>& reactant_nu_indices() const { return reactant_nu_indices_;}

		/**
		*@brief Returns the indices of product species (in the stoichiometric sense)
		*/
		const std::vector<unsigned int>& product_nu_indices() const { return product_nu_indices_;}

		/**
		*@brief Returns the sum of stoichiometric coefficients of reactants
		*/
		double sumNuReactants() const { return sumNuReactants_; }

		/**
		*@brief Returns the sum of stoichiometric coefficients of products
		*/
		double sumNuProducts() const { return sumNuProducts_; }

		/**
		*@brief Returns the reaction orders of reactants
		*/
		const std::vector<double>& reactant_lambda() const { return reactant_lambda_;}

		/**
		*@brief Returns the reaction orders of products
		*/
		const std::vector<double>& product_lambda() const { return product_lambda_;}

		/**
		*@brief Returns the indices of reactants (in the kinetic sense)
		*/
		const std::vector<unsigned int>& reactant_lambda_indices() const { return reactant_lambda_indices_;}

		/**
		*@brief Returns the indices of products (in the kinetic sense)
		*/
		const std::vector<unsigned int>& product_lambda_indices() const { return product_lambda_indices_;}

		/**
		*@brief Returns the sum of reaction orders of reactants
		*/
//		double sumLambdaReactants() const { return sumLambdaReactants_; }

		/**
		*@brief Returns the sum of reaction orders of products
		*/
//		double sumLambdaProducts() const { return sumLambdaProducts_; }
	
		/**
		*@brief Was the reaction declared as DUPLICATE (or DUP) in the input kinetic file?
		*/
		bool IsDuplicate() const { return iDuplicate_; }

		/**
		*@brief Was the reaction of type FORD?
		*/
		bool IsFORD() const { return iFord_; }

		/**
		*@brief Was the reaction of type RORD?
		*/
		bool IsRORD() const { return iRord_; }

		/**
		*@brief Is the reaction reversible?
		*/
		bool IsReversible() const { return iReversible_; }

		/**
		*@brief Arrhenius frequency factor (units: kmol, m, s)
		*/
		double A() const { return A_; }

		/**
		*@brief Arrhenius temperature exponent
		*/
		double Beta() const { return beta_; }

		/**
		*@brief Arrhenius activation temperature (i.e. activation energy divided by the gas constant R)
		*/
		double E_over_R() const { return E_/PhysicalConstants::R_J_kmol; }

		/**
		*@brief Explicitly reversible reaction: Arrhenius frequency factor (units: kmol, m, s)
		*/
		double A_reversible() const { return ARev_; }

		/**
		*@brief Explicitly reversible reaction: Arrhenius temperature exponent
		*/
		double Beta_reversible() const { return betaRev_; }

		/**
		*@brief Explicitly reversible reaction: Arrhenius activation temperature (i.e. activation energy divided by the gas constant R)
		*/
		double E_over_R_reversible() const { return ERev_/PhysicalConstants::R_J_kmol; }

		/**
		*@brief Conversion factor of pre-exponential factor from OpenSMOKE to CGS (CHEMKIN)
		*/
		double A_conversion() const;

		/**
		*@brief Conversion factor of reversible pre-exponential factor from OpenSMOKE to CGS (CHEMKIN)
		*/
		double Arev_conversion() const;

		/**
		*@brief Is the reaction reversible with explicit kinetic parametrs (REV in CHEMKIN standard)?
		*/
		bool IsExplicitlyReversible() const { return iExplicitlyReversible_; }

		/**
		*@brief Returns a string reporting the reaction in a readable format
		*/
		void GetReactionString(const std::vector<std::string>& list_species, std::string& line_reaction) const;

		/**
		*@brief Returns a string reporting the reaction in CHEMKIN format (CGS units)
		*/
		void GetReactionStringCHEMKIN(const std::vector<std::string>& list_species, std::stringstream& line_reaction, const std::vector<bool>& isReducedSpecies) const;

		/**
		*@brief Returns a string reporting the reaction in CHEMKIN format (CGS units)
		*/
		void GetReactionStringCHEMKIN(const std::vector<std::string>& list_species, std::stringstream& line_reaction) const;

		/**
		*@brief Reports the reaction data (details) on a file
		*/
		void WriteSummary(std::ofstream& fOut, const std::vector<std::string>& list_species, const unsigned int index) const;

		/**
		*@brief Reports the reaction data (short summary) on a file
		*/
		void WriteShortSummary(std::ostream& fOut, const std::vector<std::string>& list_species) const;

		/**
		*@brief Writes additional data (details) on a file for non-conventional reactions (PLOG, FIT1, CHEB, etc.)
		*/
		void WriteAdditionalDataOnASCIIFile(std::ostream& fOut) const;

		/**
		*@brief Returns the conversion factor for the forward kinetic constant, i.e. to transform 
		        the [kmol, m, s] units into the [mol, cm, s] units
		*/
		double GetForwardConversionFactor() const;

		/**
		*@brief Returns the conversion factor for the backward kinetic constant, i.e. to transform 
		        the [kmol, m, s] units into the [mol, cm, s] units
		*/
		double GetBackwardConversionFactor() const;

		/**
		*@brief Set default units
		*/
		void SetDefaultUnits();

		/**
		*@brief Returns the kinetic order of the forward reaction (only gas)
		*/
		double forward_gas_kinetic_order() const { return sumLambdaGasReactants_; }

		/**
		*@brief Returns the units for the composition
		*/
		PhysicalConstants::UNITS_REACTION_COMPOSITION composition_units() const { return composition_units_; };

		/**
		*@brief Sticking coefficients: the sum of all the stoichiometric coefficients of reactants that are surface species
		*/
		double stick_power() const { return stick_power_; }

		/**
		*@brief Sticking coefficients: the index of the only gas phase species (1-index based)
		*/
		unsigned int stick_gas_species() const { return stick_gas_species_; }

		/**
		*@brief Sticking coefficients: the index of the only gas phase species (1-index based)
		*/
		unsigned int motz_wise_correction() const { return iMotzWiseCorrection_; }

		/**
		*@brief Sticking coefficients: the index of the only gas phase species (1-index based)
		*/
		const std::vector<unsigned int>& stick_indices_site_species() const { return stick_indices_site_species_; }

		/**
		*@brief Sticking coefficients: the index of the only gas phase species (1-index based)
		*/
		const std::vector<double>& stick_exponents_site_species() const { return stick_exponents_site_species_; }

		/**
		*@brief Coverage dependent reaction: returns true if the species is a site-species
		*/
		const std::vector<bool>& coverage_dependent_species_site_type() const { return coverage_dependent_species_site_type_; }
		
		/**
		*@brief Coverage dependent reaction: the index of the species (1-index based)
		*/
		const std::vector<unsigned int>& coverage_dependent_species_index() const { return coverage_dependent_species_index_; }

		/**
		*@brief Coverage dependent reaction: the first parameter (eta)
		*/
		const std::vector<double>& coverage_dependent_eta() const { return coverage_dependent_eta_; }

		/**
		*@brief Coverage dependent reaction: the second parameter (mu)
		*/
		const std::vector<double>& coverage_dependent_mu() const { return coverage_dependent_mu_; }

		/**
		*@brief Coverage dependent reaction: the third parameter (epsilon)
		*/
		const std::vector<double>& coverage_dependent_epsilon() const { return coverage_dependent_epsilon_; }

		/**
		*@brief Langmuir-Hinshelwood reaction: the index of the species (1-index based)
		*/
		const std::vector<unsigned int>& langmuir_species_index() const { return langmuir_species_index_; }

		/**
		*@brief Langmuir-Hinshelwood reaction: the pre-exponential factor
		*/
		const std::vector<double>& langmuir_A() const { return langmuir_A_; }

		/**
		*@brief Langmuir-Hinshelwood reaction: the temperature exponent
		*/
		const std::vector<double>& langmuir_Beta() const { return langmuir_Beta_; }

		/**
		*@brief Langmuir-Hinshelwood reaction: the equilibrium enthalpy (in K)
		*/
		const std::vector<double>& langmuir_H_over_R() const { return langmuir_H_over_R_; }

		/**
		*@brief Langmuir-Hinshelwood reaction: the reaction order
		*/
		const std::vector<double>& langmuir_order() const { return langmuir_order_; }

		/**
		*@brief Langmuir-Hinshelwood reaction: the exponent for the denominator (default equal to 2)
		*/
		double langmuir_denominator_order() const { return langmuir_denominator_order_; }

		/**
		*@brief Langmuir-Hinshelwood reaction: true if the species contributes to the numerator
		*/
		const std::vector<bool>& langmuir_numerator_species() const { return langmuir_numerator_species_; }

		/**
		*@brief Langmuir-Hinshelwood reaction: returns the units of measure of composition
		*/
		PhysicalConstants::UNITS_REACTION_COMPOSITION langmuir_units() const { return langmuir_units_; }
		
		/**
		*@brief Lumped reaction: returns the name of the function to be called to calculate the reaction rate
		*/
		const std::string& name_of_lumped_function() const { return name_of_lumped_function_; }

		/**
		*@brief Returns the difference between the product sites and the reactant sites
		*/
		double delta_occupancy_sites() const { return delta_occupancy_sites_; }

		/**
		*@brief Returns the difference between the stoichiometric coefficients of gas species
		*/
		double delta_nu_gas() const { return delta_nu_gas_; }

		/**
		*@brief Returns the site phase to which the reaction belongs (zero if the reaction
		        does not involve any surface species)
		*/
		unsigned int surface_reaction_membership() const { return surface_reaction_membership_; }


		bool isUBIQEP() const { return iUBIQEP_; }
		bool isUBIQEP_Direct() const { return ubiqep_direct_; }
		unsigned int UBIQEP_Reaction_Class() const { return ubiqep_reaction_class_; }
		bool isLumped() const { return iLumped_;  }

		bool WriteUBIParametersOnFile(std::ostream& fOutput) const;
		

	protected:

		PhysicalConstants::TAG_REACTION_SURFACE			tag_reaction_;			//!< type of reaction
		PhysicalConstants::UNITS_REACTION				e_units_;				//!< units of activation energy
		PhysicalConstants::UNITS_REACTION				a_units_;				//!< units of frequency factor
		PhysicalConstants::UNITS_REACTION_COMPOSITION	composition_units_;		//!< units of composition
		double a_conversion_;													//!< conversion factor for the frequency factor
		double e_conversion_;													//!< conversion factor for the activation energy

		// Arrhenius kinetic parameters
		double A_;												//!< pre-exponential factor
		double beta_;											//!< temperature exponent
		double E_;												//!< activation energy

		// Explicit reverse Arrhenius kinetic parameters
		double ARev_;											//!< pre-exponential factor (of reversible reaction, if REV is enabled)
		double betaRev_;										//!< temperature exponent (of reversible reaction, if REV is enabled)
		double ERev_;											//!< activation energy (of reversible reaction, if REV is enabled)

		// Reversible/Irreversible reaction
		bool iReversible_;										//!< is the reaction reversible?
		bool iExplicitlyReversible_;							//!< is the reaction reversible with explicitly assigned kinetic coefficients (REV)?
	
		// Stoichiometric coefficients
		std::vector<double> reactant_nu_;						//!< stoichiometric coefficients of reactants
		std::vector<double> product_nu_;						//!< stoichiometric coefficients of products
		std::vector<unsigned int>	reactant_nu_indices_;		//!< indices of stoichiometric coefficients of reactants
		std::vector<unsigned int>	product_nu_indices_;		//!< indices of stoichiometric coefficients of products
		double sumNuReactants_;									//!< sum of stoichiometric coefficients of reactants
		double sumNuProducts_;									//!< sum of stoichiometric coefficients of products

		// Reaction orders
		bool iFord_;											//!< are the reactant reaction orders different than stoichiometric coefficients?
		bool iRord_;											//!< are the reactant reaction orders different than stoichiometric coefficients?
		std::vector<double> reactant_lambda_;					//!< reaction orders of reactants
		std::vector<double> product_lambda_;					//!< reaction orders of products
		std::vector<unsigned int>	reactant_lambda_indices_;	//!< indices of reactants (in the kinetic sense)
		std::vector<unsigned int>	product_lambda_indices_;	//!< indices of products (in the kinetic sense)
		
		double sumLambdaGasReactants_;								//!< sum of reaction orders of reactants
		double sumLambdaGasProducts_;								//!< sum of reaction orders of products
		double sumLambdaSiteReactants_;								//!< sum of reaction orders of reactants
		double sumLambdaSiteProducts_;								//!< sum of reaction orders of products
		double sumLambdaBulkReactants_;								//!< sum of reaction orders of reactants
		double sumLambdaBulkProducts_;								//!< sum of reaction orders of products

		// Conservation of sites
		double delta_occupancy_sites_;
		double delta_nu_gas_;
		bool global_non_conservation_of_sites_;
		unsigned int surface_reaction_membership_;

		// Duplicate reaction
		bool iDuplicate_;										//!< is the reaction declared as duplicate?

		// Sticking coefficients
		bool iStick_;
		unsigned int stick_gas_species_;
		double stick_power_;
		bool iMotzWiseCorrection_;
		std::vector<unsigned int> stick_indices_site_species_;
		std::vector<double> stick_exponents_site_species_;

		// Coverage dependent reactions
		bool iCoverageDependent_;
		std::vector<bool> coverage_dependent_species_site_type_;
		std::vector<unsigned int> coverage_dependent_species_index_;
		std::vector<double> coverage_dependent_eta_;
		std::vector<double> coverage_dependent_mu_;
		std::vector<double> coverage_dependent_epsilon_;


		// Langmuir-Hinshelwood
		bool iLangmuir_;
		std::vector<unsigned int> langmuir_species_index_;
		std::vector<unsigned int> langmuir_numerator_species_index_provisional_;
		std::vector<bool>	langmuir_numerator_species_;
		std::vector<double> langmuir_A_;
		std::vector<double> langmuir_Beta_;
		std::vector<double> langmuir_H_over_R_;
		std::vector<double> langmuir_order_;
		double langmuir_denominator_order_;
		PhysicalConstants::UNITS_REACTION_COMPOSITION langmuir_units_;
		
		// Lumped
		bool iLumped_;
		std::string name_of_lumped_function_;

		// UBI-QEP Reactions
		bool iUBIQEP_;
		bool ubiqep_direct_;
		unsigned int ubiqep_reaction_class_;
		PhysicalConstants::UBIQEP_TYPE ubiqep_reaction_type_;
		double ubiqep_A_;
		double ubiqep_Beta_;
		double ubiqep_BondIndex_;
		unsigned int ubiqep_index_A_;

		unsigned int ubiqep_index_B_;
		unsigned int ubiqep_index_C_;
		unsigned int ubiqep_index_D_;

		unsigned int ubiqep_index_A2_;
		unsigned int ubiqep_index_AB_;
		unsigned int ubiqep_index_Star_;
		unsigned int ubiqep_index_AStar_;
		unsigned int ubiqep_index_BStar_;
		unsigned int ubiqep_index_CStar_;
		unsigned int ubiqep_index_DStar_;
		unsigned int ubiqep_index_ABStar_;
		std::vector<unsigned int> original_reactant_nu_indices_;
		std::vector<unsigned int> original_product_nu_indices_;
		std::vector<double> original_reactant_nu_;
		std::vector<double> original_product_nu_;			

		//
		unsigned int number_of_gas_species_;
		unsigned int number_of_site_species_;
		unsigned int number_of_bulk_species_;

	private:

		/**
		*@brief Calculates the global reaction order of the reaction
		*/
		void ReactionOrders(); 

		/**
		*@brief Converts the units for the pre-exponential factor and the activation energy
		*/
		void ConvertUnits();

		/**
		*@brief Calcultaes the conversion factors for the pre-exponential factor and the activation energy
		*/
		void ConversionFactors();

		/**
		*@brief Error message
		*/
		bool FatalErrorMessage(const std::string message);
	};

}

#include "ReactionPolicy_Surface_CHEMKIN.hpp"

#endif	/* OpenSMOKE_ReactionPolicy_Surface_CHEMKIN_H*/

