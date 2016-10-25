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

#ifndef OpenSMOKE_ReactionPolicy_CHEMKIN_H
#define	OpenSMOKE_ReactionPolicy_CHEMKIN_H

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

	class ReactionPolicy_CHEMKIN 
	{
	public:

		/**
		* Default constructor
		*/
		ReactionPolicy_CHEMKIN();

		/**
		* Copy constructor
		*/
		ReactionPolicy_CHEMKIN(const ReactionPolicy_CHEMKIN& orig);
    
		/**
		* Default destructor
		*/
		//virtual ~ReactionPolicy_CHEMKIN();

		/**
		*@brief Returns the type of reaction
		*/
		PhysicalConstants::TAG_REACTION Tag() const; 

		/**
		*@brief Returns the type of reaction in a readable format
		*/
		std::string TagASCII() const;

		/**
		*@brief Reads and interprets the reaction from a vector of lines extracted from a kinetic file in CHEMKIN format
		*/
		bool ReadReactionFromStrings(const std::vector<std::string>& lines, const std::map<std::string, unsigned int>& map_of_species);

		/**
		*@brief Sets the units for the pre-exponential factor and the activation energy, according to what reported in 
		        the kinetic input file in CHEMKIN file
		*/
		bool SetUnits(const PhysicalConstants::UNITS_REACTION a_units, const PhysicalConstants::UNITS_REACTION e_units);
	
		/**
		*@brief Returns the stoichiometric coefficients of reactants
		*/
		const std::vector<double>& reactant_nu() const { return reactant_nu_;}

		/**
		*@brief Returns the stoichiometric coefficients of products
		*/
		const std::vector<double>& product_nu() const { return product_nu_;}

		/**
		*@brief Sets the stoichiometric coefficients of reactants
		*@param i index of species
		*@param value value of new stoichiometric coefficient
		*/
		void set_reactant_nu(const unsigned int i, const double value) { reactant_nu_[i] = value; }

		/**
		*@brief Sets the stoichiometric coefficients of product
		*@param i index of species
		*@param value value of new stoichiometric coefficient
		*/
		void set_product_nu(const unsigned int i, const double value)  { product_nu_[i] = value; }

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
		double sumLambdaReactants() const { return sumLambdaReactants_; }

		/**
		*@brief Returns the sum of reaction orders of products
		*/
		double sumLambdaProducts() const { return sumLambdaProducts_; }
	
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
		*@brief Pressure dependent reaction (fall-off or cabr): Arrhenius frequency factor (units: kmol, m, s)
		*/
		double A_inf() const { return AInf_; }

		/**
		*@brief Pressure dependent reaction (fall-off or cabr): Arrhenius temperature exponent
		*/
		double Beta_inf() const { return betaInf_; }

		/**
		*@brief Pressure dependent reaction (fall-off or cabr): Arrhenius activation temperature (i.e. activation energy divided by the gas constant R)
		*/
		double E_over_R_inf() const { return EInf_/PhysicalConstants::R_J_kmol; }

		/**
		*@brief Conversion factor of pre-exponential factor from OpenSMOKE to CGS (CHEMKIN)
		*/
		double A_conversion() const;

    	/**
		*@brief Conversion factor of high pressure pre-exponential factor from OpenSMOKE to CGS (CHEMKIN)
		*/    
		double A_inf_conversion() const;

		/**
		*@brief Conversion factor of reversible pre-exponential factor from OpenSMOKE to CGS (CHEMKIN)
		*/
		double Arev_conversion() const;

		/**
		*@brief Third body reaction: the indices of third-body species 
		*/
		const std::vector<unsigned int>& third_body_indices() const { return third_body_indices_; }

		/**
		*@brief Third body reaction: the efficiencies of third-body species 
		*/
		const std::vector<double>& third_body_efficiencies() const { return third_body_efficiencies_; }


		/**
		*@brief Is the reaction reversible with explicit kinetic parametrs (REV in CHEMKIN standard)?
		*/
		bool IsExplicitlyReversible() const { return iExplicitlyReversible_; }

		/**
		*@brief Is the reaction a pressure logarithmic interpolated?
		*/
		bool IsPressureLog() const { return iPlog_; }

		/**
		*@brief Is the reaction a Janev-Langer reaction?
		*/
		bool IsJanevLanger() const { return iJan_; }

		/**
		*@brief Is the reaction a fitted reaction?
		*/
		bool IsFit1() const { return iFit1_; }

		/**
		*@brief Is the reaction a Landau-Teller reaction?
		*/
		bool IsLandauTeller() const { return iLandauTeller_; }

		/**
		*@brief Is the reaction reversible with explicit kinetic parametrs (REV in CHEMKIN standard)?
		*/
		int pressureDependentSpeciesIndex() const { return pressureDependentSpeciesIndex_; }

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
		*@brief Writes the preprocessed reaction data on a file in a old-style format which can be read 
		        by a previous version of OpenSMOKE. This function has been kept only for 
		        compatibility reasons with previous versions.
		*/
		void WriteKineticsDataOnASCIIFileOldStyle(std::ofstream &fOutput) const;

		/**
		*@brief Writes additional data (details) on a file for non-conventional reactions (PLOG, FIT1, CHEB, etc.)
		*/
		void WriteAdditionalDataOnASCIIFile(std::ostream& fOut) const;

		/**
		*@brief Writes the three-body parameters on file
		*/
		void WriteThirdBodyParametersOnASCIIFile(std::ostream &fOutput) const;

		/**
		*@brief Writes additional data for the pressure dependent reactions (LOW and HIGH)
		*/
		void WritePressureDependentParametersOnASCIIFile(std::ostream &fOutput) const;

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


	protected:

		PhysicalConstants::TAG_REACTION tag_reaction_;			//!< type of reaction
		PhysicalConstants::UNITS_REACTION e_units_;				//!< units of activation energy
		PhysicalConstants::UNITS_REACTION a_units_;				//!< units of frequency factor
		double a_conversion_;									//!< conversion factor for the frequency factor
		double e_conversion_;									//!< conversion factor for the activation energy

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
		double sumLambdaReactants_;								//!< sum of reaction orders of reactants
		double sumLambdaProducts_;								//!< sum of reaction orders of products

		// Third body efficiencies
		bool iThirdBody_;										//!< are third-body coefficients supplied?
		std::vector<unsigned int>   third_body_indices_;		//!< indices of third-body species
		std::vector<double>			third_body_efficiencies_;	//!< third-body efficiencies
	
		// Low/High pressure limit
		int pressureDependentSpeciesIndex_;						//!< index of third-body species in pressure-dependent reactions (if available)
		bool iHigh_;											//!< is a Chemcally Activated Bimolecular Reaction?
		bool iLow_;												//!< is a fall-off reaction?
		double AInf_;											//!< pre-exponential factor
		double betaInf_;										//!< temperature exponent
		double EInf_;											//!< activation energy

		// Lindemann pressure dependent reaction
		bool iLindemann_;										//!< is the Lindemann form enabled?

		// SRI pressure dependent reaction
		bool iSRI_;												//!< is the SRI form enabled?
		std::vector<double> sri_;								//!< SRI coefficients

		// TROE pressure dependent reaction
		bool iTroe_;											//!< is the TROE form enabled?
		std::vector<double> troe_;								//!< TROE coefficients

		// Duplicate reaction
		bool iDuplicate_;										//!< is the reaction declared as duplicate?

		// Chebishev reaction
		bool iChebyshev_;										//!< is the reaction a Chebyshev reaction?
		std::vector<double> chebyshev_coefficients_;			//!< Chebishev coefficients
		std::vector<double> chebyshev_pressure_limits_;			//!< temperature limits
		std::vector<double> chebyshev_temperature_limits_;		//!< temperature limits

		// Power series modified Arrhenius law
		bool iFit1_;											//!< is a Power series modified Arrhenius law (FIT1)?
		std::vector<double> fit1_coefficients_;					//!< coefficients

		// Janev-Langer reaction rate
		bool iJan_;												//!< is a Janev-Langer reaction (JAN)?
		std::vector<double> janev_langer_coefficients_;			//!< coefficients

		// Landau-Teller reaction rate
		bool iLandauTeller_;									//!< is a Landau-Teller reaction (LT)?
		std::vector<double> landau_teller_coefficients_;		//!< coefficients

		// Pressure Dependence Through Logarithmic Interp.
		bool iPlog_;											//!< is a Pressure Dependence Through Logarithmic Interp. reaction (PLOG)?
		std::vector<double> plog_coefficients_;					//!< coefficients

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
		*@brief Gets the number of decimal places after comma
		*/               
                int GetDecimalPlaces(const double nu) const;
	};

}

#include "ReactionPolicy_CHEMKIN.hpp"

#endif	/* OpenSMOKE_ReactionPolicy_CHEMKIN_H*/

