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

#ifndef OpenSMOKE_KineticsMap_CHEMKIN_CHEMKIN_H
#define OpenSMOKE_KineticsMap_CHEMKIN_CHEMKIN_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include "KineticsMap.h"
#include "StoichiometricMap.h"
#include "JacobianSparsityPatternMap.h"
#include "rapidxml.hpp"

namespace OpenSMOKE
{
	enum EnergyEquationType { CONSTANT_VOLUME_SYSTEM, CONSTANT_PRESSURE_SYSTEM, PLUGFLOW_SYSTEM };

	typedef OpenSMOKEVector<PhysicalConstants::TAG_REACTION, OpenSMOKE::OneIndexPolicy > OpenSMOKEVectorReactionTag;

	//!  A class to efficiently evaluate the reaction and formation rates, to be used in production codes
	/*!
		 This class provides the tools to calculate in a very efficient way the reaction rates and the
		 formation rates. In order to ensure a good efficiency a map is created to store all the data
		 depending on the temperature. Inthis way they are recalculated only if strictly needed, i.e. only
		 if the temperature changes
	*/

	template<typename map> 
	class KineticsMap_CHEMKIN : public KineticsMap<map>
	{

	public:

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties (obsolete, TOREMOVE)
		*@param thermo the thermodynamic map
		*@param nSpecies number of species 
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN<map>& thermo, const unsigned int nSpecies, const unsigned int nPoints = 1);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file  
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN<map>& thermo, rapidxml::xml_document<>& doc, const unsigned int nPoints = 1);
                
                /**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file  
                *@param verbose activate output
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN<map>& thermo, rapidxml::xml_document<>& doc, bool verbose, const unsigned int nPoints = 1);

                /**
		*@brief Copy constructor
                *       IMPORTANT: Please consider that the thermodynamic map to which the current kinetic map (this) is linked, is
                *       still the same thermodynamic map to which the rhs map is linked! If you want to link the current 
                *       kinetic map to a different thermodynamic map you should use a different copy constructor (see below)
		*@param rhs the object to be copied in the current object
		*/
                KineticsMap_CHEMKIN( const KineticsMap_CHEMKIN& rhs );
                
                /**
		*@brief Copy constructor
                *       IMPORTANT: no checks are done about the consistency between the thermodynamic map linked to the kinetic map
                *       (rhs.thermodynamics_) and the provided map (thermo). It is up to the user to be sure about the consistency.
		*@param rhs the object to be copied in the current object
		*/
                KineticsMap_CHEMKIN( const KineticsMap_CHEMKIN& rhs, ThermodynamicsMap_CHEMKIN<map>& thermo );
                		
		/**
		*@brief Default destructor
		*/
		~KineticsMap_CHEMKIN();

		ThermodynamicsMap_CHEMKIN<map>& thermodynamics() const { return thermodynamics_; }

		/**
		*@brief Sets the verbose output
		*/
		void SetVerboseOutput(const bool verbose_output)		{ verbose_output_ = verbose_output; }

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
		void WriteKineticData(std::ostream& fOut, const unsigned int k, const double T, const double P_Pa, OpenSMOKEVectorDouble& c);

		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k);

		/**
		*@brief Calculates the formation rates for all the species in the kinetic mechanism
		*/
		void FormationRates(OpenSMOKEVectorDouble* R);

		/**
		*@brief Calculates the heat release
		*/
		double HeatRelease(const OpenSMOKEVectorDouble& R);

		/**
		*@brief Calculates the production and the destruction rates for all the species in the kinetic mechanism [kmol/m3/s]
		*/
		void ProductionAndDestructionRates(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D);
        void ProductionAndDestructionRatesGross(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D);


		/**
		*@brief Returns the forward reaction rates for all the reactions in the kinetic scheme [kmol/m3/s]
		*/
		void GetForwardReactionRates(OpenSMOKEVectorDouble* r);

		/**
		*@brief Returns the backward reaction rates for all the reactions in the kinetic scheme [kmol/m3/s]
		        If a reaction is irreversible, it returns zero
		*/
		void GetBackwardReactionRates(OpenSMOKEVectorDouble* r);

		/**
		*@brief Calculates the reaction rates for all the reactions in the kinetic scheme
		*/
		void ReactionRates(const OpenSMOKEVectorDouble& c);
		void DerivativesOfReactionRatesWithRespectToKineticParameters(const PhysicalConstants::sensitivity_type type, unsigned int jReaction, const OpenSMOKEVectorDouble& c, double& parameter);

		/**
		*@brief Calculates the reaction rates for all the reactions in the kinetic scheme 
		*/
		void ReactionRates(const OpenSMOKEVectorDouble& c, const double cTot);

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
		void KineticConstants();

		/**
		*@brief Return the net reaction rates in [kmol/m3/s]
		*/
		const OpenSMOKEVectorDouble& GetReactionRates();

		/**
		*@brief Return the net reaction rates in [kmol/m3/s]
		*/
		void GetReactionRates(OpenSMOKEVectorDouble* r);

		void RateOfProductionAnalysis(const bool iNormalize) const;

		void RateOfProductionAnalysis(std::ostream& fout) const;
		void RateOfProductionAnalysis(ROPA_Data& ropa) const;
		void RateOfProductionAnalysis(ROPA_Data& ropa, const OpenSMOKE::OpenSMOKEVectorDouble& rf, const OpenSMOKE::OpenSMOKEVectorDouble& rb) const;


		const OpenSMOKEVectorDouble& kArrheniusModified() const { return kArrheniusModified_; }
		const OpenSMOKEVectorDouble& kArrhenius() const { return kArrhenius_; }

		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const unsigned int k, const OpenSMOKEVectorDouble& c, OpenSMOKEVectorDouble* Jalfa, double& parameter);
		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& parameter);
		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& mole_fractions, OpenSMOKEVectorDouble* Jalfa, double& JT, double& Jrho, double& parameter);


		void DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* dR_over_dC);
		void DerivativesOfFormationRates(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* dR_over_domega);

		void Derivatives(const OpenSMOKEVectorDouble& c, OpenSMOKEMatrixDouble* derivatives, const bool constant_density = false);
		void Derivatives(const OpenSMOKEVectorDouble& c, const OpenSMOKEVectorDouble& omega, OpenSMOKEMatrixDouble* derivatives);

		unsigned int NumberOfReversibleReactions() const { return number_of_reversible_reactions_; }

		unsigned int NumberOfThirdBodyReactions() const { return number_of_thirdbody_reactions_; }

		unsigned int NumberOfFallOffReactions() const { return number_of_falloff_reactions_; }

		unsigned int NumberOfCABRReactions() const { return number_of_cabr_reactions_; }

		StoichiometricMap& stoichiometry() { return *stoichiometry_; }

		const OpenSMOKEVectorUnsignedInt& indices_of_thirdbody_reactions() const { return indices_of_thirdbody_reactions_; };
		const std::vector<OpenSMOKEVectorUnsignedInt>& indices_of_thirdbody_species() const { return indices_of_thirdbody_species_; };
		const OpenSMOKEVectorUnsignedInt& indices_of_falloff_reactions() const { return indices_of_falloff_reactions_; };
		const OpenSMOKEVectorUnsignedInt& indices_of_cabr_reactions() const { return indices_of_cabr_reactions_; };

		const OpenSMOKEVectorUnsignedInt& isThermodynamicallyReversible() const { return isThermodynamicallyReversible_; };
		const OpenSMOKEVectorUnsignedInt& isExplicitlyReversible() const { return isExplicitlyReversible_; };

		OpenSMOKEVectorDouble& netReactionRates() { return netReactionRates_; }

		OpenSMOKEVectorDouble& Meff() { return Meff_; }

		const std::vector<OpenSMOKEVectorDouble>&	indices_of_thirdbody_efficiencies() const { return indices_of_thirdbody_efficiencies_ ; }

		void WeakThirdBodyConcentrationEfficiencies(std::vector<unsigned int>& reaction, std::vector<unsigned int>& species);
		void WeakFallOffConcentrationEfficiencies(std::vector<unsigned int>& reaction, std::vector<unsigned int>& species);
		void WeakCABRConcentrationEfficiencies(std::vector<unsigned int>& reaction, std::vector<unsigned int>& species);
		void StrongConcentrationEfficiencies(std::vector<unsigned int>& reaction);

		double FallOffReactionsCorrection(const unsigned int local_k, const double cTot, const OpenSMOKEVectorDouble& c);

		JacobianSparsityPatternMap<KineticsMap_CHEMKIN>* jacobian_sparsity_pattern_map() { return jacobian_sparsity_pattern_map_; };

		/**
		*@brief Return the frequency factor of a single reaction [kmol, m, s]
		*@param j index of reaction (starting from zero)
		*/
		double A(const unsigned int j) { return std::exp(lnA_[j + 1]); }

		/**
		*@brief Return the temperature exponent a single reaction 
		*@param j index of reaction (starting from zero)
		*/
		double Beta(const unsigned int j) { return Beta_[j + 1]; }

		/**
		*@brief Return the activation temperature of a single reaction [K]
		*@param j index of reaction (starting from zero)
		*/
		double E_over_R(const unsigned int j) { return E_over_R_[j + 1]; }
		

	private:

		/**
		*@brief Calculates the efficiencies for threebody reactions
		*/
		void ThirdBodyReactions(const double cTot, const OpenSMOKEVectorDouble& c);

		/**
		*@brief Calculates the corrections for the falloff reactions
		*/
		void FallOffReactions(const double cTot, const OpenSMOKEVectorDouble& c);

		/**
		*@brief Calculates the corrections for the cabr reactions
		*/
		void ChemicallyActivatedBimolecularReactions(const double cTot, const OpenSMOKEVectorDouble& c);
                
                /**
		*@brief Copies the data from another kinetic map (used by copy constructors)
		*/
                void CopyFromMap( const KineticsMap_CHEMKIN& rhs );
                
                // TODO
                void FallOffReactions(const unsigned int k, const double cTot, const OpenSMOKEVectorDouble& c, double &F, double &dF_over_dA0, double &dF_over_dAInf);
		void ChemicallyActivatedBimolecularReactions(const unsigned int k, const double cTot, const OpenSMOKEVectorDouble& c, double &F, double &dF_over_dA0, double &dF_over_dAInf);




	private:

		ThermodynamicsMap_CHEMKIN<map>& thermodynamics_;		//!< reference to the thermodynamics

		OpenSMOKEVectorDouble reaction_entropy_over_R_;			//!< list of reaction entropies
		OpenSMOKEVectorDouble reaction_enthalpy_over_RT_;		//!< list of reaction enthalpies

		OpenSMOKEVectorUnsignedInt indices_of_irreversible_reactions_;				//!< indices of irreversible reactions
		OpenSMOKEVectorUnsignedInt indices_of_reversible_reactions_;				//!< indices of reversible reactions
		OpenSMOKEVectorUnsignedInt indices_of_thermodynamic_reversible_reactions_;	//!< indices of reversible (thermodynamic) reactions
		OpenSMOKEVectorUnsignedInt indices_of_explicitly_reversible_reactions_;		//!< indices of reversible (explicit) reactions
		OpenSMOKEVectorUnsignedInt indices_of_thirdbody_reactions_;					//!< indices of three-body reactions
		OpenSMOKEVectorUnsignedInt indices_of_falloff_reactions_;					//!< indices of falloff reactions
		OpenSMOKEVectorUnsignedInt indices_of_cabr_reactions_;						//!< indices of cabr reactions
		OpenSMOKEVectorUnsignedInt indices_of_chebyshev_reactions_;					//!< indices of chebyshev reactions
		OpenSMOKEVectorUnsignedInt indices_of_pressurelog_reactions_;				//!< indices of pressurelog (PLOG) reactions
		OpenSMOKEVectorUnsignedInt indices_of_fit1_reactions_;						//!< indices of FIT1 reactions
		OpenSMOKEVectorUnsignedInt indices_of_janevlanger_reactions_;				//!< indices of JAN reactions
		OpenSMOKEVectorUnsignedInt indices_of_landauteller_reactions_;				//!< indices of LT reactions

		OpenSMOKEVectorUnsignedInt indices_of_falloff_lindemann_reactions_;			//!< indices of falloff (Lindemann) reactions
		OpenSMOKEVectorUnsignedInt indices_of_cabr_lindemann_reactions_;			//!< indices of cabr (Lindemann) reactions
		OpenSMOKEVectorUnsignedInt indices_of_falloff_troe_reactions_;				//!< indices of falloff (Troe) reactions
		OpenSMOKEVectorUnsignedInt indices_of_cabr_troe_reactions_;					//!< indices of cabr (Troe) reactions
		OpenSMOKEVectorUnsignedInt indices_of_falloff_sri_reactions_;				//!< indices of falloff (SRI) reactions
		OpenSMOKEVectorUnsignedInt indices_of_cabr_sri_reactions_;					//!< indices of cabr (Troe) reactions

		unsigned int number_of_irreversible_reactions_;
		unsigned int number_of_reversible_reactions_;
		unsigned int number_of_thermodynamic_reversible_reactions_;
		unsigned int number_of_explicitly_reversible_reactions_;
		unsigned int number_of_thirdbody_reactions_;
		unsigned int number_of_falloff_reactions_;
		unsigned int number_of_cabr_reactions_;
		unsigned int number_of_chebyshev_reactions_;
		unsigned int number_of_pressurelog_reactions_;
		unsigned int number_of_fit1_reactions_;
		unsigned int number_of_janevlanger_reactions_;
		unsigned int number_of_landauteller_reactions_;

		unsigned int number_of_falloff_lindemann_reactions_;
		unsigned int number_of_cabr_lindemann_reactions_;
		unsigned int number_of_falloff_troe_reactions_;
		unsigned int number_of_cabr_troe_reactions_;
		unsigned int number_of_falloff_sri_reactions_;
		unsigned int number_of_cabr_sri_reactions_;

		bool verbose_output_;

		OpenSMOKEVectorDouble lnA_;								//!< frequency factors (log)
		OpenSMOKEVectorDouble Beta_;							//!< temperature exponents
		OpenSMOKEVectorDouble E_over_R_;						//!< activation temperatures

		OpenSMOKEVectorDouble lnA_reversible_;					//!< frequency factors (log) for explicitly reversible reactions
		OpenSMOKEVectorDouble Beta_reversible_;					//!< temperature exponents for explicitly reversible reactions
		OpenSMOKEVectorDouble E_over_R_reversible_;				//!< activation temperatures for explicitly reversible reactions

		OpenSMOKEVectorDouble Meff_;														//!< threebody efficiencies
		std::vector<OpenSMOKEVectorUnsignedInt>		indices_of_thirdbody_species_;		//!< indices of threebody species
		std::vector<OpenSMOKEVectorDouble>		indices_of_thirdbody_efficiencies_;	//!< efficiencies of threebody species

		OpenSMOKEVectorDouble lnA_falloff_inf_;			//!< frequency factors (log) for falloff reactions
		OpenSMOKEVectorDouble Beta_falloff_inf_;		//!< temperature exponents for falloff reactions
		OpenSMOKEVectorDouble E_over_R_falloff_inf_;		//!< activation temperatures for falloff reactions

		std::vector<OpenSMOKEVectorUnsignedInt>         falloff_indices_of_thirdbody_species_;
		std::vector<OpenSMOKEVectorDouble>		falloff_indices_of_thirdbody_efficiencies_;
		OpenSMOKEVectorInt                              falloff_index_of_single_thirdbody_species_;

		OpenSMOKEVectorReactionTag falloff_reaction_type_;
		OpenSMOKEVectorDouble  a_falloff_;
		OpenSMOKEVectorDouble  b_falloff_;
		OpenSMOKEVectorDouble  c_falloff_;
		OpenSMOKEVectorDouble  d_falloff_;
		OpenSMOKEVectorDouble  e_falloff_;
		OpenSMOKEVectorDouble  logFcent_falloff_;

		OpenSMOKEVectorDouble lnA_cabr_inf_;				//!< frequency factors (log) for cabr reactions
		OpenSMOKEVectorDouble Beta_cabr_inf_;				//!< temperature exponents for cabr reactions
		OpenSMOKEVectorDouble E_over_R_cabr_inf_;			//!< activation temperatures for cabr reactions

		std::vector<OpenSMOKEVectorUnsignedInt>     cabr_indices_of_thirdbody_species_;
		std::vector<OpenSMOKEVectorDouble>          cabr_indices_of_thirdbody_efficiencies_;
		OpenSMOKEVectorInt                          cabr_index_of_single_thirdbody_species_;

		OpenSMOKEVectorReactionTag cabr_reaction_type_;
		OpenSMOKEVectorDouble  a_cabr_;
		OpenSMOKEVectorDouble  b_cabr_;
		OpenSMOKEVectorDouble  c_cabr_;
		OpenSMOKEVectorDouble  d_cabr_;
		OpenSMOKEVectorDouble  e_cabr_;
		OpenSMOKEVectorDouble  logFcent_cabr_;

		OpenSMOKEVectorDouble changeOfMoles_;		//!< list of change of moles

		StoichiometricMap* stoichiometry_;		//!< pointer to the stoichiometry

		bool arrhenius_kinetic_constants_must_be_recalculated_;
		bool nonconventional_kinetic_constants_must_be_recalculated_;
		bool reaction_h_and_s_must_be_recalculated_;
        bool isJacobianSparsityMapAvailable_;

		OpenSMOKEVectorDouble reaction_s_over_R_;
		OpenSMOKEVectorDouble reaction_h_over_RT_;
		OpenSMOKEVectorDouble kArrheniusModified_;
		OpenSMOKEVectorDouble kArrhenius_;
		OpenSMOKEVectorDouble uKeq_;
		OpenSMOKEVectorDouble kArrhenius_reversible_;
		OpenSMOKEVectorDouble kArrhenius_falloff_inf_;
		OpenSMOKEVectorDouble kArrhenius_cabr_inf_;

		OpenSMOKEVectorDouble forwardReactionRates_;
		OpenSMOKEVectorDouble reverseReactionRates_;
		OpenSMOKEVectorDouble netReactionRates_;

		double Patm_over_RT_;
		double log_Patm_over_RT_;

		OpenSMOKEVectorDouble correction_falloff_;						//!< correction factors for the falloff reactions
		OpenSMOKEVectorDouble correction_cabr_;							//!< correction factors for the cabr reactions

		ChebyshevPolynomialRateExpression* chebyshev_reactions_;		//!< pointer to the list of Chebyshev reactions
		PressureLogarithmicRateExpression* pressurelog_reactions_;		//!< pointer to the list of PLOG reactions
	//	Fit1RateExpression* fit1_reactions_;					//!< pointer to the list of FIT1 reactions
	//	JanevLangerRateExpression* janevlanger_reactions_;			//!< pointer to the list of JAN reactions
	//	LandauTellerRateExpression* landauteller_reactions_;			//!< pointer to the list of LT reactions

		OpenSMOKEVectorUnsignedInt isThermodynamicallyReversible_;		//!< vector containing the local index of thermodynamically reversible reactions
		OpenSMOKEVectorUnsignedInt isExplicitlyReversible_;			//!< vector containing the local index of explicitly reversible reactions

		OpenSMOKEVector<PhysicalConstants::TAG_REACTION> type_of_reaction_;
		OpenSMOKEVectorUnsignedInt local_family_index_;


		JacobianSparsityPatternMap<KineticsMap_CHEMKIN>* jacobian_sparsity_pattern_map_;
	};


}

#include "KineticsMap_CHEMKIN.hpp"

#endif /* OpenSMOKE_KineticsMap_CHEMKIN_H */
