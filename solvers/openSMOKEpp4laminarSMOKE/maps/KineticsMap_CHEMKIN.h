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
|	License                                                           |
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
#include "TransportPropertiesMap_CHEMKIN.h"
#include "StoichiometricMap.h"
#include "JacobianSparsityPatternMap.h"
#include "rapidxml.hpp"

namespace OpenSMOKE
{
	enum EnergyEquationType { CONSTANT_VOLUME_SYSTEM, CONSTANT_PRESSURE_SYSTEM, PLUGFLOW_SYSTEM };

	typedef std::vector<PhysicalConstants::TAG_REACTION> VectorReactionTags;

	//!  A class to efficiently evaluate the reaction and formation rates, to be used in production codes
	/*!
	This class provides the tools to calculate in a very efficient way the reaction rates and the
	formation rates. In order to ensure a good efficiency a map is created to store all the data
	depending on the temperature. Inthis way they are recalculated only if strictly needed, i.e. only
	if the temperature changes
	*/

	class KineticsMap_CHEMKIN : public KineticsMap
	{

	public:

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties (obsolete, TOREMOVE)
		*@param thermo the thermodynamic map
		*@param nSpecies number of species
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN& thermo, const unsigned int nSpecies, const unsigned int nPoints = 1);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN& thermo, rapidxml::xml_document<>& doc, const unsigned int nPoints = 1);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file
		*@param verbose activate output
		*@param nPoints number op points in the map (for a scalar map the number of oints is 1)
		*/
		KineticsMap_CHEMKIN(ThermodynamicsMap_CHEMKIN& thermo, rapidxml::xml_document<>& doc, bool verbose, const unsigned int nPoints = 1);

		/**
		*@brief Copy constructor
		*       IMPORTANT: Please consider that the thermodynamic map to which the current kinetic map (this) is linked, is
		*       still the same thermodynamic map to which the rhs map is linked! If you want to link the current
		*       kinetic map to a different thermodynamic map you should use a different copy constructor (see below)
		*@param rhs the object to be copied in the current object
		*/
		KineticsMap_CHEMKIN(const KineticsMap_CHEMKIN& rhs);

		/**
		*@brief Copy constructor
		*       IMPORTANT: no checks are done about the consistency between the thermodynamic map linked to the kinetic map
		*       (rhs.thermodynamics_) and the provided map (thermo). It is up to the user to be sure about the consistency.
		*@param rhs the object to be copied in the current object
		*/
		KineticsMap_CHEMKIN(const KineticsMap_CHEMKIN& rhs, ThermodynamicsMap_CHEMKIN& thermo);

		/**
		*@brief Default destructor
		*/
		~KineticsMap_CHEMKIN();

		/**
		*@brief Return the thermodynamics map associated to the current kinetics map
		*/
		ThermodynamicsMap_CHEMKIN& thermodynamics() const { return thermodynamics_; }

		/**
		*@brief Sets the verbose output
		*/
		void SetVerboseOutput(const bool verbose_output) { verbose_output_ = verbose_output; }

		/**
		*@brief Set the temperature at which the properties have to be evaluated
		*@param T the temperature value in K
		*/
		virtual void SetTemperature(const double& T);

		/**
		*@brief Set the pressure at which the properties have to be evaluated
		*@param P the pressure value in Pa
		*/
		virtual void SetPressure(const double& P);

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
		*@brief Calculates the formation rates for all the species in the kinetic mechanism (in kmol/m3/s)
		*/
		void FormationRates(double* R);

		/**
		*@brief Calculates the heat release on the basis of formation rates
		*@param R formation rates of species (in kmol/m3/s)
		*returns the heat release rate (in J/m3/s)
		*/
		double HeatRelease(const double* R);

		/**
		*@brief Calculates the production and the destruction rates for all the species in the kinetic mechanism [kmol/m3/s]
		*       The contributions are calculated on the basis of net reaction rates
		*/
		void ProductionAndDestructionRates(double* P, double* D);

		/**
		*@brief Calculates the production and the destruction rates for all the species in the kinetic mechanism [kmol/m3/s]
		*       The contributions are calculated on the basis of separate forward and backward reaction rates
		*/
		void ProductionAndDestructionRatesGross(double* P, double* D);

		/**
		*@brief Returns the net reaction rates (in kmol/m3/s)
		*/
		const std::vector<double>& GiveMeReactionRates();

		/**
		*@brief Return the net reaction rates (in kmol/m3/s)
		*/
		void GiveMeReactionRates(double* r);

		/**
		*@brief Returns the forward reaction rates for all the reactions in the kinetic scheme [kmol/m3/s]
		*/
		void GetForwardReactionRates(double* r);

		/**
		*@brief Returns the backward reaction rates for all the reactions in the kinetic scheme [kmol/m3/s]
		If a reaction is irreversible, it returns zero
		*/
		void GetBackwardReactionRates(double* r);

		/**
		*@brief Calculates the reaction rates for all the reactions in the kinetic scheme
		*@param c concentrations of species (in kmol/m3)
		*/
		void ReactionRates(const double* c);

		/**
		*@brief Calculates the reaction rates for all the reactions in the kinetic scheme (if total concentration is available)
		*@param c concentrations of species (in kmol/m3)
		*@param cTot total concentration (in kmol/m3)
		*/
		void ReactionRates(const double* c, const double cTot);

		/**
		*@brief Returns the indices of the reversible reactions
		*/
		const std::vector<unsigned int>& IndicesOfReversibleReactions() const { return indices_of_reversible_reactions__; }

		/**
		*@brief Calculates the reaction enthalpies and entropies (to be used for the kinetic constants)
		*/
		void ReactionEnthalpiesAndEntropies();

		/**
		*@brief Calculates the kinetic constants
		*/
		void KineticConstants();

		/**
		*@brief Return the frequency factor of a single reaction [kmol, m, s]
		*@param j index of reaction (starting from zero)
		*/
		double A(const unsigned int j) { return sign_lnA__[j] * std::exp(lnA__[j]); }

		/**
		*@brief Return the temperature exponent a single reaction
		*@param j index of reaction (starting from zero)
		*/
		double Beta(const unsigned int j) { return Beta__[j]; }

		/**
		*@brief Return the activation temperature of a single reaction [K]
		*@param j index of reaction (starting from zero)
		*/
		double E_over_R(const unsigned int j) { return E_over_R__[j]; }

		/**
		*@brief Return the third-body efficiency of a species in a given reaction
		*@param j index of reaction (starting from zero)
		*@param j index of species (starting from zero)
		*@return the efficiency (if the species does not exist in the list it returns 1.)
		*/
		double ThirdBody(const unsigned int j, const unsigned int k);

		/**
		*@brief Sets the frequency factor of a single reaction [kmol, m, s]
		*@param j index of reaction (starting from zero)
		*@param A the new frequency factor
		*/
		void Set_A(const unsigned int j, const double A) { lnA__[j] = std::log(A); }

		/**
		*@brief Sets the temperature exponent a single reaction
		*@param j index of reaction (starting from zero)
		*@param Beta the new index of reaction
		*/
		void Set_Beta(const unsigned int j, const double Beta) { Beta__[j] = Beta; }

		/**
		*@brief Sets the activation temperature of a single reaction [K]
		*@param j index of reaction (starting from zero)
		*@param E_over_R the new activation temperature
		*/
		void Set_E_over_R(const unsigned int j, const double E_over_R) { E_over_R__[j] = E_over_R; }

		/**
		*@brief Sets the third body efficiency of a species in a give reaction
		*@param j index of reaction (starting from zero)
		*@param k index of species (starting from zero)
		*@param efficiency the new efficiency coefficient
		*/
		void Set_ThirdBody(const unsigned int j, const unsigned int k, const double efficiency);

		/**
		*@brief Return the frequency factor of the high pressure limit of a single reaction [kmol, m, s]
		*@param j index of reaction (starting from zero)
		*/
		double A_falloff_inf(const unsigned int j) { return  std::exp(lnA_falloff_inf__[j]); }

		/**
		*@brief Return the temperature exponent of the high pressure limit of a single reaction
		*@param j index of reaction (starting from zero)
		*/
		double Beta_falloff_inf(const unsigned int j) { return Beta_falloff_inf__[j]; }

		/**
		*@brief Return the activation temperature of the high pressure limit of a single reaction [K]
		*@param j index of reaction (starting from zero)
		*/
		double E_over_R_falloff_inf(const unsigned int j) { return E_over_R_falloff_inf__[j]; }

		/**
		*@brief Sets the frequency factor for the high pressure limit of a single reaction [kmol, m, s]
		*@param j index of reaction (starting from zero)
		*@param A the new frequency factor
		*/
		void Set_A_falloff_inf(const unsigned int j, const double A_falloff_inf) { lnA_falloff_inf__[j] = std::log(A_falloff_inf); }

		/**
		*@brief Sets the temperature exponent of the high pressure limit of a single reaction
		*@param j index of reaction (starting from zero)
		*@param Beta the new index of reaction
		*/
		void Set_Beta_falloff_inf(const unsigned int j, const double Beta_falloff_inf) { Beta_falloff_inf__[j] = Beta_falloff_inf; }

		/**
		*@brief Sets the activation temperature of the high pressure limit of a single reaction [K]
		*@param j index of reaction (starting from zero)
		*@param E_over_R the new activation temperature
		*/
		void Set_E_over_R_falloff_inf(const unsigned int j, const double E_over_R_falloff_inf) { E_over_R_falloff_inf__[j] = E_over_R_falloff_inf; }

		/**
		*@brief Returns the total number of reversible reactions
		*/
		unsigned int NumberOfReversibleReactions() const { return number_of_reversible_reactions_; }

		/**
		*@brief For each reaction returns 1 if thermodynamically reversible, 0 otherwise
		*/
		const std::vector<unsigned int>& IsThermodynamicallyReversible() const { return isThermodynamicallyReversible__; };

		/**
		*@brief For each reaction returns 1 if explicitly reversible, 0 otherwise
		*/
		const std::vector<unsigned int>& IsExplicitlyReversible() const { return isExplicitlyReversible__; };

		/**
		*@brief Returns the total number of reactions with third body effects
		*/
		unsigned int NumberOfThirdBodyReactions() const { return number_of_thirdbody_reactions_; }

		/**
		*@brief Returns the indices of reactions with third-body effects
		*/
		const std::vector<unsigned int>& IndicesOfThirdbodyReactions() const { return indices_of_thirdbody_reactions__; };

		/**
		*@brief Returns the indices of of third body species for each third body reaction
		*/
		const std::vector< std::vector<unsigned int> >& IndicesOfThirdbodySpecies() const { return indices_of_thirdbody_species__; };

		/**
		*@brief Returns the efficiencies of of third body species for each third body reaction
		*/
		const std::vector< std::vector<double> >&	IndicesOfThirdbodyEfficiencies() const { return indices_of_thirdbody_efficiencies__; }

		/**
		*@brief Return the third-body global efficiency for each third-body reaction
		*/
		std::vector<double>& M() { return Meff__; }

		/**
		*@brief Returns the total number of fall off reactions
		*/
		unsigned int NumberOfFallOffReactions() const { return number_of_falloff_reactions_; }

		/**
		*@brief Returns the indices of fall-off reactions
		*/
		const std::vector<unsigned int>& IndicesOfFalloffReactions() const { return indices_of_falloff_reactions__; };

		/**
		*@brief Returns the total number of Chemically Activated Bimolecular Reactions (CABR)
		*/
		unsigned int NumberOfCABRReactions() const { return number_of_cabr_reactions_; }

		/**
		*@brief Returns the indices of Chemically Activated Bimolecular Reactions (CABR)
		*/
		const std::vector<unsigned int>& IndicesOfCabrReactions() const { return indices_of_cabr_reactions__; };

		/**
		*@brief Returns the modified Arrhenius kinetic constants (in kmol, m, s)
		*/
		const std::vector<double>& KArrheniusModified() const { return kArrheniusModified__; }

		/**
		*@brief Returns the Arrhenius kinetic constants (in kmol, m, s)
		*/
		const std::vector<double>& KArrhenius() const { return kArrhenius__; }

		/**
		*@brief Returns the net reaction rates (in kmol/m3/s) for each reaction
		*/
		std::vector<double>& NetReactionRates() { return netReactionRates__; }

		/**
		*@brief Returns the correction due to fall-off effects for a given fall-off reaction
		*@param local_k index of reaction (among the fall-off reactions)
		*@param cTot total concentration (in kmol/m3)
		*@param concentrations of species (in kmol/m3)
		*returns the fall-off correction
		*/
		double FallOffReactionsCorrection(const unsigned int local_k, const double cTot, const double* c);


		const std::vector<unsigned int>& FallOffIndexOfSingleThirdbodySpecies() const { return falloff_index_of_single_thirdbody_species__; }
		const std::vector< std::vector<unsigned int> >& FallOffIndicesOfThirdbodySpecies() const { return falloff_indices_of_thirdbody_species__; }

		ChebyshevPolynomialRateExpression& chebyshev_reactions(const unsigned int j) const { return chebyshev_reactions_[j]; }
		unsigned int number_of_chebyshev_reactions() const { return number_of_chebyshev_reactions_; }
		const std::vector<unsigned int>& indices_of_chebyshev_reactions() const { return indices_of_chebyshev_reactions__; }

	public:

		/**
		*@brief Calculates the derivatives of formation rates with respect to concentrations
		*@param c current concentrations (in kmol/m3)
		*@param dR_over_dC the calculated derivatives (in 1/s) of formation rates with respect to concentrations
		*/
		void DerivativesOfFormationRates(const double* c, Eigen::MatrixXd* dR_over_dC);

		/**
		*@brief Calculates the derivatives of formation rates with respect to mass fractions
		*@param c current concentrations (in kmol/m3)
		*@param omega current mass fractions
		*@param dR_over_domega the calculated derivatives (in kmol/m3/s) of formation rates with respect to mass fractions
		*/
		void DerivativesOfFormationRates(const double* c, const double* omega, Eigen::MatrixXd* dR_over_domega);

		/**
		*@brief Returns the stoichiometric map
		*/
		StoichiometricMap& stoichiometry() { return *stoichiometry_; }

		/**
		*@brief Returns the sparsity pattern of Jacobian matrix
		*/
		JacobianSparsityPatternMap<KineticsMap_CHEMKIN>* jacobian_sparsity_pattern_map() { return jacobian_sparsity_pattern_map_; };

	public:

		/**
		*@brief Calculates the kinetic constants of the reverse reactions
		*/
		void FittedReverseKineticConstants(const double* x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters, const bool only_reversible);

		/**
		*@brief Calculates the kinetic constants of the reverse reactions
		*/
		void FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters);

		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k, const double* c_bath, const std::vector<double> list_of_temperatures, const double conversion_forward = 1., const double conversion_backward = 1.);

		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k, const double T, const double P_Pa, const double* c);

		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k);

		/**
		*@brief Write collision rate constants for bimolecular reactions
		*@param transport map describing the transport properties
		*@param fOut output stream where to write the results of the analysis
		*@param reaction_names list of reaction strings (for a better output)
		*@param list_of_temperatures list of temperatures (in K) at which the analysis is carried out
		*/
		void WriteCollisionRateConstantForBimolecularReactions(
			TransportPropertiesMap_CHEMKIN& transport,
			std::ostream& fOut, const std::vector<std::string>& reaction_names,
			const std::vector<double>& list_of_temperatures);


	public:	// Functions for builing the sparsity patterns

		void WeakThirdBodyConcentrationEfficiencies(std::vector<unsigned int>& reaction, std::vector<unsigned int>& species);
		void WeakFallOffConcentrationEfficiencies(std::vector<unsigned int>& reaction, std::vector<unsigned int>& species);
		void WeakCABRConcentrationEfficiencies(std::vector<unsigned int>& reaction, std::vector<unsigned int>& species);
		void StrongConcentrationEfficiencies(std::vector<unsigned int>& reaction);


	public:	// Rate of Production Analysis (ROPA) utilities

		void RateOfProductionAnalysis(const bool iNormalize) const;
		void RateOfProductionAnalysis(std::ostream& fout) const;
		void RateOfProductionAnalysis(ROPA_Data& ropa) const;
		void RateOfProductionAnalysis(ROPA_Data& ropa, const double* rf, const double* rb) const;


	public:	// Sensitivity Analysis (SA) utilities

		void DerivativesOfReactionRatesWithRespectToKineticParameters(const PhysicalConstants::sensitivity_type type, unsigned int jReaction, const double* c, double& parameter);
		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const unsigned int k, const double* c, double* Jalfa, double& parameter);
		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const double* c, const double* mole_fractions, double* Jalfa, double& JT, double& parameter);
		void SensitivityWithRespectKineticParameter(const PhysicalConstants::sensitivity_type type, const EnergyEquationType energy_type, const unsigned int k, const double* c, const double* mole_fractions, double* Jalfa, double& JT, double& Jrho, double& parameter);

	private:

		/**
		*@brief Calculates the efficiencies for threebody reactions
		*/
		void ThirdBodyReactions(const double cTot, const double* c);

		/**
		*@brief Calculates the corrections for the falloff reactions
		*/
		void FallOffReactions(const double cTot, const double* c);

		/**
		*@brief Calculates the corrections for the cabr reactions
		*/
		void ChemicallyActivatedBimolecularReactions(const double cTot, const double* c);

		/**
		*@brief Calculates the kinetic constants for extended pressure log reactions
		*/
		void ExtendedPressureLogReactions(const double cTot, const double* c);

		/**
		*@brief Calculates the kinetic constants for extended falloff reactions
		*/
		void ExtendedFallOffReactions(const double cTot, const double* c);

		/**
		*@brief Copies the data from another kinetic map (used by copy constructors)
		*/
		void CopyFromMap(const KineticsMap_CHEMKIN& rhs);

	private:

		// TODO
		void FallOffReactions(const unsigned int k, const double cTot, const double* c, double &F, double &dF_over_dA0, double &dF_over_dAInf);
		void ChemicallyActivatedBimolecularReactions(const unsigned int k, const double cTot, const double* c, double &F, double &dF_over_dA0, double &dF_over_dAInf);
		void Derivatives(const double* c, Eigen::MatrixXd* derivatives, const bool constant_density = false);
		void Derivatives(const double* c, const double* omega, Eigen::MatrixXd* derivatives);

	private:

		ThermodynamicsMap_CHEMKIN & thermodynamics_;					//!< reference to the thermodynamics

		std::vector<unsigned int> indices_of_irreversible_reactions__;			//!< indices of irreversible reactions
		std::vector<unsigned int> indices_of_reversible_reactions__;			//!< indices of reversible reactions
		std::vector<unsigned int> indices_of_thermodynamic_reversible_reactions__;	//!< indices of reversible (thermodynamic) reactions
		std::vector<unsigned int> indices_of_explicitly_reversible_reactions__;		//!< indices of reversible (explicit) reactions
		std::vector<unsigned int> indices_of_thirdbody_reactions__;			//!< indices of three-body reactions
		std::vector<unsigned int> indices_of_falloff_reactions__;			//!< indices of falloff reactions
		std::vector<unsigned int> indices_of_extendedfalloff_reactions__;		//!< indices of extended falloff reactions
		std::vector<unsigned int> indices_of_cabr_reactions__;				//!< indices of cabr reactions
		std::vector<unsigned int> indices_of_chebyshev_reactions__;			//!< indices of chebyshev reactions
		std::vector<unsigned int> indices_of_pressurelog_reactions__;			//!< indices of pressurelog (PLOG) reactions
		std::vector<unsigned int> indices_of_extendedpressurelog_reactions__;		//!< indices of extended pressurelog (EXTPLOG) reactions
		std::vector<unsigned int> indices_of_fit1_reactions__;				//!< indices of FIT1 reactions
		std::vector<unsigned int> indices_of_janevlanger_reactions__;			//!< indices of JAN reactions
		std::vector<unsigned int> indices_of_landauteller_reactions__;			//!< indices of LT reactions

		std::vector<unsigned int> indices_of_falloff_lindemann_reactions__;		//!< indices of falloff (Lindemann) reactions
		std::vector<unsigned int> indices_of_cabr_lindemann_reactions__;		//!< indices of cabr (Lindemann) reactions
		std::vector<unsigned int> indices_of_falloff_troe_reactions__;			//!< indices of falloff (Troe) reactions
		std::vector<unsigned int> indices_of_cabr_troe_reactions__;			//!< indices of cabr (Troe) reactions
		std::vector<unsigned int> indices_of_falloff_sri_reactions__;			//!< indices of falloff (SRI) reactions
		std::vector<unsigned int> indices_of_cabr_sri_reactions__;			//!< indices of cabr (Troe) reactions

		unsigned int number_of_irreversible_reactions_;
		unsigned int number_of_reversible_reactions_;
		unsigned int number_of_thermodynamic_reversible_reactions_;
		unsigned int number_of_explicitly_reversible_reactions_;
		unsigned int number_of_thirdbody_reactions_;
		unsigned int number_of_falloff_reactions_;
		unsigned int number_of_extendedfalloff_reactions_;
		unsigned int number_of_cabr_reactions_;
		unsigned int number_of_chebyshev_reactions_;
		unsigned int number_of_pressurelog_reactions_;
		unsigned int number_of_extendedpressurelog_reactions_;
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

		std::vector<double> lnA__;				//!< frequency factors (log)
		std::vector<double> Beta__;				//!< temperature exponents
		std::vector<double> E_over_R__;				//!< activation temperatures
		std::vector<int> negative_lnA__;			//!< list of reactions with negative frequency factor (1-index based)
		std::vector<int> sign_lnA__;				//!< sign of frequency factors of reactions (+1 or -1)

		std::vector<double> lnA_reversible__;			//!< frequency factors (log) for explicitly reversible reactions
		std::vector<double> Beta_reversible__;			//!< temperature exponents for explicitly reversible reactions
		std::vector<double> E_over_R_reversible__;		//!< activation temperatures for explicitly reversible reactions

		std::vector<double> Meff__;																//!< threebody efficiencies
		std::vector< std::vector<unsigned int> >		indices_of_thirdbody_species__;		//!< indices of threebody species
		std::vector< std::vector<double> >			indices_of_thirdbody_efficiencies__;	//!< efficiencies of threebody species

		std::vector<double> lnA_falloff_inf__;			//!< frequency factors (log) for falloff reactions
		std::vector<double> Beta_falloff_inf__;			//!< temperature exponents for falloff reactions
		std::vector<double> E_over_R_falloff_inf__;		//!< activation temperatures for falloff reactions

		std::vector< std::vector<unsigned int> >	falloff_indices_of_thirdbody_species__;
		std::vector< std::vector<double> >		falloff_indices_of_thirdbody_efficiencies__;
		std::vector<unsigned int>			falloff_index_of_single_thirdbody_species__;

		VectorReactionTags falloff_reaction_type__;
		std::vector<double>  a_falloff__;
		std::vector<double>  b_falloff__;
		std::vector<double>  c_falloff__;
		std::vector<double>  d_falloff__;
		std::vector<double>  e_falloff__;
		std::vector<double>  logFcent_falloff__;

		std::vector<double> lnA_cabr_inf__;		//!< frequency factors (log) for cabr reactions
		std::vector<double> Beta_cabr_inf__;		//!< temperature exponents for cabr reactions
		std::vector<double> E_over_R_cabr_inf__;	//!< activation temperatures for cabr reactions

		std::vector< std::vector<unsigned int> >	cabr_indices_of_thirdbody_species__;
		std::vector< std::vector<double> >          	cabr_indices_of_thirdbody_efficiencies__;
		std::vector<unsigned int>                   	cabr_index_of_single_thirdbody_species__;

		VectorReactionTags cabr_reaction_type__;
		std::vector<double>  a_cabr__;
		std::vector<double>  b_cabr__;
		std::vector<double>  c_cabr__;
		std::vector<double>  d_cabr__;
		std::vector<double>  e_cabr__;
		std::vector<double>  logFcent_cabr__;

		std::vector<double> changeOfMoles__;		//!< list of change of moles

		StoichiometricMap* stoichiometry_;		//!< pointer to the stoichiometry

		bool arrhenius_kinetic_constants_must_be_recalculated_;
		bool nonconventional_kinetic_constants_must_be_recalculated_;
		bool reaction_h_and_s_must_be_recalculated_;
		bool isJacobianSparsityMapAvailable_;

		std::vector<double> reaction_s_over_R__;
		std::vector<double> reaction_h_over_RT__;
		std::vector<double> kArrheniusModified__;
		std::vector<double> kArrhenius__;
		std::vector<double> uKeq__;
		std::vector<double> kArrhenius_reversible__;
		std::vector<double> kArrhenius_falloff_inf__;
		std::vector<double> kArrhenius_cabr_inf__;

		std::vector<double> forwardReactionRates__;
		std::vector<double> reverseReactionRates__;
		std::vector<double> netReactionRates__;

		double Patm_over_RT_;
		double log_Patm_over_RT_;

		std::vector<double> correction_falloff__;						//!< correction factors for the falloff reactions
		std::vector<double> correction_cabr__;							//!< correction factors for the cabr reactions

		ChebyshevPolynomialRateExpression* chebyshev_reactions_;				//!< pointer to the list of Chebyshev reactions
		PressureLogarithmicRateExpression* pressurelog_reactions_;				//!< pointer to the list of PLOG reactions
		ExtendedPressureLogarithmicRateExpression* extendedpressurelog_reactions_;		//!< pointer to the list of PLOGMX/PLOGSP reactions
		ExtendedFallOff* extendedfalloff_reactions_;						//!< pointer to the list of extended falloff reactions

		std::vector<unsigned int> isThermodynamicallyReversible__;		//!< vector containing the local index of thermodynamically reversible reactions
		std::vector<unsigned int> isExplicitlyReversible__;			//!< vector containing the local index of explicitly reversible reactions

		VectorReactionTags type_of_reaction__;
		std::vector<unsigned int> local_family_index__;

		JacobianSparsityPatternMap<KineticsMap_CHEMKIN>* jacobian_sparsity_pattern_map_;
	};

}

#include "KineticsMap_CHEMKIN.hpp"

#endif /* OpenSMOKE_KineticsMap_CHEMKIN_H */
