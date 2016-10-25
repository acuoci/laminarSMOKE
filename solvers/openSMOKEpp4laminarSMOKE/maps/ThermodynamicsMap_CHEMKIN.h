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

#ifndef OpenSMOKE_ThermodynamicsMap_CHEMKIN_CHEMKIN_H
#define OpenSMOKE_ThermodynamicsMap_CHEMKIN_CHEMKIN_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include "ThermodynamicsMap.h"
#include "rapidxml.hpp"

namespace OpenSMOKE
{
	//!  A class to efficiently evaluate the thermodynamic properties
	/*!
	This class provides the tools to calculate in a very efficient way the thermodynamic properties.
	In order to ensure a good efficiency a map is created to store all the data
	depending on the temperature. In this way they are recalculated only if strictly needed, i.e. only
	if the temperature changes.
	*/

	template<typename map>
	class ThermodynamicsMap_CHEMKIN : public ThermodynamicsMap<map>
	{

	public:

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties (obsolete TOREMOVE)
		*@param nSpecies number of species
		*@param nPoints number of points in the map (for a scalar map the number of points is 1)
		*/
		ThermodynamicsMap_CHEMKIN(const unsigned int nSpecies, const unsigned int nPoints = 1);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param doc file in XML format
		*@param nPoints number of points in the map (for a scalar map the number of points is 1)
		*/
		ThermodynamicsMap_CHEMKIN(rapidxml::xml_document<>& doc, const unsigned int nPoints = 1);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param doc file in XML format
		*@param verbose activate or deactivate output
		*@param nPoints number of points in the map (for a scalar map the number of points is 1)
		*/
		ThermodynamicsMap_CHEMKIN(rapidxml::xml_document<>& doc, bool verbose, const unsigned int nPoints = 1);

		/**
		*@brief Copy constructor
		*@param rhs the object to be copied in the current object
		*/
		ThermodynamicsMap_CHEMKIN(const ThermodynamicsMap_CHEMKIN& rhs);

		/**
		*@brief Default destructor
		*/
		~ThermodynamicsMap_CHEMKIN(void);

		/**
		*@brief Sets the verbose output
		*/
		void SetVerboseOutput(const bool verbose_output) { verbose_output_ = verbose_output; }


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
		*@brief Returns the molecular weight of the mixture from the mole fractions
		*@param x mole fractions of species
		*@param MW the molecular weight in kg/kmol
		*/
		virtual inline void MolecularWeight_From_MoleFractions(map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Returns the molecular weight of the mixture from the mass fractions
		*@param y mass fractions of species
		*@param MW the molecular weight in kg/kmol
		*/
		virtual inline void MolecularWeight_From_MassFractions(map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& y);

		/**
		*@brief Calculates the mass fractions and the molecular weight from the mole fractions
		*/
		virtual inline void MassFractions_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& y, map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the mole fractions and the molecular weight from the mass fractions
		*/
		virtual inline void MoleFractions_From_MassFractions(OpenSMOKE::OpenSMOKEVectorDouble& x, map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& y);

		/**
		*@brief Calculates the mixture averaged specific heat (equivalent to CKCPBL in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param cpmix the mixture-averaged specific heat in J/kmol/K
		*/
		virtual void cpMolar_Mixture_From_MoleFractions(map& cpmix, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the mixture averaged enthalpy (equivalent to CKHBML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param hmix the mixture-averaged enthalpy in J/kmol
		*/
		virtual void hMolar_Mixture_From_MoleFractions(map& hmix, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the mixture averaged entropy (equivalent to CKSBML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param smix the mixture-averaged entropy in J/kmol/K
		*/
		void sMolar_Mixture_From_MoleFractions(map& smix, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the mixture averaged internal energy (equivalent to CKUBML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param umix mixture-averaged internal energy in J/kmol
		*/
		virtual void uMolar_Mixture_From_MoleFractions(map& umix, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the mixture averaged Gibb's free energy (equivalent to CKGBML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param gmix the mixture-averaged Gibb's free energy in J/kmol
		*/
		virtual void gMolar_Mixture_From_MoleFractions(map& gmix, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the mixture averaged Helmholtz free energy (equivalent to CKABML in CHEMKIN software)
		*@param x the vector of mole fractions for all the species
		*@param amix the mixture-averaged Helmholtz free energy in J/kmol
		*/
		virtual void aMolar_Mixture_From_MoleFractions(map& amix, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the standard specific heats of species at constant pressure (equivalent to CKCPOR in CHEMKIN software)
		*@return the specific heats at constant pressure in J/kmol/K
		*/
		virtual void cpMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& cp_species);

		/**
		*@brief Calculates the standard enthalpies of species (equivalent to CKHORT in CHEMKIN software)
		*@return the standard enthalpies in J/kmol
		*/
		virtual void hMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& h_species);

		/**
		*@brief Calculates the standard entropies of species (equivalent to CKSOR in CHEMKIN software)
		*@return the standard entropies in J/kmol/K
		*/
		virtual void sMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& s_species);

		/**
		*@brief Calculates the standard internal energies of species (equivalent to CKUML in CHEMKIN software)
		*@return the standard entropies in J/kmol
		*/
		virtual void uMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& u_species);

		/**
		*@brief Calculates the standard Gibb's free energies of species (equivalent to CKGML in CHEMKIN software)
		*@return the standard Gibb's free energies in J/kmol
		*/
		virtual void gMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& g_species);

		/**
		*@brief Calculates the standard Helmholtz's free energies of species (equivalent to CKAML in CHEMKIN software)
		*@return the standard Helmholtz's free energies in J/kmol
		*/
		virtual void aMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& a_species);

		/**
		*@brief Calculates the mixture averaged entropies of species (no equivalent in CHEMKIN software)
		*       The mixture averaged species entropies are not the same as the standard state values
		*       and must account for the appropriate pressure abd entropy-of-mixing terms
		*@return the mixture averaged entropies in J/kmol/K
		*/
		virtual void sMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& s_species, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the mixture averaged Gibb's free energies of species (no equivalent in CHEMKIN software)
		*       The mixture averaged Gibb's free energies are not the same as the standard state values
		*       and must account for the appropriate pressure abd entropy-of-mixing terms
		*@return the mixture averaged Gibb's free energies in J/kmol
		*/
		virtual void gMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& g_species, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Calculates the mixture averaged Helmholtz's free energies of species (no equivalent in CHEMKIN software)
		*       The mixture averaged Helmholtz's free energies are not the same as the standard state values
		*       and must account for the appropriate pressure abd entropy-of-mixing terms
		*@return the mixture averaged Helmholtz's free energies in J/kmol
		*/
		virtual void aMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& a_species, const OpenSMOKE::OpenSMOKEVectorDouble& x);

		/**
		*@brief Returns the normalized species enthalpies
		*/
		virtual const OpenSMOKE::OpenSMOKEVectorDouble& species_h_over_RT() { h_over_RT(); return species_h_over_RT_; }

		/**
		*@brief Returns the normalized species entropies
		*/
		virtual const OpenSMOKE::OpenSMOKEVectorDouble& species_s_over_R() { s_over_R(); return species_s_over_R_; }

		/**
		*@brief Import the coefficients from a file in ASCII format (obsolete, TOREMOVE)
		*/
		virtual void ImportCoefficientsFromASCIIFile(std::ifstream& fInput);

		/**
		*@brief Import the coefficients from a file in XML format
		*/
		virtual void ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc);

		/**
		*@brief Import the species from a file in XML format
		*/
		virtual void ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc);

		/**
		*@brief TODO (This is just a test subroutine!)
		*/
		void Test(const int nLoops, const map& T, int* index);

		/**
		*@brief Calculates the normalized specific heats for the species
		*/
		inline void cp_over_R();

		void DerivativesOfConcentrationsWithRespectToMassFractions(const double cTot, const double MW, const OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEMatrixDouble* dc_over_omega);


		const Eigen::MatrixXd& atomic_composition() const { return this->atomic_composition_; }
		const std::vector<std::string>& elements() const { return this->elements_; }


		/**
		*@brief Calculates the temperature from the molar enthalpy
		*@param H molar enthalpy in J/kmol
		*@param P_Pa pressure in Pascal
		*@param x molar fraction
		*@param TFirstGuess first guess for the temperature
		*/
		double GetTemperatureFromEnthalpyAndMoleFractions(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double TFirstGuess);

		/**
		*@brief Copies the data from another thermodynamic map (used by copy constructors)
		*/
		void CopyFromMap(const ThermodynamicsMap_CHEMKIN& rhs);

		bool FindTheRightInterval(const double H, const OpenSMOKEVectorDouble& x, const double TFirstGuess,
			double& TA, double&TB, double& HA, double& HB,
			const double Diff_Temperature);
		/*
		double DivideTheInterval(	const double H, const OpenSMOKEVectorDouble& x,
		const double diff_temperature_to_switch_to_newtons,
		const double TA_, const double TB_, const double HA_, const double HB_);
		double DivideTheIntervalBy2(	const double H, const OpenSMOKEVectorDouble& x,
		const double diff_temperature_to_switch_to_newtons,
		const double TA_, const double TB_, const double HA_, const double HB_);


		double NewtonMethod(	const double H, const OpenSMOKEVectorDouble& x,
		const double T0);
		*/

		double Brent(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double x1, const double x2, const double tol);
		double Ridder(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double x1, const double x2, const double xacc);
		double NewtonRaphson(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double x1, const double x2, const double xacc);

		double Function(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double T);
		void Function(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double T, double& f, double& df);

		void NASA_LowT(const unsigned int i, double* coefficients) const;
		void NASA_HighT(const unsigned int i, double* coefficients) const;
		void NASA_Coefficients(const unsigned int i, double* coefficients) const;
		void NASA_Coefficients(double* coefficients) const;

	protected:

		/**
		*@brief Allocates the memory
		*/
		void MemoryAllocation();

		/**
		*@brief Sets the thermodynamic coefficients for each species in the map
		*/
		void SetCoefficients(const unsigned int i, const double* coefficients);

		/**
		*@brief Calculates the normalized specific enthalpies for the species
		*/
		inline void h_over_RT();

		/**
		*@brief Calculates the normalized entropies for the species
		*/
		inline void s_over_R();

		/**
		*@brief Calculates the free Gibbs' energies for the species
		*/
		inline void g_over_RT();

	private:

		OpenSMOKEVectorDouble species_cp_over_R_;
		OpenSMOKEVectorDouble species_h_over_RT_;
		OpenSMOKEVectorDouble species_g_over_RT_;
		OpenSMOKEVectorDouble species_s_over_R_;

		double*	Cp_LT;	/**< correlation coefficients, specific heat, low temperature region */
		double*	Cp_HT;	/**< correlation coefficients, specific heat, high temperature region */
		double*	DH_LT;	/**< correlation coefficients, enthalpy, low temperature region */
		double*	DH_HT;	/**< correlation coefficients, enthalpy, high temperature region */
		double*	DS_LT;	/**< correlation coefficients, entropy, low temperature region */
		double*	DS_HT;	/**< correlation coefficients, entropy, high temperature region */

		double* TL;			/**< low  temperature limit */
		double* TM;			/**< mean temperature limit */
		double* TH;			/**< high temperature limit */

		bool cp_must_be_recalculated_;
		bool h_must_be_recalculated_;
		bool s_must_be_recalculated_;

		bool verbose_output_;			/**< Print video info  */
	};
}

#include "ThermodynamicsMap_CHEMKIN.hpp"

#endif /* OpenSMOKE_ThermodynamicsMap_CHEMKIN_H */
