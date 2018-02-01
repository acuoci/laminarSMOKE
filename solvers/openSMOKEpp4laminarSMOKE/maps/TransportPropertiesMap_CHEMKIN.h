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

#ifndef OpenSMOKE_TransportPropertiesMap_CHEMKIN_H
#define OpenSMOKE_TransportPropertiesMap_CHEMKIN_H

#include "TransportPropertiesMap.h"
#include "rapidxml.hpp"

namespace OpenSMOKE
{
	//!  A class to efficiently evaluate the transport properties
	/*!
	This class provides the tools to calculate in a very efficient way the transport properties.
	In order to ensure a good efficiency a map is created to store all the data
	depending on the temperature. In this way they are recalculated only if strictly needed, i.e. only
	if the temperature changes.
	*/

	class TransportPropertiesMap_CHEMKIN : public TransportPropertiesMap
	{
	public:

		/**
		*@brief Creates a transport map for the evaluation of transport properties (obsolete TOREMOVE)
		*@param nSpecies number of species
		*/
		TransportPropertiesMap_CHEMKIN(const unsigned int nSpecies);

		/**
		*@brief Creates a transport map for the evaluation of transport properties
		*@param doc file in XML format
		*/
		TransportPropertiesMap_CHEMKIN(rapidxml::xml_document<>& doc);

		/**
		*@brief Copy constructor
		*@param rhs the object to be copied in the current object
		*/
		TransportPropertiesMap_CHEMKIN(const TransportPropertiesMap_CHEMKIN& rhs);

		/**
		*@brief Default destructor
		*/
		~TransportPropertiesMap_CHEMKIN();

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
		*@brief Sets the transport coefficient for the requested species
		*/
		virtual void SetCoefficients(const unsigned int i, const double* coefficients);

		/**
		*@brief Import the coefficients from a file in ASCII format (obsolete, TOREMOVE)
		*/
		virtual void ImportCoefficientsFromASCIIFile(std::istream& fInput);

		/**
		*@brief Import the coefficients from a file in ASCII format (obsolete, TOREMOVE)
		*@param fInput input stream
		*/
		virtual void ImportLennardJonesCoefficientsFromASCIIFile(std::istream& fInput);

		/**
		*@brief Import the coefficients from a file in XML format
		*/
		virtual void ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc);

		/**
		*@brief Import the species bundling coefficients from a file in XML format
		*/
		virtual void ImportSpeciesBundlingFromXMLFile(rapidxml::xml_document<>& doc, const double epsilon);

		/**
		*@brief Import the species from a file in XML format
		*/
		virtual void ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc);

		/**
		*@brief Import the viscosity model from a file in XML format
		*/
		virtual void ImportViscosityModelFromXMLFile(rapidxml::xml_document<>& doc);

		/**
		*@brief Returns the reduced collision integral Omega11* according to
		        the regression fit reported in Chen et al., Comb. and Flame 186, p. 208 (2017)
		*@param TStar reduced temperature (dimensionless) calculated as TStar=kb*T/eps
		*/
		double Omega11(const double TStar);

		/**
		*@brief Returns the reduced collision integral Omega11* according to
		        a 2D interpolation table
		*@param TStar reduced temperature (dimensionless) calculated as TStar=kb*T/eps
		*@param DStar reduced dipole moment (see CHEMKIN Theory Guide)
		*/
		double Omega11(const double TStar, const double DStar);

		/**
		*@brief Returns the mass (in kg) of species (0-based index)
		*/
		const std::vector<double>& mu() const { return mu_; }

		/**
		*@brief Returns the collision diameters (in m) of species (0-based index)
		*/
		const std::vector<double>& sigma() const { return sigma_; }

		/**
		*@brief Returns the scaled well depths (in K) of species (0-based index)
		*/
		const std::vector<double>& epsilon_over_kb() const { return epsilon_over_kb_; }

	protected:

		/**
		*@brief TODO (This is just a test subroutine!)
		*/
		void Test(const int nLoops, const double& T, int* index);

	public:

		/**
		*@brief Combines the species thermal conductivities to calculate the mixture thermal conductivity
		*/
		virtual double lambdaMix(const double* x);

		/**
		*@brief Combines the species dynamic viscosities to calculate the mixture dynamic viscosity
		*/
		virtual double etaMix(const double* x);

		/**
		*@brief Combines the species mass diffusion coefficients to calculate the mixture mass diffusion coefficients
		*/
		virtual void gammaMix(double* gammamix, const double* x);

		/**
		*@brief Combines the species mass diffusion coefficients to calculate the mixture mass diffusion coefficients using the bundling algorithm
		*/
		virtual void bundling_gammaMix(double* gammamix, const double* x);

		/**
		*@brief Combines the species thermal diffusion coefficients to calculate the mixture thermal diffusion coefficients
		*/
		virtual void tetaMix(double* tetamix, const double* x);

		/**
		*@brief Combines the planck mean absorption coefficients of relevant species, according to their mole fractions
		Returns the Planck mean absorption coefficient of the mixture in [1/m]
		*/
		virtual double kPlanckMix(const double* moleFractions);

		/**
		*@brief Calculates the thermal conductivities for all the species
		*/
		inline virtual void lambda();

		/**
		*@brief Calculates the dynamic viscosities for all the species
		*/
		inline virtual void eta();

		/**
		*@brief Calculates the mass diffusion coefficients for all the species
		*/
		inline virtual void gamma();

		/**
		*@brief Calculates the mass diffusion coefficients for all the species using the bundling algorithm
		*/
		inline virtual void bundling_gamma();

		/**
		*@brief Calculates the thermal diffusion coefficients for all the species
		*/
		inline virtual void teta();

		/**
		*@brief Return the vector of species (1-index based) for which the Soret effect is active
		*/
		const std::vector<unsigned int>& iThermalDiffusionRatios() const { return iThermalDiffusionRatios_; }

		/**
		*@brief Calculates the collision rate constant for bimolecular reactions [m3/kmol/s]
		*/
		double kCollision(const unsigned int i, const unsigned int k, const double T);

	private:

		/**
		*@brief Allocates the memory
		*/
		void MemoryAllocation();

		/**
		*@brief Precalculates useful data
		*/
		void CompleteInitialization();

		/**
		*@brief Copies the data from another transport map (used by copy constructors)
		*/
		void CopyFromMap(const TransportPropertiesMap_CHEMKIN& rhs);

	private:

		PhysicalConstants::OpenSMOKE_GasMixture_Viscosity_Model viscosity_model;

		double* M;				//!< molecular weights
		double* fittingLambda;			//!< fitting coefficients for the thermal conductivities
		double* fittingEta;			//!< fitting coefficients for the dynamic viscosities
		double* fittingTeta;			//!< fitting coefficients for the thermal diffusion coefficients
		double* fittingGamma;			//!< fitting coefficients for the mass diffusion coefficients

		double* MWRatio1over4;			//!< auxiliary vector for the dynamic viscosity calculation
		double* phi_eta_sup;			//!< auxiliary vector for the dynamic viscosity calculation
		double* phi_eta_inf;			//!< auxiliary vector for the dynamic viscosity calculation

		double* sqrtEta;			//!< auxiliary vector for the dynamic viscosity calculation
		double* usqrtEta;			//!< auxiliary vector for the dynamic viscosity calculation
		Eigen::VectorXd sumK;			//!< auxiliary vector for the dynamic viscosity calculation

		double* sqrtMWRatio_inf;		//!< auxiliary vector for the dynamic viscosity calculation
		double* sqrtMWRatio_sup;		//!< auxiliary vector for the dynamic viscosity calculation

		double* sqrtMW;				//!< auxiliary vector for the dynamic viscosity calculation

		double* sum_diffusion_coefficients;	//!< auxiliary vector used for the calculation of the mass diffusion coefficients
		double* x_corrected;			//!< auxiliary vector used for the calculation of the mass diffusion coefficients

		unsigned int count_species_thermal_diffusion_ratios_;   //!< number of species for which the thermal diffusion coefficients are evaluated 
		std::vector<unsigned int> iThermalDiffusionRatios_;	//!< indices of species tracked for the thermal diffusion coefficients

		static const double threshold_;
		double sum_threshold_;

		bool temperature_lambda_must_be_recalculated_;
		bool temperature_eta_must_be_recalculated_;
		bool temperature_gamma_must_be_recalculated_;
		bool temperature_teta_must_be_recalculated_;
		bool pressure_gamma_must_be_recalculated_;

		// Indices of relevant species (1-based)
		unsigned int index_H2O_;
		unsigned int index_CO_;
		unsigned int index_CH4_;
		unsigned int index_CO2_;

		// Bundling
		unsigned int bundling_number_groups_;
		std::vector<unsigned int> bundling_reference_species_;
		std::vector< std::vector<unsigned int> > bundling_groups_;
		std::vector<unsigned int> bundling_species_group_;
		std::vector<double> bundling_sum_diffusion_coefficients_;
		std::vector<double> bundling_sum_x_groups_;

		double* bundling_fittingGammaSelfDiffusion_;
		double* bundling_fittingGamma_;
		double* bundling_gammaSpecies_;
		double* bundling_gammaSpeciesSelfDiffusion_;

		// Lennard-Jones parameters
		bool is_lennard_jones_available_ = false;	//!< true if the Lennard-Jones coefficients are available
		std::vector<double> mu_;					//!< species masses (in kg)
		std::vector<double> sigma_;					//!< species collision diameters (in m)
		std::vector<double> epsilon_over_kb_;		//!< species scaled well depths (in K)

	};
}

#include "TransportPropertiesMap_CHEMKIN.hpp"

#endif // OpenSMOKE_TransportPropertiesMap_CHEMKIN_H