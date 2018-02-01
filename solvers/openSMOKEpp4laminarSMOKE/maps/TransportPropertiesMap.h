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

#ifndef OpenSMOKE_TransportPropertiesMap_H
#define OpenSMOKE_TransportPropertiesMap_H

namespace OpenSMOKE
{
	//!  A struct to store common data for transport properties maps
	/*!
	This struct provides the tools o store common data for transport properties maps.
	*/

	struct TransportPropertiesMapBaseClass
	{
	protected:

		unsigned int nspecies_;

		Eigen::VectorXd lambdaSpecies_;
		Eigen::VectorXd etaSpecies_;
		Eigen::VectorXd gammaSpecies_;
		Eigen::VectorXd tetaSpecies_;
	};


	//!  A virtual class to provide a common interface to transport properties maps
	/*!
	This virtual class provides a common interface to transport properties maps
	*/

	class TransportPropertiesMap : public TransportPropertiesMapBaseClass
	{

	public:

		/**
		* Sets the temperature (in K)
		*/
		virtual void SetTemperature(const double& T) = 0;

		/**
		* Sets the pressure (in Pa)
		*/
		virtual void SetPressure(const double& P) = 0;

		/**
		*@brief Calculates the thermal conductivity of a mixture from the mole fractions
		*/
		double ThermalConductivity(const double* moleFractions);

		/**
		*@brief Calculates the dynamic viscosity of a mixture from the mole fractions
		*/
		double DynamicViscosity(const double* moleFractions);

		/**
		*@brief Calculates the mass diffusion coefficients (mixture averaged formulation) of a mixture from the mole fractions
		*/
		void MassDiffusionCoefficients(double* gammamix, const double* moleFractions, const bool bundling = false);

		/**
		*@brief Calculates the thermal diffusion coefficients (mixture averaged formulation) of a mixture from the mole fractions
		*/
		void ThermalDiffusionRatios(double* tetamix, const double* moleFractions);

		/**
		*@brief Sets the coefficients of transport properties for each species
		*/
		virtual void SetCoefficients(const unsigned int i, const double* coefficients) = 0;

		/**
		*@brief Imports the transport coefficients for all the species from an external file
		*/
		virtual void ImportCoefficientsFromASCIIFile(std::istream& fInput) = 0;

		/**
		*@brief Returns true if the species bundling for calculation of mass diffusion coefficients is activated
		*/
		bool is_species_bundling() const { return species_bundling_; }


	protected:

		/**
		*@brief TODO
		*/
		virtual void Test(const int nLoops, const double& T, int* index) = 0;

	protected:

		/**
		*@brief Combines the species thermal conductivities to calculate the mixture thermal conductivity
		*/
		virtual double lambdaMix(const double* moleFractions) = 0;

		/**
		*@brief Combines the species dynamic viscosities to calculate the mixture dynamic viscosity
		*/
		virtual double etaMix(const double* moleFractions) = 0;

		/**
		*@brief Combines the species mass diffusion coefficients to calculate the mixture mass diffusion coefficients
		*/
		virtual void gammaMix(double* gammamix, const double* moleFractions) = 0;

		/**
		*@brief Combines the species mass diffusion coefficients to calculate the mixture mass diffusion coefficients using the bundling algorithm
		*/
		virtual void bundling_gammaMix(double* gammamix, const double* moleFractions) = 0;

		/**
		*@brief Combines the species thermal diffusion coefficients to calculate the mixture thermal diffusion coefficients
		*/
		virtual void tetaMix(double* tetamix, const double* moleFractions) = 0;

		/**
		*@brief Combines the planck mean absorption coefficients of relevant species, according to their mole fractions
		*/
		virtual double kPlanckMix(const double* moleFractions) = 0;

		/**
		*@brief Calculates the thermal conductivities for all the species
		*/
		virtual void lambda() = 0;

		/**
		*@brief Calculates the dynamic viscosities for all the species
		*/
		virtual void eta() = 0;

		/**
		*@brief Calculates the mass diffusion coefficients for all the species
		*/
		virtual void gamma() = 0;

		/**
		*@brief Calculates the mass diffusion coefficients for all the species using the bundling algorithm
		*/
		virtual void bundling_gamma() = 0;

		/**
		*@brief Calculates the thermal diffusion coefficients for all the species
		*/
		virtual void teta() = 0;

	protected:

		double T_;			//!< temperature [K]
		double P_;			//!< pressure [Pa]

		double T_old_;			//!< temperature [K] (previous value)
		double P_old_;			//!< pressure [Pa] (previous value)

		bool species_bundling_;		//!< bundling of species for calculation of mass diffusion coefficients
	};
}

#include "TransportPropertiesMap.hpp"

#endif // OpenSMOKE_TransportPropertiesMap_Hpp