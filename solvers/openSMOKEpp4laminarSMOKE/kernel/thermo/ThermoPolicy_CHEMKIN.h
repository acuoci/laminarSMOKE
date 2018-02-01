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

#ifndef OpenSMOKE_ThermoPolicy_CHEMKIN_H
#define	OpenSMOKE_ThermoPolicy_CHEMKIN_H

#include "math/OpenSMOKEVector.h"
#include "AtomicComposition.h"
#include "math/PhysicalConstants.h"
#include "math/Conversions.h"

namespace OpenSMOKE
{
	//!  A class to manage the thermodynamic properties of a species in CHEMKIN format
	/*!
		 This class reads and analyzes the thermodynamic properties in CHEMKIN format
		 (atomic composition, phase, temperature limits for thermodynamic properties and
		 thermodynamic coefficients)
	*/

	class ThermoPolicy_CHEMKIN 
	{
	
	public:

		/**
		*@brief Default constructor
		*/
		ThermoPolicy_CHEMKIN();

		/**
		*@brief Default copy constructor
		*/
		ThermoPolicy_CHEMKIN(const ThermoPolicy_CHEMKIN& orig);
		
		/**
		*@brief Default destructor
		*/		
		virtual ~ThermoPolicy_CHEMKIN();
    
		/**
		*@brief Name of the species as reported in the CHEMKIN thermodynamic file
		*/
		std::string name_thermo() const { return name_thermo_; }
		
		/**
		*@brief Phase of the species
		*/
		char phase() const { return phase_; }

		/**
		*@brief Returns the specific heat at constant pressure of the species [J/kg/K]
		*/
		double cp(const double T) const;

		/**
		*@brief Returns the specific heat at constant volume of the species [J/kg/K]
		*/
		double cv(const double T) const;

		/**
		*@brief Returns the mass enthalpy the species [J/kg]
		*/
		double enthalpy(const double T) const;

		/**
		*@brief Returns the mass entropy of the species [J/kg/K]
		*/
		double entropy(const double T) const;

		/**
		*@brief Returns the Gibbs energy of the species [J/kg]
		*/
		double gibbs_energy(const double T) const;
    
		/**
		*@brief Reads the thermodynamic data from several lines, usually extracted from the CHEMKIN thermodynamic file
		*/
		bool ReadThermodynamicsFromStrings(const std::vector<std::string> &line, const OpenSMOKE::OpenSMOKEVectorDouble& temperature_limits_default_);

		/**
		*@brief Reports some useful information
		*/
		void ThermodynamicsStatus(std::ostream& fOut, const std::vector<double>& T) const;

		/**
		*@brief Reports some useful information
		*/
		void ThermodynamicsStatus(std::ostream& fOut) const;
	
		/**
		*@brief Writes the thermodynamic data on a binary file
		*/
		bool WriteOnBinaryFile(std::ofstream& fOut);

		/**
		*@brief Reads the thermodynamic data from a binary file
		*/
		bool ReadFromBinaryFile(std::ifstream& fIn);

		/**
		*@brief Writes the thermodynamic data on a file (ASCII)
		*/
		bool WriteOnASCIIFile(std::ostream& fOutput) const;

		/**
		*@brief Writes the thermodynamic data on a file (ASCII) in a old-style
		*/
		bool WriteOnASCIIFileOldStyle(std::ofstream& fOutput) const;

		/**
		*@brief Returns the molecular weight
		*/
		double MolecularWeight() const { return atomic_composition_.mw(); }

		/**
		*@brief Returns the low (first) limit temperature (look at the CHEMKIN manual)
		*/
		double Tlow() const { return thermodynamics_parameters_[1]; }

		/**
		*@brief Returns the medium (second) limit temperature (look at the CHEMKIN manual)
		*/
		double Tmedium() const { return thermodynamics_parameters_[2]; }

		/**
		*@brief Returns the high (third) limit temperature (look at the CHEMKIN manual)
		*/
		double Thigh() const { return thermodynamics_parameters_[3]; }

		/**
		*@brief Copies the thermodynamic data to the rhs species
		*/
		void CopyThermodynamics(ThermoPolicy_CHEMKIN& rhs) const;

		/**
		*@brief Copies the atomic composition to the atomic object
		*/
		void AtomicComposition(OpenSMOKE::AtomicComposition *atomic) const;

		/**
		*@brief Returns the thermodynamic properties
		*/
		const OpenSMOKEVectorDouble& thermodynamics_parameters() const { return thermodynamics_parameters_;}

		/**
		*@brief Returns the ratio between the R constnat and the molecular weight [J/kg/K]
		*/
		double RGAS_over_MW() const { return PhysicalConstants::R_J_mol/(MolecularWeight()*1.e-3); }

		/**
		*@brief Checks the consistency of thermodynamic properties and corrects the thermodynamic coefficients
		        if small errors are found
		*/		
		int CheckThermodynamicConsistency(std::ostream& fout);

		/**
		*@brief Reformultes the thermodynamic properties in order to ensure their consistency at the
		        intermediate temperature
		*/	
		void ReformulationOfThermodynamics(std::ostream& fout, const unsigned int policy_intermediate_temperature, const double intermediate_temperature, const double max_temperature);

		/**
		*@brief Reformultes the thermodynamic properties in order to ensure their consistency at the
		intermediate temperature
		*/
		void ReformulationOfThermodynamicsFixedIntermediateTemperature(std::ostream& fout, const unsigned int policy_intermediate_temperature, const double intermediate_temperature, const double max_temperature);

		/**
		*@brief Reports information about the consistency of the specific heat at the intermediate temperature
		*/
		int ReportStatusSpecificHeatCoefficient(std::ostream& fout);

		/**
		*@brief Reports information about the consistency of the specific enthalpy at the intermediate temperature
		*/
		int ReportStatusEnthalpy(std::ostream& fout);

		/**
		*@brief Reports information about the consistency of the specific entropy at the intermediate temperature
		*/
		int ReportStatusEntropy(std::ostream& fout);

		/**
		*@brief Checks if anomalous behaviours of the thermodynamic properties exist (basically local maxima)
		*/
		int CheckForAnomalies(std::ostream& fout);

	protected:
    
		char phase_;										//!< phase of the species
		std::string name_thermo_;							//!< name of the species
		OpenSMOKE::AtomicComposition atomic_composition_;	//!< atomic composition of the species
		OpenSMOKEVectorDouble thermodynamics_parameters_;	//!< thermodynamic parameters
	};
}

#include "ThermoPolicy_CHEMKIN.hpp"

#endif	/* OpenSMOKE_ThermoPolicy_CHEMKIN_H */

