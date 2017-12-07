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

#ifndef OpenSMOKE_PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport_H
#define	OpenSMOKE_PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport_H

#include <Eigen/Dense>
#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{

	//!  A class to preprocess thermodynamic properties provided in CHEMKIN format
	/*!
		 This class preprocess the thermodynamic properties provided
		 in CHEMKIN format.
	*/

	template<typename Species>
	class PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport 
	{
	public:

		/**
		* Default constructor
		*/
		PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport();

		/**
		* Copy constructor
		*/
		PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport(const PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport& orig);
    
		/**
		* Default destructor
		*/
		virtual ~PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport();
    
		/**
		*@brief This function returns the vector containg the list of species
		*@return the vector containing the list of species
		*/
		std::vector<Species>& species()				{ return species_; }

		/**
		*@brief This function returns the vector containg the list of species
		*@return the vector containing the list of species
		*/
		const std::vector<Species>& species() const { return species_; }

		/** 
		* This function allocates the memory and prepares the data to be used for
		  the internal preprocessing
		*/
		bool Setup();

		/**
		* This function write the thermodynamic properties in a readeable format.
		  This function can be called only after the thermodynamic properties have 
		  been preprocessed
		*/
		bool WriteThermodynamicDataOnASCIIFileOldStyle(std::ofstream &fOutput) const;

		/**
		* This function write the transport properties in a readeable format
		  This function can be called only after the transport properties have 
		  been preprocessed
		*/
		bool WriteTransportDataOnASCIIFileOldStyle(std::ofstream &fOutput) const;

		/**
		* This function writes the thermodynamic properties in a readeable format
		  to be read by the Thermodynamic Maps (obsolete, TOREMOVE)
		*/
		bool WriteThermodynamicsDataOnASCIIFile(std::ofstream &fOutput) const;

		/**
		* This function writes the thermodynamic properties in a readeable XML format
		  to be read by the Thermodynamic Maps
		*/
		bool WriteThermodynamicsDataOnXMLFile(std::stringstream& xml_string) const;

		/**
		* This function writes the list of species in a readeable XMLformat
		  to be read by the Thermodynamic Maps
		*/
		bool WriteSpeciesOnXMLFile(std::stringstream &xml_string) const;

		/**
		* This function write the thermodynamic coefficients in a file in a readable format
		*/
		bool WriteThermodynamicCoefficientsOnASCIIFile(const std::string file_name) const;

		/**
		* This function write the thermodynamic data for each molecule in a readable format
		*/
		bool WriteThermodynamicTablesOnASCIIFile(const std::string file_name) const;

		/**
		*@brief Checks the thermodynamic onsistency
		*/
		bool CheckThermodynamicConsistency(std::ostream& flog);

		/**
		* This function reformulates the thermodynamic properties in order to ensure their consistency.
		  In particular, the new formulation ensure the continuity of the first, second and third derivates
		  and chooses the intermediate temperature in order to minimize the fitting error with respect to
		  the original thermodynamic data. This function can be used only if the original thermodynamic
		  data are correct, i.e. no macroscopic error exist. This means that only small inaccuracies
		  can be accepted.
		*/
		bool ReformulationOfThermodynamics(const std::string file_name, const std::string original_file_name);

		/**
		* This function reformulates the thermodynamic properties in order to ensure their consistency.
		In particular, the new formulation ensure the continuity of the first, second and third derivates. 
		This function can be used only if the original thermodynamic data are correct, 
		i.e. no macroscopic error exist. This means that only small inaccuracies
		can be accepted.
		*/
		bool ReformulationOfThermodynamicsFixedIntermediateTemperature(const std::string file_name, const std::string original_file_name);

		/**
		* This function analyzes the thermodnamic data for all the species and checks the continuity
		  of the derivatives at the intermediate temperature. A report is written on a file in ascii format.
		*/
		bool StatusOfThermodynamics(const std::string file_name);

		/**
		* This function writes the thermodynamic properties in ascii file to be post-processed and 
		  analyzed in MATLAB.
		*/
		bool WriteThermodynamicsForMATLAB(const std::string file_name);

		/**
		* This function writes the transport properties in a readeable format
		  to be read by the TransportProperties Maps (obsolete, TOREMOVE)
		*/
		bool WriteTransportDataOnASCIIFile(std::ofstream &fOutput) const;

		/**
		* This function writes the transport properties in a readeable XML format
		  to be read by the TransportProperties Maps
		*/
		bool WriteTransportDataOnXMLFile(std::stringstream &xml_string) const;

		/**
		* This function writes the transport properties (Lennard-Jones parameters)
		  in a readable format (tabulation)
		*/
		bool WriteTransportTableOnASCIIFile(std::ostream& fOutput) const;

	private:

		/**
		* Allocation of memory
		*/
		void Allocate();

	protected:

		/**
		*@brief This function returns the vector containg the species
		*@return the vector containing the species
		*/
		std::vector<Species> species_;

		/**
		*@brief This function returns the vector containg the list of species names
		*@return the vector containing the list of species names
		*/
		std::vector<std::string> names_;

	protected:

		/** @name Evaluation of specific heats of species
		 *  This group provides the functions to calculate the specific heats of single species
		 */
		///@{
		/** calculates the specific heat at constant pressure */
		void SpeciesCp(const double T);
		/** calculates the specific heat at constant volume */
		void SpeciesCv(void);
		///@}

	protected:

		unsigned int NC;	//!< number of species

		OpenSMOKE::OpenSMOKEVectorDouble MW;			//!< molecular weight [kg/kmol]
		OpenSMOKE::OpenSMOKEVectorDouble uMW;			//!< reciprocal of molecular weight [kmol/kg]

		// Thermodynamic properties
		OpenSMOKE::OpenSMOKEVectorDouble Cp;			//!< specific heat at constant pressure
		OpenSMOKE::OpenSMOKEVectorDouble Cv;			//!< specific heat at constant volume
	};

}

#include "PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport.hpp"

#endif	/* OpenSMOKE_PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport_H */

