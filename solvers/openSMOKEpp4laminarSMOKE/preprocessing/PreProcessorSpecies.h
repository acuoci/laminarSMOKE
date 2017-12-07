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

#ifndef OpenSMOKE_PreProcessorSpecies_H
#define	OpenSMOKE_PreProcessorSpecies_H

#include "kernel/thermo/AtomicElementMap.h"
#include "kernel/thermo/AtomicCompositionTable.h"

namespace OpenSMOKE
{
	//!  Abstract class to preprocess thermodynamic and transport properties
	/*!
		 This class provides a common interface to different user-defined
		 pre-processors for thermodynamic and transport properties.
	*/

	template<typename Species>
	class PreProcessorSpecies : public Species
	{
	public:

		/**
		* Default constructor
		*/
		PreProcessorSpecies();

		/**
		* Copy constructor
		*/
		PreProcessorSpecies(const PreProcessorSpecies& orig);
    
		/**
		* Default destructor
		*/
		virtual ~PreProcessorSpecies();
    
		/**
		* Constructor
		  The constructor needs as input variables the readers for thermodynamic and transport 
		  properties and the kinetic preprocessor. The constructor only checks if all the species
		  contained in the kinetic mechanism are also available in the thermodynamic and transport
		  readers. If not, a fatal error message will be shown, together with a list of species
		  for which the properties are missing.
		*/
		template<typename Thermo, typename Transport, typename Kinetics>
		PreProcessorSpecies(const Thermo&  thermo, const Transport&  transport, const Kinetics& kinetics, std::ostream& flog);

		/**
		* Constructor
		The constructor needs as input variables the readers for thermodynamic and transport
		properties and the kinetic preprocessor. The constructor only checks if all the species
		contained in the kinetic mechanism are also available in the thermodynamic and transport
		readers. If not, a fatal error message will be shown, together with a list of species
		for which the properties are missing.
		*/
		template<typename Thermo, typename Transport >
		PreProcessorSpecies(const Thermo&  thermo, const Transport& transport, std::vector<std::string> kinetics, std::ostream& flog);

		/**
		* Constructor
		  The constructor needs as input variables the reader for thermodynamic properties only.
		  In this way, the preprocessor can be used only for thermodynamic analyses..
		*/
		template<typename Thermo>
		PreProcessorSpecies(const Thermo&  thermo, std::ostream& flog);

		/**
		* Constructor
		  The constructor needs as input variables the readers for thermodynamic 
		  properties and the kinetic preprocessor. The constructor only checks if all the species
		  contained in the kinetic mechanism are also available in the thermodynamic reader.
		  If not, a fatal error message will be shown, together with a list of species
		  for which the properties are missing.
		*/
		template<typename Thermo, typename Kinetics>
		PreProcessorSpecies(const Thermo&  thermo, const Kinetics& kinetics, std::ostream& flog);

		/**
		* Constructor
		  The constructor needs as input variables the readers for thermodynamic 
		  properties and the kinetic preprocessor. The constructor only checks if all the species
		  contained in the kinetic mechanism are also available in the thermodynamic reader.
		  If not, a fatal error message will be shown, together with a list of species
		  for which the properties are missing.
		*/
		template<typename Thermo, typename Kinetics >
		PreProcessorSpecies(const Thermo&  thermo, const std::vector<std::string> additional_species, const Kinetics& kinetics, std::ostream& flog);


		/**
		* Constructor
		  The constructor needs as input variables the readers for thermodynamic 
		  properties and the kinetic preprocessor. The constructor only checks if all the species
		  contained in the thermodynamic reader are also present in the transport reader.
		  If not, a fatal error message will be shown, together with a list of species
		  for which the properties are missing.
		*/
		template<typename Thermo, typename Transport>
		PreProcessorSpecies(std::ostream& flog, const Thermo&  thermo, const Transport&  transport);

		/**
		* Writes the preprocessed thermodynamic and transport properties in a old-style
		  format which can be read by a previous version of OpenSMOKE. This function has 
		  been kept only for compatibility reasons with previous versions.
		*/
		void WriteASCIIFileOldStyle(const std::string file_name) const;

		/**
		* Writes the preprocessed thermodynamic and transport properties in a 
		  format which can be read by the OpenSMOKE Maps objects (obsolete, TOREMOVE).
		*/
		void WriteASCIIFile(const std::string file_name) const;

		/**
		* Writes the preprocessed thermodynamic and transport properties in a 
		  XML format which can be read by the OpenSMOKE Maps objects
		*/
		void WriteXMLFile(std::stringstream& xml_string) const;

		/**
		*@brief This function returns the matrix of atomic composition
		*@return the matrix of atomic composition
		*/
		const AtomicCompositionTable& AtomicTable() const { return atomicTable_; }

		/**
		*@brief This function returns the matrix of atomic composition
		*@return the matrix of atomic composition
		*/
		AtomicCompositionTable& AtomicTable() { return atomicTable_; }

		/**
		*@brief Check the thermodynamic consistency of the properties contained in the thermodynamic database
		*/
		bool ThermodynamicConsistency();

		/**
		*@brief Writes the atomic composition of each species in a readable format (tabulated)
		*/
		bool WriteElementTableOnASCIIFile(std::ostream& fOutput) const;

		/**
		*@brief Writes the atomic composition of each species in a readable format (tabulated) 
		        and ordered according to the elemental composition
		*/
		bool WriteReorderedElementTableOnASCIIFile(std::ostream& fOutput) const;

		/**
		*@brief Species bundling accordingly to the specified maximum error
		*/
		bool SpeciesBundling(std::stringstream& xml_string, const double epsilon);
		

	protected:

		AtomicCompositionTable atomicTable_;	//!< matrix of atomic composition 

		std::ostream& flog_;	//!< log file where useful information is written during the pre-processing of a kinetic mechanism
	};
}

#include "PreProcessorSpecies.hpp"

#endif	/* OpenSMOKE_PreProcessorSpecies_H */

