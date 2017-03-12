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

#ifndef OpenSMOKE_ThermoReaderPolicy_CHEMKIN_H
#define	OpenSMOKE_ThermoReaderPolicy_CHEMKIN_H

#include "AtomicComposition.h"
#include "kernel/kinetics/KineticsUtilityFunctions.h"

namespace OpenSMOKE
{
	//!  A class for reading a thermodynamic database (from file) in CHEMKIN format
	/*!
		 The purpose of this class is to read a file containing the thermodynamic database in CHEMKIN format
	*/

	template<typename Species>
	class ThermoReaderPolicy_CHEMKIN 
	{

	public:

		typedef std::map<std::string, Species, OpenSMOKE_Utilities::ciLessLibC> map_species;

		/**
		*@brief Default constructor
		*/
		ThermoReaderPolicy_CHEMKIN();
		
		/**
		*@brief Default copy constructor
		*/		
		ThermoReaderPolicy_CHEMKIN(const ThermoReaderPolicy_CHEMKIN& orig);
		
		/**
		*@brief Default destructor
		*/		
		virtual ~ThermoReaderPolicy_CHEMKIN();

		/**
		*@brief Returns the map of species included in the thermodynamic database
		*/	
		map_species& species() { return species_; }

		/**
		*@brief Returns the map of species included in the thermodynamic database
		*/	
		const map_species& species() const { return species_; }

	protected:

		map_species species_;	//!<  map of species included in the thermodynamic database

		/**
		*@brief Read the thermodynamic database in CHEMKIN format
		*/	
		bool ReadFromASCIIFile(const std::string file_name);
	};
}

#include "ThermoReaderPolicy_CHEMKIN.hpp"

#endif	/* OpenSMOKE_ThermoReaderPolicy_CHEMKIN_H */

