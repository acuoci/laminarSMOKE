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

#ifndef OpenSMOKE_Species_H
#define	OpenSMOKE_Species_H

#include <string>
#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	//!  A class providing a common interface for storing and manage information about a single species
	/*!
		 This class basically provides a common interface for storing and manage basic information
		 about a single chemical species.
	*/

	template<typename ThermoPolicy, typename TransportPolicy>
	class Species : public ThermoPolicy, public TransportPolicy
	{
	public:

		/**
		*@brief Default constructor
		*/
		Species();

		/**
		*@brief Default copy constructor
		*/
		Species(const Species& orig);

		/**
		*@brief Default destructor
		*/
		virtual ~Species();
	
		/**
		*@brief Constructor from thermodynamic and transport data
		*/    
		Species(const ThermoPolicy& thermo, const TransportPolicy& transport);

		/**
		*@brief Constructor from thermodynamic (TODO)
		*/    
		Species(const ThermoPolicy& thermo);
	
		/**
		*@brief Assignment
		*/	
		Species& operator=(const Species& rhs);
	
		/**
		*@brief Overloading () operator
		*/	
		void operator() (const ThermoPolicy& thermo, const TransportPolicy& transport);

		/**
		*@brief Overloading () operator
		*/	
		void operator() (const ThermoPolicy& thermo);

		/**
		*@brief Returns the name of the species
		*/	
		std::string name() const { return name_; };

	private:

		std::string name_;				//!< name of the species
		std::string description_;		//!< short description (if any)
		double mw_;						//!< molecular weight [kg/kmol]
	};

}

#include "Species.hpp"

#endif	/* OpenSMOKE_Species_H */

