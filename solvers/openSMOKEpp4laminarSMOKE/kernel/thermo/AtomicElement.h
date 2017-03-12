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

#ifndef OpenSMOKE_AtomicElement_H
#define OpenSMOKE_AtomicElement_H

namespace OpenSMOKE
{
	//!  A class describing asingle atomic element
	/*!
		 This class provides the information for each atomic element (name,
		 molecular weight, short description).
	*/

	class AtomicElement
	{
	
	public:

		/**
		*@brief Constructor: default
		*/
		AtomicElement();

		/**
		*  Constructor: description, name and molecular weight
		*/
		AtomicElement(const std::string description, const std::string name, const double mw);

		/**
		*  Constructor: name and molecular weight (description is missing)
		*/
		AtomicElement(const std::string name, const double mw);

		/**
		*  Constructor: copy from an existing element
		*/
		AtomicElement( const AtomicElement& other );

		/**
		*  Assignment from an existing element
		*/
		AtomicElement& operator=( const AtomicElement& rhs );
	
		/**
		*  Returns the molecular weight of the atomic element
		*/
		double mw()	const				{return mw_;}

		/**
		*  Returns the name of the atomic element
		*/
		std::string name()	const 		{return name_;}

		/**
		*  Returns a brief description of the atomic element
		*/
		std::string description() const	{return description_;}

	private:
		
		double mw_;						//!< molecular weight
		std::string name_;				//!< name
		std::string description_;		//!< description

	};

}

#include "AtomicElement.hpp"

#endif // OpenSMOKE_AtomicElement_H