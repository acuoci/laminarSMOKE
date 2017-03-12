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

#ifndef OpenSMOKE_AtomicElementMap_H
#define OpenSMOKE_AtomicElementMap_H

#include "AtomicElement.h"

namespace OpenSMOKE
{
	//!  A table of atomic weights for all the elements 
	/*!
		 This class provides the atomic weights for (almost) all the elements in the 
		 periodic table. In the current version 103 elements are available.
	*/

	class AtomicElementMap : public std::map<std::string, double>
	{
	
	public:

		/**
		* Total number of hard-coded atomic elements
		*/
		static const int nElements = 103;

		/**
		* Static table of the weights of all known elements
		*/
		static const AtomicElement AtomicWeights[nElements];

		/**
		* Constructor: default
		*/
		AtomicElementMap();

	};

	/**
	* Global data: Atomic weights table for every element in the periodic table
	*/
	extern AtomicElementMap AtomicWeights;

}

#include "AtomicElementMap.hpp"

#endif // OpenSMOKE_AtomicElementMap_H