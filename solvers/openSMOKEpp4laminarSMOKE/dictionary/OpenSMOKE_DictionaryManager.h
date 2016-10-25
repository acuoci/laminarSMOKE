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

#ifndef OpenSMOKE_DictionaryManager_H
#define	OpenSMOKE_DictionaryManager_H

#include <string.h>
#include <vector>
#include <iostream>
#include "InputFileDictionary.h"
#include "OpenSMOKE_Dictionary.h"

namespace OpenSMOKE
{	
	//!  A class which manages a list of different dictionaries
	/*!
			This class manages a set of different dictionaries. The user can freely chose the dictionary contained
			in this list (which is stored as a map).
	*/

	class OpenSMOKE_DictionaryManager
	{
	public:

		/**
		* Reads all the dictionaries defined in the requested file
		*/
		void ReadDictionariesFromFile(const std::string file_name);

		/**
		* Return the requested dictionary
		* @param name name of the requested dictionary
		*/
		OpenSMOKE_Dictionary& operator() (const std::string& name);

	private:

		std::map<std::string, OpenSMOKE_Dictionary> map_of_dictionaries_;	//!< map containing the dictionaries
		
		/**
		* Error message utility
		*/		
		void ErrorMessage(const std::string message) const;
	};
}

#include "OpenSMOKE_DictionaryManager.hpp"

#endif // OpenSMOKE_DictionaryManager_H


