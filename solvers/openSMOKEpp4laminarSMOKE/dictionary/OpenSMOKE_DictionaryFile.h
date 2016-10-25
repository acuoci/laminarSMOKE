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


#ifndef OpenSMOKE_DictionaryFile_H
#define	OpenSMOKE_DictionaryFile_H

#include <string.h>
#include <vector>

#include "InputFileDictionary.h"
#include "OpenSMOKE_Dictionary.h"

namespace OpenSMOKE
{	
	//!  A class which manages the dictionary written on a file
	/*!
			This class mainly provides the tools to extract data from the ascii file containing the 
			dictionary definition.
	*/

	class OpenSMOKE_DictionaryFile
	{
		public:

			/**
			* Sets the name of the file from which the dictionary is imported
			*/
			void SetFileName(const std::string file_name) { file_name_ = file_name; }

			/**
			* Sets the index corresponding to the first line of the dictionary with respect to the 
			  ascii file from which the dictionary is imported
			*/
			void SetFirstLine(const unsigned int index_first_line) { index_first_line_ = index_first_line; }

			/**
			* Adds a line to the dictionary
			*/
			void AddLine(const std::string line) { lines_.push_back(line); }

			/**
			* Analyzes all the lines in order to check if the dictionary is correctly written (i.e. if
			  syntax errors are present)
			*/
			void Analyze();

			/**
			* All the relevant information is moved to the OpenSMOKE_Dictionary object
			*/
			void Transfer(OpenSMOKE_Dictionary& dictionary);

			/**
			* Writes the dictionary (useful for debugging purposes)
			*/
			void Write(std::ostream& fout) const;
		
		private:

			std::string name_;									//!< name of the dictionary
			std::string file_name_;								//!< name of the file from which the dictionary is imported
			unsigned int index_first_line_;						//!< index of the first line
			std::vector<std::string> lines_;					//!< vector containing all the lines
			std::vector<std::string> good_lines_;				//!< vector containing all and only the useful lines
			std::vector<std::string> extended_lines_;			//!< TODO
			std::vector<unsigned int> indices_of_good_lines_;	//!< indices of all and onlythe useful lines
			std::vector<unsigned int> i_start_;					//!< local/global indeces correspondence

			/**
			* Error message utility
			*/
			void ErrorMessage(const std::string message) const;
	};
}

#include "OpenSMOKE_DictionaryFile.hpp"

#endif // OpenSMOKE_DictionaryFile_H


