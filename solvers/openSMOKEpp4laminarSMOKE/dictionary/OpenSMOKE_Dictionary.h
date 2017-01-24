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

#ifndef OpenSMOKE_Dictionary_H
#define	OpenSMOKE_Dictionary_H

#include <string.h>
#include <vector>
#include <iostream>
#include "OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{	
	//!  A class containing the dictionary read from a file
	/*!
			This class mainly provides the tools to extract data from a dictionary. The dictionary requires a grammar
			(which can be hard-coded or defined in ascii files), which defines a set of rules which have to be respected
			by the dictionary. The dictionary is a collection of keywords and options. Several functions automatically check
			the internal consistency of the dictionary and the agreement with the grammar rules.
	*/

	class OpenSMOKE_Dictionary
	{
	public:

		/**
		* Writes a summary of the dictionary, i.e. the list of all the keywords
		*/
		void Summary(std::ostream& fout) const;

		/**
		* Sets the name of the dictionary
		*/
		void SetName(const std::string name) { name_ = name; } 
		
		/**
		* Sets the name of the file from which the dictionary was imported (useful for debugging purposes)
		*/	
		void SetFileName(const std::string file_name) { file_name_ = file_name; } 

		/**
		* Sets the dictionary, i.e. the list of keyword, together with the options. Moreover, in order to
		  better find possible errors by the user, the class tracks also the numbers of the lines from which
		  the dictionary is imported
		*/	
		void SetDictionary(const std::vector<std::string>& keywords, const std::vector<std::string>& options, 
			               const std::vector<unsigned int>& starting_lines, const std::vector<unsigned int>& ending_lines);

		/**
		* Sets the grammar by importing it from a file (very unusual)
		*/
		void SetGrammar(const std::string file_name);

		/**
		* Sets the grammar (which was hard-coded)
		*/
		void SetGrammar(OpenSMOKE_DictionaryGrammar& grammar);

		/**
		* Checks if the keyword exists
		*/
		bool CheckOption(const std::string name_of_keyword);

		/**
		* Returns the option value as int 
		*/
		void ReadInt(const std::string option, int& value);

		/**
		* Returns the option value as double 
		*/
		void ReadDouble(const std::string option, double& value);

		/**
		* Returns the option value as char
		*/
		void ReadChar(const std::string option, char& value);

		/**
		* Returns the option value as string 
		*/
		void ReadString(const std::string option, std::string& value);

		/**
		* Returns the option value as bool 
		*/
		void ReadBool(const std::string option, bool& value);

		/**
		* Returns the option value as boost::filesystem::path
		*/
		void ReadPath(const std::string name_of_keyword, boost::filesystem::path& value);

		/**
		* Returns the option value as string (actually the returned string contains a list of substrings)
		*/
		void ReadSequence(const std::string name_of_keyword, std::string& value);

		/**
		* Returns the option value as double 
		*/
		void ReadMeasure(const std::string option, double& value, std::string& units);

		/**
		* Returns the option value as string (actually the returned string contains the name of a dictionary)
		*/
		void ReadDictionary(const std::string name_of_keyword, std::string& value);
		
		/**
		* Returns the option values as a vector of int
		*/
		void ReadOption(const std::string option, std::vector<int>& value);

		/**
		* Returns the option values as a vector of double
		*/
		void ReadOption(const std::string option, std::vector<double>& value);

		/**
		* Returns the option values as a vector of char
		*/
		void ReadOption(const std::string option, std::vector<char>& value);

		/**
		* Returns the option values as a vector of string
		*/
		void ReadOption(const std::string option, std::vector<std::string>& value);

		/**
		* Returns the option values as a vector of bool
		*/
		void ReadOption(const bool option, std::vector<bool>& value);

		/**
		* Returns the option values as a vector of strings and a vector of double
		*/	
		void ReadOption(const std::string name_of_keyword, std::vector<std::string>& strings, std::vector<double>& values);

                
       		/**
		* Returns the option values as a vector of double and a vector of strings
		*/	
		void ReadOption(const std::string name_of_keyword, std::vector<double>& values, std::vector<std::string>& strings);


		/**
		* Returns the name of the dictionary
		*/
		std::string name() const { return name_; }

	private:
			
		std::string name_;							//!< name of the dictionary
		std::string file_name_;						//!< name of the file from which the dictionary is imported
		std::vector<std::string> keywords_;			//!< vector containing all the keywords
		std::vector<std::string> options_;			//!< vector containing all the options associated to the kywords
		std::vector<unsigned int> starting_lines_;	//!< vector containing the values of the starting lines for each keyword
		std::vector<unsigned int> ending_lines_;	//!< vector containing the values of the ending lines for each keyword

		OpenSMOKE_DictionaryGrammar grammar_;		//!< grammar defining the rules on which the dictionary is based

	private:

		/**
		* Checks the internal consistency of the dictionary and the agreement with the grammar rules
		*/
		void Checks();

		/**
		* Checks the internal consistency of the required option
		*/
		int CheckOption(const std::string name_of_keyword, const OpenSMOKE_DictionaryKeyWordTypes expected_type);

		/**
		* Error message utility
		*/
		void ErrorMessage(const std::string message) const;
	};
}

#include "OpenSMOKE_Dictionary.hpp"

#endif // OpenSMOKE_Dictionary_H


