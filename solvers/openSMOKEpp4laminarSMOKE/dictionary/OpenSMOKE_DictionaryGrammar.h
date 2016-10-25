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

#ifndef OpenSMOKE_DictionaryGrammar_H
#define	OpenSMOKE_DictionaryGrammar_H

#include <boost/filesystem.hpp>
#include "OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	//!  A class which defines the grammar rules to be used by the OpenSMOKE dictionary
	/*!
			This class defines the grammar rules which have to be respected by the OpenSMOKE dictionary.
	*/

	class OpenSMOKE_DictionaryGrammar {
		
		friend class OpenSMOKE_Dictionary;

	public:

		/**
		* Initializes the grammar from a file (quite unusual)
		*/
		void ReadFromFile(const std::string file_name);

		/**
		* Writes the grammar in a readable format
		*/
		void ShortSummary(std::ostream& fout) const;

		/**
		* Default destructor
		*/
		virtual ~OpenSMOKE_DictionaryGrammar() {};
    
	private:
    
		std::string name_;								//!< name of grammar
		boost::filesystem::path* file_name_;			//!< name of the file containing the grammar (if any)

	protected: 

		/**
		* Initializes the grammar from a file (quite unusual)
		*/
		void UserDefined();

		/**
		* Checks the keywords in a dictionary to see if some of them are undefined in the current grammar
		*/
		bool CheckUndefinedKeyWords(std::vector<std::string>& list_of_keywords);

		/**
		* Checks the keywords in a dictionary to see if all the compulsory keywords are present
		*/
		bool CheckCompulsoryKeyWords(std::vector<std::string>& list_of_keywords);

		/**
		* Checks the keywords in a dictionary to see if all the required keywords are present
		*/
		bool CheckRequiredKeyWords(std::vector<std::string>& list_of_keywords);

		/**
		* Checks the keywords in a dictionary to see if there are conflicts between them
		*/
		bool CheckConflictingKeyWords(std::vector<std::string>& list_of_keywords);

		/**
		* Checks the keywords type
		*/
		bool CheckType(const std::string keyword_name, const OpenSMOKE_DictionaryKeyWordTypes expected_type);

		unsigned int number_of_keywords_;						//!< number of keywords defined in the grammar
		std::vector<OpenSMOKE_DictionaryKeyWord> keywords_;		//!< list of keywords defined in the grammar
		std::vector<std::string> list_of_keywords_names_;		//!< list of keywords names defined in the grammar

		/**
		* Checks the internal consistency of the grammar
		*/
		void Check();
		
		/**
		* Additional checks
		*/
		void CheckExistence();

		/**
		* Adds a new keyword to the grammar
		*/	
		void AddKeyWord(const OpenSMOKE_DictionaryKeyWord& keyword);

		/**
		* This function must be over-written by the derived Grammar which are hard-coded by the user. 
		  If the grammar is imported from a file this function is never called.
		*/
		virtual void DefineRules();

		/**
		* Error message utility
		*/
		void ErrorMessage(const std::string message) const;
		
	};
}

#include "OpenSMOKE_DictionaryGrammar.hpp"

#endif	/* OpenSMOKE_DictionaryGrammar_H */

