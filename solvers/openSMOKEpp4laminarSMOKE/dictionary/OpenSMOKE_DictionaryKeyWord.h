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

#ifndef OpenSMOKE_DictionaryKeyWord_H
#define	OpenSMOKE_DictionaryKeyWord_H

#include <boost/filesystem.hpp>

namespace OpenSMOKE
{
	//!  A class which defines the keywords to be used in the OpenSMOKE dictionaries
	/*!
			This class defines the keywords to be used in the OpenSMOKE dictionaries
	*/

	enum OpenSMOKE_DictionaryKeyWordTypes { NONE, SINGLE_INT, SINGLE_DOUBLE, SINGLE_STRING, SINGLE_CHAR, SINGLE_BOOL, SINGLE_MEASURE, SINGLE_PATH, SINGLE_DICTIONARY,
											VECTOR_INT, VECTOR_DOUBLE, VECTOR_STRING, VECTOR_CHAR, VECTOR_BOOL, VECTOR_MEASURE,
											VECTOR_STRING_DOUBLE, SEQUENCE_STRING };

	class OpenSMOKE_DictionaryKeyWord {
	public:

		/**
		* Default constructor
		*/
		OpenSMOKE_DictionaryKeyWord() {};

		/**
		* Constructor
		* @param name name of the keyword (this name will be directly used in the ascii file)
		* @param type type of the keyword
		* @param short_comment short comment which will be displayed in the summary
		* @param is_compulsory if the keyword is compulsory
		*/
		OpenSMOKE_DictionaryKeyWord(const std::string& name, const OpenSMOKE_DictionaryKeyWordTypes type, 
									const std::string& short_comment, const bool is_compulsory);

		/**
		* Constructor
		* @param name name of the keyword (this name will be directly used in the ascii file)
		* @param type type of the keyword
		* @param short_comment short comment which will be displayed in the summary
		* @param is_compulsory if the keyword is compulsory
		* @param compulsory_alternatives list of alternative keywords to the current compulsory keyword
		* @param required_keywords list of keyword which are needed by the current keyword
		* @param conflicting_keywords list of conflicting keywords
		*/
		OpenSMOKE_DictionaryKeyWord(const std::string& name, const OpenSMOKE_DictionaryKeyWordTypes type, 
									const std::string& short_comment, const bool is_compulsory, 
									const std::string& compulsory_alternatives, const std::string& required_keywords, 
									const std::string& conflicting_keywords);

		/**
		* Constructor
		* Initializes the keyword from a block of lines extracted from a file
		* @param lines
		*/
		OpenSMOKE_DictionaryKeyWord(std::vector<std::string>& lines);

		/**
		*@brief Default destructor
		*/
		virtual ~OpenSMOKE_DictionaryKeyWord() {};

		/**
		* Shows a short summary of the keyword
		*/
		void ShortSummary(std::ostream& fout) const;

		/**
		* Returns the name of the keyword
		*/
		const std::string& name() const { return name_; }

		/**
		* Returns the type of the keyword
		*/
		OpenSMOKE_DictionaryKeyWordTypes type() const { return type_; }

		/**
		* Returns the type of the keyword in ASCII format
		*/
		const std::string type_ascii() const { return type_ascii_; }

		/**
		* Returns true is the keyword is compulsory
		*/
		bool is_compulsory() const { return is_compulsory_; }

		/**
		* Returns the list of alternatives to the current compulsory keyword
		*/
		const std::vector<std::string>& compulsory_alternatives() const { return compulsory_alternatives_; }

		/**
		* Returns the the list of keywords needed by the current keyword
		*/
		const std::vector<std::string>& required_keywords() const { return required_keywords_; }

		/**
		* Returns the the list of keywords conflicting with the current keyword
		*/
		const std::vector<std::string>& conflicting_keywords() const { return conflicting_keywords_; }

	private:
    
		std::string name_;									//!< name of the keyword
		std::string type_ascii_;							//!< type of keyword
		OpenSMOKE_DictionaryKeyWordTypes type_;				//!< type of keyword
		std::string comment_short_;							//!< short comment
		bool is_compulsory_;								//!< is the key word compulsory?
		std::vector<std::string> compulsory_alternatives_;	//!< list of compulsory alternatives (if any)
		std::vector<std::string> required_keywords_;		//!< list of required keywords (if any)
		std::vector<std::string> conflicting_keywords_;		//!< list of excluded keywords (if any)
		
	private:

		/**
		* Set the keyword
		*/
		void SetKeyword(const std::string& line);
		
		/**
		* Set the keyword type
		*/
		void SetType(const std::string& line);

		/**
		* Set the keyword type (ASCII format)
		*/
		void SetTypeASCII();

		/**
		* Set the comment
		*/
		void SetComment(const std::string& line);

		/**
		* Set if the keyword is compulsory
		*/
		void SetIsCompulsory(const std::string& line);

		/**
		* Set the list of alternative keywords
		*/
		void SetCompulsoryAlternatives(const std::string& line);

		/**
		* Set the list of keywords which are needed by the current keyword
		*/
		void SetRequiredKeyWords(const std::string& line);

		/**
		* Set the list of keywords which conflicts whit the current keyword
		*/
		void SetConflictingKeywords(const std::string& line);

		/**
		* Error message utility
		*/
		void ErrorMessage(const std::string message) const;
	};
}

#include "OpenSMOKE_DictionaryKeyWord.hpp"

#endif	/* OpenSMOKE_DictionaryKeyWord_H */

