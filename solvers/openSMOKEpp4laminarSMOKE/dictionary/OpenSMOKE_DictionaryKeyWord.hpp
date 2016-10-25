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

#include "math/OpenSMOKEUtilities.h"
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

namespace OpenSMOKE
{
	void OpenSMOKE_DictionaryKeyWord::ErrorMessage(const std::string message) const
	{
		std::cout << "Keyword  " << name_ << std::endl;
		std::cout << "Fatal error:     " << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}
	
	OpenSMOKE_DictionaryKeyWord::OpenSMOKE_DictionaryKeyWord(const std::string& name, const OpenSMOKE_DictionaryKeyWordTypes type, 
															 const std::string& short_comment, const bool is_compulsory) :
		name_(name), 
		type_(type),
		comment_short_(short_comment),
		is_compulsory_(is_compulsory)
	{
		SetTypeASCII();
	}
	
	OpenSMOKE_DictionaryKeyWord::OpenSMOKE_DictionaryKeyWord(const std::string& name, const OpenSMOKE_DictionaryKeyWordTypes type, 
															 const std::string& short_comment, const bool is_compulsory, 
															 const std::string& compulsory_alternatives, const std::string& required_keywords, 
															 const std::string& conflicting_keywords) :
		name_(name), 
		type_(type),
		comment_short_(short_comment),
		is_compulsory_(is_compulsory)
	{
		
		std::string line1 = "compulsory_alternatives: " + compulsory_alternatives;
		SetCompulsoryAlternatives(line1);

		std::string line2 = "required_keywords: " + required_keywords;
		SetRequiredKeyWords(line2);

		std::string line3 = "conflicting_keywords: " + conflicting_keywords;
		SetConflictingKeywords(line3);

		SetTypeASCII();
	}

	void OpenSMOKE_DictionaryKeyWord::SetTypeASCII()
	{
		     if (type_ == NONE )					type_ascii_ = "none";
		else if (type_ == SINGLE_INT )				type_ascii_ = "single_int";
		else if (type_ == SINGLE_DOUBLE )			type_ascii_ = "single_double";
		else if (type_ == SINGLE_CHAR )				type_ascii_ = "single_char";
		else if (type_ == SINGLE_STRING )			type_ascii_ = "single_string";
		else if (type_ == SINGLE_BOOL )				type_ascii_ = "single_bool";
		else if (type_ == SINGLE_MEASURE )			type_ascii_ = "single_measure";
		else if (type_ == SINGLE_PATH )				type_ascii_ = "single_path";
		else if (type_ == SINGLE_DICTIONARY )		type_ascii_ = "single_dictionary";
		else if (type_ == VECTOR_INT )				type_ascii_ = "vector_int";
		else if (type_ == VECTOR_DOUBLE )			type_ascii_ = "vector_double";
		else if (type_ == VECTOR_CHAR )				type_ascii_ = "vector_char";
		else if (type_ == VECTOR_STRING )			type_ascii_ = "vector_string";
		else if (type_ == VECTOR_BOOL )				type_ascii_ = "vector_bool";
		else if (type_ == VECTOR_MEASURE )			type_ascii_ = "vector_measure";
		else if (type_ == VECTOR_STRING_DOUBLE )	type_ascii_ = "vector_string_double";
		else if (type_ == SEQUENCE_STRING )			type_ascii_ = "sequence_string";
		else ErrorMessage("The keyword type is not specified correctly.");
	}

	void OpenSMOKE_DictionaryKeyWord::SetKeyword(const std::string& line)
	{
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");

		tokenizer_blank tokens(line, sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		tokenizer_blank::iterator tok_blank = tokens.begin();

		if (n !=2 || *tok_blank!="keyword:" )
				ErrorMessage("The keyword name is not specified correctly.");
		++tok_blank;
		name_ = *tok_blank;
	}

	void OpenSMOKE_DictionaryKeyWord::SetType(const std::string& line)
	{
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");

		tokenizer_blank tokens(line, sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		tokenizer_blank::iterator tok_blank = tokens.begin();

		if (n !=2 || *tok_blank!="type:" )
				ErrorMessage("The keyword type is not specified correctly.");
		++tok_blank;
		if (*tok_blank == "none")						type_ = NONE; 
		else if (*tok_blank == "single_int")			type_ = SINGLE_INT;
		else if (*tok_blank == "single_double")			type_ = SINGLE_DOUBLE;
		else if (*tok_blank == "single_string")			type_ = SINGLE_STRING;
		else if (*tok_blank == "single_char")			type_ = SINGLE_CHAR;
		else if (*tok_blank == "single_bool")			type_ = SINGLE_BOOL;
		else if (*tok_blank == "single_measure")		type_ = SINGLE_MEASURE;
		else if (*tok_blank == "single_path")			type_ = SINGLE_PATH;
		else if (*tok_blank == "single_dictionary")		type_ = SINGLE_DICTIONARY;
		else if (*tok_blank == "vector_int")			type_ = VECTOR_INT;
		else if (*tok_blank == "vector_double")			type_ = VECTOR_DOUBLE;
		else if (*tok_blank == "vector_string")			type_ = VECTOR_STRING;
		else if (*tok_blank == "vector_char")			type_ = VECTOR_CHAR;
		else if (*tok_blank == "vector_bool")			type_ = VECTOR_BOOL;
		else if (*tok_blank == "vector_measure")		type_ = VECTOR_MEASURE;
		else if (*tok_blank == "vector_string_double")	type_ = VECTOR_STRING_DOUBLE;
		else if (*tok_blank == "sequence_string")	    type_ = SEQUENCE_STRING;
		else ErrorMessage("The keyword type is not specified correctly.");

		type_ascii_ = *tok_blank;
	}

	void OpenSMOKE_DictionaryKeyWord::SetComment(const std::string& line)
	{
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");

		tokenizer_blank tokens(line, sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		tokenizer_blank::iterator tok_blank = tokens.begin();

		if (*tok_blank!="short_comment:" )
				ErrorMessage("The keyword short comment is not specified correctly.");
		for (;;)
		{
			++tok_blank;
			if ( tok_blank == tokens.end() )	break;
			comment_short_ += *tok_blank + " ";
		}
	}

	void OpenSMOKE_DictionaryKeyWord::SetIsCompulsory(const std::string& line)
	{
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");

		tokenizer_blank tokens(line, sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		tokenizer_blank::iterator tok_blank = tokens.begin();

		if (n !=2 || *tok_blank!="is_compulsory:" )
				ErrorMessage("The keyword is_compulsory is not specified correctly.");
		++tok_blank;
		if (*tok_blank == "true")	is_compulsory_ = true;
		else if (*tok_blank == "false")	is_compulsory_ = false;
		else ErrorMessage("The keyword is_compulsory is not specified correctly.");
	}

	void OpenSMOKE_DictionaryKeyWord::SetCompulsoryAlternatives(const std::string& line)
	{
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");

		tokenizer_blank tokens(line, sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		tokenizer_blank::iterator tok_blank = tokens.begin();

		if (*tok_blank!="compulsory_alternatives:" )
				ErrorMessage("The keyword compulsory_alternatives are not specified correctly.");
			
		for (;;)
		{
			++tok_blank;
			if ( tok_blank == tokens.end() )	break;
			compulsory_alternatives_.push_back(*tok_blank);
		}

		if (is_compulsory_ == false)
		{
			if (compulsory_alternatives_.size() != 0)
				if (compulsory_alternatives_[0] != "none")
					ErrorMessage("The keyword compulsory_alternatives are not specified correctly.");
		}

		if (compulsory_alternatives_.size() == 1)
		{
			if (compulsory_alternatives_[0] == "none")
				compulsory_alternatives_.resize(0);
		}
	}

	void OpenSMOKE_DictionaryKeyWord::SetRequiredKeyWords(const std::string& line)
	{
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");

		tokenizer_blank tokens(line, sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		tokenizer_blank::iterator tok_blank = tokens.begin();

		if (*tok_blank!="required_keywords:" )
				ErrorMessage("The keyword required_keywords are not specified correctly.");
			
		for (;;)
		{
			++tok_blank;
			if ( tok_blank == tokens.end() )	break;
			required_keywords_.push_back(*tok_blank);
		}

		if (required_keywords_.size() != 0)
		{
			if (required_keywords_[0] == "none")
				required_keywords_.resize(0);
		}
	}

	void OpenSMOKE_DictionaryKeyWord::SetConflictingKeywords(const std::string& line)
	{
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");

		tokenizer_blank tokens(line, sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		tokenizer_blank::iterator tok_blank = tokens.begin();

		if (*tok_blank!="conflicting_keywords:" )
				ErrorMessage("The keyword conflicting_keywords are not specified correctly.");
			
		for (;;)
		{
			++tok_blank;
			if ( tok_blank == tokens.end() )	break;
			conflicting_keywords_.push_back(*tok_blank);
		}

		if (conflicting_keywords_.size() != 0)
		{
			if (conflicting_keywords_[0] == "none")
				conflicting_keywords_.resize(0);
		}
	}

	OpenSMOKE_DictionaryKeyWord::OpenSMOKE_DictionaryKeyWord(std::vector<std::string>& lines)
	{
		if (lines.size() != 7)
			ErrorMessage("Wrong numer of lines defining the keyword");

		for (unsigned int i=0;i<lines.size();i++)
			boost::algorithm::trim(lines[i]);

		SetKeyword(lines[0]);
		SetType(lines[1]);
		SetComment(lines[2]);
		SetIsCompulsory(lines[3]);
		SetCompulsoryAlternatives(lines[4]);
		SetRequiredKeyWords(lines[5]);
		SetConflictingKeywords(lines[6]);
	}

	void OpenSMOKE_DictionaryKeyWord::ShortSummary(std::ostream& fout) const
	{
		fout << name_ << " (" << type_ascii_ << ")" << std::endl;
		fout << " * " << comment_short_ << std::endl;
	}

}
