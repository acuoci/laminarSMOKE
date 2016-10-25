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

namespace OpenSMOKE
{	
	void OpenSMOKE_Dictionary::ErrorMessage(const std::string message) const
	{
		std::cout << "Dictionary " << name_ << " defined in file: " << file_name_ << std::endl;
		std::cout << "Fatal error:     " << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}

	void OpenSMOKE_Dictionary::SetDictionary(const std::vector<std::string>& keywords, const std::vector<std::string>& options, 
											 const std::vector<unsigned int>& starting_lines, const std::vector<unsigned int>& ending_lines)
	{
		keywords_ = keywords;
		options_ = options;
		starting_lines_ = starting_lines;
		ending_lines_ = ending_lines;
	}

	void OpenSMOKE_Dictionary::Summary(std::ostream& fout) const
	{
		std::cout << "Dictionary " << name_  << " defined in file " << file_name_ << std::endl;
		for (unsigned int i=0;i<keywords_.size();i++)
			fout << keywords_[i] << "*" << options_[i] << "*" << starting_lines_[i] << "-" << ending_lines_[i] << std::endl;
	}

	void OpenSMOKE_Dictionary::Checks()
	{
		// Check for undefined keywords
		bool undefined_keywords = grammar_.CheckUndefinedKeyWords(keywords_);
		if (undefined_keywords == false)
			ErrorMessage("Error in the dictionary (undefined keywords). See the messages reported above.");

		// Check for compulsory keywords
		bool compulsory_keywords = grammar_.CheckCompulsoryKeyWords(keywords_);
		if (compulsory_keywords == false)
			ErrorMessage("Error in the dictionary (compulsory keywords). See the messages reported above.");

		// Check for required keywords
		bool required_keywords = grammar_.CheckRequiredKeyWords(keywords_);
		if (required_keywords == false)
			ErrorMessage("Error in the dictionary (required keywords). See the messages reported above.");
		
		// Check for conflicting keywords
		bool conflicting_keywords = grammar_.CheckConflictingKeyWords(keywords_);
		if (conflicting_keywords == false)
			ErrorMessage("Error in the dictionary (conflicting keywords). See the messages reported above.");
	}

	void OpenSMOKE_Dictionary::SetGrammar(OpenSMOKE_DictionaryGrammar& grammar)
	{
		grammar.UserDefined();
		grammar_ = grammar;
		Checks();
	}

	void OpenSMOKE_Dictionary::SetGrammar(const std::string file_name)
	{
		grammar_.ReadFromFile(file_name);
		Checks();
	}

	bool OpenSMOKE_Dictionary::CheckOption(const std::string name_of_keyword)
	{
		const std::vector<std::string>::iterator it = find (keywords_.begin(), keywords_.end(), name_of_keyword);
		
		if (it == keywords_.end() )
			return false;
		else
			return true;
	}

	int OpenSMOKE_Dictionary::CheckOption(const std::string name_of_keyword, const OpenSMOKE_DictionaryKeyWordTypes expected_type)
	{
		const std::vector<std::string>::iterator it = find (keywords_.begin(), keywords_.end(), name_of_keyword);
		
		if (it == keywords_.end() )
		{
			std::string message = "The required keyword (" + name_of_keyword + ") is not present in the dictionary.";
			ErrorMessage(message);
			return 0;
		}
		else
		{
			size_t index = it-keywords_.begin();

			if ( grammar_.CheckType(name_of_keyword, expected_type) == false)
			{
				std::stringstream line; line << starting_lines_[index];
				std::string message = "Error in the keyword at line : " + line.str() + "\n";
				message += "The required keyword (" + name_of_keyword + ") is of different type";
				std::cout << "Expected type: " << expected_type << std::endl;
				ErrorMessage(message);
			}

			return boost::lexical_cast<int>(index);
		}
	}

	void OpenSMOKE_Dictionary::ReadDouble(const std::string name_of_keyword, double& value)
	{
		size_t index = CheckOption(name_of_keyword, SINGLE_DOUBLE);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) != 1)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a number (double).";
			ErrorMessage(message);
		}

		try
		{
			const std::string number = *tokens.begin();
			value = boost::lexical_cast<double>(number.c_str());
		}
		catch(boost::bad_lexical_cast &)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "Failure in the numerical conversion.";
			ErrorMessage(message);
		}
	}

	void OpenSMOKE_Dictionary::ReadMeasure(const std::string name_of_keyword, double& value, std::string& units)
	{
		size_t index = CheckOption(name_of_keyword, SINGLE_MEASURE);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) != 2)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a measure (double + string).";
			ErrorMessage(message);
		}

		try
		{
			const std::string number = *tokens.begin();
			value = boost::lexical_cast<double>(number.c_str());
		}
		catch(boost::bad_lexical_cast &)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "Failure in the numerical conversion.";
			ErrorMessage(message);
		}
		tokenizer_blank::iterator tok_iter = tokens.begin();
		++tok_iter;
		units = *tok_iter;
	}

	void OpenSMOKE_Dictionary::ReadInt(const std::string name_of_keyword, int& value)
	{
		size_t index = CheckOption(name_of_keyword, SINGLE_INT);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) != 1)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a number (int).";
			ErrorMessage(message);
		}

		try
		{
			const std::string number = *tokens.begin();
			value = boost::lexical_cast<int>(number.c_str());
		}
		catch(boost::bad_lexical_cast &)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "Failure in the numerical conversion.";
			ErrorMessage(message);
		}
	}

	void OpenSMOKE_Dictionary::ReadPath(const std::string name_of_keyword, boost::filesystem::path& value)
	{
		size_t index = CheckOption(name_of_keyword, SINGLE_PATH);

		try
		{
			std::string pathcomplete = options_[index];
			value = boost::filesystem::path(pathcomplete);                        
		}
		catch(boost::bad_lexical_cast &)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "Failure in reading the file name";
			ErrorMessage(message);
		}
	}

	void OpenSMOKE_Dictionary::ReadDictionary(const std::string name_of_keyword, std::string& value)
	{
		size_t index = CheckOption(name_of_keyword, SINGLE_DICTIONARY);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) != 1)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a dictionary (string)";
			ErrorMessage(message);
		}

		value = *tokens.begin();
	}

	void OpenSMOKE_Dictionary::ReadSequence(const std::string name_of_keyword, std::string& value)
	{
		size_t index = CheckOption(name_of_keyword, SEQUENCE_STRING);
		
		value = options_[index];
	}

	void OpenSMOKE_Dictionary::ReadBool(const std::string name_of_keyword, bool& value)
	{
		size_t index = CheckOption(name_of_keyword, SINGLE_BOOL);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) != 1)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a number (int).";
			ErrorMessage(message);
		}

		const std::string number = *tokens.begin();
		if (number == "true" || number == "on") value = true;
		else if (number == "false" || number == "off") value = false;
		else
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "Only boolean types can be accepted.";
			ErrorMessage(message);
		}
	}

	void OpenSMOKE_Dictionary::ReadChar(const std::string name_of_keyword, char& value)
	{
		size_t index = CheckOption(name_of_keyword, SINGLE_CHAR);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) != 1)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a char";
			ErrorMessage(message);
		}

		try
		{
			const std::string character = *tokens.begin();
			value = boost::lexical_cast<char>(character.c_str());
		}
		catch(boost::bad_lexical_cast &)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "Failure in the character conversion.";
			ErrorMessage(message);
		}
	}

	void OpenSMOKE_Dictionary::ReadString(const std::string name_of_keyword, std::string& value)
	{
		size_t index = CheckOption(name_of_keyword, SINGLE_STRING);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) != 1)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a string";
			ErrorMessage(message);
		}

		value = *tokens.begin();
	}

	void OpenSMOKE_Dictionary::ReadOption(const std::string name_of_keyword, std::vector<std::string>& values)
	{
		size_t index = CheckOption(name_of_keyword, VECTOR_STRING);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) == 0)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a list of strings";
			ErrorMessage(message);
		}

		values.resize(std::distance (tokens.begin(), tokens.end()));
		unsigned int count = 0;
		for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
			values[count++] = *tok_iter;
	}

	void OpenSMOKE_Dictionary::ReadOption(const std::string name_of_keyword, std::vector<double>& values)
	{
		size_t index = CheckOption(name_of_keyword, VECTOR_DOUBLE);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) == 0)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a list of numbers (double).";
			ErrorMessage(message);
		}

		values.resize(std::distance (tokens.begin(), tokens.end()));
		unsigned int count = 0;
		for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
		{
			try
			{
				const std::string number = *tok_iter;
				values[count++]  = boost::lexical_cast<double>(number.c_str());
			}
			catch(boost::bad_lexical_cast &)
			{
				std::stringstream line; line << starting_lines_[index];
				std::string message = "Error in the keyword at line : " + line.str() + "\n";
				message += "Failure in the numerical conversion.";
				ErrorMessage(message);
			}
		}
	}

	void OpenSMOKE_Dictionary::ReadOption(const std::string name_of_keyword, std::vector<int>& values)
	{
		size_t index = CheckOption(name_of_keyword, VECTOR_INT);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) == 0)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a list of numbers (int).";
			ErrorMessage(message);
		}

		values.resize(std::distance (tokens.begin(), tokens.end()));
		unsigned int count = 0;
		for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
		{
			try
			{
				const std::string number = *tok_iter;
				values[count++]  = boost::lexical_cast<int>(number.c_str());
			}
			catch(boost::bad_lexical_cast &)
			{
				std::stringstream line; line << starting_lines_[index];
				std::string message = "Error in the keyword at line : " + line.str() + "\n";
				message += "Failure in the numerical conversion.";
				ErrorMessage(message);
			}
		}
	}

	void OpenSMOKE_Dictionary::ReadOption(const std::string name_of_keyword, std::vector<char>& values)
	{
		size_t index = CheckOption(name_of_keyword, VECTOR_CHAR);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		if ( std::distance (tokens.begin(), tokens.end()) == 0)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a list of chars.";
			ErrorMessage(message);
		}

		values.resize(std::distance (tokens.begin(), tokens.end()));
		unsigned int count = 0;
		for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
		{
			try
			{
				const std::string number = *tok_iter;
				values[count++]  = boost::lexical_cast<char>(number.c_str());
			}
			catch(boost::bad_lexical_cast &)
			{
				std::stringstream line; line << starting_lines_[index];
				std::string message = "Error in the keyword at line : " + line.str() + "\n";
				message += "Failure in the character conversion.";
				ErrorMessage(message);
			}
		}
	}

	void OpenSMOKE_Dictionary::ReadOption(const std::string name_of_keyword, std::vector<std::string>& names, std::vector<double>& values)
	{
		size_t index = CheckOption(name_of_keyword, VECTOR_STRING_DOUBLE);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		
		if (n == 0 || n%2 !=0)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a list of string and numbers (double)";
			ErrorMessage(message);
		}


		names.resize(n/2);
		values.resize(n/2);
		
		unsigned int count = 0;
		for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
		{
			names[count] = *tok_iter;
			++tok_iter;

			try
			{
				const std::string number = *tok_iter;
				values[count]  = boost::lexical_cast<double>(number.c_str());
			}
			catch(boost::bad_lexical_cast &)
			{
				std::stringstream line; line << starting_lines_[index];
				std::string message = "Error in the keyword at line : " + line.str() + "\n";
				message += "Failure in the numerical conversion.";
				ErrorMessage(message);
			}

			count++;
		}
	}
        
        void OpenSMOKE_Dictionary::ReadOption(const std::string name_of_keyword, std::vector<double>& values, std::vector<std::string>& names)
	{
		size_t index = CheckOption(name_of_keyword, VECTOR_MEASURE);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens(options_[index], sep_blank);
		const std::size_t n = std::distance(tokens.begin(), tokens.end());
		
		if (n == 0 || n%2 !=0)
		{
			std::stringstream line; line << starting_lines_[index];
			std::string message = "Error in the keyword at line : " + line.str() + "\n";
			message += "The required keyword (" + name_of_keyword + ") requires a list of numbers (double) and strings";
			ErrorMessage(message);
		}


		names.resize(n/2);
		values.resize(n/2);
		
		unsigned int count = 0;
		for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
		{
			try
			{
				const std::string number = *tok_iter;
				values[count]  = boost::lexical_cast<double>(number.c_str());
			}
			catch(boost::bad_lexical_cast &)
			{
				std::stringstream line; line << starting_lines_[index];
				std::string message = "Error in the keyword at line : " + line.str() + "\n";
				message += "Failure in the numerical conversion.";
				ErrorMessage(message);
			}
                        
                        ++tok_iter;
                        
                        names[count] = *tok_iter;

			count++;
		}
	}

}
