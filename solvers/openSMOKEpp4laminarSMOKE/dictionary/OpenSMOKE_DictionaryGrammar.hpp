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
	void OpenSMOKE_DictionaryGrammar::ErrorMessage(const std::string message) const
	{
		std::cout << "Grammar defined in file " << file_name_->leaf() << std::endl;
		std::cout << "Fatal error:     " << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}

	void OpenSMOKE_DictionaryGrammar::DefineRules()
	{
		ErrorMessage("No grammar rules were defined by the user.");
	}

	void OpenSMOKE_DictionaryGrammar::UserDefined()
	{
		DefineRules();

		// Number of keywords
		number_of_keywords_ = boost::lexical_cast<unsigned int>(keywords_.size());
		
		// List of keywords names
		list_of_keywords_names_.resize(number_of_keywords_);
		for (unsigned int i=0;i<number_of_keywords_;i++)
			list_of_keywords_names_[i] = keywords_[i].name();

		// Check the grammar
		Check();
	}

	void OpenSMOKE_DictionaryGrammar::AddKeyWord(const OpenSMOKE_DictionaryKeyWord& keyword)
	{
		keywords_.push_back(keyword);
	}

	void OpenSMOKE_DictionaryGrammar::ReadFromFile(const std::string file_name)
	{
		file_name_ = new boost::filesystem::path(file_name);
        
		std::ifstream fInput(file_name.c_str(), std::ios::in);
		CheckIfFileIsOpen(fInput, file_name);

		// Reads the lines
		std::vector<std::string> lines;
		while ( fInput.good() )
		{
			std::string line;
			std::getline(fInput,line);
			lines.push_back(line);
		}
		fInput.close();
		
		// Count number of keywords
		std::vector<unsigned int> line_keywords;
		for(unsigned int i=0;i<lines.size();i++)
		{
			typedef boost::find_iterator<std::string::iterator> find_iterator_dictionaries;
			for(find_iterator_dictionaries It=boost::make_find_iterator(lines[i], boost::first_finder("keyword:", boost::is_iequal())); It!=find_iterator_dictionaries();++It)
				line_keywords.push_back(i+1);
		}
		number_of_keywords_ = boost::lexical_cast<unsigned int>(line_keywords.size());

		// Reads the grammar rules
		keywords_.resize(number_of_keywords_);
		for(unsigned int i=0;i<number_of_keywords_;i++)
		{
			std::vector<std::string> block_of_lines;
			for(unsigned int j=line_keywords[i]-1;j<line_keywords[i]-1+7;j++)
				block_of_lines.push_back(lines[j]);
			OpenSMOKE_DictionaryKeyWord tmp(block_of_lines);
			keywords_[i] = tmp;
		}

		// List of keywords names
		list_of_keywords_names_.resize(number_of_keywords_);
		for (unsigned int i=0;i<number_of_keywords_;i++)
			list_of_keywords_names_[i] = keywords_[i].name();

		// Check the grammar
		Check();
	}

	void OpenSMOKE_DictionaryGrammar::Check()
	{
		// Check if one option is specified more than once
		for(unsigned int i=0;i<number_of_keywords_;i++)
			for(unsigned int j=i+1;j<number_of_keywords_;j++)
				if ( keywords_[i].name() == keywords_[j].name() )
					FatalErrorMessage("The following keyword is defined more than once: " + keywords_[i].name() );

		// Check the existence of additional keywords
		CheckExistence();

		// Check the relationships between compulsory keywords
		for(unsigned int j=0;j<number_of_keywords_;j++)
		{
			if (keywords_[j].is_compulsory() == true)
			{	
				if (keywords_[j].compulsory_alternatives().size() != 0)
				{
					for(unsigned int k=0;k<keywords_[j].compulsory_alternatives().size();k++)
						for(unsigned int i=0;i<number_of_keywords_;i++)
							if (keywords_[i].name() == keywords_[j].compulsory_alternatives()[k])
							{
								if (keywords_[i].is_compulsory() == false)
									ErrorMessage("Conflicting compulsory options between the following keywords: " + keywords_[j].name() + " and " + keywords_[i].name());
									
								bool found = false; 
								for(unsigned int kk=0;kk<keywords_[i].compulsory_alternatives().size();kk++)
									if (keywords_[i].compulsory_alternatives()[kk] == keywords_[j].name())
									{
										found = true;
										break;
									}
								if (found==false)
									ErrorMessage("Unconsistent compulsory options between the following keywords: " + keywords_[j].name() + " and " + keywords_[i].name());
							}
				}
			}
		}	

		// Check the relationships between conflicting keywords
		for(unsigned int j=0;j<number_of_keywords_;j++)
		{
			if (keywords_[j].conflicting_keywords().size() != 0)
			{
				for(unsigned int k=0;k<keywords_[j].conflicting_keywords().size();k++)
					for(unsigned int i=0;i<number_of_keywords_;i++)
						if (keywords_[i].name() == keywords_[j].conflicting_keywords()[k])
						{									
							bool found = false; 
							for(unsigned int kk=0;kk<keywords_[i].conflicting_keywords().size();kk++)
								if (keywords_[i].conflicting_keywords()[kk] == keywords_[j].name())
								{
									found = true;
									break;
								}
							if (found==false)
								ErrorMessage("Unconsistent conflicting options between the following keywords: " + keywords_[j].name() + " and " + keywords_[i].name());
						}
			}
		}	

	}

	void OpenSMOKE_DictionaryGrammar::ShortSummary(std::ostream& fout) const
	{
		fout << "Grammar defined in file: " << file_name_->leaf() << std::endl;
		fout << "-----------------------------------------------------------------------------------------------------" << std::endl;
		fout << std::endl;
		for(unsigned int i=0;i<number_of_keywords_;i++)
		{
			keywords_[i].ShortSummary(fout);
			fout << std::endl;
		}
		fout << "-----------------------------------------------------------------------------------------------------" << std::endl;
	}

	bool OpenSMOKE_DictionaryGrammar::CheckUndefinedKeyWords(std::vector<std::string>& list_of_keywords)
	{
		bool global_error = true;
		for(unsigned int i=0;i<list_of_keywords.size();i++)
		{
			bool is_defined = false;
			for(unsigned int j=0;j<number_of_keywords_;j++)
				if (list_of_keywords[i] == keywords_[j].name())
				{
					is_defined = true;
					break;
				}
			
			if (is_defined == false)
			{
				global_error = false;
				std::cout << "The " << list_of_keywords[i] << " type is not recognized as a keyword." << std::endl;
			}
		}
		return global_error;
	}

	bool OpenSMOKE_DictionaryGrammar::CheckCompulsoryKeyWords(std::vector<std::string>& list_of_keywords)
	{
		bool global_error = true;
		for(unsigned int j=0;j<number_of_keywords_;j++)
		{
			if (keywords_[j].is_compulsory() == true)
			{	
				if (keywords_[j].compulsory_alternatives().size() == 0)
				{
					bool is_defined = false;
					for(unsigned int i=0;i<list_of_keywords.size();i++)
						if (list_of_keywords[i] == keywords_[j].name())
						{
							is_defined = true;
							break;
						}

					if (is_defined == false)
					{
						global_error = false;
						std::cout << "The compulsory keyword " << keywords_[j].name() << " is not defined." << std::endl;
					}
				}
				else
				{
					unsigned int number_of_alternatives_found = 0;
					for(unsigned int i=0;i<list_of_keywords.size();i++)
						if (list_of_keywords[i] == keywords_[j].name())
						{
							number_of_alternatives_found++;
							break;
						}
					
					for(unsigned int k=0;k<keywords_[j].compulsory_alternatives().size();k++)
						for(unsigned int i=0;i<list_of_keywords.size();i++)
								if (list_of_keywords[i] == keywords_[j].compulsory_alternatives()[k])
								{
									number_of_alternatives_found++;
									break;
								}

					if (number_of_alternatives_found == 0)	
					{
						global_error = false;
						std::cout << "Neither the compulsory keyword " << keywords_[j].name() << ", neither its alternatives, are defined." << std::endl;
						return global_error;
					}
					else if (number_of_alternatives_found > 1)	
					{
						global_error = false;
						std::cout << "The compulsory keyword " << keywords_[j].name() << ", or one of its alternatives, are defined more than once." << std::endl;
						return global_error;
					}
				}
			}
		}
			
		return global_error;
	}

	bool OpenSMOKE_DictionaryGrammar::CheckRequiredKeyWords(std::vector<std::string>& list_of_keywords)
	{
		bool global_error = true;
		for(unsigned int j=0;j<number_of_keywords_;j++)
		{
			// Check only the keywords requiring additional keywords
			if (keywords_[j].required_keywords().size() != 0)
			{
				// Check if this keyword is used in the dictionary
				std::vector<std::string>::iterator index = find (list_of_keywords.begin(), list_of_keywords.end(), keywords_[j].name());
				if (index != list_of_keywords.end() )
				{
					// Loop over all the required keywords
					for(unsigned int k=0;k<keywords_[j].required_keywords().size();k++)
					{
						bool is_defined = false;
						std::vector<std::string>::iterator it = find (list_of_keywords.begin(), list_of_keywords.end(), keywords_[j].required_keywords()[k]);
						if (it != list_of_keywords.end() )
							is_defined = true;

						if (is_defined == false)
						{
							global_error = false;
							std::cout << "The " << keywords_[j].name() << " keyword requires the " << keywords_[j].required_keywords()[k] << " keyword, which is not present in the dictionary." << std::endl;
						}
					}
				}
			}
		}
		return global_error;
	}

	bool OpenSMOKE_DictionaryGrammar::CheckConflictingKeyWords(std::vector<std::string>& list_of_keywords)
	{
		bool global_error = true;
		for(unsigned int j=0;j<number_of_keywords_;j++)
		{
			// Check only the keywords with conflicting_keywords option enabled on
			if (keywords_[j].conflicting_keywords().size() != 0)
			{
				// Check if this keyword is used in the dictionary
				std::vector<std::string>::iterator index = find (list_of_keywords.begin(), list_of_keywords.end(), keywords_[j].name());
				if (index != list_of_keywords.end() )
				{
					// Loop over all the conflicting keywords
					for(unsigned int k=0;k<keywords_[j].conflicting_keywords().size();k++)
					{
						bool is_found = false;
						std::vector<std::string>::iterator it = find (list_of_keywords.begin(), list_of_keywords.end(), keywords_[j].conflicting_keywords()[k]);
						if (it != list_of_keywords.end() )
							is_found = true;

						if (is_found == true)
						{
							global_error = false;
							std::cout << "The " << keywords_[j].name() << " and the " << keywords_[j].conflicting_keywords()[k] << " keywords are used together, but they are mutually exclusive." << std::endl;
						}
					}
				}
			}
		}

		return global_error;
	}

	void OpenSMOKE_DictionaryGrammar::CheckExistence()
	{
		for(unsigned int j=0;j<number_of_keywords_;j++)
		{
			// Compulsory alternatives
			if (keywords_[j].compulsory_alternatives().size() != 0)
			{
				for(unsigned int k=0;k<keywords_[j].compulsory_alternatives().size();k++)
				{
					bool is_found = false;
					std::vector<std::string>::iterator it = find (list_of_keywords_names_.begin(), list_of_keywords_names_.end(), keywords_[j].compulsory_alternatives()[k]);
					if (it != list_of_keywords_names_.end() )
						is_found = true;
					
					if (is_found ==false)
						ErrorMessage("The undefined " + keywords_[j].compulsory_alternatives()[k] + " keyword is present inside the definition of the " + keywords_[j].name() + " keyword");
				}
			}

			// Required keywords
			if (keywords_[j].required_keywords().size() != 0)
			{
				for(unsigned int k=0;k<keywords_[j].required_keywords().size();k++)
				{
					bool is_found = false;
					std::vector<std::string>::iterator it = find (list_of_keywords_names_.begin(), list_of_keywords_names_.end(), keywords_[j].required_keywords()[k]);
					if (it != list_of_keywords_names_.end() )
						is_found = true;
					
					if (is_found ==false)
						ErrorMessage("The undefined " + keywords_[j].required_keywords()[k] + " keyword is present inside the definition of the " + keywords_[j].name() + " keyword");
				}
			}

			// Conflicting keywords
			if (keywords_[j].conflicting_keywords().size() != 0)
			{
				for(unsigned int k=0;k<keywords_[j].conflicting_keywords().size();k++)
				{
					bool is_found = false;
					std::vector<std::string>::iterator it = find (list_of_keywords_names_.begin(), list_of_keywords_names_.end(), keywords_[j].conflicting_keywords()[k]);
					if (it != list_of_keywords_names_.end() )
						is_found = true;
					
					if (is_found ==false)
						ErrorMessage("The undefined " + keywords_[j].conflicting_keywords()[k] + " keyword is present inside the definition of the " + keywords_[j].name() + " keyword");
				}
			}
		}
	}

	bool OpenSMOKE_DictionaryGrammar::CheckType(const std::string keyword_name, const OpenSMOKE_DictionaryKeyWordTypes expected_type)
	{
		std::vector<std::string>::iterator it = find (list_of_keywords_names_.begin(), list_of_keywords_names_.end(), keyword_name);
		if (it == list_of_keywords_names_.end() )
			return false;
		else
		{
			size_t index = it-list_of_keywords_names_.begin();
			return (keywords_[index].type() == expected_type);
		}
	}
}
