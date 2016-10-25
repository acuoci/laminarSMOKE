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

#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <sstream>

namespace OpenSMOKE
{
	void OpenSMOKE_DictionaryFile::ErrorMessage(const std::string message) const
	{
		std::cout << "Dictionary " << name_  << " starting at line " << index_first_line_ << " of file " << file_name_ << std::endl;
		std::cout << "Fatal error:     " << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}

	void OpenSMOKE_DictionaryFile::Write(std::ostream& fout) const
	{
		fout << "Dictionary: " << name_ << std::endl;
		for(unsigned int i=0;i<good_lines_.size();i++)
			fout << indices_of_good_lines_[i] << "*" << good_lines_[i] << "*" << std::endl;
		fout << "Dictionary(ext): " << name_ << std::endl;
		for(unsigned int i=0;i<extended_lines_.size();i++)
			fout << "*" << extended_lines_[i] << "*" << std::endl;
	}

	void OpenSMOKE_DictionaryFile::Analyze()
	{
		// Recognize the name of the dictionary
		{
			boost::algorithm::trim(lines_[0]);
			typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
			boost::char_separator<char> sep_blank(" ");
			tokenizer_blank tokens(lines_[0], sep_blank);
			const std::size_t n = std::distance(tokens.begin(), tokens.end());
		
			if (n !=2)
				ErrorMessage("The name of the dictionary is not specified correctly.");

			tokenizer_blank::iterator tok_blank = tokens.begin();
			if (*tok_blank != "Dictionary")
				ErrorMessage("The name of the dictionary is not specified correctly.");
			++tok_blank;
			name_ = *tok_blank;
		}

		// Control the syntax
		std::vector<unsigned int> i_forward;
		std::vector<unsigned int> i_backward;
		{
			unsigned int n_forward = 0;
			unsigned int n_backward = 0;
			for(unsigned int i=0;i<lines_.size();i++)
			{
				typedef boost::find_iterator<std::string::iterator> find_iterator_forward;
				for(find_iterator_forward It=boost::make_find_iterator(lines_[i], boost::first_finder("{", boost::is_iequal())); It!=find_iterator_forward();++It)
				{
					n_forward++;
					i_forward.push_back(i);
				}

				typedef boost::find_iterator<std::string::iterator> find_iterator_backward;
				for(find_iterator_backward It=boost::make_find_iterator(lines_[i], boost::first_finder("}", boost::is_iequal())); It!=find_iterator_backward();++It)
				{
					n_backward++;
					i_backward.push_back(i);
				}
			}

			if (n_forward == n_backward && n_forward > 1)
				ErrorMessage("It seems that sub-dictionaries are defined inside the main dictionary. This is not possible.");

			if (n_forward != 1 || n_backward != 1)
				ErrorMessage("Wrong number of curly brackets {}");
		
			if (i_forward[0] >= i_backward[0])
				ErrorMessage("Wrong position of curly brackets {}");
		}
		
		// Check the curly braces
		{
			boost::algorithm::trim(lines_[i_forward[0]]);
			typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
			boost::char_separator<char> sep_blank(" ");
			tokenizer_blank tokens(lines_[i_forward[0]], sep_blank);
			const std::size_t n = std::distance(tokens.begin(), tokens.end());
		
			if (n != 1 && *tokens.begin() != "{")
				ErrorMessage("The left hand curly bracket { is not specified correctly.");
		}
		// Check the curly braces
		{
			boost::algorithm::trim(lines_[i_backward[0]]);
			typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
			boost::char_separator<char> sep_blank(" ");
			tokenizer_blank tokens(lines_[i_backward[0]], sep_blank);
			const std::size_t n = std::distance(tokens.begin(), tokens.end());
		
			if (n != 1 && *tokens.begin() != "}")
				ErrorMessage("The right hand curly bracket } is not specified correctly.");
		}

		// Select lines
		{
			for(unsigned int i=i_forward[0]+1;i<i_backward[0];i++)
			{
				if (lines_[i].find_first_not_of (' ') != lines_[i].npos)   // has only spaces?
				{
					good_lines_.push_back(lines_[i]);
					indices_of_good_lines_.push_back(index_first_line_+i);
				}
			}

			boost::algorithm::trim(good_lines_[0]);
			if (good_lines_[0].at(0) != '@')
			{
				std::stringstream line_index; line_index << indices_of_good_lines_[0];
				std::string message = "The line " + line_index.str() + " is not written correctly. Remember that a keyword must be preceeded by the @ character.";
				ErrorMessage(message);
			}
		}

		// Analyze and compact the lines
		{
			unsigned int n_start = 0;
			for(unsigned int i=0;i<good_lines_.size();i++)
			{
				typedef boost::find_iterator<std::string::iterator> find_iterator_start;
				for(find_iterator_start It=boost::make_find_iterator(good_lines_[i], boost::first_finder("@", boost::is_iequal())); It!=find_iterator_start();++It)
				{
					n_start++;
					i_start_.push_back(i);
				}
			}
			i_start_.push_back(boost::lexical_cast<int>(good_lines_.size()));

			extended_lines_.resize(n_start);
			for(unsigned int i=0;i<n_start;i++)
				for (unsigned int j=i_start_[i];j<i_start_[i+1];j++)
					extended_lines_[i] += good_lines_[j] + " "; 
			
			for(unsigned int i=0;i<n_start;i++)
			{
				boost::algorithm::trim(extended_lines_[i]);

				// Check if closure char is present
				if (extended_lines_[i].find_first_of (';') == extended_lines_[i].npos )
				{
					std::stringstream line1; line1 << indices_of_good_lines_[i_start_[i]];
					std::stringstream line2; line2 << indices_of_good_lines_[i_start_[i+1]-1];
					
					std::string message;
					if (i_start_[i]==i_start_[i+1]-1)
					{
						message = "The line " + line1.str() + " is not written correctly: no ; character.";
					}
					else
						message = "The lines " + line1.str() + "-" + line2.str() + " are not written correctly: no ; character.";
					ErrorMessage(message);
				}

				// Check if the line starts with a @ char
				if (extended_lines_[i].find_first_of ('@') != 0 )
				{
					std::stringstream line1; line1 << indices_of_good_lines_[i_start_[i]];
					std::stringstream line2; line2 << indices_of_good_lines_[i_start_[i+1]-1];
					
					std::string message;
					if (i_start_[i]==i_start_[i+1]-1)
						message = "The line " + line1.str() + " is not written correctly: A keyword must be preceeded by the @ character.";
					else
						message = "The lines " + line1.str() + "-" + line2.str() + " are not written correctly: A keyword must be preceeded by the @ character.";
					ErrorMessage(message);
				}

				// Check if only one ; char is present at the end
				if (extended_lines_[i].find_first_of (';') != (extended_lines_[i].size()-1) )
				{
					std::stringstream line1; line1 << indices_of_good_lines_[i_start_[i]];
					std::stringstream line2; line2 << indices_of_good_lines_[i_start_[i+1]-1];
			
					std::string message;
					if (i_start_[i]==i_start_[i+1]-1)
						message = "The line " + line1.str() + " is not written correctly: Too many ; characters.";
					else
						message = "The lines " + line1.str() + "-" + line2.str() + " are not written correctly: Too many ; characters.";
					ErrorMessage(message);
				}

				// Check if only one @ char is present as the first character
				if (extended_lines_[i].find_last_of ('@') != 0 )
				{
					std::stringstream line1; line1 << indices_of_good_lines_[i_start_[i]];
					std::stringstream line2; line2 << indices_of_good_lines_[i_start_[i+1]-1];
					
					std::string message;
					if (i_start_[i]==i_start_[i+1]-1)
						message = "The line " + line1.str() + " is not written correctly: Too many @ characters.";
					else
						message = "The lines " + line1.str() + "-" + line2.str() + " are not written correctly: Too many @ characters.";
					ErrorMessage(message);
				}
			}
		}
	}

	void OpenSMOKE_DictionaryFile::Transfer(OpenSMOKE_Dictionary& dictionary)
	{
		std::vector<std::string> keywords(extended_lines_.size());
		std::vector<std::string> options(extended_lines_.size());
		std::vector<unsigned int> startingline(extended_lines_.size());
		std::vector<unsigned int> endingline(extended_lines_.size());
		for(unsigned int i=0;i<extended_lines_.size();i++)
		{
			typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
			boost::char_separator<char> sep_blank(" ");
			tokenizer_blank tokens(extended_lines_[i], sep_blank);
				
			unsigned int count = 0; 
			for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
			{
				if (count == 0) keywords[i] = *tok_iter;
				else            options[i] += *tok_iter + " ";
				count++;
			}

			startingline[i] = indices_of_good_lines_[i_start_[i]];
			endingline[i]   = indices_of_good_lines_[i_start_[i+1]-1];
		}

		for(unsigned int i=0;i<extended_lines_.size();i++)
		{
			boost::algorithm::trim(keywords[i]);
				
			if (options[i].size() == 0)
				keywords[i].erase(keywords[i].size()-1);
			else
			{
				boost::algorithm::trim(options[i]);
				options[i].erase(options[i].size()-1);
			}
			if (options[i].size() != 0)
				boost::algorithm::trim(options[i]);
		}

		// Check for duplicates
		for(unsigned int i=0;i<keywords.size();i++)
			for(unsigned int j=i+1;j<keywords.size();j++)
				if (keywords[i] == keywords[j])
				{
					std::stringstream line1; line1 << startingline[i];
					std::stringstream line2; line2 << startingline[j];
					std::string message = "The same keyword (" + keywords[i] + ") is specified twice at lines " + line1.str() + " and " + line2.str();
					ErrorMessage(message);
				}

		// Set the dictionary
		dictionary.SetName(name_);
		dictionary.SetFileName(file_name_);
		dictionary.SetDictionary(keywords, options, startingline, endingline);
	}

}