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

namespace OpenSMOKE
{
	InputFileDictionary::InputFileDictionary() {
	}

	InputFileDictionary::InputFileDictionary(const std::string file_name) 
	{
			file_name_ = new boost::filesystem::path(file_name);
        
			std::ifstream myfile(file_name.c_str(), std::ios::in);
			CheckIfFileIsOpen(myfile, file_name);

			int count=1;
			std::string line;
			while ( myfile.good() )
			{
					std::getline(myfile,line);

					size_t found=line.find("//");
					if (found!=line.npos)
							line.erase(found);
					
					// replacing tabs with spaces
					boost::replace_all(line, "\t", " ");
					if (line.find_first_not_of (' ') == line.npos)   // has only spaces?
					{
						indices_of_blank_lines_.push_back(count);
					}
					else
					{	
							good_lines_.push_back(line);
							indices_of_good_lines_.push_back(count);
					}
					clean_lines_.push_back(line);
					count++;
			}
			myfile.close();

			number_of_blank_lines_ = boost::lexical_cast<int>(indices_of_blank_lines_.size());
			number_of_good_lines_ = boost::lexical_cast<int>(indices_of_good_lines_.size());
			number_of_lines_ = number_of_blank_lines_ + number_of_good_lines_;
	}

	InputFileDictionary::InputFileDictionary(const InputFileDictionary& orig) {
	}

	InputFileDictionary::~InputFileDictionary() {
	}

	void InputFileDictionary::Status(std::ostream &fOut) const
	{ 
		fOut << "Name:        " << file_name_->leaf() << std::endl;
		fOut << "Path:        " << file_name_->parent_path() << std::endl;
		fOut << "Size:        " << boost::filesystem::file_size(*file_name_)/1000. << " kB" << std::endl;
		fOut << "Lines:       " << number_of_lines_ << std::endl;
		fOut << "Blank lines: " << number_of_blank_lines_ << std::endl;
	}

}
