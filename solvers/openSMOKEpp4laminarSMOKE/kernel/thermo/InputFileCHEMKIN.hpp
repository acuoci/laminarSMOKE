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
	InputFileCHEMKIN::InputFileCHEMKIN() {
	}

	InputFileCHEMKIN::InputFileCHEMKIN(const std::string file_name) 
	{
            #if defined __linux || defined __APPLE__
                        
				//Check the file format through the awk command
				const std::string check_string = "if awk  '/\\r$/{exit 0;} 1{exit 1;}' " + file_name + 
					"; then exit 0; else  exit 1; fi; ";
				const int ff_result = system(check_string.c_str()) / 256;

				if (ff_result == 0)
				{
							// Check if the dos2unix application exists
							{
								const std::string exec_test = "which dos2unix >/dev/null";
								const int result_test = system(exec_test.c_str()) / 256;
					
								if (result_test == 0)
								{
									const std::string exec_dos2unix = "dos2unix " + file_name + " 2>/dev/null";
									const int result_dos2unix = system(exec_dos2unix.c_str());
								}
								else
								{
									const std::string exec_sed = "perl -pi -e 's/\r\n|\n|\r/\n/g' " + file_name;
									const int result_sed = system(exec_sed.c_str());
								}
							}
				}
                      
            #endif
			
            file_name_ = new boost::filesystem::path(file_name);
        
			std::ifstream myfile(file_name.c_str(), std::ios::in);
			CheckIfFileIsOpen(myfile, file_name);

			int count=1;
			std::string line;
			while ( myfile.good() )
			{
					std::getline(myfile,line);

					size_t found_comment = line.find("!#");
					std::string comment = "";
					if (found_comment != line.npos)
						comment = line.substr(found_comment);

					size_t found=line.find_first_of("!");
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
							strong_comments_.push_back(comment);
					}
					count++;
			}
			myfile.close();

        
			number_of_blank_lines_ = boost::lexical_cast<int>(indices_of_blank_lines_.size());
			number_of_good_lines_ = boost::lexical_cast<int>(indices_of_good_lines_.size());
			number_of_lines_ = number_of_blank_lines_ + number_of_good_lines_;
	}

	void InputFileCHEMKIN::ConvertGoodLinesIntoBlankLines(const std::vector<unsigned int> lines_to_remove)
	{
		std::vector<unsigned int> indices = lines_to_remove;
		std::sort(indices.begin(), indices.end());
		std::reverse(indices.begin(), indices.end());

		number_of_blank_lines_ += indices.size();
		number_of_good_lines_ -= indices.size();

		for (unsigned int j = 0; j < indices.size(); j++)
		{
			indices_of_blank_lines_.push_back(indices_of_good_lines_[indices[j]]);
			indices_of_good_lines_.erase(indices_of_good_lines_.begin() + (indices[j]-1));
		}

		for (unsigned int j = 0; j<indices.size(); j++)
			good_lines_.erase(good_lines_.begin() + (indices[j]-1));	
	}

	InputFileCHEMKIN::InputFileCHEMKIN(const InputFileCHEMKIN& orig) {
	}

	InputFileCHEMKIN::~InputFileCHEMKIN() {
	}

	void InputFileCHEMKIN::Status(std::ostream &fOut) const
	{ 
		fOut << "Name:        " << file_name_->leaf() << std::endl;
		fOut << "Path:        " << file_name_->parent_path() << std::endl;
		fOut << "Size:        " << boost::filesystem::file_size(*file_name_)/1000. << " kB" << std::endl;
		fOut << "Lines:       " << number_of_lines_ << std::endl;
		fOut << "Blank lines: " << number_of_blank_lines_ << std::endl;
	}

}
