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

#include "OpenSMOKE_DictionaryFile.h"

namespace OpenSMOKE
{	
	void OpenSMOKE_DictionaryManager::ErrorMessage(const std::string message) const
	{
		std::cout << "OpenSMOKE_DictionaryManager Class" << std::endl;
		std::cout << "Fatal error:     " << message << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}

	void OpenSMOKE_DictionaryManager::ReadDictionariesFromFile(const std::string file_name)
	{
		OpenSMOKE::InputFileDictionary dict_file(file_name);

		const int number_of_block_lines = boost::lexical_cast<int>(dict_file.clean_lines().size());
		
		// Count number of dictionaries
		unsigned int n_dictionaries = 0;
		std::vector<unsigned int> line_dictionaries;
		for(int i=0;i<number_of_block_lines;i++)
		{
			typedef boost::find_iterator<std::string::iterator> find_iterator_dictionaries;
			for(find_iterator_dictionaries It=boost::make_find_iterator(dict_file.clean_lines()[i], boost::first_finder("Dictionary", boost::is_iequal())); It!=find_iterator_dictionaries();++It)
			{
				n_dictionaries++;
				line_dictionaries.push_back(i+1);
			}
		}
		line_dictionaries.push_back(number_of_block_lines+1);

		if (n_dictionaries == 0)
		{
			std::cout << "No dictionaries are specified in the current file..." << std::endl;
			std::cout << "Press enter to exit..." << std::endl;
			getchar();
			exit(OPENSMOKE_FATAL_ERROR_EXIT);
		}

		for (unsigned int i = 0; i<n_dictionaries; i++)
		{
			OpenSMOKE_DictionaryFile dictionary_file;
			
			dictionary_file.SetFileName(file_name);
			dictionary_file.SetFirstLine(line_dictionaries[i]);
			for (unsigned int j = line_dictionaries[i] - 1; j<line_dictionaries[i + 1] - 1; j++)
				dictionary_file.AddLine(dict_file.clean_lines()[j]);
			dictionary_file.Analyze();

			OpenSMOKE_Dictionary dictionary;
			dictionary_file.Transfer(dictionary);

			if (map_of_dictionaries_.insert(std::pair<std::string, OpenSMOKE_Dictionary>(dictionary.name(), dictionary)).second == false)
				std::cout << "Insertion of new dictionary failed. Key was present." << std::endl;
		}
	}

	OpenSMOKE_Dictionary& OpenSMOKE_DictionaryManager::operator() (const std::string& name)
	{
		std::map<std::string, OpenSMOKE_Dictionary>::iterator dict;
		dict = map_of_dictionaries_.find(name);

		if(dict == map_of_dictionaries_.end())
			ErrorMessage("The " + name + " dictionary was not defined.");
		
		return dict->second;
	}
}
