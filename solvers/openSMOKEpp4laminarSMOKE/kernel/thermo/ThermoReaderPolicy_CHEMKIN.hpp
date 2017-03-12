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

#include "kernel/thermo/InputFileCHEMKIN.h"

namespace OpenSMOKE
{
	template<typename Species>
	ThermoReaderPolicy_CHEMKIN<Species>::ThermoReaderPolicy_CHEMKIN() {
	}

	template<typename Species>
	ThermoReaderPolicy_CHEMKIN<Species>::ThermoReaderPolicy_CHEMKIN(const ThermoReaderPolicy_CHEMKIN& orig) {
	}

	template<typename Species>
	ThermoReaderPolicy_CHEMKIN<Species>::~ThermoReaderPolicy_CHEMKIN() {
	}

	template<typename Species>
	bool ThermoReaderPolicy_CHEMKIN<Species>::ReadFromASCIIFile(const std::string file_name)
	{
		std::cout << " * Reading thermodynamic file... " << std::endl;
		InputFileCHEMKIN myThermo(file_name);
		//myThermo.Status(std::cout);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_only_blanks;
		boost::char_separator<char> sep(" ");

                OpenSMOKE::OpenSMOKEVectorDouble temperature_limits_;
		{
			tokenizer_only_blanks tokens(myThermo.good_lines()[0], sep);
			const std::size_t count = std::distance( tokens.begin(), tokens.end() );
                        
			if (count > 2)
			{
				std::cout << "Expected THERMO || THERMO ALL. Found: " << myThermo.good_lines()[0] << std::endl;
				return false;        
			}

			if (*tokens.begin() != "THERMO" && *tokens.begin() != "thermo")
			{
				std::cout << "Expected THERMO || thermo. Found: " << *tokens.begin() << std::endl;
				return false;
			}

			if (count == 2)
			{
                            unsigned int i = 1;
                            for (tokenizer_only_blanks::iterator tok_iter=tokens.begin();tok_iter != tokens.end(); ++tok_iter)
                            {   
                                if (i==2)
                                {
                                    if( *tok_iter!="ALL"  && *tok_iter!="all" )
                                    {
                                        std::cout << "Expected ALL || all. Found: " << *tok_iter << std::endl;
                                        return false;
                                    }
                                }
                                i++;
                            }
			}
		}

		{       
                        ChangeDimensions(3, &temperature_limits_, true);
        
			tokenizer_only_blanks tokens(myThermo.good_lines()[1], sep);
			const std::size_t count = std::distance( tokens.begin(), tokens.end() );
			if (count != 3)
			{
				std::cout << "Wrong number of temperature ranges." << std::endl;
				return false;        
			}
			try
			{
				tokenizer_only_blanks::iterator tok_iter = tokens.begin();
				temperature_limits_[1] = boost::lexical_cast<double>(*tok_iter++);  // low
				temperature_limits_[2] = boost::lexical_cast<double>(*tok_iter++);  // medium
                                temperature_limits_[3] = boost::lexical_cast<double>(*tok_iter++);  // high
			}
			catch(boost::bad_lexical_cast &)
			{
				std::cout << "Numerical conversion failure of temperature ranges." << std::endl;
				return false;
			}
			if (temperature_limits_[1] >= temperature_limits_[2] || temperature_limits_[2] >= temperature_limits_[3])
			{
				std::cout << "Wrong temperature ranges." << std::endl;
				return false;
			}
		} 

		std::cout << " * Analyzing thermodynamic file... " << std::endl;
		{
			int index_of_end_file_ = 0;
			std::vector<int> indices_of_starting_species_;
			int offsets[] = {15};
			boost::offset_separator f(offsets, offsets+1, false);
			for (unsigned int j=2;j<myThermo.good_lines().size();j++)
			{                            
				boost::tokenizer<boost::offset_separator> tokens(myThermo.good_lines()[j],f);
				std::string word=*tokens.begin();
				boost::algorithm::trim(word);
				try
				{
					boost::lexical_cast<double>(word);
				}
				catch(boost::bad_lexical_cast &)
				{
					tokenizer_only_blanks tok(word, sep);

					if (*tok.begin()=="END" || *tok.begin()=="end" )
					{
						if (std::distance( tok.begin(), tok.end() ) == 1)
						{
							 index_of_end_file_ = j;
							 break;
						}
						else
						{
							std::cout << "Expected END. Found: " << myThermo.good_lines()[j] << std::endl;
							return false;
						}
					}
                
					if (*tok.begin()!="TEMP" && *tok.begin()!="temp")
					{
						if (myThermo.good_lines()[j].size()>=80)
							if (myThermo.good_lines()[j].at(79) == '1')
								indices_of_starting_species_.push_back(j);
					}
				}
			}

			if (index_of_end_file_ != 0)
				indices_of_starting_species_.push_back(index_of_end_file_);
			else
			{
				indices_of_starting_species_.push_back(myThermo.number_of_good_lines());
			}
            
			std::cout << "     number of species available: " << indices_of_starting_species_.size()-1 << std::endl;
        
			for (unsigned int i=0;i<indices_of_starting_species_.size()-1;i++)
			{
				std::vector<std::string> block(indices_of_starting_species_[i+1]-indices_of_starting_species_[i]);
				for(int j=0;j<indices_of_starting_species_[i+1]-indices_of_starting_species_[i];j++)
					block[j] = myThermo.good_lines()[indices_of_starting_species_[i]+j];

				try
				{
						Species species;
						bool success = species.ReadThermodynamicsFromStrings(block, temperature_limits_);
						if (success == false)
								throw indices_of_starting_species_[i];
						species_.insert(make_pair(species.name_thermo(), species));
				}
				catch(int k)
				{
					std::cout << "Reading thermodynamic properties: error in species starting at line " << myThermo.indices_of_good_lines()[k] << std::endl;
					return false;
				}
			}
		}

		return true;
	}
}
