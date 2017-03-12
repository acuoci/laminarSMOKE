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
	template<typename Species>
	TransportReaderPolicy_CHEMKIN<Species>::TransportReaderPolicy_CHEMKIN () {
	}

	template<typename Species>
	TransportReaderPolicy_CHEMKIN<Species>::TransportReaderPolicy_CHEMKIN(const TransportReaderPolicy_CHEMKIN& orig) {
	}

	template<typename Species>
	TransportReaderPolicy_CHEMKIN<Species>::~TransportReaderPolicy_CHEMKIN() {
	}

	template<typename Species>
	bool TransportReaderPolicy_CHEMKIN<Species>::ReadFromASCIIFile(const std::string file_name)
	{
		std::cout << " * Reading transport file... " << std::endl;

		InputFileCHEMKIN myTransport(file_name);
	
		for (unsigned int j=0;j<myTransport.good_lines().size();j++)
		{
			size_t found1=myTransport.good_lines()[j].find("END");
			if (found1!=std::string::npos )	break;
			size_t found2=myTransport.good_lines()[j].find("end");
			if (found2!=std::string::npos )	break;

			try
			{
					Species species;
					std::vector<std::string> lines(1);
					lines[0] = myTransport.good_lines()[j];
					bool success = species.ReadTransportFromStrings(lines);
					if (success == false)
						throw myTransport.indices_of_good_lines()[j];
					species_.insert(make_pair(species.name_transport(), species));
			}
			catch(int k)
			{
				std::cout << "Reading transport properties: error in line " << k << std::endl;
				return false;
			}
		}

		return true;
	}

}