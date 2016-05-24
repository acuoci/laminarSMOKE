/*-----------------------------------------------------------------------*\
|                                                                         |
|                    ╔═══╦═╗╔═╦═══╦╗╔═╦═══╗                               |
|                    ║╔═╗║║╚╝║║╔═╗║║║╔╣╔══╝                               | 
|   ╔╗╔══╦╗╔╦╦═╗╔══╦═╣╚══╣╔╗╔╗║║ ║║╚╝╝║╚══╗                               |
|   ║║║╔╗║╚╝╠╣╔╗╣╔╗║╔╩══╗║║║║║║║ ║║╔╗║║╔══╝                               |
|   ║╚╣╔╗║║║║║║║║╔╗║║║╚═╝║║║║║║╚═╝║║║╚╣╚══╗                               |
|   ╚═╩╝╚╩╩╩╩╩╝╚╩╝╚╩╝╚═══╩╝╚╝╚╩═══╩╝╚═╩═══╝                               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci                                                     |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of laminarSMOKE solver.                             |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2016, 2015, 2014 A. Cuoci                                |
|   laminarSMOKE is free software: you can redistribute it and/or modify  |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   laminarSMOKE is distributed in the hope that it will be useful,       |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with laminarSMOKE. If not, see <http://www.gnu.org/licenses/>.  |
|                                                                         |
\*-----------------------------------------------------------------------*/

if (calculateRatesAcrossBoundaries == true)
{
	std::vector<double> mass_flow_rates;
	std::vector<std::string> boundary_names;
 
	const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
	{
		scalar mass_flow_rate = sum(phi.boundaryField()[patchi]); 
		mass_flow_rates.push_back(mass_flow_rate);
		boundary_names.push_back(patches[patchi].name());
	}

	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "                                 Mass flow rates                                          " << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	for(int i=0;i<mass_flow_rates.size();i++)
	{
		std::cout << std::setw(25) << std::left << boundary_names[i];
		std::cout << std::setw(18) << std::right << mass_flow_rates[i];
		std::cout << std::endl;
	}
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(25) << std::left << "Net";
	std::cout << std::setw(18) << std::right << std::accumulate(mass_flow_rates.begin(), mass_flow_rates.end(), 0.);
	std::cout << std::endl;
}

