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
	std::vector<double> enthalpy_rates;
	std::vector<double> enthalpy_rates_inlet;	
	std::vector<std::string> boundary_names;

	const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
	{
		boundary_names.push_back(patches[patchi].name());
		scalar enthalpy_rate = sum(phi.boundaryField()[patchi]*h.boundaryField()[patchi]); 
		enthalpy_rates.push_back(enthalpy_rate);
		enthalpy_rates_inlet.push_back( std::min( enthalpy_rate, 0.) );
	}

	std::cout << std::endl;
        std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "                                       Enthalpy                                           " << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	for(int i=0;i<enthalpy_rates.size();i++)
	{
		std::cout << std::setw(25) << std::left << boundary_names[i];
		std::cout << std::setw(20) << std::right << enthalpy_rates[i];
		std::cout << std::endl;
	}
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;

	std::cout << std::setw(25) << std::left << "Net";
	std::cout << std::setw(20) << std::right << std::accumulate(enthalpy_rates.begin(), enthalpy_rates.end(), 0.);
	std::cout << std::endl;

	std::cout << std::setw(25) << std::left << "Error(%)";
	std::cout << std::setw(20) << std::right << std::accumulate(enthalpy_rates.begin(), enthalpy_rates.end(), 0.) /
						    ( std::accumulate(enthalpy_rates_inlet.begin(), enthalpy_rates_inlet.end(), 0.) + 1.e-32) * 100.;
	std::cout << std::endl;
}
