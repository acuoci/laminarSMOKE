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
	std::vector<double> omega_C_rates;
	std::vector<double> omega_H_rates;
	std::vector<double> omega_O_rates;

	std::vector<double> omega_C_rates_inlet;
	std::vector<double> omega_H_rates_inlet;
	std::vector<double> omega_O_rates_inlet;

	std::vector<std::string> boundary_names;

	const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
	{
		boundary_names.push_back(patches[patchi].name());
		scalar omega_O_rate = sum(phi.boundaryField()[patchi]*omega_O.boundaryField()[patchi]); 
		scalar omega_C_rate = sum(phi.boundaryField()[patchi]*omega_C.boundaryField()[patchi]); 
		scalar omega_H_rate = sum(phi.boundaryField()[patchi]*omega_H.boundaryField()[patchi]); 
		
		omega_O_rates.push_back( omega_O_rate);
		omega_C_rates.push_back( omega_C_rate);
		omega_H_rates.push_back( omega_H_rate);

		omega_O_rates_inlet.push_back( std::min( omega_O_rate, 0.) );
		omega_C_rates_inlet.push_back( std::min( omega_C_rate, 0.) );
		omega_H_rates_inlet.push_back( std::min( omega_H_rate, 0.) );

	}

	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(25) << std::left  << "";
	std::cout << std::setw(20) << std::right << "H atom";
	std::cout << std::setw(20) << std::right << "O atom";
	std::cout << std::setw(20) << std::right << "C atom";
	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	for(int i=0;i<boundary_names.size();i++)
	{
		std::cout << std::setw(25) << std::left << boundary_names[i];
		std::cout << std::setw(20) << std::right << omega_H_rates[i];
		std::cout << std::setw(20) << std::right << omega_O_rates[i];
		std::cout << std::setw(20) << std::right << omega_C_rates[i];
		std::cout << std::endl;
	}
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(25) << std::left << "Net";
	std::cout << std::setw(20) << std::right << std::accumulate(omega_H_rates.begin(), omega_H_rates.end(), 0.);
	std::cout << std::setw(20) << std::right << std::accumulate(omega_O_rates.begin(), omega_O_rates.end(), 0.);
	std::cout << std::setw(20) << std::right << std::accumulate(omega_C_rates.begin(), omega_C_rates.end(), 0.);
	std::cout << std::endl;

	std::cout << std::setw(25) << std::left << "Error(%)";
	std::cout << std::setw(20) << std::right << std::accumulate(omega_H_rates.begin(), omega_H_rates.end(), 0.) /
						    ( std::accumulate(omega_H_rates_inlet.begin(), omega_H_rates_inlet.end(), 0.) + 1.e-32) * 100.;
	std::cout << std::setw(20) << std::right << std::accumulate(omega_O_rates.begin(), omega_O_rates.end(), 0.) /
						    ( std::accumulate(omega_O_rates_inlet.begin(), omega_O_rates_inlet.end(), 0.) + 1.e-32) * 100.;
	std::cout << std::setw(20) << std::right << std::accumulate(omega_C_rates.begin(), omega_C_rates.end(), 0.) /
						    ( std::accumulate(omega_C_rates_inlet.begin(), omega_C_rates_inlet.end(), 0.) + 1.e-32) * 100.;
	std::cout << std::endl;

}
