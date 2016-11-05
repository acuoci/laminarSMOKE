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

if (exportDisks == true)
{
	label ns = Y.size();

	OpenSMOKE::OpenSMOKEVectorDouble y(ns);
	OpenSMOKE::OpenSMOKEVectorDouble x(ns);

	const scalarField& TCells = T.internalField();
	
	for (unsigned int i=1;i<100;i++)
	{
		std::stringstream id;	id << i;

		std::string north_name = "North." + id.str();
		std::string south_name = "South." + id.str();
		std::string east_name = "East." + id.str();
		std::string west_name = "West." + id.str();

		bool disk_exist = false;
		label patchID[4];
		patchID[0] = mesh.boundaryMesh().findPatchID(north_name);
		patchID[1] = mesh.boundaryMesh().findPatchID(south_name);
		patchID[2] = mesh.boundaryMesh().findPatchID(west_name);
		patchID[3] = mesh.boundaryMesh().findPatchID(east_name);
		if (patchID[0]!=-1 && patchID[1]!=-1 && patchID[2]!=-1 && patchID[3]!=-1)
			disk_exist = true;

		if (disk_exist == true)
		{
			std::string filename = "Disk." + id.str() + ".xml";
			std::ofstream fDisk(filename.c_str(), std::ios::out);
			fDisk.setf(std::ios::scientific);

			fDisk << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			fDisk << "<opensmoke version=\"0.1a\">" << std::endl;

			fDisk << "<NumberOfSpecies>" << std::endl;
			fDisk << ns << std::endl;
			fDisk << "</NumberOfSpecies>" << std::endl;

			fDisk << "<NamesOfSpecies>" << std::endl;
			for(unsigned int i=0;i<ns;i++)
				fDisk << thermodynamicsMapXML->NamesOfSpecies()[i] << std::endl;
			fDisk << "</NamesOfSpecies>" << std::endl;			


			for(int k=0;k<4;k++)
			{
				label patchi = patchID[k];

				std::string side_string = "";
				if (k==0) side_string = "North";
				if (k==1) side_string = "South";
				if (k==2) side_string = "West";
				if (k==3) side_string = "East";

				fDisk << "<Data side=\"" << side_string << "\">" << std::endl;
				
				int number_points = 0;
				forAll(mesh.boundaryMesh()[patchi].faceCentres(), facei)
					number_points++;
				
				fDisk << "<NumberOfPoints>" << std::endl;
				fDisk << number_points << std::endl;
				fDisk << "</NumberOfPoints>" << std::endl;
			
				fDisk << "<Points>" << std::endl;
				forAll(mesh.boundaryMesh()[patchi].faceCentres(), facei)
				{
					fDisk 	<< mesh.boundaryMesh()[patchi].faceCentres()[facei].x() << " " 
						<< mesh.boundaryMesh()[patchi].faceCentres()[facei].y() << " "
						<< mesh.boundaryMesh()[patchi].faceCentres()[facei].z() << " "
						<< std::endl;
				}
				fDisk << "</Points>" << std::endl;

				fDisk << "<MassFractions>" << std::endl;
				forAll(mesh.boundaryMesh()[patchi].faceCentres(), facei)
				{
					for(unsigned int i=1;i<=ns;i++)
						fDisk << Y[i-1].boundaryField()[patchi][facei] << " ";
					fDisk << std::endl;
				}
				fDisk << "</MassFractions>" << std::endl;

				fDisk << "<Temperature>" << std::endl;
				forAll(mesh.boundaryMesh()[patchi].faceCentres(), facei)
				{
					fDisk << T.boundaryField()[patchi][facei] << std::endl;
				}
				fDisk << "</Temperature>" << std::endl;

				fDisk << "</Data>" << std::endl;
			}
			
			fDisk << "</opensmoke>" << std::endl;
			fDisk.close();
		}
	}	
}
