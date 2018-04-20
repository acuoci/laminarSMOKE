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
	const std::string alignement = "yx";
	const bool invert_north_south = false;
	const bool invert_east_west = true;

	label ns = Y.size();

	OpenSMOKE::OpenSMOKEVectorDouble y(ns);
	OpenSMOKE::OpenSMOKEVectorDouble x(ns);

	const scalarField& TCells = T.internalField();

	std::ofstream fDiskTopology("DiskTopology.xml", std::ios::out);
	fDiskTopology.setf(std::ios::scientific);

	fDiskTopology << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
	fDiskTopology << "<opensmoke version=\"0.1a\">" << std::endl;

	unsigned int ndisks = 0;
	for (unsigned int i=1;i<100;i++)
	{
		std::stringstream id;	id << i;

		std::string north_name_u = "North." + id.str();
		std::string south_name_u = "South." + id.str();
		std::string north_name_l = "north." + id.str();
		std::string south_name_l = "south." + id.str();

		std::string east_name_u = "East." + id.str();
		std::string west_name_u = "West." + id.str();
		std::string east_name_l = "east." + id.str();
		std::string west_name_l = "west." + id.str();

		if (invert_north_south == true)
		{
			north_name_u = "South." + id.str();
			south_name_u = "North." + id.str();
			north_name_l = "south." + id.str();
			south_name_l = "north." + id.str();
		}

		if (invert_east_west == true)
		{
			east_name_u = "West." + id.str();
			west_name_u = "East." + id.str();
			east_name_l = "west." + id.str();
			west_name_l = "east." + id.str();
		}

		bool disk_exist = false;
		label patchID[4];
		
		patchID[0] = mesh.boundaryMesh().findPatchID(north_name_u);
		if (patchID[0]==-1)
			patchID[0] = mesh.boundaryMesh().findPatchID(north_name_l);

		patchID[1] = mesh.boundaryMesh().findPatchID(south_name_u);
		if (patchID[1]==-1)
			patchID[1] = mesh.boundaryMesh().findPatchID(south_name_l);

		patchID[2] = mesh.boundaryMesh().findPatchID(west_name_u);
		if (patchID[2]==-1)
			patchID[2] = mesh.boundaryMesh().findPatchID(west_name_l);

		patchID[3] = mesh.boundaryMesh().findPatchID(east_name_u);
		if (patchID[3]==-1)
			patchID[3] = mesh.boundaryMesh().findPatchID(east_name_l);

		if (patchID[0]!=-1 && patchID[1]!=-1 && patchID[3]!=-1)
			disk_exist = true;

		if (disk_exist == true)
		{
			ndisks++;
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

			unsigned int ncells = 0;
			double volume = 0.;
			double x_west = 0.;
			double x_east = 0.;
			double y_north = 0.;
			double y_south = 0.;
			fDiskTopology << "<Disk." << i << ">" << std::endl;
			fDiskTopology << "<Cells>" << std::endl;

			for(int k=0;k<4;k++)
			{
				label patchi = patchID[k];

				std::string side_string = "";
				if (k==0) side_string = "North";
				if (k==1) side_string = "South";
				if (k==2) side_string = "West";
				if (k==3) side_string = "East";

				if (patchi != -1)
				{
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
						if (alignement == "yx")
						fDisk 	<< mesh.boundaryMesh()[patchi].faceCentres()[facei].x() << " " 
							<< mesh.boundaryMesh()[patchi].faceCentres()[facei].y() << " "
							<< 5.e-4 << " "
							<< std::endl;

						if (alignement == "xy")
						fDisk 	<< mesh.boundaryMesh()[patchi].faceCentres()[facei].y() << " " 
							<< mesh.boundaryMesh()[patchi].faceCentres()[facei].x() << " "
							<< 5.e-4 << " "
							<< std::endl;

						if (k==0 && alignement == "yx")	y_north = mesh.boundaryMesh()[patchi].faceCentres()[facei].x();
						if (k==0 && alignement == "xy")	y_north = mesh.boundaryMesh()[patchi].faceCentres()[facei].y();
						if (k==1 && alignement == "yx")	y_south = mesh.boundaryMesh()[patchi].faceCentres()[facei].x();
						if (k==1 && alignement == "xy")	y_south = mesh.boundaryMesh()[patchi].faceCentres()[facei].y();

						if (k==2 && alignement == "yx")	x_west = mesh.boundaryMesh()[patchi].faceCentres()[facei].y();
						if (k==2 && alignement == "xy")	x_west = mesh.boundaryMesh()[patchi].faceCentres()[facei].x();
						if (k==3 && alignement == "yx")	x_east = mesh.boundaryMesh()[patchi].faceCentres()[facei].y();
						if (k==3 && alignement == "xy")	x_east = mesh.boundaryMesh()[patchi].faceCentres()[facei].x();

					}
					fDisk << "</Points>" << std::endl;

					fDisk << "<Areas>" << std::endl;
					label startFace = mesh.boundaryMesh()[patchi].start();
					label nFaces    = mesh.boundaryMesh()[patchi].size();
					for (label facei = startFace; facei < startFace + nFaces; facei++)
						fDisk 	<< mag(mesh.Sf()[facei]) << std::endl;
					fDisk << "</Areas>" << std::endl;

					fDisk << "<Cells>" << std::endl;
					forAll(mesh.boundaryMesh()[patchi], facei)
					{
						const label index = mesh.boundaryMesh()[patchi].faceCells()[facei];
						fDisk 	<< index << " " << mesh.V()[index] << std::endl;
					
						ncells++;
						volume += mesh.V()[index];
						fDiskTopology 	<< index << " " << mesh.V()[index] << std::endl;
					}
					fDisk << "</Cells>" << std::endl;

					fDisk << "<MassFractions>" << std::endl;
					forAll(mesh.boundaryMesh()[patchi].faceCentres(), facei)
					{
						double sum = 0.;
						for(unsigned int i=1;i<=ns;i++)
							sum += Y[i-1].boundaryField()[patchi][facei];
						for(unsigned int i=1;i<=ns;i++)
							fDisk << Y[i-1].boundaryField()[patchi][facei]/sum << " ";
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
				else
				{
					fDisk << "<Data side=\"" << side_string << "\">" << std::endl;				

					int number_points = 0;			
					fDisk << "<NumberOfPoints>" << std::endl;
					fDisk << number_points << std::endl;
					fDisk << "</NumberOfPoints>" << std::endl;
					fDisk << "</Data>" << std::endl;
				}
			}

			fDiskTopology << "</Cells>" << std::endl;
			fDiskTopology << "<NumberOfPoints>" << std::endl;
			fDiskTopology << ncells << std::endl;
			fDiskTopology << "</NumberOfPoints>" << std::endl;
			fDiskTopology << "<Volume>" << std::endl;
			fDiskTopology << volume << std::endl;
			fDiskTopology << "</Volume>" << std::endl;
			fDiskTopology << "</Disk." << i << ">" << std::endl;
			
			fDisk << "</opensmoke>" << std::endl;
			fDisk.close();

			if (x_west < 0 || x_east < 0|| y_south < 0|| y_north < 0)
			{
				Info << "Error in exportDisks: the mesh must be entirely in the positive xy space" << endl;
				abort();
			}

			if (x_west > x_east)
			{
				Info << "West: " << x_west << " East: " << x_east << endl;
				Info << "Error in exportDisks: the east faces must be always the external faces" << endl;
				Info << "Please switch east/west faces" << endl;
				abort();
			}
			if (y_south > y_north)
			{
				Info << "South: " << y_south << " North: " << y_north << endl;
				Info << "Error in exportDisks: the north faces must be always the top faces" << endl;
				Info << "Please switch top/south faces" << endl;
				abort();
			}

			if ((y_north-y_south) > x_east)
			{
				Info << "Thickness: " << y_north-y_south << " Radius: " << x_east << endl;
				Info << "Error in exportDisks: the disk thickness is larger than the external radius" << endl;
				Info << "Please switch x/y directions" << endl;
				abort();
			}
		}
	}	

	fDiskTopology << "<Disks>" << std::endl;
	fDiskTopology << ndisks << std::endl;
	fDiskTopology << "</Disks>" << std::endl;
	fDiskTopology << "</opensmoke>" << std::endl;
	fDiskTopology.close();
}
