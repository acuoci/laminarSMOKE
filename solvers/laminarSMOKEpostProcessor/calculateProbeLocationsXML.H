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

// XML files in specific locations
if (pnts_xml.size() > 0)
{
	Info << "Writing XML files in singl points..." << endl;

	OpenSMOKE::OpenSMOKEVectorDouble y(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML->NumberOfSpecies());

	for (unsigned int i=0;i<pnts_xml.size();i++)
	{
		const double xx = pnts_xml[i][0];
		const double yy = pnts_xml[i][1];
		const double zz = pnts_xml[i][2];

		Info << " * point " << i+1 << " @ " << xx << " " << yy << " " << zz << endl;

		vector pos(xx,yy,zz);
		label cellI = mesh.findCell(pos);

		if (cellI > -1)
		{
			// Temperature and pressure
			autoPtr<interpolation<scalar> > Tinterp = interpolation<scalar>::New("cellPoint", T);
			scalar Tint = Tinterp->interpolate(pos, cellI);
			autoPtr<interpolation<scalar> > pinterp = interpolation<scalar>::New("cellPoint", p);
			scalar pint = pinterp->interpolate(pos, cellI);

			Info << " * cell: " << cellI << " T: " << Tint << " P: " << pint << endl; 

			// Extract the mass fractions
			for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			{
				autoPtr<interpolation<scalar> > Yinterp = interpolation<scalar>::New("cellPoint", Y[i]);
				y[i+1] = Yinterp->interpolate(pos, cellI);
			}

			const double sum = y.SumElements();
			for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				y[i+1] /= sum;

			// Molecular weight
			double mw;
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

			// Density
			const double rho = pint*mw/PhysicalConstants::R_J_kmol/Tint;


			// Write XML file
			{
				std::stringstream tag;
				tag << i;
				std::string name_file = "/Output_point_" + tag.str() + ".xml";

				// File name
				fileName outputFile(runTime.path()/runTime.timeName()/name_file.c_str());
				std::ofstream fXML(outputFile.c_str(), std::ios::out);
				fXML.setf(std::ios::scientific);

				unsigned int n_additional = 5;

				fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
				fXML << "<opensmoke version=\"0.1a\">" << std::endl;
		
				fXML << "<Type> SinglePoint </Type>" << std::endl;
		
				fXML << "<additional>" << std::endl;
				fXML << n_additional << std::endl;
				fXML << "dummy [-] 2" << std::endl;
				fXML << "temperature [K] 3" << std::endl;
				fXML << "pressure [Pa] 4" << std::endl;
				fXML << "mol-weight [kg/kmol] 5" << std::endl;
				fXML << "density [kg/m3] 6" << std::endl;
				fXML << "</additional>" << std::endl;
		
				fXML << "<t-p-mw>" << std::endl;
				fXML << "1 2 3" << std::endl;
				fXML << "</t-p-mw>" << std::endl;
		
				fXML << "<mass-fractions>" << std::endl;
				fXML << thermodynamicsMapXML->NumberOfSpecies() << std::endl;
				for(unsigned int i=1;i<=thermodynamicsMapXML->NumberOfSpecies();i++)
					fXML << thermodynamicsMapXML->NamesOfSpecies()[i-1] << " " << thermodynamicsMapXML->MW()[i] << " " << n_additional+i << std::endl;
				fXML << "</mass-fractions>" << std::endl;

				fXML << "<profiles>" << std::endl;
				for(unsigned int i=1;i<=2;i++)
				{
					fXML << i << " ";
					fXML << Tint << " ";
					fXML << pint << " ";
					fXML << mw << " ";
					fXML << rho << " ";

					for(unsigned int j=1;j<=thermodynamicsMapXML->NumberOfSpecies();j++)
						fXML << y[j] << " ";
					fXML << std::endl;
				}
				fXML << "</profiles>" << std::endl;

				fXML << "<profiles-size>" << std::endl;
				fXML << 2 << " " << thermodynamicsMapXML->NumberOfSpecies()+n_additional << std::endl;
				fXML << "</profiles-size>" << std::endl;
				fXML << "</opensmoke>" << std::endl;
				fXML.close();
			}
		}
	}
}
