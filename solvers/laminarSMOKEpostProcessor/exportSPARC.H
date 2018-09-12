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

// Export SPARC
if (exportSPARC == true)
{
	Info<< "Reading field zMix..." << endl;
	volScalarField zMix
	(
	    IOobject
	    (
		"zMix",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	    mesh
	);


	Info<< "Reading field clusterIndex..." << endl;
	volScalarField clusterIndex
	(
	    IOobject
	    (
		"clusterIndex",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	    mesh
	);

        Info<< "Reading field clusterSpecies..." << endl;
	volScalarField clusterSpecies
	(
	    IOobject
	    (
		"clusterSpecies",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	    mesh
	);

        Info<< "Reading field clusterReactions..." << endl;
	volScalarField clusterReactions
	(
	    IOobject
	    (
		"clusterReactions",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	    mesh
	);

	std::string filename = "output." + runTime.timeName();
	std::ofstream fOutput(filename.c_str(), std::ios::out);
	fOutput.setf(std::ios::scientific);
	fOutput << "Data from OpenFOAM simulation" << std::endl;

	const scalarField& TCells = T.internalField();
	const scalarField& zMixCells = zMix.internalField();
	const scalarField& clusterIndexCells = clusterIndex.internalField();
	const scalarField& clusterSpeciesCells = clusterSpecies.internalField();
	const scalarField& clusterReactionsCells = clusterReactions.internalField();

	const double T_threshold = 200.;
	
	forAll(TCells, celli)
	{	
		if (TCells[celli] > T_threshold)
		{
			fOutput << std::setw(15) << std::scientific << std::setprecision(6) << runTime.timeName();

			fOutput << std::setw(15) << std::scientific << std::setprecision(6) << mesh.C()[celli][0];
			fOutput << std::setw(15) << std::scientific << std::setprecision(6) << mesh.C()[celli][1];
			fOutput << std::setw(15) << std::scientific << std::setprecision(6) << mesh.C()[celli][2];

			fOutput << std::setw(15) << std::fixed << clusterIndexCells[celli];
			fOutput << std::setw(15) << std::fixed << clusterSpeciesCells[celli];
			fOutput << std::setw(15) << std::fixed << clusterReactionsCells[celli];

			fOutput << std::setw(15) << std::scientific << std::setprecision(6) << zMixCells[celli];


			fOutput << std::setw(15) << std::scientific << std::setprecision(6) << TCells[celli];
			for(unsigned int i=0;i<Y.size();i++)
				fOutput << std::setw(15) << std::scientific << std::setprecision(6) << Y[i].internalField()[celli];

			fOutput << std::endl;
		}
	}

	fOutput.close();		
}
