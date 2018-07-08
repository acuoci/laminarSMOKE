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
|   Authors: A. Cuoci, G. D'Alessio, A. Parente                           |
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
|   Copyright(C) 2018 A. Cuoci                                            |
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

#if SPARC==1

Info<< "Creating Cluster fields" << endl;
volScalarField clusterIndex
(
	IOobject
	(
		"clusterIndex",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("dummy", dimensionSet(0, 0, 0, 0, 0), 0.),
	#if OPENFOAM_VERSION >= 40
	extrapolatedCalculatedFvPatchScalarField::typeName
	#else
	zeroGradientFvPatchScalarField::typeName
	#endif
);
volScalarField clusterWeight
(
	IOobject
	(
		"clusterWeight",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("dummy", dimensionSet(0, 0, 0, 0, 0), 0.),
	#if OPENFOAM_VERSION >= 40
	extrapolatedCalculatedFvPatchScalarField::typeName
	#else
	zeroGradientFvPatchScalarField::typeName
	#endif
);
volScalarField clusterSpecies
(
	IOobject
	(
		"clusterSpecies",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("dummy", dimensionSet(0, 0, 0, 0, 0), 0.),
	#if OPENFOAM_VERSION >= 40
	extrapolatedCalculatedFvPatchScalarField::typeName
	#else
	zeroGradientFvPatchScalarField::typeName
	#endif
);
volScalarField clusterReactions
(
        IOobject
        (
                "clusterReactions",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dummy", dimensionSet(0, 0, 0, 0, 0), 0.),
        #if OPENFOAM_VERSION >= 40
        extrapolatedCalculatedFvPatchScalarField::typeName
        #else
        zeroGradientFvPatchScalarField::typeName
        #endif
);

#endif

