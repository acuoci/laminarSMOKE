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

//- Allocating memory for mean molecular weight [kg/mol]
Info<< "Creating field MWmix\n" << endl;
volScalarField MWmix
(
    IOobject
    (
        "gas_MW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputGasProperties
    ),
    mesh,
    dimensionSet(1, 0, 0, 0, -1)
);

//- Allocating memory for total concentration [kmol/m3]
Info<< "Creating field cTot\n" << endl;
volScalarField cTot
(
    IOobject
    (
        "gas_cTot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputGasProperties
    ),
    mesh,
    dimensionSet(0, -3, 0, 0, 1)
);

//- Allocating memory for dynamic viscosity [kg/m/s]
Info<< "Creating field mu\n" << endl;
volScalarField mu
(
    IOobject
    (
        "gas_mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputGasProperties
    ),
    mesh,
    dimensionSet(1, -1, -1, 0, 0)
);

//- Allocating memory for thermal conductivity [W/m/K]
Info<< "Creating field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "gas_lambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputGasProperties
    ),
    mesh,
    dimensionSet(1, 1, -3, -1, 0)
);

//- Allocating memory for radiative heat transfer term [W/m3]
volScalarField Qrad
(
    IOobject
    (
        "rad_Qrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputGasRadiation
    ),
    mesh,
    dimensionedScalar("rad_Qrad", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

//- Allocating memory for dummy field useful only to tune the radiative heat transfer model
volScalarField dQ
(
	IOobject
	(
	    "dQ",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

//- Allocating memory for reaction heat [W/m3]
volScalarField Q
(
    IOobject
    (
        "gas_Q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gas_Q", dimensionSet(1, -1, -3, 0, 0), 0.0)
); 

//- Allocating memory for mean specific heat [J/kg/K]
Info<< "Creating field Cp\n" << endl;
volScalarField cp
(
    IOobject
    (
        "gas_Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputGasProperties
    ),
    mesh,
    dimensionSet(0, 2, -2, -1, 0)
);

//- Allocating memory for mean specific heat [J/kg/K]
Info<< "Creating field Cv\n" << endl;
volScalarField cv
(
    IOobject
    (
        "gas_Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputGasProperties
    ),
    mesh,
    dimensionSet(0, 2, -2, -1, 0)
);

//- Allocating memory for total enthalpy [m2/s2]
Info<< "Creating field enthalpy\n" << endl;
volScalarField h
(
    IOobject
    (
        "gas_h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputGasProperties
    ),
    mesh,
    dimensionSet(0, 2, -2, 0, 0)
);

Info<< "Creating sumDiffusionCorrections\n" << endl;
volScalarField sumDiffusionCorrections
(
    IOobject
    (
        "gas_sumCorrections",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.)
);

Info<< "Creating massDiffusionInEnergyEquation\n" << endl;
volScalarField massDiffusionInEnergyEquation
(
    IOobject
    (
		"gas_massDiffusionInEnergyEquation",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.)
);
