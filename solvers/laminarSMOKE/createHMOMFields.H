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

// Avogadro's Number (6.022e23 #/mol)
dimensionedScalar NAV = dimensionedScalar("NAV", dimensionSet(0., 0., 0., 0., -1.), 6.022e23);

// Surface and volume of primary particles
dimensionedScalar S0 = dimensionedScalar("S0", dimensionSet(0., 2., 0., 0., 0.), 1.);
dimensionedScalar V0 = dimensionedScalar("V0", dimensionSet(0., 3., 0., 0., 0.), 1.);


Info<< "HMOM_Norm_N0 [mol/m3]\n" << endl;
volScalarField HMOM_Norm_N0
(
    IOobject
    (
        "HMOM_Norm_N0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "HMOM_Norm_M00 [mol/m3]\n" << endl;
volScalarField HMOM_Norm_M00
(
    IOobject
    (
        "HMOM_Norm_M00",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "HMOM_Norm_M10 [mol/m3]\n" << endl;
volScalarField HMOM_Norm_M10
(
    IOobject
    (
        "HMOM_Norm_M10",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "HMOM_Norm_M01 [mol/m3]\n" << endl;
volScalarField HMOM_Norm_M01
(
    IOobject
    (
        "HMOM_Norm_M01",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "HMOM_N0\n" << endl;
volScalarField HMOM_N0
(
    IOobject
    (
        "HMOM_N0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_N0", dimensionSet(0., -3., 0., 0., 0.), 0.)
);

Info<< "HMOM_M00\n" << endl;
volScalarField HMOM_M00
(
    IOobject
    (
        "HMOM_M00",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_M00", dimensionSet(0., -3., 0., 0., 0.), 0.)
);

Info<< "HMOM_M10\n" << endl;
volScalarField HMOM_M10
(
    IOobject
    (
        "HMOM_M10",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_M10", dimensionSet(0., 0., 0., 0., 0.), 0.)
);

Info<< "HMOM_M10\n" << endl;
volScalarField HMOM_M01
(
    IOobject
    (
        "HMOM_M01",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_M01", dimensionSet(0., -1., 0., 0., 0.), 0.)
);

Info<< "HMOM_fv\n" << endl;
volScalarField HMOM_fv
(
    IOobject
    (
        "HMOM_fv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_fv", dimensionSet(0., 0., 0., 0., 0.), 0.)
);

Info<< "HMOM_N\n" << endl;
volScalarField HMOM_N
(
    IOobject
    (
        "HMOM_N",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_N", dimensionSet(0., -3., 0., 0., 0.), 0.)
);

Info<< "HMOM_S\n" << endl;
volScalarField HMOM_S
(
    IOobject
    (
        "HMOM_S",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_S", dimensionSet(0., 2., 0., 0., 0.), 0.)
);

Info<< "HMOM_V\n" << endl;
volScalarField HMOM_V
(
    IOobject
    (
        "HMOM_V",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_V", dimensionSet(0., 3., 0., 0., 0.), 0.)
);

