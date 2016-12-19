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
const double AvogadroNumber = 6.022e23;
dimensionedScalar NAV = dimensionedScalar("NAV", dimensionSet(0., 0., 0., 0., -1.), AvogadroNumber);

// HMOM Constants
double pi = 3.141592653589793;
double hmom_one  = 1.;
double hmom_zero = 0.;
double hmom_minus_one = -1.;
double hmom_minus_two = -2.;
double hmom_three = 3.;
double hmom_Avcoll = 1.-2./1.8;
double hmom_Ascoll = 3./1.8-1.;
double hmom_K_coll = 6./std::pow(36.*pi, 1./1.8);
double hmom_K_sphe = std::pow(36.*pi, 1./3.);
double WC = 12e-3;

// Init soot module
if (hmom_analysis == true)
{
	int zero = 0;
	int one = 1;

	__mod_soot_init_MOD_init_soot_hmom (	&hmom_pah_nc, &hmom_bivariate, &hmom_number_moments, &hmom_volume_surface, &hmom_collision_model,
						&one,  &one, &one, &one, &one, &one	);
}

const double VPAH=(WC/hmom_soot_density/AvogadroNumber) * double(hmom_pah_nc);
const double VDIM=2.*VPAH;
const double SDIM=hmom_K_sphe*Foam::pow(VDIM, 2./3.);
const double VNUCL=2.*VDIM; 
const double SNUCL=hmom_K_sphe*Foam::pow(VNUCL, 2./3.);  

// Surface and volume of primary particles
dimensionedScalar S0 = dimensionedScalar("S0", dimensionSet(0., 2., 0., 0., 0.), SNUCL);
dimensionedScalar V0 = dimensionedScalar("V0", dimensionSet(0., 3., 0., 0., 0.), VNUCL);


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


Info<< "HMOM_N0 [#/m3]\n" << endl;
volScalarField HMOM_N0
(
    IOobject
    (
        "HMOM_N0",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_N0", dimensionSet(0., -3., 0., 0., 0.), 0.)
);

Info<< "HMOM_M00 [1/m3]\n" << endl;
volScalarField HMOM_M00
(
    IOobject
    (
        "HMOM_M00",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_M00", dimensionSet(0., -3., 0., 0., 0.), 0.)
);

Info<< "HMOM_M10 [-]\n" << endl;
volScalarField HMOM_M10
(
    IOobject
    (
        "HMOM_M10",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_M10", dimensionSet(0., 0., 0., 0., 0.), 0.)
);

Info<< "HMOM_M01 [1/m]" << endl;
volScalarField HMOM_M01
(
    IOobject
    (
        "HMOM_M01",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_M01", dimensionSet(0., -1., 0., 0., 0.), 0.)
);

Info<< "HMOM_fv [-]\n" << endl;
volScalarField HMOM_fv
(
    IOobject
    (
        "HMOM_fv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_fv", dimensionSet(0., 0., 0., 0., 0.), 0.)
);

Info<< "HMOM_N [#/m3]\n" << endl;
volScalarField HMOM_N
(
    IOobject
    (
        "HMOM_N",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_N", dimensionSet(0., -3., 0., 0., 0.), 0.)
);

Info<< "HMOM_dp [m]\n" << endl;
volScalarField HMOM_dp
(
    IOobject
    (
        "HMOM_dp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_dp", dimensionSet(0., 1., 0., 0., 0.), 0.)
);

Info<< "HMOM_dc [m]\n" << endl;
volScalarField HMOM_dc
(
    IOobject
    (
        "HMOM_dc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_dc", dimensionSet(0., 1., 0., 0., 0.), 0.)
);

Info<< "HMOM_np [-]\n" << endl;
volScalarField HMOM_np
(
    IOobject
    (
        "HMOM_np",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_np", dimensionSet(0., 0., 0., 0., 0.), 0.)
);

Info<< "HMOM_S\n" << endl;
volScalarField HMOM_S
(
    IOobject
    (
        "HMOM_S",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("HMOM_V", dimensionSet(0., 3., 0., 0., 0.), 0.)
);

