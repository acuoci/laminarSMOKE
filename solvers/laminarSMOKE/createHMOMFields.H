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

OpenSMOKE::HMOM hmom;

hmom.SetSurfaceDensity(hmom_surface_density);
hmom.SetNucleation(hmom_nucleation);
hmom.SetSurfaceGrowth(hmom_surface_growth);
hmom.SetOxidation(hmom_oxidation);
hmom.SetCondensation(hmom_condensation);
hmom.SetCoagulation(hmom_coagulation);
hmom.SetCoagulationContinous(hmom_coagulation_continous);
hmom.SetThermophoreticModel(hmom_thermophoretic_model);
hmom.SetPAH(hmom_pah_species_names);

hmom.SetNumberCarbonPAH(hmom_pah_nc);
hmom.SetFractalDiameterModel(hmom_fractal_diameter_model);
hmom.SetCollisionDiameterModel(hmom_collision_diameter_model);
hmom.SetSchmidtNumber(hmom_schmidt_number);
hmom.SetStickingCoefficient(hmom_sticking_coefficient);
hmom.SetSootDensity(hmom_soot_density);

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

Info<< "soot_fv_large [-]\n" << endl;
volScalarField soot_fv_large
(
    IOobject
    (
        "soot_fv_large",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("soot_fv_large", dimensionSet(0., 0., 0., 0., 0.), 0.)
);

soot_fv_large = HMOM_fv;
