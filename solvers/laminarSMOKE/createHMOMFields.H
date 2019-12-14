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
hmom.SetSurfaceDensityCorrectionCoefficient(hmom_surface_density_correction_coefficient);
hmom.SetSurfaceDensityCorrectionCoefficientA1(hmom_surface_density_coeff_a1);
hmom.SetSurfaceDensityCorrectionCoefficientA2(hmom_surface_density_coeff_a2);
hmom.SetSurfaceDensityCorrectionCoefficientB1(hmom_surface_density_coeff_b1);
hmom.SetSurfaceDensityCorrectionCoefficientB2(hmom_surface_density_coeff_b2);

hmom.SetA1f(hmom_A1f);
hmom.SetA1b(hmom_A1b);
hmom.SetA2f(hmom_A2f);
hmom.SetA2b(hmom_A2b);
hmom.SetA3f(hmom_A3f);
hmom.SetA3b(hmom_A3b);
hmom.SetA4(hmom_A4);
hmom.SetA5(hmom_A5);

hmom.SetE1f(hmom_E1f);
hmom.SetE1b(hmom_E1b);
hmom.SetE2f(hmom_E2f);
hmom.SetE2b(hmom_E2b);
hmom.SetE3f(hmom_E3f);
hmom.SetE3b(hmom_E3b);
hmom.SetE4(hmom_E4);
hmom.SetE5(hmom_E5);

hmom.Setn1f(hmom_n1f);
hmom.Setn1b(hmom_n1b);
hmom.Setn2f(hmom_n2f);
hmom.Setn2b(hmom_n2b);
hmom.Setn3f(hmom_n3f);
hmom.Setn3b(hmom_n3b);
hmom.Setn4(hmom_n4);
hmom.Setn5(hmom_n5);

hmom.SetEfficiency6(hmom_eff6);

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

hmom.Summary();


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
