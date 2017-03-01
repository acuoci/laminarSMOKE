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

// Post processing: soot (HMOM)
if (postProcessingHMOM == true)
{
	Switch hmom_analysis = false;
	Switch hmom_pah_gas_consumption = true;
	int hmom_index_H = -1;
	int hmom_index_OH = -1;
	int hmom_index_H2O = -1;
	int hmom_index_H2 = -1;
	int hmom_index_O2 = -1;
	int hmom_index_C2H2 = -1;
	std::vector<int> hmom_pah_species_indices;
	std::vector<std::string> hmom_pah_species_names;

	int hmom_pah_nc = 16;
	int hmom_number_moments = 4;
	int hmom_collision_diameter_model = 2;
	int hmom_fractal_diameter_model = 1;

	int hmom_nucleation = 1;
	int hmom_surface_growth = 1;
	int hmom_oxidation = 1;
	int hmom_condensation = 1;
	int hmom_coagulation = 1;
	int hmom_coagulation_continous = 1;

	scalar hmom_minimum_normalized_value = 1.e-16;
	scalar hmom_minimum_fv_for_integration = 1.e-16;


	const dictionary& hmomDictionary = solverOptionsDictionary.subDict("HMOM");
	{
		hmom_analysis = Switch(hmomDictionary.lookup(word("hmom")));

		if ( hmom_analysis == true)
		{
			hmom_index_H = thermodynamicsMapXML->IndexOfSpecies("H")-1;
			hmom_index_OH = thermodynamicsMapXML->IndexOfSpecies("OH")-1;
			hmom_index_H2O = thermodynamicsMapXML->IndexOfSpecies("H2O")-1;
			hmom_index_H2 = thermodynamicsMapXML->IndexOfSpecies("H2")-1;
			hmom_index_O2 = thermodynamicsMapXML->IndexOfSpecies("O2")-1;
			hmom_index_C2H2 = thermodynamicsMapXML->IndexOfSpecies("C2H2")-1;
		
			List<word>  listSpecies(hmomDictionary.lookup("pah_species"));
			hmom_pah_species_indices.resize(listSpecies.size());
			hmom_pah_species_names.resize(listSpecies.size());
			for (int i=0;i<listSpecies.size();i++)
			{
				hmom_pah_species_indices[i] = thermodynamicsMapXML->IndexOfSpecies(listSpecies[i])-1;
				hmom_pah_species_names[i] = listSpecies[i];
			}

			hmom_pah_nc = readLabel(hmomDictionary.lookup("pah_nc"));
			hmom_collision_diameter_model = readLabel(hmomDictionary.lookup("collision_diameter_model"));
			hmom_fractal_diameter_model = readLabel(hmomDictionary.lookup("fractal_diameter_model"));

			hmom_nucleation = readLabel(hmomDictionary.lookup("nucleation"));
			hmom_surface_growth = readLabel(hmomDictionary.lookup("surface_growth"));
			hmom_oxidation = readLabel(hmomDictionary.lookup("oxidation"));
			hmom_condensation = readLabel(hmomDictionary.lookup("condensation"));
			hmom_coagulation = readLabel(hmomDictionary.lookup("coagulation"));
			hmom_coagulation_continous = readLabel(hmomDictionary.lookup("coagulation_continous"));

			hmom_pah_gas_consumption = Switch(hmomDictionary.lookup(word("pah_gas_consumption")));
		}
	}

	OpenSMOKE::HMOM hmom;

	hmom.SetNucleation(hmom_nucleation);
	hmom.SetSurfaceGrowth(hmom_surface_growth);
	hmom.SetOxidation(hmom_oxidation);
	hmom.SetCondensation(hmom_condensation);
	hmom.SetCoagulation(hmom_coagulation);
	hmom.SetCoagulationContinous(hmom_coagulation_continous);
	hmom.SetPAH(hmom_pah_species_names);

	hmom.SetNumberCarbonPAH(hmom_pah_nc);
	hmom.SetFractalDiameterModel(hmom_fractal_diameter_model);
	hmom.SetCollisionDiameterModel(hmom_collision_diameter_model);

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

	Info<< "HMOM_Source_M00 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_M00
	(
	    IOobject
	    (
		"HMOM_Source_M00",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_M00", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_M10 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_M10
	(
	    IOobject
	    (
		"HMOM_Source_M10",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_M10", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_M01 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_M01
	(
	    IOobject
	    (
		"HMOM_Source_M01",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_M01", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_N0 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_N0
	(
	    IOobject
	    (
		"HMOM_Source_N0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_N0", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Nucleation_M00 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Nucleation_M00
	(
	    IOobject
	    (
		"HMOM_Source_Nucleation_M00",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Nucleation_M00", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Nucleation_M10 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Nucleation_M10
	(
	    IOobject
	    (
		"HMOM_Source_Nucleation_M10",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Nucleation_M10", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Nucleation_M01 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Nucleation_M01
	(
	    IOobject
	    (
		"HMOM_Source_Nucleation_M01",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Nucleation_M01", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Nucleation_N0 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Nucleation_N0
	(
	    IOobject
	    (
		"HMOM_Source_Nucleation_N0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Nucleation_N0", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Growth_M00 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Growth_M00
	(
	    IOobject
	    (
		"HMOM_Source_Growth_M00",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Growth_M00", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Growth_M10 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Growth_M10
	(
	    IOobject
	    (
		"HMOM_Source_Growth_M10",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Growth_M10", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Growth_M01 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Growth_M01
	(
	    IOobject
	    (
		"HMOM_Source_Growth_M01",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Growth_M01", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Growth_N0 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Growth_N0
	(
	    IOobject
	    (
		"HMOM_Source_Growth_N0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Growth_N0", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Oxidation_M00 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Oxidation_M00
	(
	    IOobject
	    (
		"HMOM_Source_Oxidation_M00",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Oxidation_M00", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Oxidation_M10 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Oxidation_M10
	(
	    IOobject
	    (
		"HMOM_Source_Oxidation_M10",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Oxidation_M10", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Oxidation_M01 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Oxidation_M01
	(
	    IOobject
	    (
		"HMOM_Source_Oxidation_M01",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Oxidation_M01", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Oxidation_N0 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Oxidation_N0
	(
	    IOobject
	    (
		"HMOM_Source_Oxidation_N0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Oxidation_N0", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Condensation_M00 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Condensation_M00
	(
	    IOobject
	    (
		"HMOM_Source_Condensation_M00",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Condensation_M00", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Condensation_M10 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Condensation_M10
	(
	    IOobject
	    (
		"HMOM_Source_Condensation_M10",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Condensation_M10", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Condensation_M01 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Condensation_M01
	(
	    IOobject
	    (
		"HMOM_Source_Condensation_M01",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Condensation_M01", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Condensation_N0 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Condensation_N0
	(
	    IOobject
	    (
		"HMOM_Source_Condensation_N0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Condensation_N0", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Coagulation_M00 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Coagulation_M00
	(
	    IOobject
	    (
		"HMOM_Source_Coagulation_M00",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Coagulation_M00", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Coagulation_M10 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Coagulation_M10
	(
	    IOobject
	    (
		"HMOM_Source_Coagulation_M10",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Coagulation_M10", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Coagulation_M01 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Coagulation_M01
	(
	    IOobject
	    (
		"HMOM_Source_Coagulation_M01",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Coagulation_M01", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info<< "HMOM_Source_Coagulation_N0 [mol/m3/s]\n" << endl;
	volScalarField HMOM_Source_Coagulation_N0
	(
	    IOobject
	    (
		"HMOM_Source_Coagulation_N0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("HMOM_Source_Coagulation_N0", dimensionSet(0., -3., -1., 0., 1.), 0.)
	);

	Info << "Post processing soot (hmom)..." << endl;

	label ns = Y.size();

	OpenSMOKE::OpenSMOKEVectorDouble y(ns);
	OpenSMOKE::OpenSMOKEVectorDouble x(ns);
	OpenSMOKE::OpenSMOKEVectorDouble c(thermodynamicsMapXML->NumberOfSpecies());
	Eigen::VectorXd y_eigen(ns);
	Eigen::VectorXd x_eigen(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();

	const scalarField& HMOM_Norm_N0Cells  = HMOM_Norm_N0.internalField();
	const scalarField& HMOM_Norm_M00Cells = HMOM_Norm_M00.internalField();
	const scalarField& HMOM_Norm_M01Cells = HMOM_Norm_M01.internalField();
	const scalarField& HMOM_Norm_M10Cells = HMOM_Norm_M10.internalField();

	scalarField& HMOM_fvCells = HMOM_fv.internalField(); 
	scalarField& HMOM_NCells  = HMOM_N.internalField(); 
	scalarField& HMOM_dpCells = HMOM_dp.internalField(); 
	scalarField& HMOM_dcCells = HMOM_dc.internalField(); 
	scalarField& HMOM_npCells = HMOM_np.internalField(); 

	scalarField& HMOM_Source_M00Cells = HMOM_Source_M00.internalField(); 
	scalarField& HMOM_Source_M10Cells = HMOM_Source_M10.internalField(); 
	scalarField& HMOM_Source_M01Cells = HMOM_Source_M01.internalField(); 
	scalarField& HMOM_Source_N0Cells  = HMOM_Source_N0.internalField(); 
	
	scalarField& HMOM_Source_Nucleation_M00Cells = HMOM_Source_Nucleation_M00.internalField(); 
	scalarField& HMOM_Source_Nucleation_M10Cells = HMOM_Source_Nucleation_M10.internalField(); 
	scalarField& HMOM_Source_Nucleation_M01Cells = HMOM_Source_Nucleation_M01.internalField(); 
	scalarField& HMOM_Source_Nucleation_N0Cells  = HMOM_Source_Nucleation_N0.internalField(); 

	scalarField& HMOM_Source_Growth_M00Cells = HMOM_Source_Growth_M00.internalField(); 
	scalarField& HMOM_Source_Growth_M10Cells = HMOM_Source_Growth_M10.internalField(); 
	scalarField& HMOM_Source_Growth_M01Cells = HMOM_Source_Growth_M01.internalField(); 
	scalarField& HMOM_Source_Growth_N0Cells  = HMOM_Source_Growth_N0.internalField(); 

	scalarField& HMOM_Source_Oxidation_M00Cells = HMOM_Source_Oxidation_M00.internalField(); 
	scalarField& HMOM_Source_Oxidation_M10Cells = HMOM_Source_Oxidation_M10.internalField(); 
	scalarField& HMOM_Source_Oxidation_M01Cells = HMOM_Source_Oxidation_M01.internalField(); 
	scalarField& HMOM_Source_Oxidation_N0Cells  = HMOM_Source_Oxidation_N0.internalField();

	scalarField& HMOM_Source_Condensation_M00Cells = HMOM_Source_Condensation_M00.internalField(); 
	scalarField& HMOM_Source_Condensation_M10Cells = HMOM_Source_Condensation_M10.internalField(); 
	scalarField& HMOM_Source_Condensation_M01Cells = HMOM_Source_Condensation_M01.internalField(); 
	scalarField& HMOM_Source_Condensation_N0Cells  = HMOM_Source_Condensation_N0.internalField();

	scalarField& HMOM_Source_Coagulation_M00Cells = HMOM_Source_Coagulation_M00.internalField(); 
	scalarField& HMOM_Source_Coagulation_M10Cells = HMOM_Source_Coagulation_M10.internalField(); 
	scalarField& HMOM_Source_Coagulation_M01Cells = HMOM_Source_Coagulation_M01.internalField(); 
	scalarField& HMOM_Source_Coagulation_N0Cells  = HMOM_Source_Coagulation_N0.internalField();

	Info << " * internal fields..." << endl;
	forAll(TCells, celli)
	{	
		double mw;

		// Extract the mean mass fractions
		for(unsigned int i=0;i<ns;i++)
			y[i+1] = Y[i].internalField()[celli];
		const double sum = y.SumElements();
		for(unsigned int i=0;i<ns;i++)
			y[i+1] /= sum;
		

		// Molecular weight
		thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

		// Eigen vectors
		for(unsigned int i=0;i<ns;i++)
		{
			y_eigen(i) = y[i+1];
			x_eigen(i) = x[i+1];
		}

		// Density
		const double rhoGas = pCells[celli]*mw/PhysicalConstants::R_J_kmol/TCells[celli];

		// Concentrations
		const double cTot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
		Product(cTot,x, &c);

		// HMOM
		if (hmom_analysis == true)
		{
			// Viscosity
			const double viscosity = 1.716e-5*Foam::pow(TCells[celli]/273.15,1.5)*(273.15+110.4)/(TCells[celli]+110.4);		// [kg/m/s]

			// PAH concentration
			double conc_PAH = 0.;
			for (unsigned int j = 0; j<hmom.pah_species().size(); j++)
				conc_PAH += cTot*x[hmom_pah_species_indices[j] + 1];

			// Setup
			hmom.SetNormalizedMoments(HMOM_Norm_M00Cells[celli], HMOM_Norm_M10Cells[celli], HMOM_Norm_M01Cells[celli], HMOM_Norm_N0Cells[celli]);
			hmom.SetTemperatureAndPressure(TCells[celli], pCells[celli]);
			hmom.SetMassFractions(y[hmom_index_OH+1], y[hmom_index_H+1]);
			hmom.SetConcentrations("kmol/m3", c[hmom_index_OH+1], c[hmom_index_H+1], c[hmom_index_H2O+1], c[hmom_index_H2+1], c[hmom_index_C2H2+1], c[hmom_index_O2+1], conc_PAH);
			hmom.SetViscosity(viscosity);

			// Calculates
			hmom.CalculateSourceMoments();

			// Calculating soot volume fraction [-]
			HMOM_fvCells[celli] = hmom.SootVolumeFraction();

			// Calculating soot particle number density [#/m3]
			HMOM_NCells[celli] = hmom.SootParticleNumberDensity();

			// Calculating particle diameters [nm]
			HMOM_dpCells[celli] = hmom.SootParticleDiameter() * 1.e9;

			// Calculating collision diameters [nm]
			HMOM_dcCells[celli] = hmom.SootCollisionParticleDiameter() * 1.e9;

			// Calculating number of primary particles [-]
			HMOM_npCells[celli] = hmom.SootNumberOfPrimaryParticles();

			// Overall source terms [mol/m3/s]
			HMOM_Source_M00Cells[celli] = hmom.sources()[0];
			HMOM_Source_M10Cells[celli] = hmom.sources()[1];
			HMOM_Source_M01Cells[celli] = hmom.sources()[2];
			HMOM_Source_N0Cells[celli]  = hmom.sources()[3];

			// Nucleation source terms [mol/m3/s]
			HMOM_Source_Nucleation_M00Cells[celli] = hmom.sources_nucleation()[0];
			HMOM_Source_Nucleation_M10Cells[celli] = hmom.sources_nucleation()[1];
			HMOM_Source_Nucleation_M01Cells[celli] = hmom.sources_nucleation()[2];
			HMOM_Source_Nucleation_N0Cells[celli]  = hmom.sources_nucleation()[3];

			// Growth source terms [mol/m3/s]
			HMOM_Source_Growth_M00Cells[celli] = hmom.sources_growth()[0];
			HMOM_Source_Growth_M10Cells[celli] = hmom.sources_growth()[1];
			HMOM_Source_Growth_M01Cells[celli] = hmom.sources_growth()[2];
			HMOM_Source_Growth_N0Cells[celli]  = hmom.sources_growth()[3];

			// Oxidation source terms [mol/m3/s]
			HMOM_Source_Oxidation_M00Cells[celli] = hmom.sources_oxidation()[0];
			HMOM_Source_Oxidation_M10Cells[celli] = hmom.sources_oxidation()[1];
			HMOM_Source_Oxidation_M01Cells[celli] = hmom.sources_oxidation()[2];
			HMOM_Source_Oxidation_N0Cells[celli]  = hmom.sources_oxidation()[3];

			// Condensation source terms [mol/m3/s]
			HMOM_Source_Condensation_M00Cells[celli] = hmom.sources_condensation()[0];
			HMOM_Source_Condensation_M10Cells[celli] = hmom.sources_condensation()[1];
			HMOM_Source_Condensation_M01Cells[celli] = hmom.sources_condensation()[2];
			HMOM_Source_Condensation_N0Cells[celli]  = hmom.sources_condensation()[3];

			// Coagulation source terms [mol/m3/s]
			HMOM_Source_Coagulation_M00Cells[celli] = hmom.sources_coagulation_overall()[0];
			HMOM_Source_Coagulation_M10Cells[celli] = hmom.sources_coagulation_overall()[1];
			HMOM_Source_Coagulation_M01Cells[celli] = hmom.sources_coagulation_overall()[2];
			HMOM_Source_Coagulation_N0Cells[celli]  = hmom.sources_coagulation_overall()[3];
		}		
	}

	Info << " * boundaries..." << endl;
	forAll(T.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT  = T.boundaryField()[patchi];
		const fvPatchScalarField& pp  = p.boundaryField()[patchi];

		const fvPatchScalarField& pHMOM_Norm_N0  = HMOM_Norm_N0.boundaryField()[patchi];
		const fvPatchScalarField& pHMOM_Norm_M00 = HMOM_Norm_M00.boundaryField()[patchi];
		const fvPatchScalarField& pHMOM_Norm_M01 = HMOM_Norm_M01.boundaryField()[patchi];
		const fvPatchScalarField& pHMOM_Norm_M10 = HMOM_Norm_M10.boundaryField()[patchi];

		fvPatchScalarField& pHMOM_fv = HMOM_fv.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_N  = HMOM_N.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_dp = HMOM_dp.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_dc = HMOM_dc.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_np = HMOM_np.boundaryField()[patchi];

		fvPatchScalarField& pHMOM_Source_M00 = HMOM_Source_M00.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_M10 = HMOM_Source_M10.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_M01 = HMOM_Source_M01.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_N0  = HMOM_Source_N0.boundaryField()[patchi]; 

		fvPatchScalarField& pHMOM_Source_Nucleation_M00 = HMOM_Source_Nucleation_M00.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Nucleation_M10 = HMOM_Source_Nucleation_M10.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Nucleation_M01 = HMOM_Source_Nucleation_M01.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Nucleation_N0  = HMOM_Source_Nucleation_N0.boundaryField()[patchi]; 

		fvPatchScalarField& pHMOM_Source_Growth_M00 = HMOM_Source_Growth_M00.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Growth_M10 = HMOM_Source_Growth_M10.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Growth_M01 = HMOM_Source_Growth_M01.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Growth_N0  = HMOM_Source_Growth_N0.boundaryField()[patchi]; 

		fvPatchScalarField& pHMOM_Source_Oxidation_M00 = HMOM_Source_Oxidation_M00.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Oxidation_M10 = HMOM_Source_Oxidation_M10.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Oxidation_M01 = HMOM_Source_Oxidation_M01.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Oxidation_N0  = HMOM_Source_Oxidation_N0.boundaryField()[patchi]; 

		fvPatchScalarField& pHMOM_Source_Condensation_M00 = HMOM_Source_Condensation_M00.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Condensation_M10 = HMOM_Source_Condensation_M10.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Condensation_M01 = HMOM_Source_Condensation_M01.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Condensation_N0  = HMOM_Source_Condensation_N0.boundaryField()[patchi]; 

		fvPatchScalarField& pHMOM_Source_Coagulation_M00 = HMOM_Source_Coagulation_M00.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Coagulation_M10 = HMOM_Source_Coagulation_M10.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Coagulation_M01 = HMOM_Source_Coagulation_M01.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_Source_Coagulation_N0  = HMOM_Source_Coagulation_N0.boundaryField()[patchi]; 
	
		forAll(pT, facei)
		{
			double mw;

			// Extract the mass fractions
			for(unsigned int i=0;i<ns;i++)
				y[i+1] = Y[i].boundaryField()[patchi][facei];
			const double sum = y.SumElements();
			for(unsigned int i=0;i<ns;i++)
				y[i+1] /= sum;

			// Molecular weight
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

			// Eigen vectors
			for(unsigned int i=0;i<ns;i++)
			{
				y_eigen(i) = y[i+1];
				x_eigen(i) = x[i+1];
			}

			// Density
			const double rhoGas = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];

			// Concentrations
			const double cTot = pp[facei]/PhysicalConstants::R_J_kmol/pT[facei];
			Product(cTot,x, &c);

			// HMOM
			if (hmom_analysis == true)
			{
				// Viscosity
				const double viscosity = 1.716e-5*Foam::pow(pT[facei]/273.15,1.5)*(273.15+110.4)/(pT[facei]+110.4);		// [kg/m/s]

				// PAH concentration
				double conc_PAH = 0.;
				for (unsigned int j = 0; j<hmom.pah_species().size(); j++)
					conc_PAH += cTot*x[hmom_pah_species_indices[j] + 1];

				// Setup
				hmom.SetNormalizedMoments(pHMOM_Norm_M00[facei], pHMOM_Norm_M10[facei], pHMOM_Norm_M01[facei], pHMOM_Norm_N0[facei]);
				hmom.SetTemperatureAndPressure(pT[facei], pp[facei]);
				hmom.SetMassFractions(y[hmom_index_OH+1], y[hmom_index_H+1]);
				hmom.SetConcentrations("kmol/m3", c[hmom_index_OH+1], c[hmom_index_H+1], c[hmom_index_H2O+1], c[hmom_index_H2+1], c[hmom_index_C2H2+1], c[hmom_index_O2+1], conc_PAH);
				hmom.SetViscosity(viscosity);

				// Calculates
				hmom.CalculateSourceMoments();

				// Calculating soot volume fraction [-]
				pHMOM_fv[facei] = hmom.SootVolumeFraction();

				// Calculating soot particle number density [#/m3]
				pHMOM_N[facei] = hmom.SootParticleNumberDensity();

				// Calculating particle diameters [nm]
				pHMOM_dp[facei] = hmom.SootParticleDiameter() * 1.e9;

				// Calculating collision diameters [nm]
				pHMOM_dc[facei] = hmom.SootCollisionParticleDiameter() * 1.e9;
	
				// Calculating number of primary particles [-]
				pHMOM_np[facei] = hmom.SootNumberOfPrimaryParticles();

				// Overall source terms [mol/m3/s]
				pHMOM_Source_M00[facei] = hmom.sources()[0];
				pHMOM_Source_M10[facei] = hmom.sources()[1];
				pHMOM_Source_M01[facei] = hmom.sources()[2];
				pHMOM_Source_N0[facei]  = hmom.sources()[3];

				// Nucleation source terms [mol/m3/s]
				pHMOM_Source_Nucleation_M00[facei] = hmom.sources_nucleation()[0];
				pHMOM_Source_Nucleation_M10[facei] = hmom.sources_nucleation()[1];
				pHMOM_Source_Nucleation_M01[facei] = hmom.sources_nucleation()[2];
				pHMOM_Source_Nucleation_N0[facei]  = hmom.sources_nucleation()[3];

				// Growth source terms [mol/m3/s]
				pHMOM_Source_Growth_M00[facei] = hmom.sources_growth()[0];
				pHMOM_Source_Growth_M10[facei] = hmom.sources_growth()[1];
				pHMOM_Source_Growth_M01[facei] = hmom.sources_growth()[2];
				pHMOM_Source_Growth_N0[facei]  = hmom.sources_growth()[3];

				// Oxidation source terms [mol/m3/s]
				pHMOM_Source_Oxidation_M00[facei] = hmom.sources_oxidation()[0];
				pHMOM_Source_Oxidation_M10[facei] = hmom.sources_oxidation()[1];
				pHMOM_Source_Oxidation_M01[facei] = hmom.sources_oxidation()[2];
				pHMOM_Source_Oxidation_N0[facei]  = hmom.sources_oxidation()[3];

				// Condensation source terms [mol/m3/s]
				pHMOM_Source_Condensation_M00[facei] = hmom.sources_condensation()[0];
				pHMOM_Source_Condensation_M10[facei] = hmom.sources_condensation()[1];
				pHMOM_Source_Condensation_M01[facei] = hmom.sources_condensation()[2];
				pHMOM_Source_Condensation_N0[facei]  = hmom.sources_condensation()[3];

				// Coagulation source terms [mol/m3/s]
				pHMOM_Source_Coagulation_M00[facei] = hmom.sources_coagulation_overall()[0];
				pHMOM_Source_Coagulation_M10[facei] = hmom.sources_coagulation_overall()[1];
				pHMOM_Source_Coagulation_M01[facei] = hmom.sources_coagulation_overall()[2];
				pHMOM_Source_Coagulation_N0[facei]  = hmom.sources_coagulation_overall()[3];
			}
		}
	}

	Info << " * writing on files..." << endl;
	{
		// HMOM
		HMOM_fv.write();
		HMOM_N.write();
		HMOM_dp.write();
		HMOM_dc.write();
 		HMOM_np.write();

		// Overall source terms [mol/m3/s]
		HMOM_Source_M00.write();
		HMOM_Source_M10.write();
		HMOM_Source_M01.write();
		HMOM_Source_N0.write();

		// Nucleation source terms [mol/m3/s]
		HMOM_Source_Nucleation_M00.write();
		HMOM_Source_Nucleation_M10.write();
		HMOM_Source_Nucleation_M01.write();
		HMOM_Source_Nucleation_N0.write();

		// Growth source terms [mol/m3/s]
		HMOM_Source_Growth_M00.write();
		HMOM_Source_Growth_M10.write();
		HMOM_Source_Growth_M01.write();
		HMOM_Source_Growth_N0.write();

		// Oxidation source terms [mol/m3/s]
		HMOM_Source_Oxidation_M00.write();
		HMOM_Source_Oxidation_M10.write();
		HMOM_Source_Oxidation_M01.write();
		HMOM_Source_Oxidation_N0.write();

		// Condensation source terms [mol/m3/s]
		HMOM_Source_Condensation_M00.write();
		HMOM_Source_Condensation_M10.write();
		HMOM_Source_Condensation_M01.write();
		HMOM_Source_Condensation_N0.write();

		// Coagulation source terms [mol/m3/s]
		HMOM_Source_Coagulation_M00.write();
		HMOM_Source_Coagulation_M10.write();
		HMOM_Source_Coagulation_M01.write();
		HMOM_Source_Coagulation_N0.write();
	}
}
