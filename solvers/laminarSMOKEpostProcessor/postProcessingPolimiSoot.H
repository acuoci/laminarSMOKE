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

// Post processing: soot (Polimi)
if (postProcessingPolimiSoot == true)
{
	volScalarField soot_fv_large
	(
		IOobject
		(
			"soot_fv_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_fv_small
	(
		IOobject
		(
			"soot_fv_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_rho_large
	(
		IOobject
		(
			"soot_rho_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(1, -3, 0, 0, 0)
	);
	volScalarField soot_rho_small
	(
		IOobject
		(
			"soot_rho_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(1, -3, 0, 0, 0)
	);
	volScalarField soot_N_large
	(
		IOobject
		(
			"soot_N_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, -3, 0, 0, 0)
	);
	volScalarField soot_N_small
	(
		IOobject
		(
			"soot_N_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, -3, 0, 0, 0)
	);
	volScalarField soot_omega_large
	(
		IOobject
		(
			"soot_omega_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_omega_small
	(
		IOobject
		(
			"soot_omega_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_x_large
	(
		IOobject
		(
			"soot_x_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_x_small
	(
		IOobject
		(
			"soot_x_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_h_over_c_large
	(
		IOobject
		(
			"soot_h_over_c_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_h_over_c_small
	(
		IOobject
		(
			"soot_h_over_c_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_o_over_c_large
	(
		IOobject
		(
			"soot_o_over_c_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_o_over_c_small
	(
		IOobject
		(
			"soot_o_over_c_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_o_over_h_large
	(
		IOobject
		(
			"soot_o_over_h_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);
	volScalarField soot_o_over_h_small
	(
		IOobject
		(
			"soot_o_over_h_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);

        volScalarField soot_omega_precursors
	(
		IOobject
		(
			"soot_omega_precursors",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);

        volScalarField soot_x_precursors
	(
		IOobject
		(
			"soot_x_precursors",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionSet(0, 0, 0, 0, 0)
	);

	volScalarField soot_R_large
	(
		IOobject
		(
			"soot_R_large",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionSet(1, -3, -1, 0, 0)
	);

	volScalarField soot_R_small
	(
		IOobject
		(
			"soot_R_small",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionSet(1, -3, -1, 0, 0)
	);

 	volScalarField pah_R_1_2
	(
		IOobject
		(
			"pah::R_1_2",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionSet(1, -3, -1, 0, 0)
	);

	volScalarField pah_R_3_4
	(
		IOobject
		(
			"pah::R_3_4",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionSet(1, -3, -1, 0, 0)
	);

	volScalarField pah_R_more_4
	(
		IOobject
		(
			"pah::R_more_4",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionSet(1, -3, -1, 0, 0)
	);

	// --------------------------------------------------------------------------
	// Creating soot class fields
	// --------------------------------------------------------------------------
	PtrList<volScalarField> sootClasses;
	std::vector<double> sootClassesIntegrals;
	if (calculateSootClasses == true)
	{
		Info << "Post processing soot classes..." << endl;

		sootClasses.resize(soot_classes_reader.number_of_classes());
		for (int i=0;i<soot_classes_reader.number_of_classes();i++)
		{
			sootClasses.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"soot__" + soot_classes_reader.class_name(i),
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					mesh,
					dimensionSet(0, -3, -1, 0, 1)
				)
			);
		}

		sootClassesIntegrals.resize(soot_classes_reader.number_of_classes());
		for (int i=0;i<soot_classes_reader.number_of_classes();i++)
			sootClassesIntegrals[i] = 0.;
	}

	Info << "Post processing soot..." << endl;

	label ns = Y.size();

	OpenSMOKE::OpenSMOKEVectorDouble y(ns);
	OpenSMOKE::OpenSMOKEVectorDouble x(ns);
	OpenSMOKE::OpenSMOKEVectorDouble c(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble R(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble r(kineticsMapXML->NumberOfReactions());
	Eigen::VectorXd y_eigen(ns);
	Eigen::VectorXd x_eigen(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();

	#if OPENFOAM_VERSION >= 40
	scalarField& soot_fv_largeCells = soot_fv_large.ref();
	scalarField& soot_fv_smallCells = soot_fv_small.ref();
	scalarField& soot_rho_largeCells = soot_rho_large.ref();
	scalarField& soot_rho_smallCells = soot_rho_small.ref();
	scalarField& soot_N_largeCells = soot_N_large.ref();
	scalarField& soot_N_smallCells = soot_N_small.ref();
	scalarField& soot_x_largeCells = soot_x_large.ref();
	scalarField& soot_x_smallCells = soot_x_small.ref();
	scalarField& soot_omega_largeCells = soot_omega_large.ref();
	scalarField& soot_omega_smallCells = soot_omega_small.ref();
	scalarField& soot_h_over_c_largeCells = soot_h_over_c_large.ref();
	scalarField& soot_h_over_c_smallCells = soot_h_over_c_small.ref();
	scalarField& soot_o_over_c_largeCells = soot_o_over_c_large.ref();
	scalarField& soot_o_over_c_smallCells = soot_o_over_c_small.ref();
	scalarField& soot_o_over_h_largeCells = soot_o_over_h_large.ref();
	scalarField& soot_o_over_h_smallCells = soot_o_over_h_small.ref();
	scalarField& soot_omega_precursorsCells = soot_omega_precursors.ref();
	scalarField& soot_x_precursorsCells = soot_x_precursors.ref();

	scalarField& soot_R_largeCells = soot_R_large.ref();
	scalarField& soot_R_smallCells = soot_R_small.ref();
	scalarField& pah_R_1_2Cells = pah_R_1_2.ref();
	scalarField& pah_R_3_4Cells = pah_R_3_4.ref();
	scalarField& pah_R_more_4Cells = pah_R_more_4.ref();
	#else
	scalarField& soot_fv_largeCells = soot_fv_large.internalField();
	scalarField& soot_fv_smallCells = soot_fv_small.internalField();
	scalarField& soot_rho_largeCells = soot_rho_large.internalField();
	scalarField& soot_rho_smallCells = soot_rho_small.internalField();
	scalarField& soot_N_largeCells = soot_N_large.internalField();
	scalarField& soot_N_smallCells = soot_N_small.internalField();
	scalarField& soot_x_largeCells = soot_x_large.internalField();
	scalarField& soot_x_smallCells = soot_x_small.internalField();
	scalarField& soot_omega_largeCells = soot_omega_large.internalField();
	scalarField& soot_omega_smallCells = soot_omega_small.internalField();
	scalarField& soot_h_over_c_largeCells = soot_h_over_c_large.internalField();
	scalarField& soot_h_over_c_smallCells = soot_h_over_c_small.internalField();
	scalarField& soot_o_over_c_largeCells = soot_o_over_c_large.internalField();
	scalarField& soot_o_over_c_smallCells = soot_o_over_c_small.internalField();
	scalarField& soot_o_over_h_largeCells = soot_o_over_h_large.internalField();
	scalarField& soot_o_over_h_smallCells = soot_o_over_h_small.internalField();
	scalarField& soot_omega_precursorsCells = soot_omega_precursors.internalField();
	scalarField& soot_x_precursorsCells = soot_x_precursors.internalField();

	scalarField& soot_R_largeCells = soot_R_large.internalField();
	scalarField& soot_R_smallCells = soot_R_small.internalField();
	scalarField& pah_R_1_2Cells = pah_R_1_2.internalField();
	scalarField& pah_R_3_4Cells = pah_R_3_4.internalField();
	scalarField& pah_R_more_4Cells = pah_R_more_4.internalField();
	#endif

	// Total volume
	double Vtot = 0.;
	forAll(TCells, celli)
		Vtot += mesh.V()[celli];

	// Integral quantities
	double sootFvLargeIntegral 	= 0;
	double sootFvSmallIntegral 	= 0;
	double sootRhoLargeIntegral 	= 0;
	double sootRhoSmallIntegral 	= 0;
	double sootNLargeIntegral 	= 0;
	double sootNSmallIntegral 	= 0;
	double 	sootOmegaLargeIntegral 	= 0;
	double sootOmegaSmallIntegral 	= 0;
	double sootRLargeIntegral  	= 0;
	double sootRSmallIntegral  	= 0;
	double sootR12Integral  	= 0;
	double sootR34Integral 	 	= 0;
	double sootRmore4Integral  	= 0;

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
		thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), mw, y.GetHandle());

		// Eigen vectors
		for(unsigned int i=0;i<ns;i++)
		{
			y_eigen(i) = y[i+1];
			x_eigen(i) = x[i+1];
		}

		// Density
		const double rhoGas = pCells[celli]*mw/PhysicalConstants::R_J_kmol/TCells[celli];

		// Soot analysis
		sootAnalyzer->Analysis(TCells[celli], pCells[celli], rhoGas, y_eigen, x_eigen);
		soot_fv_largeCells[celli] = sootAnalyzer->fv_large();
		soot_fv_smallCells[celli] = sootAnalyzer->fv_small();
		soot_rho_largeCells[celli] = sootAnalyzer->rho_large();
		soot_rho_smallCells[celli] = sootAnalyzer->rho_small();
		soot_N_largeCells[celli] = sootAnalyzer->N_large();
		soot_N_smallCells[celli] = sootAnalyzer->N_small();
		soot_omega_largeCells[celli] = sootAnalyzer->omega_large();
		soot_omega_smallCells[celli] = sootAnalyzer->omega_small();
		soot_x_largeCells[celli] = sootAnalyzer->x_large();
		soot_x_smallCells[celli] = sootAnalyzer->x_small();
		soot_h_over_c_largeCells[celli] = sootAnalyzer->h_over_c_large();
		soot_h_over_c_smallCells[celli] = sootAnalyzer->h_over_c_small();
		soot_o_over_c_largeCells[celli] = sootAnalyzer->o_over_c_large();
		soot_o_over_c_smallCells[celli] = sootAnalyzer->o_over_c_small();
		soot_o_over_h_largeCells[celli] = sootAnalyzer->o_over_h_large();
		soot_o_over_h_smallCells[celli] = sootAnalyzer->o_over_h_small();

		soot_omega_precursorsCells[celli] = 0.;
		soot_x_precursorsCells[celli] = 0.;
		for(unsigned int k=0;k<soot_precursors_indices.size();k++)
		{
			int index = soot_precursors_indices[k];
			soot_omega_precursorsCells[celli] += y_eigen(index);
			soot_x_precursorsCells[celli] += x_eigen(index);
		}

		// Concentrations
		const double cTot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
		Product(cTot,x, &c);

		// Kinetics
		kineticsMapXML->SetTemperature(TCells[celli]);
		kineticsMapXML->SetPressure(pCells[celli]);
		kineticsMapXML->KineticConstants();
		kineticsMapXML->ReactionRates(c.GetHandle());
		kineticsMapXML->FormationRates(R.GetHandle());
		kineticsMapXML->GiveMeReactionRates(r.GetHandle());
		for(unsigned int i=1;i<=ns;i++)
			R[i] *= thermodynamicsMapXML->MW(i-1);

		soot_R_largeCells[celli] = largeBinMassFormationRates(sootAnalyzer, R);
		soot_R_smallCells[celli] = smallBinMassFormationRates(sootAnalyzer, R);
		pah_R_1_2Cells[celli]    = pah12RingsFormationRates(sootAnalyzer, R);
		pah_R_3_4Cells[celli]    = pah34RingsFormationRates(sootAnalyzer, R);
		pah_R_more_4Cells[celli] = pahMoreThan4RingsFormationRates(sootAnalyzer, R);	


		sootFvLargeIntegral 	+= soot_fv_largeCells[celli]*mesh.V()[celli];
		sootFvSmallIntegral 	+= soot_fv_smallCells[celli]*mesh.V()[celli];
		sootRhoLargeIntegral 	+= soot_rho_largeCells[celli]*mesh.V()[celli];
		sootRhoSmallIntegral 	+= soot_rho_smallCells[celli]*mesh.V()[celli];
		sootNLargeIntegral 	+= soot_N_largeCells[celli]*mesh.V()[celli];
		sootNSmallIntegral 	+= soot_N_largeCells[celli]*mesh.V()[celli];
		sootOmegaLargeIntegral  += soot_omega_largeCells[celli]*mesh.V()[celli];
		sootOmegaSmallIntegral  += soot_omega_largeCells[celli]*mesh.V()[celli];
		sootRLargeIntegral 	+= soot_R_largeCells[celli]*mesh.V()[celli];
		sootRSmallIntegral 	+= soot_R_smallCells[celli]*mesh.V()[celli];
		sootR12Integral 	+= pah_R_1_2Cells[celli]*mesh.V()[celli];
		sootR34Integral 	+= pah_R_3_4Cells[celli]*mesh.V()[celli];
		sootRmore4Integral 	+= pah_R_more_4Cells[celli]*mesh.V()[celli];

		if (calculateSootClasses == true)
		{
			for (int i=0;i<soot_classes_reader.number_of_classes();i++)
			{
				double sum = 0.;
				for (int j=0;j<soot_classes_reader.number_reactions_per_class(i);j++)
				{
					sum += r[soot_classes_reader.reaction_indices(i)[j]];
				}
			
				#if OPENFOAM_VERSION >= 40
				sootClasses[i].ref()[celli] = sum;
				#else
				sootClasses[i].internalField()[celli] = sum;
				#endif	
				sootClassesIntegrals[i] += sum*mesh.V()[celli];	
			}
		}
	}

	Info << " * boundaries..." << endl;
	forAll(T.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT  = T.boundaryField()[patchi];
		const fvPatchScalarField& pp  = p.boundaryField()[patchi];
		
		#if OPENFOAM_VERSION >= 40
		fvPatchScalarField& psoot_fv_large = soot_fv_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_fv_small = soot_fv_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_rho_large = soot_rho_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_rho_small = soot_rho_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_N_large = soot_N_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_N_small = soot_N_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_omega_large = soot_omega_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_omega_small = soot_omega_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_x_large = soot_x_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_x_small = soot_x_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_h_over_c_large = soot_h_over_c_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_h_over_c_small = soot_h_over_c_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_o_over_c_large = soot_o_over_c_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_o_over_c_small = soot_o_over_c_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_o_over_h_large = soot_o_over_h_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_o_over_h_small = soot_o_over_h_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_x_precursors = soot_x_precursors.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_omega_precursors = soot_omega_precursors.boundaryFieldRef()[patchi];

		fvPatchScalarField& psoot_R_large = soot_R_large.boundaryFieldRef()[patchi];
		fvPatchScalarField& psoot_R_small = soot_R_small.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppah_R_1_2 = pah_R_1_2.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppah_R_3_4 = pah_R_3_4.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppah_R_more_4 = pah_R_more_4.boundaryFieldRef()[patchi];
		#else
		fvPatchScalarField& psoot_fv_large = soot_fv_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_fv_small = soot_fv_small.boundaryField()[patchi];
		fvPatchScalarField& psoot_rho_large = soot_rho_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_rho_small = soot_rho_small.boundaryField()[patchi];
		fvPatchScalarField& psoot_N_large = soot_N_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_N_small = soot_N_small.boundaryField()[patchi];
		fvPatchScalarField& psoot_omega_large = soot_omega_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_omega_small = soot_omega_small.boundaryField()[patchi];
		fvPatchScalarField& psoot_x_large = soot_x_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_x_small = soot_x_small.boundaryField()[patchi];
		fvPatchScalarField& psoot_h_over_c_large = soot_h_over_c_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_h_over_c_small = soot_h_over_c_small.boundaryField()[patchi];
		fvPatchScalarField& psoot_o_over_c_large = soot_o_over_c_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_o_over_c_small = soot_o_over_c_small.boundaryField()[patchi];
		fvPatchScalarField& psoot_o_over_h_large = soot_o_over_h_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_o_over_h_small = soot_o_over_h_small.boundaryField()[patchi];
		fvPatchScalarField& psoot_x_precursors = soot_x_precursors.boundaryField()[patchi];
		fvPatchScalarField& psoot_omega_precursors = soot_omega_precursors.boundaryField()[patchi];

		fvPatchScalarField& psoot_R_large = soot_R_large.boundaryField()[patchi];
		fvPatchScalarField& psoot_R_small = soot_R_small.boundaryField()[patchi];
		fvPatchScalarField& ppah_R_1_2 = pah_R_1_2.boundaryField()[patchi];
		fvPatchScalarField& ppah_R_3_4 = pah_R_3_4.boundaryField()[patchi];
		fvPatchScalarField& ppah_R_more_4 = pah_R_more_4.boundaryField()[patchi];
		#endif

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
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), mw, y.GetHandle());

			// Eigen vectors
			for(unsigned int i=0;i<ns;i++)
			{
				y_eigen(i) = y[i+1];
				x_eigen(i) = x[i+1];
			}

			// Density
			const double rhoGas = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];

			// Soot analysis
			sootAnalyzer->Analysis(pT[facei], pp[facei], rhoGas, y_eigen, x_eigen);	
			psoot_fv_large[facei] = sootAnalyzer->fv_large();
			psoot_fv_small[facei] = sootAnalyzer->fv_small();
			psoot_rho_large[facei] = sootAnalyzer->rho_large();
			psoot_rho_small[facei] = sootAnalyzer->rho_small();
			psoot_N_large[facei] = sootAnalyzer->N_large();
			psoot_N_small[facei] = sootAnalyzer->N_small();
			psoot_omega_large[facei] = sootAnalyzer->omega_large();
			psoot_omega_small[facei] = sootAnalyzer->omega_small();
			psoot_x_large[facei] = sootAnalyzer->x_large();
			psoot_x_small[facei] = sootAnalyzer->x_small();
			psoot_h_over_c_large[facei] = sootAnalyzer->h_over_c_large();
			psoot_h_over_c_small[facei] = sootAnalyzer->h_over_c_small();
			psoot_o_over_c_large[facei] = sootAnalyzer->o_over_c_large();
			psoot_o_over_c_small[facei] = sootAnalyzer->o_over_c_small();
			psoot_o_over_h_large[facei] = sootAnalyzer->o_over_h_large();
			psoot_o_over_h_small[facei] = sootAnalyzer->o_over_h_small();

			// Soot precursors
        		psoot_omega_precursors[facei] = 0.;
			psoot_x_precursors[facei] = 0.;
			for(unsigned int k=0;k<soot_precursors_indices.size();k++)
			{
				int index = soot_precursors_indices[k];
				psoot_omega_precursors[facei] += y_eigen(index);
				psoot_x_precursors[facei] += x_eigen(index);
			}

			// Concentrations
			const double cTot = pp[facei]/PhysicalConstants::R_J_kmol/pT[facei];
			Product(cTot,x, &c);

			// Kinetics
			kineticsMapXML->SetTemperature(pT[facei]);
			kineticsMapXML->SetPressure(pp[facei]);
			kineticsMapXML->KineticConstants();
			kineticsMapXML->ReactionRates(c.GetHandle());
			kineticsMapXML->FormationRates(R.GetHandle());
			kineticsMapXML->GiveMeReactionRates(r.GetHandle());
			for(unsigned int i=1;i<=ns;i++)
				R[i] *= thermodynamicsMapXML->MW(i-1);

			psoot_R_large[facei] = largeBinMassFormationRates(sootAnalyzer, R);
			psoot_R_small[facei] = smallBinMassFormationRates(sootAnalyzer, R);
			ppah_R_1_2[facei]    = pah12RingsFormationRates(sootAnalyzer, R);
			ppah_R_3_4[facei]    = pah34RingsFormationRates(sootAnalyzer, R);
			ppah_R_more_4[facei] = pahMoreThan4RingsFormationRates(sootAnalyzer, R);

			if (calculateSootClasses == true)
			{
				for (int i=0;i<soot_classes_reader.number_of_classes();i++)
				{
					double sum = 0.;
					for (int j=0;j<soot_classes_reader.number_reactions_per_class(i);j++)
					{
						sum += r[soot_classes_reader.reaction_indices(i)[j]];
					}

					#if OPENFOAM_VERSION >= 40
					sootClasses[i].boundaryFieldRef()[patchi][facei] = sum;
					#else
					sootClasses[i].boundaryField()[patchi][facei] = sum;
					#endif				
				}
			}
		}
	}

	Info << " * writing on files..." << endl;
	{
		soot_fv_large.write();	
		soot_fv_small.write();
		soot_rho_large.write();	
		soot_rho_small.write();	
		soot_N_large.write();	
		soot_N_small.write();	
		soot_omega_large.write();	
		soot_omega_small.write();
		soot_x_large.write();	
		soot_x_small.write();
		soot_h_over_c_large.write();
		soot_h_over_c_small.write();
		soot_o_over_c_large.write();
		soot_o_over_c_small.write();
		soot_o_over_h_large.write();
		soot_o_over_h_small.write();	
		soot_omega_precursors.write();
		soot_x_precursors.write();

		soot_R_large.write();
		soot_R_small.write();
		pah_R_1_2.write();
		pah_R_3_4.write();
		pah_R_more_4.write();

		fSootFvLarge << runTime.timeName() 	<< " " << sootFvLargeIntegral 	<< " " << sootFvLargeIntegral/Vtot << endl;
		fSootFvSmall << runTime.timeName() 	<< " " << sootFvSmallIntegral 	<< " " << sootFvSmallIntegral/Vtot << endl;
		fSootRhoLarge << runTime.timeName() 	<< " " << sootRhoLargeIntegral << " " << sootRhoLargeIntegral/Vtot << endl;
		fSootRhoSmall << runTime.timeName() 	<< " " << sootRhoSmallIntegral << " " << sootRhoSmallIntegral/Vtot << endl;
		fSootNLarge << runTime.timeName() 	<< " " << sootNLargeIntegral 	<< " " << sootNLargeIntegral/Vtot << endl;
		fSootNSmall << runTime.timeName() 	<< " " << sootNSmallIntegral 	<< " " << sootNSmallIntegral/Vtot << endl;
		fSootOmegaLarge << runTime.timeName() 	<< " " << sootOmegaLargeIntegral << " " << sootOmegaLargeIntegral/Vtot << endl;
		fSootOmegaSmall << runTime.timeName() 	<< " " << sootOmegaSmallIntegral << " " << sootOmegaSmallIntegral/Vtot << endl;
		fSootRLarge << runTime.timeName() 	<< " " << sootRLargeIntegral 	<< " " << sootRLargeIntegral/Vtot << endl;
		fSootRSmall << runTime.timeName() 	<< " " << sootRSmallIntegral 	<< " " << sootRSmallIntegral/Vtot << endl;
		fSootR12 << runTime.timeName() 		<< " " << sootR12Integral 	<< " " << sootR12Integral/Vtot << endl;
		fSootR34 << runTime.timeName() 		<< " " << sootR34Integral 	<< " " << sootR34Integral/Vtot << endl;
		fSootRmore4 << runTime.timeName()	<< " " << sootRmore4Integral 	<< " " << sootRmore4Integral/Vtot << endl;

		if (calculateSootClasses == true)
		{
			Info << "* writing on files (soot classes)..." << endl;

			for (int i=0;i<soot_classes_reader.number_of_classes();i++)
				sootClasses[i].write();

			for (int i=0;i<soot_classes_reader.number_of_classes();i++)
				fsootClassesIntegrals[i] << runTime.timeName() << " " << sootClassesIntegrals[i] << " " << sootClassesIntegrals[i]/Vtot << endl;
		}
	}

	// Distribution along a boundary
	if (polimiSootBoundaries.size() > 0)
	{
		Info << "Post processing PSDF along boundaries..." << endl;

		for (int i=0;i<polimiSootBoundaries.size();i++)
		{
			Info << " * boundary:" << polimiSootBoundaries[i] << endl;

			label patchi = mesh.boundaryMesh().findPatchID(polimiSootBoundaries[i]); 

			const fvPatchScalarField& pT  = T.boundaryField()[patchi];
			const fvPatchScalarField& pp  = p.boundaryField()[patchi];
			const vectorField& faceCentres = mesh.Cf().boundaryField()[patchi];

			const std::string name_file = "/soot_PSDF_" + polimiSootBoundaries[i] + ".out";
			fileName outputFile(runTime.path()/runTime.timeName()/name_file.c_str());
			std::ofstream fOutput(outputFile.c_str(), std::ios::out);
			fOutput.setf(std::ios::scientific);
		
			sootAnalyzer->WriteDistributionLabel(fOutput);

			forAll(pT, facei)
			{
				scalar face_x = faceCentres[facei].x();
				scalar face_y = faceCentres[facei].y();
				scalar face_z = faceCentres[facei].z();

				double mw;

				// Extract the mass fractions
				for(unsigned int i=0;i<ns;i++)
					y[i+1] = Y[i].boundaryField()[patchi][facei];
				const double sum = y.SumElements();
				for(unsigned int i=0;i<ns;i++)
					y[i+1] /= sum;
	
				// Molecular weight
				thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), mw, y.GetHandle());

				// Eigen vectors
				for(unsigned int i=0;i<ns;i++)
				{
					y_eigen(i) = y[i+1];
					x_eigen(i) = x[i+1];
				}

				// Density
				const double rhoGas = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];

				// Soot analysis
				sootAnalyzer->Analysis(pT[facei], pp[facei], rhoGas, y_eigen, x_eigen);	
				sootAnalyzer->Distribution();
				sootAnalyzer->WriteDistribution(fOutput, 0., face_x, 0., 0., pT[facei]);
			}
		
			fOutput.close();
		}
	}
	
	// Distribution in a specific point
	if (pnts_soot_psdf.size() > 0)
	{
		Info << "Post processing PSDF in single points..." << endl;

		for (unsigned int i=0;i<pnts_soot_psdf.size();i++)
		{
			const double xx = pnts_soot_psdf[i][0];
			const double yy = pnts_soot_psdf[i][1];
			const double zz = pnts_soot_psdf[i][2];

			Info << " * point " << i+1 << " @ " << xx << " " << yy << " " << zz << endl;

			vector pos(xx,yy,zz);
			label cellI = mesh.findCell(pos);

			if (cellI > -1)
			{
				std::stringstream tag;
				tag << i;
				std::string name_file = "/soot_PSDF_point_" + tag.str();

				// File name
				fileName outputFile(runTime.path()/runTime.timeName()/name_file.c_str());
				std::ofstream fOutput(outputFile.c_str(), std::ios::out);
				fOutput.setf(std::ios::scientific);
				sootAnalyzer->WriteDistributionLabel(fOutput);

				// Temperature and pressure
				autoPtr<interpolation<scalar> > Tinterp = interpolation<scalar>::New("cellPoint", T);
				scalar Tint = Tinterp->interpolate(pos, cellI);
				autoPtr<interpolation<scalar> > pinterp = interpolation<scalar>::New("cellPoint", p);
				scalar pint = pinterp->interpolate(pos, cellI);

				Info << " * cell: " << cellI << " T: " << Tint << " P: " << pint << endl; 


				double mw;

				// Extract the mass fractions
				for(unsigned int i=0;i<ns;i++)
				{
					autoPtr<interpolation<scalar> > Yinterp = interpolation<scalar>::New("cellPoint", Y[i]);
					y[i+1] = Yinterp->interpolate(pos, cellI);
				}

				const double sum = y.SumElements();
				for(unsigned int i=0;i<ns;i++)
					y[i+1] /= sum;

				// Molecular weight
				thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), mw, y.GetHandle());

				// Eigen vectors
				for(unsigned int i=0;i<ns;i++)
				{
					y_eigen(i) = y[i+1];
					x_eigen(i) = x[i+1];
				}

				// Density
				const double rhoGas = pint*mw/PhysicalConstants::R_J_kmol/Tint;

				// Soot analysis
				sootAnalyzer->Analysis(Tint, pint, rhoGas, y_eigen, x_eigen);	
				sootAnalyzer->Distribution();
				sootAnalyzer->WriteDistribution(fOutput, 0., xx, yy, zz, Tint);
				fOutput.close();
			}
		}
	}
}
