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

Info << endl;
Info<< " * Reading chemical species mole fractions or mass fractions..." << endl;

PtrList<volScalarField> Y(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> X(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> Dmix(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> CpSpecies(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> FormationRates(outputFormationRatesIndices.size());
PtrList<volScalarField> RR(thermodynamicsMapXML->NumberOfSpecies());

#if STEADYSTATE == 1
PtrList<volScalarField> sourceImplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
PtrList<volScalarField> sourceExplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
#endif

PtrList<volScalarField> Dsoret(transportMapXML->iThermalDiffusionRatios().size());
if (soretEffect == true)
{
	for(unsigned int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
	{ 
		unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
		indexSoret[index-1] = i;
		soretEffectList[index-1] = true;

		std::cout << "Soret effect: species " 	<< i+1 << "/" << transportMapXML->iThermalDiffusionRatios().size() 
							<< " " << thermodynamicsMapXML->NamesOfSpecies()[index-1] << " (" << index << ")" << std::endl;

		Dsoret.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"gas_Dsoret_"+ thermodynamicsMapXML->NamesOfSpecies()[index-1],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar("Dsoret", dimensionSet(0, 2, -1, 0, 0), 1e-10)
			)
		);
	}
}

// Loop over all the species in the kinetic mechanism
for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{
	Info << "   Reading species: " << thermodynamicsMapXML->NamesOfSpecies()[i] << " > ";

	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	// Sets of diffusivities and specific heats of species
	{
		Dmix.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"gas_Dmix_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("Dmix", dimensionSet(0, 2, -1, 0, 0), 1e-10)
			)
		);

		CpSpecies.set
                (
                        i,
                        new volScalarField
                        (
                                IOobject
                                (
                                "thermo_Cp_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
                                mesh.time().timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                                ),
                                mesh,
                                dimensionSet(0, 2, -2, -1, 0)
                        )
                );

		X.set
                (
                        i,
                        new volScalarField
                        (
                                IOobject
                                (
                                "thermo_X_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
                                mesh.time().timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                                ),
                                mesh,
                                dimensionSet(0, 0, 0, 0, 0)
                        )
                );

		#if STEADYSTATE == 0
		RR.set
                (
                        i,
                        new volScalarField
                        (
                                IOobject
                                (
                                "thermo_RR_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
                                mesh.time().timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                                ),
                                mesh,
                                dimensionedScalar("RR", dimensionSet(1, -3, -1, 0, 0), 0.0)
                        )
                );
		#endif
	}

	// check if field exists and can be read
	if (header.headerOk())
	{
		Info << thermodynamicsMapXML->NamesOfSpecies()[i] << endl;
		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);
	}
	else
	{
		Info << "Ydefault" << endl;

		volScalarField Ydefault
		(
			IOobject
			(
				"Ydefault",
				mesh.time().timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);

		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				Ydefault
			)
		);
	}
}

#if STEADYSTATE == 0

Info<< "Creating field RT\n" << endl;
volScalarField RT
(
	IOobject
	(
		"RT",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("RT", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

#endif

#if STEADYSTATE == 1

for (int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{
	sourceImplicit.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"SI_" + thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("SI", dimensionSet(1, -3, -1, 0, 0), 0.0)
		)
	);

	sourceExplicit.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"SE_" + thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("SE", dimensionSet(1, -3, -1, 0, 0), 0.0)
		)
	);
}

sourceImplicit.set
(
	thermodynamicsMapXML->NumberOfSpecies(),
	new volScalarField
	(
		IOobject
		(
			"SI_T",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("SIT", dimensionSet(1, -1, -3, -1, 0), 0.0)
	)
);

sourceExplicit.set
(
	thermodynamicsMapXML->NumberOfSpecies(),
	new volScalarField
	(
		IOobject
		(
			"SE_T",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("SET", dimensionSet(1, -1, -3, 0, 0), 0.0)
	)
);

#endif

if (outputFormationRatesIndices.size() != 0)
{
	for (int i=0;i<outputFormationRatesIndices.size();i++)
	{
		FormationRates.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"R_" + thermodynamicsMapXML->NamesOfSpecies()[outputFormationRatesIndices(i)],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("R", dimensionSet(1, -3, -1, 0, 0), 0.0)
			)
		);
	}
}

PtrList<volScalarField> sootFields(18);
if (sootOnTheFlyFill == true)
{
	sootFields.set
	(
		0,
		new volScalarField
		(
			IOobject
			(
				"soot_fv_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("soot_fv_large", dimensionSet(0, 0, 0, 0, 0), 0.)
		)
	);

	sootFields.set
	(
		1,
		new volScalarField
		(
			IOobject
			(
				"soot_fv_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("soot_fv_small", dimensionSet(0, 0, 0, 0, 0), 0.)
		)
	);

	sootFields.set
	(
		2,
		new volScalarField
		(
			IOobject
			(
				"soot_rho_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, 0, 0, 0)
		)
	);

	sootFields.set
	(
		3,
		new volScalarField
		(
			IOobject
			(
				"soot_rho_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, 0, 0, 0)
		)
	);

	sootFields.set
	(
		4,
		new volScalarField
		(
			IOobject
			(
				"soot_N_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, -3, 0, 0, 0)
		)
	);

	sootFields.set
	(
		5,
		new volScalarField
		(
			IOobject
			(
				"soot_N_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, -3, 0, 0, 0)
		)
	);

	sootFields.set
	(
		6,
		new volScalarField
		(
			IOobject
			(
				"soot_x_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		7,
		new volScalarField
		(
			IOobject
			(
				"soot_x_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		8,
		new volScalarField
		(
			IOobject
			(
				"soot_omega_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		9,
		new volScalarField
		(
			IOobject
			(
				"soot_omega_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		10,
		new volScalarField
		(
			IOobject
			(
				"pah_omega_1_2",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		11,
		new volScalarField
		(
			IOobject
			(
				"pah_omega_3_4",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		12,
		new volScalarField
		(
			IOobject
			(
				"pah_omega_more_4",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

        sootFields.set
	(
		13,
		new volScalarField
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
		)
	);

        sootFields.set
	(
		14,
		new volScalarField
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
		)
	);

        sootFields.set
	(
		15,
		new volScalarField
		(
			IOobject
			(
				"pah_R_1_2",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, -1, 0, 0)
		)
	);

        sootFields.set
	(
		16,
		new volScalarField
		(
			IOobject
			(
				"pah_R_3_4",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, -1, 0, 0)
		)
	);

        sootFields.set
	(
		17,
		new volScalarField
		(
			IOobject
			(
				"pah_R_more_4",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, -1, 0, 0)
		)
	);
}

volScalarField cpuChemistry
(
        IOobject
        (
                "cpuChemistry",
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


#include "createFluxes.H"
