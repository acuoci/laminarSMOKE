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
PtrList<volScalarField> CpSpecies(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> FormationRates(outputFormationRatesIndices.size());
PtrList<volScalarField> RR(thermodynamicsMapXML->NumberOfSpecies());

#if STEADYSTATE == 1
PtrList<volScalarField> sourceImplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
PtrList<volScalarField> sourceExplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
#endif


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
	#if DEVVERSION == 1
	if (header.typeHeaderOk<volScalarField>(true))
	#else
	if (header.headerOk())
	#endif
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
