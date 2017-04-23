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

volScalarField omega_O
(
	IOobject
	(
		"omega_O",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
    	dimensionedScalar("omega_O", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField omega_H
(
	IOobject
	(
		"omega_H",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
    	dimensionedScalar("omega_H", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField omega_C
(
	IOobject
	(
		"omega_C",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
    	dimensionedScalar("omega_O", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

{
	unsigned int ns = Y.size();
	OpenSMOKE::OpenSMOKEVectorDouble massFractions(ns);
	OpenSMOKE::OpenSMOKEVectorDouble moleFractions(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	#if OPENFOAM_VERSION >= 40
	scalarField& rhoCells = rho.ref();
	scalarField& omega_CCells = omega_C.ref();
	scalarField& omega_HCells = omega_H.ref();
	scalarField& omega_OCells = omega_O.ref();
	#else
	scalarField& rhoCells = rho.internalField();
	scalarField& omega_CCells = omega_C.internalField();
	scalarField& omega_HCells = omega_H.internalField();
	scalarField& omega_OCells = omega_O.internalField();
	#endif


	const unsigned int jC = thermodynamicsMapXML->IndexOfElementWithoutError("C");
	const unsigned int jO = thermodynamicsMapXML->IndexOfElementWithoutError("O");
	const unsigned int jH = thermodynamicsMapXML->IndexOfElementWithoutError("H");

	forAll(TCells, celli)
	{
		// Mass fractions
		for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			massFractions[i+1] = Y[i][celli];

		// Molecular weight
		double MWmix;
		thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(), MWmix, massFractions.GetHandle());

		if (jC>0)
		{
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				omega_CCells[celli] += thermodynamicsMapXML->atomic_composition()(i,jC-1)*moleFractions[i+1]*OpenSMOKE::AtomicWeights["C"]/MWmix;
		}
		if (jO>0)
		{
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				omega_OCells[celli] += thermodynamicsMapXML->atomic_composition()(i,jO-1)*moleFractions[i+1]*OpenSMOKE::AtomicWeights["O"]/MWmix;
		}
		if (jH>0)
		{
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				omega_HCells[celli] += thermodynamicsMapXML->atomic_composition()(i,jH-1)*moleFractions[i+1]*OpenSMOKE::AtomicWeights["H"]/MWmix;
		}	
	}

	const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
	{
		#if OPENFOAM_VERSION >= 40
		fvPatchScalarField& pomega_O  = omega_O.boundaryFieldRef()[patchi];
		fvPatchScalarField& pomega_H  = omega_H.boundaryFieldRef()[patchi];
		fvPatchScalarField& pomega_C  = omega_C.boundaryFieldRef()[patchi];
		#else
		fvPatchScalarField& pomega_O  = omega_O.boundaryField()[patchi];
		fvPatchScalarField& pomega_H  = omega_H.boundaryField()[patchi];
		fvPatchScalarField& pomega_C  = omega_C.boundaryField()[patchi];
		#endif
		
		forAll(pomega_H, facei)
		{
			pomega_H[facei] = 0.;
			pomega_O[facei] = 0.;
			pomega_C[facei] = 0.;
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			// Molecular weight
			double MWmix;
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(), MWmix, massFractions.GetHandle());

			if (jC>0)
			{
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					pomega_C[facei] += thermodynamicsMapXML->atomic_composition()(i,jC-1)*moleFractions[i+1]*OpenSMOKE::AtomicWeights["C"]/MWmix;
			}
			if (jO>0)
			{
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					pomega_O[facei] += thermodynamicsMapXML->atomic_composition()(i,jO-1)*moleFractions[i+1]*OpenSMOKE::AtomicWeights["O"]/MWmix;
			}
			if (jH>0)
			{
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					pomega_H[facei] += thermodynamicsMapXML->atomic_composition()(i,jH-1)*moleFractions[i+1]*OpenSMOKE::AtomicWeights["H"]/MWmix;
			}	
		}
	}
}
