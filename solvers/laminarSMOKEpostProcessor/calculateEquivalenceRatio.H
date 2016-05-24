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

if (calculateEquivalenceRatio == true)
{
	volScalarField equivalenceratio
	(
		IOobject
		(
			"gas_phi",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
	    	dimensionedScalar("gas_phi", dimensionSet(0, 0, 0, 0, 0), 0.0)
	);

	unsigned int ns = Y.size();
	OpenSMOKE::OpenSMOKEVectorDouble y(ns);
	OpenSMOKE::OpenSMOKEVectorDouble x(ns);
	std::vector<double> moles(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	      scalarField& equivalenceratioCells = equivalenceratio.internalField();

	const std::vector<std::string> fuel_names(2); fuel_names[0] = "C2H2"; fuel_names[1] = "N2";
	const std::vector<std::string> ox_names(2); ox_names[0] = "O2"; ox_names[1] = "N2";
	const std::vector<double> fuel_moles(2); fuel_moles[0] = 0.25; fuel_moles[1] = 0.75;
	const std::vector<double> ox_moles(2); ox_moles[0] = 0.21; ox_moles[1] = 0.79;


	// Stoichiometric composition
	std::vector<double> moles_st(ns);
	moles_st = thermodynamicsMapXML->GetMoleFractionsFromEquivalenceRatio(1.,fuel_names, fuel_moles, ox_names,ox_moles );

	
	
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

		// Vector
		for(unsigned int i=0;i<ns;i++)
			moles[i] = x[i+1];

		// Equivalence ratio
		equivalenceratioCells[celli] = thermodynamicsMapXML->GetLocalEquivalenceRatio(moles, moles_st, fuel_names);
	}


	forAll(T.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT  = T.boundaryField()[patchi];
		const fvPatchScalarField& pp  = p.boundaryField()[patchi];
		      fvPatchScalarField& pequivalenceratio  = equivalenceratio.boundaryField()[patchi];


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

			// Vector
			for(unsigned int i=0;i<ns;i++)
			moles[i] = x[i+1];

			// Equivalence ratio
			pequivalenceratio[facei] = thermodynamicsMapXML->GetLocalEquivalenceRatio(moles, moles_st, fuel_names);;
		}
	}

	equivalenceratio.write();
}
