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

if (reconstructMixtureFraction == true)
{
	//- Mixture fraction
	Info<< "Reconstructing mixture fraction\n" << endl;
	volScalarField zMix_rec
	(
		IOobject
		(
			"zMix_rec",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
	    	dimensionedScalar("zMix_rec", dimensionSet(0, 0, 0, 0, 0), 0.0)
	);

	volScalarField phi_rec
	(
		IOobject
		(
			"phi_rec",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
	    	dimensionedScalar("phi_rec", dimensionSet(0, 0, 0, 0, 0), 0.0)
	);

	double r = 1.;

	// Calculations
	{
		label ns = Y.size();
		std::vector<double> y_stl(ns);

		// Stoichiometric composition
		{		
			OpenSMOKE::OpenSMOKEVectorDouble y(ns);
			OpenSMOKE::OpenSMOKEVectorDouble x(ns);
			std::vector<double> mole_fuel(mass_fuel.size());
			std::vector<double> mole_oxidizer(mass_oxidizer.size());

			// Moles fuel side
			{
				y = 0.;
				for(unsigned int j=0;j<fuel_names.size();j++)
					y[thermodynamicsMapXML->IndexOfSpecies(fuel_names[j])] = mass_fuel[j];
				const double sum = y.SumElements();
				for(unsigned int i=0;i<ns;i++)
					y[i+1] /= sum;
		
				// Molecular weight
				double mw;
				thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), mw, y.GetHandle());

				for(unsigned int j=0;j<fuel_names.size();j++)
					mole_fuel[j] = x[thermodynamicsMapXML->IndexOfSpecies(fuel_names[j])];
			}

			// Moles oxidizer side
			{
				y = 0.;
				for(unsigned int j=0;j<oxidizer_names.size();j++)
					y[thermodynamicsMapXML->IndexOfSpecies(oxidizer_names[j])] = mass_oxidizer[j];
				const double sum = y.SumElements();
				for(unsigned int i=0;i<ns;i++)
					y[i+1] /= sum;
		
				// Molecular weight
				double mw;
				thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), mw, y.GetHandle());

				for(unsigned int j=0;j<oxidizer_names.size();j++)
					mole_oxidizer[j] = x[thermodynamicsMapXML->IndexOfSpecies(oxidizer_names[j])];
			}

			// Stoichiometric composition
			std::vector<double> moles_st(ns);
			std::vector<double> masses_st(ns);
			std::fill(moles_st.begin(), moles_st.end(), 0.);
			std::fill(masses_st.begin(), masses_st.end(), 0.);

			moles_st = thermodynamicsMapXML->GetMoleFractionsFromEquivalenceRatio( 1., fuel_names, mole_fuel, oxidizer_names, mole_oxidizer );
			
			{
				x = 0.;
				for(unsigned int j=0;j<ns;j++)
					x[j+1] = moles_st[j];
		
				double mw;
				thermodynamicsMapXML->MassFractions_From_MoleFractions(y.GetHandle(), mw, x.GetHandle());

				for(unsigned int j=0;j<ns;j++)
					masses_st[j] = y[j+1];
			}

			const double zst = thermodynamicsMapXML->GetMixtureFractionFromMassFractions(masses_st, fuel_names, mass_fuel, oxidizer_names, mass_oxidizer);
			r = (1.-zst)/zst;

			std::cout << " * Stoichiometric z:         " << zst << std::endl;
			std::cout << " * Air-to-fuel ratio (mass): " << r << std::endl;
		}

		
		// Calculating the mixture fraction
	

		const scalarField& TCells = T.internalField();
		      scalarField& zMix_recCells = zMix_rec.internalField();
		      scalarField& phi_recCells = phi_rec.internalField();

		forAll(TCells, celli)
		{	
			// Extract the mean mass fractions
			for(unsigned int i=0;i<ns;i++)
				y_stl[i] = Y[i].internalField()[celli];

			// Reconstruct mixture fraction
			zMix_recCells[celli] = std::min ( thermodynamicsMapXML->GetMixtureFractionFromMassFractions(y_stl, fuel_names, mass_fuel,oxidizer_names, mass_oxidizer), 1.);
			phi_recCells[celli]  = std::min ( zMix_recCells[celli]*r / ( (1.-zMix_recCells[celli]) + 1.e-16), 100.);
		}

		forAll(T.boundaryField(), patchi)
		{
			fvPatchScalarField& pT = T.boundaryField()[patchi];
			fvPatchScalarField& pzMix_rec = zMix_rec.boundaryField()[patchi];
			fvPatchScalarField& pphi_rec = phi_rec.boundaryField()[patchi];
	
			forAll(pT, facei)
			{
				// Extract the mass fractions
				for(unsigned int i=0;i<ns;i++)
					y_stl[i] = Y[i].boundaryField()[patchi][facei];

				// Reconstruct mixture fraction
				pzMix_rec[facei] = std::min ( thermodynamicsMapXML->GetMixtureFractionFromMassFractions(y_stl, fuel_names, mass_fuel,oxidizer_names, mass_oxidizer), 1.);
				pphi_rec[facei] = std::min ( pzMix_rec[facei]*r / ( (1.-pzMix_rec[facei]) + 1.e-16), 100.);
			}
		}

		zMix_rec.write();
		phi_rec.write();			
	}
}
