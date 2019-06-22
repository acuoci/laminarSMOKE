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

// Radiative loss
{
	if (runTime.outputTime() == true)
	{
		radiation->Qloss(T,Qrad);
	}
}

// Formation rates
{
	if (runTime.outputTime() == true && outputFormationRatesIndices.size()>0)
	{
		OpenSMOKE::OpenSMOKEVectorDouble c(thermodynamicsMapXML->NumberOfSpecies());
		OpenSMOKE::OpenSMOKEVectorDouble R(thermodynamicsMapXML->NumberOfSpecies());

		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField(); 

		forAll(TCells, celli)
		{
			// Mole fractions
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].internalField()[celli];
			double dummy;
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),dummy,massFractions.GetHandle());

			// Concentrations
			const double cTot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			Product(cTot, moleFractions, &c);

			// Kinetics
			kineticsMapXML->SetTemperature(TCells[celli]);
			kineticsMapXML->SetPressure(pCells[celli]);
			kineticsMapXML->KineticConstants();
			kineticsMapXML->ReactionRates(c.GetHandle());
			kineticsMapXML->FormationRates(R.GetHandle());

			for (int j=0;j<outputFormationRatesIndices.size();j++)
			{
				const int index = outputFormationRatesIndices(j)+1;
				#if OPENFOAM_VERSION >= 40
					FormationRates[j].ref()[celli] = thermodynamicsMapXML->MW(index-1)*R[index];
				#else
					FormationRates[j].internalField()[celli] = thermodynamicsMapXML->MW(index-1)*R[index];
				#endif
			}
		}


		// Boundaries
		forAll(T.boundaryField(), patchi)
		{
			#if OPENFOAM_VERSION >= 40
			fvPatchScalarField& pp = p.boundaryFieldRef()[patchi];
			fvPatchScalarField& pT = T.boundaryFieldRef()[patchi];
			#else
			fvPatchScalarField& pp = p.boundaryField()[patchi];
			fvPatchScalarField& pT = T.boundaryField()[patchi];
			#endif
        
			forAll(pT, facei)
			{
				// Mole fractions
				thermodynamicsMapXML->SetPressure(pp[facei]);
				thermodynamicsMapXML->SetTemperature(pT[facei]);	
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					massFractions[i+1] = Y[i].boundaryField()[patchi][facei];
				double dummy;
				thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),dummy,massFractions.GetHandle());
				
				// Concentrations
				const double cTot = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);
				Product(cTot, moleFractions, &c);

				// Kinetics
				kineticsMapXML->SetTemperature(pT[facei]);
				kineticsMapXML->SetPressure(pp[facei]);
				kineticsMapXML->KineticConstants();
				kineticsMapXML->ReactionRates(c.GetHandle());
				kineticsMapXML->FormationRates(R.GetHandle());

				for (int j=0;j<outputFormationRatesIndices.size();j++)
				{
					const int index = outputFormationRatesIndices(j)+1;

					#if OPENFOAM_VERSION >= 40
					FormationRates[j].boundaryFieldRef()[patchi][facei] = thermodynamicsMapXML->MW(index-1)*R[index];
					#else
					FormationRates[j].boundaryField()[patchi][facei] = thermodynamicsMapXML->MW(index-1)*R[index];
					#endif
				}
			}
		}
	}
}
