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

{
	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info<< "Properties evaluation... " ;

        OpenSMOKE::OpenSMOKEVectorDouble CpVector(thermodynamicsMapXML->NumberOfSpecies());
	Eigen::VectorXd massFractionsEigen(thermodynamicsMapXML->NumberOfSpecies());

	// Internal fields
	{
		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField(); 

		#if OPENFOAM_VERSION >= 40
		scalarField& cTotCells = cTot.ref();
		scalarField& psiCells =  psi.ref();
		scalarField& muCells  =  mu.ref();
		scalarField& lambdaCells = lambda.ref();
		scalarField& alphaCells = alpha.ref();
		scalarField& cpCells = cp.ref();
		scalarField& cvCells = cv.ref();
		scalarField& hCells = h.ref();
		scalarField& MWmixCells = MWmix.ref();
		#else
		scalarField& cTotCells = cTot.internalField();
		scalarField& psiCells =  psi.internalField();
		scalarField& muCells  =  mu.internalField();
		scalarField& lambdaCells = lambda.internalField();
		scalarField& alphaCells = alpha.internalField();
		scalarField& cpCells = cp.internalField();
		scalarField& cvCells = cv.internalField();
		scalarField& hCells = h.internalField();
		scalarField& MWmixCells = MWmix.internalField();
		#endif

		forAll(TCells, celli)
		{
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			
			transportMapXML->SetPressure(pCells[celli]);
			transportMapXML->SetTemperature(TCells[celli]);
				
	
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].internalField()[celli];
				
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),MWmixCells[celli],massFractions.GetHandle());

			#if OPENFOAM_VERSION >= 40
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].ref()[celli] = moleFractions[i+1];
			#else
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].internalField()[celli] = moleFractions[i+1];
			#endif

			cTotCells[celli] = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			psiCells[celli]  = cTotCells[celli]*MWmixCells[celli]/pCells[celli];
			hCells[celli] = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());		// [J/kmol]
			hCells[celli] /= MWmixCells[celli];	// [J/kg]

			muCells[celli] = transportMapXML->DynamicViscosity(moleFractions.GetHandle());

            		if (energyEquation == true)
            		{
				lambdaCells[celli] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
               			
				cpCells[celli] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());	//[J/kmol/K]
                		cvCells[celli] = (cpCells[celli]-PhysicalConstants::R_J_kmol)/MWmixCells[celli];
				cpCells[celli] = cpCells[celli]/MWmixCells[celli];

				alphaCells[celli] = lambdaCells[celli]/(pCells[celli]*psiCells[celli])/cpCells[celli];
			
				{
					thermodynamicsMapXML->cpMolar_Species(CpVector.GetHandle());

					#if OPENFOAM_VERSION >= 40
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                               CpSpecies[i].ref()[celli] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
					#else
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                               CpSpecies[i].internalField()[celli] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
					#endif	
				}
			}
		}
    }

	// Boundaries
	forAll(T.boundaryField(), patchi)
	{
		#if OPENFOAM_VERSION >= 40
		fvPatchScalarField& pp = p.boundaryFieldRef()[patchi];
		fvPatchScalarField& pT = T.boundaryFieldRef()[patchi];
		fvPatchScalarField& pcTot = cTot.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppsi  = psi.boundaryFieldRef()[patchi];
		fvPatchScalarField& pmu   = mu.boundaryFieldRef()[patchi];
		fvPatchScalarField& plambda = lambda.boundaryFieldRef()[patchi];
		fvPatchScalarField& palpha = alpha.boundaryFieldRef()[patchi];
		fvPatchScalarField& pcp = cp.boundaryFieldRef()[patchi];
		fvPatchScalarField& pcv = cv.boundaryFieldRef()[patchi];
		fvPatchScalarField& ph = h.boundaryFieldRef()[patchi];
		fvPatchScalarField& pMWmix = MWmix.boundaryFieldRef()[patchi];
		#else
		fvPatchScalarField& pp = p.boundaryField()[patchi];
		fvPatchScalarField& pT = T.boundaryField()[patchi];
		fvPatchScalarField& pcTot = cTot.boundaryField()[patchi];
		fvPatchScalarField& ppsi  = psi.boundaryField()[patchi];
		fvPatchScalarField& pmu   = mu.boundaryField()[patchi];
		fvPatchScalarField& plambda = lambda.boundaryField()[patchi];
		fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
		fvPatchScalarField& pcp = cp.boundaryField()[patchi];
		fvPatchScalarField& pcv = cv.boundaryField()[patchi];
		fvPatchScalarField& ph = h.boundaryField()[patchi];
		fvPatchScalarField& pMWmix = MWmix.boundaryField()[patchi];
		#endif
        
		forAll(pT, facei)
		{
			thermodynamicsMapXML->SetPressure(pp[facei]);
			thermodynamicsMapXML->SetTemperature(pT[facei]);

			transportMapXML->SetPressure(pp[facei]);
			transportMapXML->SetTemperature(pT[facei]);
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),pMWmix[facei],massFractions.GetHandle());

			#if OPENFOAM_VERSION >= 40
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].boundaryFieldRef()[patchi][facei] = moleFractions[i+1];
			#else
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].boundaryField()[patchi][facei] = moleFractions[i+1];
			#endif

			pcTot[facei] = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);
			ppsi[facei]  = pcTot[facei]*pMWmix[facei]/pp[facei];
			ph[facei] = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());		// [J/kmol]
			ph[facei] /= pMWmix[facei];	// [J/kg]
			
			pmu[facei] = transportMapXML->DynamicViscosity(moleFractions.GetHandle());
						
           		if (energyEquation == true)
            		{
				plambda[facei] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
				
				pcp[facei] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle()); //[J/Kmol/K]
				pcv[facei] = (pcp[facei]-PhysicalConstants::R_J_kmol)/pMWmix[facei];
				pcp[facei] = pcp[facei]/pMWmix[facei];

				palpha[facei] =  plambda[facei]/(pp[facei]*ppsi[facei])/pcp[facei];
			
				{
					thermodynamicsMapXML->cpMolar_Species(CpVector.GetHandle());
			
					#if OPENFOAM_VERSION >= 40
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                                CpSpecies[i].boundaryFieldRef()[patchi][facei] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
					#else
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                                CpSpecies[i].boundaryField()[patchi][facei] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
					#endif
				}
			}
		}
	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}

