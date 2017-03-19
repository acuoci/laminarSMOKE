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
	OpenSMOKE::OpenSMOKEVectorDouble Dmixvector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble tetamixvector(thermodynamicsMapXML->NumberOfSpecies());
	Eigen::VectorXd massFractionsEigen(thermodynamicsMapXML->NumberOfSpecies());

	// Internal fields
	{
		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField(); 

		scalarField& cTotCells = cTot.internalField();
		scalarField& psiCells =  psi.internalField();
		scalarField& muCells  =  mu.internalField();
		scalarField& lambdaCells = lambda.internalField();
		scalarField& cpCells = cp.internalField();
		scalarField& cvCells = cv.internalField();
		scalarField& hCells = h.internalField();
		scalarField& MWmixCells = MWmix.internalField();

		forAll(TCells, celli)
		{
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			
			transportMapXML->SetPressure(pCells[celli]);
			transportMapXML->SetTemperature(TCells[celli]);
				
	
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].internalField()[celli];
				
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),MWmixCells[celli],massFractions.GetHandle());

			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].internalField()[celli] = moleFractions[i+1];

			cTotCells[celli] = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			psiCells[celli]  = cTotCells[celli]*MWmixCells[celli]/pCells[celli];
			hCells[celli] = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());		// [J/kmol]
			hCells[celli] /= MWmixCells[celli];															// [J/kg]
			
			//if (viscositySpecies == true)
			{
				muCells[celli] = transportMapXML->DynamicViscosity(moleFractions.GetHandle());
			}

            		if (energyEquation == true || diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
            		{
				lambdaCells[celli] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
               			cpCells[celli] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/kmol/K]
                		cvCells[celli] = (cpCells[celli]-PhysicalConstants::R_J_kmol)/MWmixCells[celli];
				cpCells[celli] = cpCells[celli]/MWmixCells[celli];
			
				if (iMassDiffusionInEnergyEquation == true)
				{
					thermodynamicsMapXML->cpMolar_Species(CpVector.GetHandle());
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                               CpSpecies[i].internalField()[celli] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
				}
			}
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(), moleFractions.GetHandle());
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].internalField()[celli] = Dmixvector[i+1];
			}
			else if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const double coefficient = lambdaCells[celli]/(pCells[celli]*psiCells[celli])/cpCells[celli];
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].internalField()[celli] = coefficient/LewisNumbers(i);
			}

			if (physicalSootDiffusivity == true)
			{
				const double DmixReference = Dmix[physicalSootDiffusivityReferenceIndex].internalField()[celli];
				for(int i=0;i<physicalSootDiffusivityCorrectionIndex.size();i++)
					Dmix[physicalSootDiffusivityCorrectionIndex[i]].internalField()[celli] = DmixReference*physicalSootDiffusivityCorrection[i];
			}

			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector.GetHandle(), moleFractions.GetHandle());		
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
					Dsoret[i].internalField()[celli] =Dmix[index-1].internalField()[celli]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/MWmixCells[celli];
				} 
			}
		}
    }

	// Boundaries
	forAll(T.boundaryField(), patchi)
	{
		fvPatchScalarField& pp = p.boundaryField()[patchi];
		fvPatchScalarField& pT = T.boundaryField()[patchi];
        
		fvPatchScalarField& pcTot = cTot.boundaryField()[patchi];
		fvPatchScalarField& ppsi  = psi.boundaryField()[patchi];
		fvPatchScalarField& pmu   = mu.boundaryField()[patchi];
		fvPatchScalarField& plambda = lambda.boundaryField()[patchi];
		fvPatchScalarField& pcp = cp.boundaryField()[patchi];
		fvPatchScalarField& pcv = cv.boundaryField()[patchi];
		fvPatchScalarField& ph = h.boundaryField()[patchi];
		fvPatchScalarField& pMWmix = MWmix.boundaryField()[patchi];
        
		forAll(pT, facei)
		{
			thermodynamicsMapXML->SetPressure(pp[facei]);
			thermodynamicsMapXML->SetTemperature(pT[facei]);

			transportMapXML->SetPressure(pp[facei]);
			transportMapXML->SetTemperature(pT[facei]);
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),pMWmix[facei],massFractions.GetHandle());

			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].boundaryField()[patchi][facei] = moleFractions[i+1];

			pcTot[facei] = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);
			ppsi[facei]  = pcTot[facei]*pMWmix[facei]/pp[facei];
			ph[facei] = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());		// [J/kmol]
			ph[facei] /= pMWmix[facei];																// [J/kg]
			
			//if (viscositySpecies == true)
			{
				pmu[facei] = transportMapXML->DynamicViscosity(moleFractions.GetHandle());
			}
						
           		if (energyEquation == true || diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
            		{
				plambda[facei] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
				
				pcp[facei] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/Kmol/K]
				pcv[facei] = (pcp[facei]-PhysicalConstants::R_J_kmol)/pMWmix[facei];
				pcp[facei] = pcp[facei]/pMWmix[facei];
			
				if (iMassDiffusionInEnergyEquation == true)
				{
					thermodynamicsMapXML->cpMolar_Species(CpVector.GetHandle());
                                        for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                                CpSpecies[i].boundaryField()[patchi][facei] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
				}
			} 
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(),moleFractions.GetHandle());
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryField()[patchi][facei] = Dmixvector[i+1];
			}
			else if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const double coefficient =  plambda[facei]/(pp[facei]*ppsi[facei])/pcp[facei];
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryField()[patchi][facei] = coefficient/LewisNumbers(i);
			}

			if (physicalSootDiffusivity == true)
			{
				const double DmixReference = Dmix[physicalSootDiffusivityReferenceIndex].boundaryField()[patchi][facei];
				for(int i=0;i<physicalSootDiffusivityCorrectionIndex.size();i++)
					Dmix[physicalSootDiffusivityCorrectionIndex[i]].boundaryField()[patchi][facei] = DmixReference*physicalSootDiffusivityCorrection[i];
			}

			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector.GetHandle(), moleFractions.GetHandle());
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
					Dsoret[i].boundaryField()[patchi][facei] = Dmix[index-1].boundaryField()[patchi][facei]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/pMWmix[facei];
				}
			}
		}
	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
