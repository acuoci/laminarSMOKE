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
		
		#if OPENFOAM_VERSION >= 40
		scalarField& cTotCells = cTot.ref();
		scalarField& psiCells =  psi.ref();
		scalarField& muCells  =  mu.ref();
		scalarField& lambdaCells = lambda.ref();
		scalarField& cpCells = cp.ref();
		scalarField& cvCells = cv.ref();
		scalarField& hCells = h.ref();
		scalarField& MWmixCells = MWmix.ref();

		scalarField& HMOM_Norm_N0Cells  = HMOM_Norm_N0.ref();
		scalarField& HMOM_Norm_M00Cells = HMOM_Norm_M00.ref();
		scalarField& HMOM_Norm_M01Cells = HMOM_Norm_M01.ref();
		scalarField& HMOM_Norm_M10Cells = HMOM_Norm_M10.ref();

		scalarField& HMOM_fvCells = HMOM_fv.ref(); 
		scalarField& HMOM_NCells  = HMOM_N.ref(); 
		scalarField& HMOM_dpCells = HMOM_dp.ref();
		scalarField& HMOM_dcCells = HMOM_dc.ref(); 
		scalarField& HMOM_npCells = HMOM_np.ref(); 

		#else
		scalarField& cTotCells = cTot.internalField();
		scalarField& psiCells =  psi.internalField();
		scalarField& muCells  =  mu.internalField();
		scalarField& lambdaCells = lambda.internalField();
		scalarField& cpCells = cp.internalField();
		scalarField& cvCells = cv.internalField();
		scalarField& hCells = h.internalField();
		scalarField& MWmixCells = MWmix.internalField();

		scalarField& HMOM_Norm_N0Cells  = HMOM_Norm_N0.internalField();
		scalarField& HMOM_Norm_M00Cells = HMOM_Norm_M00.internalField();
		scalarField& HMOM_Norm_M01Cells = HMOM_Norm_M01.internalField();
		scalarField& HMOM_Norm_M10Cells = HMOM_Norm_M10.internalField();

		scalarField& HMOM_fvCells = HMOM_fv.internalField(); 
		scalarField& HMOM_NCells  = HMOM_N.internalField(); 
		scalarField& HMOM_dpCells = HMOM_dp.internalField(); 
		scalarField& HMOM_dcCells = HMOM_dc.internalField(); 
		scalarField& HMOM_npCells = HMOM_np.internalField(); 
		#endif


		


		forAll(TCells, celli)
		{
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			
			transportMapXML->SetPressure(pCells[celli]);
			transportMapXML->SetTemperature(TCells[celli]);
				
	
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].internalField()[celli];
				
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions,MWmixCells[celli],massFractions);

			#if OPENFOAM_VERSION >= 40
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].ref()[celli] = moleFractions[i+1];
			#else
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].internalField()[celli] = moleFractions[i+1];
			#endif

			cTotCells[celli] = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			psiCells[celli]  = cTotCells[celli]*MWmixCells[celli]/pCells[celli];
			thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(hCells[celli], moleFractions);		// [J/kmol]
			hCells[celli] /= MWmixCells[celli];															// [J/kg]
			
			//if (viscositySpecies == true)
			{
				transportMapXML->DynamicViscosity(muCells[celli],moleFractions);
			}

            		if (energyEquation == true || diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
            		{
				transportMapXML->ThermalConductivity(lambdaCells[celli],moleFractions);
               			thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(cpCells[celli], moleFractions);			//[J/kmol/K]
                		cvCells[celli] = (cpCells[celli]-PhysicalConstants::R_J_kmol)/MWmixCells[celli];
				cpCells[celli] = cpCells[celli]/MWmixCells[celli];
			
				if (iMassDiffusionInEnergyEquation == true)
				{
					thermodynamicsMapXML->cpMolar_Species(CpVector);

					#if OPENFOAM_VERSION >= 40
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                               CpSpecies[i].ref()[celli] = CpVector[i+1] / thermodynamicsMapXML->MW()[i+1];
					#else
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                               CpSpecies[i].internalField()[celli] = CpVector[i+1] / thermodynamicsMapXML->MW()[i+1];
					#endif	
				}
			}
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector, moleFractions);

				#if OPENFOAM_VERSION >= 40
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].ref()[celli] = Dmixvector[i+1];
				#else
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].internalField()[celli] = Dmixvector[i+1];
				#endif	
			}
			else if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const double coefficient = lambdaCells[celli]/(pCells[celli]*psiCells[celli])/cpCells[celli];

				#if OPENFOAM_VERSION >= 40
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].ref()[celli] = coefficient/LewisNumbers(i);
				#else
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].internalField()[celli] = coefficient/LewisNumbers(i);
				#endif	
			}

			if (physicalSootDiffusivity == true)
			{
				const double DmixReference = Dmix[physicalSootDiffusivityReferenceIndex].internalField()[celli];

				#if OPENFOAM_VERSION >= 40
				for(int i=0;i<physicalSootDiffusivityCorrectionIndex.size();i++)
					Dmix[physicalSootDiffusivityCorrectionIndex[i]].ref()[celli] = DmixReference*physicalSootDiffusivityCorrection[i];
				#else
				for(int i=0;i<physicalSootDiffusivityCorrectionIndex.size();i++)
					Dmix[physicalSootDiffusivityCorrectionIndex[i]].internalField()[celli] = 										DmixReference*physicalSootDiffusivityCorrection[i];
				#endif	
			}

			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector, moleFractions);		
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];

					#if OPENFOAM_VERSION >= 40
					Dsoret[i].ref()[celli] =Dmix[index-1].internalField()[celli]*tetamixvector[index]*thermodynamicsMapXML->MW()[index]/MWmixCells[celli];
					#else
					Dsoret[i].internalField()[celli] =Dmix[index-1].internalField()[celli]*tetamixvector[index]*thermodynamicsMapXML->MW()[index]/MWmixCells[celli];
					#endif
				} 
			}

			// HMOM
			if (hmom_analysis == true)
			{
				// Setup
				hmom.SetNormalizedMoments(HMOM_Norm_M00Cells[celli], HMOM_Norm_M10Cells[celli], HMOM_Norm_M01Cells[celli], HMOM_Norm_N0Cells[celli]);

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
		fvPatchScalarField& pcp = cp.boundaryFieldRef()[patchi];
		fvPatchScalarField& pcv = cv.boundaryFieldRef()[patchi];
		fvPatchScalarField& ph = h.boundaryFieldRef()[patchi];
		fvPatchScalarField& pMWmix = MWmix.boundaryFieldRef()[patchi];

		fvPatchScalarField& pHMOM_Norm_N0  = HMOM_Norm_N0.boundaryFieldRef()[patchi];
		fvPatchScalarField& pHMOM_Norm_M00 = HMOM_Norm_M00.boundaryFieldRef()[patchi];
		fvPatchScalarField& pHMOM_Norm_M01 = HMOM_Norm_M01.boundaryFieldRef()[patchi];
		fvPatchScalarField& pHMOM_Norm_M10 = HMOM_Norm_M10.boundaryFieldRef()[patchi];

		fvPatchScalarField& pHMOM_fv = HMOM_fv.boundaryFieldRef()[patchi];
		fvPatchScalarField& pHMOM_N  = HMOM_N.boundaryFieldRef()[patchi];
		fvPatchScalarField& pHMOM_dp = HMOM_dp.boundaryFieldRef()[patchi];
		fvPatchScalarField& pHMOM_dc = HMOM_dc.boundaryFieldRef()[patchi];
		fvPatchScalarField& pHMOM_np = HMOM_np.boundaryFieldRef()[patchi];
		#else
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

		fvPatchScalarField& pHMOM_Norm_N0  = HMOM_Norm_N0.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_Norm_M00 = HMOM_Norm_M00.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_Norm_M01 = HMOM_Norm_M01.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_Norm_M10 = HMOM_Norm_M10.boundaryField()[patchi];

		fvPatchScalarField& pHMOM_fv = HMOM_fv.boundaryField()[patchi]; 
		fvPatchScalarField& pHMOM_N  = HMOM_N.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_dp = HMOM_dp.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_dc = HMOM_dc.boundaryField()[patchi];
		fvPatchScalarField& pHMOM_np = HMOM_np.boundaryField()[patchi];
		#endif



        
		forAll(pT, facei)
		{
			thermodynamicsMapXML->SetPressure(pp[facei]);
			thermodynamicsMapXML->SetTemperature(pT[facei]);

			transportMapXML->SetPressure(pp[facei]);
			transportMapXML->SetTemperature(pT[facei]);
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions,pMWmix[facei],massFractions);

			#if OPENFOAM_VERSION >= 40
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].boundaryFieldRef()[patchi][facei] = moleFractions[i+1];
			#else
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].boundaryField()[patchi][facei] = moleFractions[i+1];
			#endif

			pcTot[facei] = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);
			ppsi[facei]  = pcTot[facei]*pMWmix[facei]/pp[facei];
			thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(ph[facei], moleFractions);	// [J/kmol]
			ph[facei] /= pMWmix[facei];								// [J/kg]
			
			//if (viscositySpecies == true)
			{
				transportMapXML->DynamicViscosity(pmu[facei],moleFractions);
			}
						
            		if (energyEquation == true || diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
            		{
				transportMapXML->ThermalConductivity(plambda[facei],moleFractions);
				
				thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(pcp[facei], moleFractions);			//[J/Kmol/K]
				pcv[facei] = (pcp[facei]-PhysicalConstants::R_J_kmol)/pMWmix[facei];
				pcp[facei] = pcp[facei]/pMWmix[facei];
			
				if (iMassDiffusionInEnergyEquation == true)
				{
					thermodynamicsMapXML->cpMolar_Species(CpVector);

					#if OPENFOAM_VERSION >= 40
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                                CpSpecies[i].boundaryFieldRef()[patchi][facei] = CpVector[i+1] / thermodynamicsMapXML->MW()[i+1];
					#else
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                                CpSpecies[i].boundaryField()[patchi][facei] = CpVector[i+1] / thermodynamicsMapXML->MW()[i+1];
					#endif
				}
			} 
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector,moleFractions);

				#if OPENFOAM_VERSION >= 40
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryFieldRef()[patchi][facei] = Dmixvector[i+1];
				#else
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryField()[patchi][facei] = Dmixvector[i+1];
				#endif
			}
			else if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const double coefficient =  plambda[facei]/(pp[facei]*ppsi[facei])/pcp[facei];

				#if OPENFOAM_VERSION >= 40
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryFieldRef()[patchi][facei] = coefficient/LewisNumbers(i);
				#else
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryField()[patchi][facei] = coefficient/LewisNumbers(i);
				#endif
			}

			if (physicalSootDiffusivity == true)
			{
				const double DmixReference = Dmix[physicalSootDiffusivityReferenceIndex].boundaryField()[patchi][facei];

				#if OPENFOAM_VERSION >= 40
				for(int i=0;i<physicalSootDiffusivityCorrectionIndex.size();i++)
					Dmix[physicalSootDiffusivityCorrectionIndex[i]].boundaryFieldRef()[patchi][facei] = DmixReference*physicalSootDiffusivityCorrection[i];
				#else
				for(int i=0;i<physicalSootDiffusivityCorrectionIndex.size();i++)
					Dmix[physicalSootDiffusivityCorrectionIndex[i]].boundaryField()[patchi][facei] = DmixReference*physicalSootDiffusivityCorrection[i];
				#endif
			}

			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector, moleFractions);
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];

					#if OPENFOAM_VERSION >= 40
					Dsoret[i].boundaryFieldRef()[patchi][facei] = Dmix[index-1].boundaryField()[patchi][facei]*tetamixvector[index]*thermodynamicsMapXML->MW()[index]/pMWmix[facei];
					#else
					Dsoret[i].boundaryField()[patchi][facei] = Dmix[index-1].boundaryField()[patchi][facei]*tetamixvector[index]*thermodynamicsMapXML->MW()[index]/pMWmix[facei];
					#endif
				}
			}

			// HMOM
			if (hmom_analysis == true)
			{
				// Setup
				hmom.SetNormalizedMoments(pHMOM_Norm_M00[facei], pHMOM_Norm_M10[facei], pHMOM_Norm_M01[facei], pHMOM_Norm_N0[facei]);

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
			}
		}
	}

	// This is needed by the radiation model
	soot_fv_large = HMOM_fv;
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
