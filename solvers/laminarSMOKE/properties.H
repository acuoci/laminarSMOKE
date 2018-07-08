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

if(virtual_chemistry == false)
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
		#else
		scalarField& cTotCells = cTot.internalField();
		scalarField& psiCells =  psi.internalField();
		scalarField& muCells  =  mu.internalField();
		scalarField& lambdaCells = lambda.internalField();
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

					#if OPENFOAM_VERSION >= 40
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                               CpSpecies[i].ref()[celli] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
					#else
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                               CpSpecies[i].internalField()[celli] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
					#endif	
				}
			}

			if (simplifiedTransportProperties == true)
			{
				const double mu0   = 1.8405e-5;	// [kg/m/s]
				const double T0    = 300.;	// [K]
				const double Beta0 = 0.6759;
				const double Pr0   = 0.7;

				muCells[celli] = mu0*std::pow(TCells[celli]/T0, Beta0);
				lambdaCells[celli] = muCells[celli]*cpCells[celli]/Pr0;
			}
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(), moleFractions.GetHandle());

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
					Dmix[physicalSootDiffusivityCorrectionIndex[i]].internalField()[celli] = DmixReference*physicalSootDiffusivityCorrection[i];
				#endif	
			}

			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector.GetHandle(), moleFractions.GetHandle());		
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
					
					#if OPENFOAM_VERSION >= 40
						Dsoret[i].ref()[celli] =Dmix[index-1].internalField()[celli]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/MWmixCells[celli];
					#else
						Dsoret[i].internalField()[celli] =Dmix[index-1].internalField()[celli]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/MWmixCells[celli];
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
			
					#if OPENFOAM_VERSION >= 40
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                                CpSpecies[i].boundaryFieldRef()[patchi][facei] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
					#else
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                                CpSpecies[i].boundaryField()[patchi][facei] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
					#endif
				}
			} 

			if (simplifiedTransportProperties == true)
			{
				const double mu0   = 1.8405e-5;	// [kg/m/s]
				const double T0    = 300.;	// [K]
				const double Beta0 = 0.6759;
				const double Pr0   = 0.7;

				pmu[facei] = mu0*std::pow(pT[facei]/T0, Beta0);	
				plambda[facei] = pmu[facei]*pcp[facei]/Pr0;
			}
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(),moleFractions.GetHandle());

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
				transportMapXML->ThermalDiffusionRatios(tetamixvector.GetHandle(), moleFractions.GetHandle());
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
					
					#if OPENFOAM_VERSION >= 40
					Dsoret[i].boundaryFieldRef()[patchi][facei] = Dmix[index-1].boundaryField()[patchi][facei]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/pMWmix[facei];
					#else
					Dsoret[i].boundaryField()[patchi][facei] = Dmix[index-1].boundaryField()[patchi][facei]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/pMWmix[facei];
					#endif
				}
			}
		}
	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}

else
{
	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info<< "Properties evaluation... " ;

	OpenSMOKE::OpenSMOKEVectorDouble Dmixvector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble CpVector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble hVector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble Omega(thermodynamicsMapXML->NumberOfSpecies());

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
		#else
		scalarField& cTotCells = cTot.internalField();
		scalarField& psiCells =  psi.internalField();
		scalarField& muCells  =  mu.internalField();
		scalarField& lambdaCells = lambda.internalField();
		scalarField& cpCells = cp.internalField();
		scalarField& cvCells = cv.internalField();
		scalarField& hCells = h.internalField();
		scalarField& MWmixCells = MWmix.internalField();
		#endif

		forAll(TCells, celli)
		{
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					massFractions[i+1] = Y[i].internalField()[celli];

			if (virtual_chemistry_table_check == true)
				#include "utilities/virtualchemistry/VirtualChemistryTableCheck.H"
	
			// Concentration [kmol/m3]
			cTotCells[celli] = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];

			// Molecular weight [kg/kmol]
			MWmixCells[celli] = virtualChemistryTable->MWMix(massFractions.GetHandle());

			// Density/Pressure ratio [kg/m3/Pa]
			psiCells[celli]  = cTotCells[celli]*MWmixCells[celli]/pCells[celli];

			// Constant pressure specific heat [J/kg/K]
			cpCells[celli] = virtualChemistryTable->CpMix(TCells[celli],pCells[celli],massFractions.GetHandle());
			cvCells[celli] = cpCells[celli]-PhysicalConstants::R_J_kmol/MWmixCells[celli];

			// Enthalpy [J/kg]
			hCells[celli] = virtualChemistryTable->HMix(TCells[celli],pCells[celli],massFractions.GetHandle());

			// Dynamic viscosity [kg/m/s]
			muCells[celli] = virtualChemistryTable->DynamicViscosity(TCells[celli]);
			
			// Thermal conductivity [W/m/K]
			lambdaCells[celli] = virtualChemistryTable->ThermalConductivity(muCells[celli], cpCells[celli]);

			// Diffusion coefficients [m2/s]
			{
				virtualChemistryTable->MassDiffusionCoefficients(lambdaCells[celli], pCells[celli]*psiCells[celli], cpCells[celli], Dmixvector.GetHandle());

				#if OPENFOAM_VERSION >= 40
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].ref()[celli] = Dmixvector[i+1];
				#else
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].internalField()[celli] = Dmixvector[i+1];
				#endif
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
		#endif
        
		forAll(pT, facei)
		{
			// Mass fractions
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			// Concentration [kmol/m3]
			pcTot[facei] = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);

			// Molecular weight [kg/kmol]
			pMWmix[facei] = virtualChemistryTable->MWMix(massFractions.GetHandle());

			// Density/Pressure ratio [kg/m3/Pa]
			ppsi[facei]  = pcTot[facei]*pMWmix[facei]/pp[facei];

			// Constant pressure specific heat [J/kg/K]
			pcp[facei] = virtualChemistryTable->CpMix(pT[facei],pp[facei],massFractions.GetHandle());
			pcv[facei] = pcp[facei]-PhysicalConstants::R_J_kmol/pMWmix[facei];

			// Enthalpy [J/kg]
			ph[facei] = virtualChemistryTable->HMix(pT[facei],pp[facei],massFractions.GetHandle());

			// Dynamic viscosity [kg/m/s]
			pmu[facei] = virtualChemistryTable->DynamicViscosity(pT[facei]);
			
			// Thermal conductivity [W/m/K]
			plambda[facei] = virtualChemistryTable->ThermalConductivity(pmu[facei], pcp[facei]);

			// Diffusion coefficients [m2/s]
			{
				virtualChemistryTable->MassDiffusionCoefficients(plambda[facei], pp[facei]*ppsi[facei], pcp[facei], Dmixvector.GetHandle());

				#if OPENFOAM_VERSION >= 40
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryFieldRef()[patchi][facei] = Dmixvector[i+1];
				#else
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryField()[patchi][facei] = Dmixvector[i+1];
				#endif	
			}
		}
	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}

#include "extensions/sparc/SPARC_classification.H"
