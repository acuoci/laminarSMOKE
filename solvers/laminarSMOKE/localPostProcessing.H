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

// Denormalize moments
{
	HMOM_N0  = HMOM_Norm_N0*NAV;
	HMOM_M00 = HMOM_Norm_M00*NAV + dimensionedScalar("eps", dimensionSet(0., -3., 0., 0., 0.), 1e-16);
	HMOM_M10 = HMOM_Norm_M10*NAV*V0;
	HMOM_M01 = HMOM_Norm_M01*NAV*S0;

	HMOM_N  = HMOM_M00;
	HMOM_fv = HMOM_M10;
	HMOM_S  = HMOM_M01/HMOM_M00;
	HMOM_V  = HMOM_M10/HMOM_M00;
}

if (iPolimiSoot == true)
{
	// Soot Analysis
	if (sootOnTheFlyFill == true && runTime.outputTime() == true)
	{
		Info<< "On the fly soot post processing... ";
		
		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

		label ns = Y.size();

		OpenSMOKE::OpenSMOKEVectorDouble y(ns);
		OpenSMOKE::OpenSMOKEVectorDouble x(ns);
		Eigen::VectorXd y_eigen(ns);
		Eigen::VectorXd x_eigen(ns);
		OpenSMOKE::OpenSMOKEVectorDouble c(thermodynamicsMapXML->NumberOfSpecies());
		OpenSMOKE::OpenSMOKEVectorDouble R(thermodynamicsMapXML->NumberOfSpecies());

		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();

		scalarField& soot_fv_largeCells = sootFields[0].internalField();
		scalarField& soot_fv_smallCells = sootFields[1].internalField();
		scalarField& soot_rho_largeCells = sootFields[2].internalField();
		scalarField& soot_rho_smallCells = sootFields[3].internalField();
		scalarField& soot_N_largeCells = sootFields[4].internalField();
		scalarField& soot_N_smallCells = sootFields[5].internalField();
		scalarField& soot_x_largeCells = sootFields[6].internalField();
		scalarField& soot_x_smallCells = sootFields[7].internalField();
		scalarField& soot_omega_largeCells = sootFields[8].internalField();
		scalarField& soot_omega_smallCells = sootFields[9].internalField();
		scalarField& pah_omega_1_2Cells = sootFields[10].internalField();
		scalarField& pah_omega_3_4Cells = sootFields[11].internalField();
		scalarField& pah_omega_more_4Cells = sootFields[12].internalField();

		scalarField& R_soot_largeCells = sootFields[13].internalField();
		scalarField& R_soot_smallCells = sootFields[14].internalField();		
		scalarField& R_pah_1_2Cells = sootFields[15].internalField();
		scalarField& R_pah_3_4Cells = sootFields[16].internalField();
		scalarField& R_pah_more_4Cells = sootFields[17].internalField();

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

			// Eigen vectors
			for(unsigned int i=0;i<ns;i++)
			{
				y_eigen(i) = y[i+1];
				x_eigen(i) = x[i+1];
			}

			// Density
			const double rhoGas = pCells[celli]*mw/PhysicalConstants::R_J_kmol/TCells[celli];

			// Soot analysis
			sootAnalyzer->Analysis(TCells[celli], pCells[celli], rhoGas, y_eigen, x_eigen);

			soot_fv_largeCells[celli] = sootAnalyzer->fv_large();
			soot_fv_smallCells[celli] = sootAnalyzer->fv_small();
			soot_rho_largeCells[celli] = sootAnalyzer->rho_large();
			soot_rho_smallCells[celli] = sootAnalyzer->rho_small();
			soot_N_largeCells[celli] = sootAnalyzer->N_large();
			soot_N_smallCells[celli] = sootAnalyzer->N_small();
			soot_omega_largeCells[celli] = sootAnalyzer->omega_large();
			soot_omega_smallCells[celli] = sootAnalyzer->omega_small();
			soot_x_largeCells[celli] = sootAnalyzer->x_large();
			soot_x_smallCells[celli] = sootAnalyzer->x_small();
			pah_omega_1_2Cells[celli] = sootAnalyzer->omega_pah_1_2_rings();
			pah_omega_3_4Cells[celli] = sootAnalyzer->omega_pah_3_4_rings();
			pah_omega_more_4Cells[celli] = sootAnalyzer->omega_pah_more_than_4_rings();


			// Concentrations
			const double cTot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			Product(cTot,x, &c);

			// Kinetics
			kineticsMapXML->SetTemperature(TCells[celli]);
			kineticsMapXML->SetPressure(pCells[celli]);
			kineticsMapXML->KineticConstants();
			kineticsMapXML->ReactionRates(c);
			kineticsMapXML->FormationRates(&R);
			for(unsigned int i=1;i<=ns;i++)
				R[i] *= thermodynamicsMapXML->MW()[i];

			R_soot_largeCells[celli] = largeBinMassFormationRates(sootAnalyzer, R);
			R_soot_smallCells[celli] = smallBinMassFormationRates(sootAnalyzer, R);
			R_pah_1_2Cells[celli]    = pah12RingsFormationRates(sootAnalyzer, R);
			R_pah_3_4Cells[celli]    = pah34RingsFormationRates(sootAnalyzer, R);
			R_pah_more_4Cells[celli] = pahMoreThan4RingsFormationRates(sootAnalyzer, R);
		}

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT  = T.boundaryField()[patchi];
			const fvPatchScalarField& pp  = p.boundaryField()[patchi];
		
			fvPatchScalarField& psoot_fv_large = sootFields[0].boundaryField()[patchi];
			fvPatchScalarField& psoot_fv_small = sootFields[1].boundaryField()[patchi];
			fvPatchScalarField& psoot_rho_large = sootFields[2].boundaryField()[patchi];
			fvPatchScalarField& psoot_rho_small = sootFields[3].boundaryField()[patchi];
			fvPatchScalarField& psoot_N_large = sootFields[4].boundaryField()[patchi];
			fvPatchScalarField& psoot_N_small = sootFields[5].boundaryField()[patchi];
			fvPatchScalarField& psoot_x_large = sootFields[6].boundaryField()[patchi];
			fvPatchScalarField& psoot_x_small = sootFields[7].boundaryField()[patchi];
			fvPatchScalarField& psoot_omega_large = sootFields[8].boundaryField()[patchi];
			fvPatchScalarField& psoot_omega_small = sootFields[9].boundaryField()[patchi];
			fvPatchScalarField& ppah_omega_1_2 = sootFields[10].boundaryField()[patchi];
			fvPatchScalarField& ppah_omega_3_4 = sootFields[11].boundaryField()[patchi];
			fvPatchScalarField& ppah_omega_more_4 = sootFields[12].boundaryField()[patchi];

			fvPatchScalarField& pR_soot_large = sootFields[13].boundaryField()[patchi];
			fvPatchScalarField& pR_soot_small = sootFields[14].boundaryField()[patchi];
			fvPatchScalarField& pR_pah_1_2 = sootFields[15].boundaryField()[patchi];
			fvPatchScalarField& pR_pah_3_4 = sootFields[16].boundaryField()[patchi];
			fvPatchScalarField& pR_pah_more_4 = sootFields[17].boundaryField()[patchi];

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

				// Eigen vectors
				for(unsigned int i=0;i<ns;i++)
				{
					y_eigen(i) = y[i+1];
					x_eigen(i) = x[i+1];
				}

				// Density
				const double rhoGas = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];

				// Soot analysis
				sootAnalyzer->Analysis(pT[facei], pp[facei], rhoGas, y_eigen, x_eigen);	
				psoot_fv_large[facei] = sootAnalyzer->fv_large();
				psoot_fv_small[facei] = sootAnalyzer->fv_small();
				psoot_rho_large[facei] = sootAnalyzer->rho_large();
				psoot_rho_small[facei] = sootAnalyzer->rho_small();
				psoot_N_large[facei] = sootAnalyzer->N_large();
				psoot_N_small[facei] = sootAnalyzer->N_small();
				psoot_omega_large[facei] = sootAnalyzer->omega_large();
				psoot_omega_small[facei] = sootAnalyzer->omega_small();
				psoot_x_large[facei] = sootAnalyzer->x_large();
				psoot_x_small[facei] = sootAnalyzer->x_small();
				ppah_omega_1_2[facei] = sootAnalyzer->omega_pah_1_2_rings();
				ppah_omega_3_4[facei] = sootAnalyzer->omega_pah_3_4_rings();
				ppah_omega_more_4[facei] = sootAnalyzer->omega_pah_more_than_4_rings();

				// Concentrations
				const double cTot = pp[facei]/PhysicalConstants::R_J_kmol/pT[facei];
				Product(cTot,x, &c);

				// Kinetics
				kineticsMapXML->SetTemperature(pT[facei]);
				kineticsMapXML->SetPressure(pp[facei]);
				kineticsMapXML->KineticConstants();
				kineticsMapXML->ReactionRates(c);
				kineticsMapXML->FormationRates(&R);
				for(unsigned int i=1;i<=ns;i++)
					R[i] *= thermodynamicsMapXML->MW()[i];

				pR_soot_large[facei] = largeBinMassFormationRates(sootAnalyzer, R);
				pR_soot_small[facei] = smallBinMassFormationRates(sootAnalyzer, R);
				pR_pah_1_2[facei]    = pah12RingsFormationRates(sootAnalyzer, R);
				pR_pah_3_4[facei]    = pah34RingsFormationRates(sootAnalyzer, R);
				pR_pah_more_4[facei] = pahMoreThan4RingsFormationRates(sootAnalyzer, R);
			}
		}

		// Distribution in a specific point
/*
		if (pnts_soot_psdf.size() > 0)
		{
			Info << "Post processing PSDF..." << endl;

			for (unsigned int i=0;i<pnts_soot_psdf.size();i++)
			{
				const double xx = pnts_soot_psdf[i][0];
				const double yy = pnts_soot_psdf[i][1];
				const double zz = pnts_soot_psdf[i][2];

				Info << "Post processing point " << i+1 << " @ " << xx << " " << yy << " " << zz << endl;

				vector pos(xx,yy,zz);
				label cellI = mesh.findCell(pos);

				if (cellI > -1)
				{
					std::stringstream tag;
					tag << i;
					std::string name_file = "/soot::PSDF_" + tag.str();
					std::string name_path = "PSDFs/" + runTime.timeName();
					// File name
					//fileName outputFile(runTime.path()/runTime.timeName()/name_file.c_str());
					std::string command1 = "mkdir PSDFs";
					std::string command2 = "mkdir " + name_path;
					Info << command1 << " " << name_file << " " << name_path << endl;
					system(command1.c_str());
					system(command2.c_str());
					fileName outputFile(name_path/name_file.c_str());
					std::ofstream fOutput(outputFile.c_str(), std::ios::out);
					fOutput.setf(std::ios::scientific);
					sootAnalyzer->WriteDistributionLabel(fOutput);

					// Temperature and pressure
					autoPtr<interpolation<scalar> > Tinterp = interpolation<scalar>::New("cellPoint", T);
					scalar Tint = Tinterp->interpolate(pos, cellI);
					autoPtr<interpolation<scalar> > pinterp = interpolation<scalar>::New("cellPoint", p);
					scalar pint = pinterp->interpolate(pos, cellI);

					Info << "Cell: " << cellI << " T: " << Tint << " P: " << pint << endl; 

					double mw;

					// Extract the mass fractions
					for(unsigned int i=0;i<ns;i++)
					{
						autoPtr<interpolation<scalar> > Yinterp = interpolation<scalar>::New("cellPoint", Y[i]);
						y[i+1] = Yinterp->interpolate(pos, cellI);
					}

					const double sum = y.SumElements();
					for(unsigned int i=0;i<ns;i++)
						y[i+1] /= sum;

					// Molecular weight
					thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

					// Eigen vectors
					for(unsigned int i=0;i<ns;i++)
					{
						y_eigen(i) = y[i+1];
						x_eigen(i) = x[i+1];
					}

					// Density
					const double rhoGas = pint*mw/PhysicalConstants::R_J_kmol/Tint;

					// Soot analysis
					sootAnalyzer->Analysis(Tint, pint, rhoGas, y_eigen, x_eigen);	
					sootAnalyzer->Distribution();
					sootAnalyzer->WriteDistribution(fOutput, xx, yy, zz, Tint);
					fOutput.close();					
				}
			}
		}
*/
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
	}

	if (sootIntegralCounter == sootIntegralUpdate)
	{
		const double large_bin_mass = largeBinMass(sootAnalyzer, rho, Y);
		const double small_bin_mass = smallBinMass(sootAnalyzer, rho, Y);
		const double large_bin_volume = largeBinVolume(sootAnalyzer, rho, Y);
		const double small_bin_volume = smallBinVolume(sootAnalyzer, rho, Y);
		const double large_bin_particles = largeBinNumberParticles(sootAnalyzer, rho, Y);
		const double small_bin_particles = smallBinNumberParticles(sootAnalyzer, rho, Y);
		const double soot_dimer_mass = sootDimerMass(sootAnalyzer, rho, Y);
		const double pah_1_2_rings_mass = pah_1_2_rings_Mass(sootAnalyzer, rho, Y);
		const double pah_3_4_rings_mass = pah_3_4_rings_Mass(sootAnalyzer, rho, Y);
		const double pah_more_than_4_rings_mass = pah_more_than_4_rings_Mass(sootAnalyzer, rho, Y);

		// Only if the processor is the master, write to the file
		if (Pstream::master())
		{
	    		fSootIntegrals() << std::setw(20) << runTime.value();
			fSootIntegrals() << std::setw(20) << large_bin_mass;
			fSootIntegrals() << std::setw(20) << large_bin_volume;
			fSootIntegrals() << std::setw(20) << large_bin_particles;		
			fSootIntegrals() << std::setw(20) << small_bin_mass;
			fSootIntegrals() << std::setw(20) << small_bin_volume;
			fSootIntegrals() << std::setw(20) << small_bin_particles;
			fSootIntegrals() << std::setw(20) << soot_dimer_mass;		
			fSootIntegrals() << std::setw(20) << pah_1_2_rings_mass;
			fSootIntegrals() << std::setw(20) << pah_3_4_rings_mass;
			fSootIntegrals() << std::setw(20) << pah_more_than_4_rings_mass;
			fSootIntegrals() << std::endl;
		}
		
		Info << "Soot mass:            " << large_bin_mass << endl;
		Info << "Soot volume:          " << large_bin_volume << endl;
		Info << "Soot particles:       " << large_bin_particles << endl;		
		Info << "Precursors mass:      " << small_bin_mass << endl;
		Info << "Precursors volume:    " << small_bin_volume << endl;
		Info << "Precursors particles: " << small_bin_particles << endl;
		Info << "Dimer mass:           " << soot_dimer_mass << endl;		
		Info << "PAH(340) mass:        " << pah_1_2_rings_mass << endl;
		Info << "PAH(400) mass:        " << pah_3_4_rings_mass << endl;
		Info << "PAH(500) mass:        " << pah_more_than_4_rings_mass << endl;

		sootIntegralCounter = 0;
	}

	sootIntegralCounter++;
}

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
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions,dummy,massFractions);

			// Concentrations
			const double cTot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			Product(cTot, moleFractions, &c);

			// Kinetics
			kineticsMapXML->SetTemperature(TCells[celli]);
			kineticsMapXML->SetPressure(pCells[celli]);
			kineticsMapXML->KineticConstants();
			kineticsMapXML->ReactionRates(c);
			kineticsMapXML->FormationRates(&R);

			for (int j=0;j<outputFormationRatesIndices.size();j++)
			{
				const int index = outputFormationRatesIndices(j)+1;
				FormationRates[j].internalField()[celli] = thermodynamicsMapXML->MW()[index]*R[index];
			}
		}


		// Boundaries
		forAll(T.boundaryField(), patchi)
		{
			fvPatchScalarField& pp = p.boundaryField()[patchi];
			fvPatchScalarField& pT = T.boundaryField()[patchi];
        
			forAll(pT, facei)
			{
				// Mole fractions
				thermodynamicsMapXML->SetPressure(pp[facei]);
				thermodynamicsMapXML->SetTemperature(pT[facei]);	
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					massFractions[i+1] = Y[i].boundaryField()[patchi][facei];
				double dummy;
				thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions,dummy,massFractions);
				
				// Concentrations
				const double cTot = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);
				Product(cTot, moleFractions, &c);

				// Kinetics
				kineticsMapXML->SetTemperature(pT[facei]);
				kineticsMapXML->SetPressure(pp[facei]);
				kineticsMapXML->KineticConstants();
				kineticsMapXML->ReactionRates(c);
				kineticsMapXML->FormationRates(&R);

				for (int j=0;j<outputFormationRatesIndices.size();j++)
				{
					const int index = outputFormationRatesIndices(j)+1;
					FormationRates[j].boundaryField()[patchi][facei] = thermodynamicsMapXML->MW()[index]*R[index];
				}
			}
		}
	}
}
