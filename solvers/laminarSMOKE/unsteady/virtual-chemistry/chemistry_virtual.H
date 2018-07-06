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
	//- Initial conditions
	#if OPENFOAM_VERSION >= 40
	scalarField& TCells = T.ref();
	scalarField& QCells = Q.ref();
	#else
	scalarField& TCells = T.internalField();
	scalarField& QCells = Q.internalField();
	#endif

	const scalarField& rhoCells = rho.internalField();
	const scalarField& cpCells = cp.internalField();
	const scalarField& vCells = mesh.V();
	const scalarField& DeltaTCells = DeltaT.internalField();

	// Select time step
	if (homogeneousReactions == true)
	{
		#if OPENFOAM_VERSION >= 40
		if (LTS)
		{
			const volScalarField& rDeltaT = trDeltaT.ref();
			dimensionedScalar maxIntegrationTime("maxIntegrationTime", dimensionSet(0,0,1,0,0,0,0), scalar(0.01)); 
			DeltaT = min(1.0/rDeltaT, maxIntegrationTime);
		}
		else
		#endif
		{
			DeltaT = mesh.time().deltaT();
		}
	}
	

	if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
		const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies()+1;
		
		// Min and max values
		Eigen::VectorXd yMin(NEQ); for(unsigned int i=0;i<NEQ;i++) yMin(i) = 0.;  yMin(NC) = 200.;
		Eigen::VectorXd yMax(NEQ); for(unsigned int i=0;i<NEQ;i++) yMax(i) = 1.;  yMax(NC) = 6000.;
		Eigen::VectorXd y0(NEQ);
		Eigen::VectorXd yf(NEQ);

		Info <<" * Solving homogeneous virtual chemistry (OpenSMOKE solver)... "<<endl;
		{			
			unsigned int counter = 0;
			
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
				//- Solving for celli:	
				if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry)
				{
					{
						for(unsigned int i=0;i<NC;i++)
							y0(i) = Y[i].internalField()[celli];
						y0(NC) = TCells[celli];

						// Check and normalize the composition
						{
							double sum = 0.;
							for(unsigned int i=0;i<NC;i++)
							{
								if (y0(i) < 0.)	y0(i) = 0.;
								sum += y0(i);
							}
							for(unsigned int i=0;i<NC;i++)
								y0(i) /= sum;
						}

						if (constPressureBatchReactor == true)
						{
							// Set reactor
							batchReactorHomogeneousConstantPressureVirtualChemistry.SetReactor(thermodynamicPressure);
							batchReactorHomogeneousConstantPressureVirtualChemistry.SetEnergyEquation(energyEquation);
						
							// Set initial conditions
							odeSolverConstantPressureVirtualChemistry().SetInitialConditions(t0, y0);

							// Additional ODE solver options
							//if (celli == 0)
							{
								// Set linear algebra options
								odeSolverConstantPressureVirtualChemistry().SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
								odeSolverConstantPressureVirtualChemistry().SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

								// Set relative and absolute tolerances
								odeSolverConstantPressureVirtualChemistry().SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
								odeSolverConstantPressureVirtualChemistry().SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

								// Set minimum and maximum values
								odeSolverConstantPressureVirtualChemistry().SetMinimumValues(yMin);
								odeSolverConstantPressureVirtualChemistry().SetMaximumValues(yMax);
							}
						
							// Solve
							OdeSMOKE::OdeStatus status = odeSolverConstantPressureVirtualChemistry().Solve(t0+DeltaTCells[celli]);
							odeSolverConstantPressureVirtualChemistry().Solution(yf);

							if (status == -6)	// Time step too small
							{
								Info << "Constant pressure reactor (virtual chemistry): " << celli << endl;
								Info << " * T: " << TCells[celli] << endl;
								for(unsigned int i=0;i<NC;i++)
								 	Info << " * " << thermodynamicsMapXML->NamesOfSpecies()[i] << ": " << y0(i) << endl;
							}

							QCells[celli] = batchReactorHomogeneousConstantPressureVirtualChemistry.QR();
						}
						/*
						else
						{
							// Set reactor pressure
							batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], thermodynamicPressure, rhoCells[celli]);
							batchReactorHomogeneousConstantVolume.SetEnergyEquation(energyEquation);
						
							// Set initial conditions
							odeSolverConstantVolume().SetInitialConditions(t0, y0);

							// Additional ODE solver options
							//if (celli == 0)
							{
								// Set linear algebra options
								odeSolverConstantVolume().SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
								odeSolverConstantVolume().SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

								// Set relative and absolute tolerances
								odeSolverConstantVolume().SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
								odeSolverConstantVolume().SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

								// Set minimum and maximum values
								odeSolverConstantVolume().SetMinimumValues(yMin);
								odeSolverConstantVolume().SetMaximumValues(yMax);
							}
						
							// Solve
							OdeSMOKE::OdeStatus status = odeSolverConstantVolume().Solve(t0+DeltaTCells[celli]);
							odeSolverConstantVolume().Solution(yf);

							if (status == -6)	// Time step too small
							{
								Info << "Constant volume reactor: " << celli << endl;
								Info << " * T: " << TCells[celli] << endl;
								for(unsigned int i=0;i<NC;i++)
								 	Info << " * " << thermodynamicsMapXML->NamesOfSpecies()[i] << ": " << y0(i) << endl;
							}

							QCells[celli] = batchReactorHomogeneousConstantVolume.QR();
						}
						*/
					}
				}
				else
				{
					for(unsigned int i=0;i<NC;i++)
						yf(i) = Y[i].internalField()[celli];
					yf(NC) = TCells[celli];
				}

				// Check mass fractions
				normalizeMassFractions(yf, celli, massFractionsTol);

				if (strangAlgorithm != STRANG_COMPACT)
				{
					// Assign mass fractions
					#if OPENFOAM_VERSION >= 40
					for(int i=0;i<NC;i++)
						Y[i].ref()[celli] = yf(i);
					#else
					for(int i=0;i<NC;i++)
						Y[i].internalField()[celli] = yf(i);
					#endif

					//- Allocating final values: temperature
					if (energyEquation == true)
						TCells[celli] = yf(NC);
				}
			/*	else
				{
					const double deltat = tf-t0;

					if (deltat>1e-14)
					{
						thermodynamicsMapXML->SetPressure(thermodynamicPressure);
						thermodynamicsMapXML->SetTemperature(yf(NC));

						double mwmix;
						double cpmix;
						for(int i=1;i<=NC;i++)
							massFractions[i] = yf(i-1);
						thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),mwmix,massFractions.GetHandle());
						cpmix = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/Kmol/K]
						cpmix /= mwmix;
						const double rhomix = thermodynamicPressure*mwmix/PhysicalConstants::R_J_kmol/yf(NC);
				
						// Assign source mass fractions
						#if OPENFOAM_VERSION >= 40
						for(int i=0;i<NC;i++)
							RR[i].ref()[celli] = rhomix*(yf(i)-Y[i].internalField()[celli])/deltat;
						#else
						for(int i=0;i<NC;i++)
							RR[i].internalField()[celli] = rhomix*(yf(i)-Y[i].internalField()[celli])/deltat;
						#endif

						//- Allocating source temperature
						if (energyEquation == true)
							RT[celli] = rhomix*cpmix*(yf(NC)-TCells[celli])/deltat;
					}
					else
					{
						// Assign source mass fractions
						#if OPENFOAM_VERSION >= 40
						for(int i=0;i<NC;i++)
							RR[i].ref()[celli] = 0.;
						#else
						for(int i=0;i<NC;i++)
							RR[i].internalField()[celli] = 0.;
						#endif

						//- Allocating source temperature
						if (energyEquation == true)
							RT[celli] = 0.;
					}
				}
				*/
				if (counter%(int(0.20*mesh.nCells())+1) == 0)
					Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

				counter++;

				// Output
				if (runTime.outputTime())
				{
					if (constPressureBatchReactor == true)
					{
						if (outputFormationRatesIndices.size() != 0)
						{
							#if OPENFOAM_VERSION >= 40
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].ref()[celli] = batchReactorHomogeneousConstantPressureVirtualChemistry.R()[outputFormationRatesIndices[i]+1];
							#else
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantPressureVirtualChemistry.R()[outputFormationRatesIndices[i]+1];
							#endif
						}
					}
					/*
					else
					{
						if (outputFormationRatesIndices.size() != 0)
						{
							#if OPENFOAM_VERSION >= 40
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].ref()[celli] = batchReactorHomogeneousConstantVolume.R()[outputFormationRatesIndices[i]+1] *
                                           	                                                           thermodynamicsMapXML->MW(outputFormationRatesIndices[i]);
							#else
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantVolume.R()[outputFormationRatesIndices[i]+1] *
                                           	                                                           thermodynamicsMapXML->MW(outputFormationRatesIndices[i]);
							#endif

						}
					}
					*/
				}
			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			Info << "   Homogeneous virtual chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
		}
	}
/*
	else if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CHEMEQ2)
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
		const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies();
		Eigen::VectorXd y0(NEQ);
		Eigen::VectorXd yf(NEQ);

		OpenSMOKE::OpenSMOKE_CHEMEQ2<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CHEMEQ2>   odeSolverConstantPressure(chemeq2SolverConstantPressure);
		odeSolverConstantPressure.SetDimensions(NEQ);
		OpenSMOKE::OpenSMOKE_CHEMEQ2<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CHEMEQ2>   odeSolverConstantVolume(chemeq2SolverConstantVolume);
		odeSolverConstantVolume.SetDimensions(NEQ);

		Info << " * Solving homogeneous chemistry (CHEMEQ2 solver)... " << endl;
		{			
			unsigned int counter = 0;
			
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
				//- Solving for celli:	
				if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry)
				{
					{
						for(unsigned int i=0;i<NC;i++)
							y0(i) = Y[i].internalField()[celli];
					
						if (constPressureBatchReactor == true)
						{
							// Set reactor
							batchReactorHomogeneousConstantPressure.SetReactor(thermodynamicPressure);
							batchReactorHomogeneousConstantPressure.SetTemperature(TCells[celli], y0.data());
							batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);
						
							// Set initial conditions
							odeSolverConstantPressure.SetInitialValues(t0, y0.data());

							// Additional ODE solver options
							if (celli == 0)
							{
								odeSolverConstantPressure.SetEpsilonMinimum(chemeq2_epsilonMinimum);
								odeSolverConstantPressure.SetEpsilonMaximum(chemeq2_epsilonMaximum );
								odeSolverConstantPressure.SetMaximumNumberSubIterations( chemeq2_subIterations);
								odeSolverConstantPressure.SetMinimumStep(chemeq2_dtMinimum);
								odeSolverConstantPressure.SetMinimumValue(chemeq2_yMinimum);
							}
						
							// Solve
							odeSolverConstantPressure.Solve(t0+DeltaTCells[celli]);

							// Recover solution
							odeSolverConstantPressure.Solution(yf.data());
							if (energyEquation == true)
								TCells[celli] = batchReactorHomogeneousConstantPressure.GetTemperature();

							QCells[celli] = batchReactorHomogeneousConstantPressure.QR();
						}
						else
						{
							// Set reactor pressure
							batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], thermodynamicPressure, rhoCells[celli]);
							batchReactorHomogeneousConstantVolume.SetTemperature(TCells[celli], y0.data());
							batchReactorHomogeneousConstantVolume.SetEnergyEquation(energyEquation);

							// Set initial conditions
							odeSolverConstantVolume.SetInitialValues(t0, y0.data());

							// Additional ODE solver options
							if (celli == 0)
							{
								odeSolverConstantVolume.SetEpsilonMinimum(chemeq2_epsilonMinimum);
								odeSolverConstantVolume.SetEpsilonMaximum(chemeq2_epsilonMaximum );
								odeSolverConstantVolume.SetMaximumNumberSubIterations( chemeq2_subIterations);
								odeSolverConstantVolume.SetMinimumStep(chemeq2_dtMinimum);
								odeSolverConstantVolume.SetMinimumValue(chemeq2_yMinimum);
							}
						
							// Solve
							odeSolverConstantVolume.Solve(t0+DeltaTCells[celli]);

							// Recover solution
							odeSolverConstantVolume.Solution(yf.data());
							if (energyEquation == true && strangAlgorithm != STRANG_COMPACT)
								TCells[celli] = batchReactorHomogeneousConstantVolume.GetTemperature();
	
							QCells[celli] = batchReactorHomogeneousConstantVolume.QR();
						}
					}
			
					// Check mass fractions
					normalizeMassFractions(yf.data(), celli, massFractionsTol, NC);


					if (strangAlgorithm != STRANG_COMPACT)
					{
						// Assign mass fractions
						#if OPENFOAM_VERSION >= 40
						for(int i=0;i<NC;i++)
							Y[i].ref()[celli] = yf(i);
						#else
						for(int i=0;i<NC;i++)
							Y[i].internalField()[celli] = yf(i);
						#endif
					}
					else
					{
						const double deltat = tf-t0;

						if (deltat>1e-14)
						{
							thermodynamicsMapXML->SetPressure(thermodynamicPressure);
							thermodynamicsMapXML->SetTemperature(yf(NC));

							double mwmix;
							double cpmix;
							for(int i=1;i<=NC;i++)
								massFractions[i] = yf(i-1);
							thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),mwmix,massFractions.GetHandle());
							cpmix = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/Kmol/K]
							cpmix /= mwmix;
							const double rhomix = thermodynamicPressure*mwmix/PhysicalConstants::R_J_kmol/yf(NC);
				
							// Assign source mass fractions
							#if OPENFOAM_VERSION >= 40
							for(int i=0;i<NC;i++)
								RR[i].ref()[celli] = rhomix*(yf(i)-Y[i].internalField()[celli])/deltat;
							#else
							for(int i=0;i<NC;i++)
								RR[i].internalField()[celli] = rhomix*(yf(i)-Y[i].internalField()[celli])/deltat;
							#endif
							
							//- Allocating source temperature
							if (energyEquation == true)
								RT[celli] = rhomix*cpmix*(batchReactorHomogeneousConstantVolume.GetTemperature()-TCells[celli])/deltat;
						}
						else
						{
							// Assign source mass fractions
							#if OPENFOAM_VERSION >= 40
							for(int i=0;i<NC;i++)
								RR[i].ref()[celli] = 0.;
							#else
							for(int i=0;i<NC;i++)
								RR[i].internalField()[celli] = 0.;
							#endif

							//- Allocating source temperature
							if (energyEquation == true)
								RT[celli] = 0.;
						}
					}
				}

				if (counter%(int(0.20*mesh.nCells())+1) == 0)
					Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

				counter++;	
			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
		}
	}
	else if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
		const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies()+1;
		
		// Min and max values
		OpenSMOKE::OpenSMOKEVectorDouble y0(NEQ);
		OpenSMOKE::OpenSMOKEVectorDouble yf(NEQ);

		Info <<" * Solving homogeneous chemistry (external ODE solver)... "<<endl;
		{			
			unsigned int counter = 0;
			
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
				//- Solving for celli:	
				if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry)
				{				
					{
						for(unsigned int i=0;i<NC;i++)
							y0[i+1] = Y[i].internalField()[celli];
						y0[NC+1] = TCells[celli];
					
						if (constPressureBatchReactor == true)
						{
							// Set reactor pressure
							batchReactorHomogeneousConstantPressure.SetReactor(thermodynamicPressure);
							batchReactorHomogeneousConstantPressure.SetCheckMassFractions(true);
							batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);
						
							if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
							{
								#if OPENSMOKE_USE_DVODE == 1
								typedef OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DVODE> Solver;
								SolveOpenSourceSolvers<Solver>(dvodeSolverConstantPressure, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);
								#endif
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODE> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodeSolverConstantPressure, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);					
								#endif
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODA> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodaSolverConstantPressure, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
							{
								#if OPENSMOKE_USE_DASPK == 1
								typedef OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DASPK> Solver;
								SolveOpenSourceSolvers<Solver>(daspkSolverConstantPressure, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);								
								#endif							
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
							{
								#if OPENSMOKE_USE_RADAU == 1
								typedef OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_RADAU5> Solver;
								SolveOpenSourceSolvers<Solver>(radauSolverConstantPressure, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);						
								#endif							
							}			
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
							{
								#if OPENSMOKE_USE_MEBDF == 1
								typedef OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_MEBDF> Solver;
								SolveOpenSourceSolvers<Solver>(mebdfSolverConstantPressure, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}	
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
							{
								#if OPENSMOKE_USE_SUNDIALS == 1
								typedef OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE> Solver;
								SolveOpenSourceSolvers<Solver>(cvodeSolverConstantPressure, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}	

							QCells[celli] = batchReactorHomogeneousConstantPressure.QR();																						
						}
						else
						{
							// Set reactor data
							batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], thermodynamicPressure, rhoCells[celli]);
							batchReactorHomogeneousConstantVolume.SetCheckMassFractions(true);
							batchReactorHomogeneousConstantVolume.SetEnergyEquation(energyEquation);
						
							if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
							{
								#if OPENSMOKE_USE_DVODE == 1
								typedef OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DVODE> Solver;
								SolveOpenSourceSolvers<Solver>(dvodeSolverConstantVolume, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);
								#endif
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODE> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodeSolverConstantVolume, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);					
								#endif
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODA> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodaSolverConstantVolume, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
							{
								#if OPENSMOKE_USE_DASPK == 1
								typedef OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DASPK> Solver;
								SolveOpenSourceSolvers<Solver>(daspkSolverConstantVolume, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);								
								#endif							
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
							{
								#if OPENSMOKE_USE_RADAU == 1
								typedef OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_RADAU5> Solver;
								SolveOpenSourceSolvers<Solver>(radauSolverConstantVolume, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);						
								#endif							
							}			
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
							{
								#if OPENSMOKE_USE_MEBDF == 1
								typedef OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_MEBDF> Solver;
								SolveOpenSourceSolvers<Solver>(mebdfSolverConstantVolume, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}	
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
							{
								#if OPENSMOKE_USE_SUNDIALS == 1
								typedef OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CVODE> Solver;
								SolveOpenSourceSolvers<Solver>(cvodeSolverConstantVolume, t0, t0+DeltaTCells[celli], y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}		

							QCells[celli] = batchReactorHomogeneousConstantVolume.QR();
						}
					}
				
					// Check mass fractions
					normalizeMassFractions(yf, celli, massFractionsTol);

					if (strangAlgorithm != STRANG_COMPACT)
					{
						// Assign mass fractions
						#if OPENFOAM_VERSION >= 40
						for(int i=0;i<NC;i++)
							Y[i].ref()[celli] = yf[i+1];
						#else
						for(int i=0;i<NC;i++)
							Y[i].internalField()[celli] = yf[i+1];
						#endif

						//- Allocating final values: temperature
						if (energyEquation == true)
							TCells[celli] = yf[NC+1];
					}
					else
					{
						const double deltat = tf-t0;

						if (deltat>1e-14)
						{
							thermodynamicsMapXML->SetPressure(thermodynamicPressure);
							thermodynamicsMapXML->SetTemperature(yf[NC+1]);

							double mwmix;
							double cpmix;
							for(int i=1;i<=NC;i++)
								massFractions[i] = yf[i];
							thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),mwmix,massFractions.GetHandle());
							cpmix = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/Kmol/K]
							cpmix /= mwmix;
							const double rhomix = thermodynamicPressure*mwmix/PhysicalConstants::R_J_kmol/yf[NC+1];
				
							// Assign source mass fractions
							#if OPENFOAM_VERSION >= 40
							for(int i=0;i<NC;i++)
								RR[i].ref()[celli] = rhomix*(yf[i+1]-Y[i].internalField()[celli])/deltat;
							#else
							for(int i=0;i<NC;i++)
								RR[i].internalField()[celli] = rhomix*(yf[i+1]-Y[i].internalField()[celli])/deltat;
							#endif
							

							//- Allocating source temperature
							if (energyEquation == true)
								RT[celli] = rhomix*cpmix*(yf[NC+1]-TCells[celli])/deltat;
						}
						else
						{
							// Assign source mass fractions
							#if OPENFOAM_VERSION >= 40
							for(int i=0;i<NC;i++)
								RR[i].ref()[celli] = 0.;
							#else
							for(int i=0;i<NC;i++)
								RR[i].internalField()[celli] = 0.;
							#endif

							//- Allocating source temperature
							if (energyEquation == true)
								RT[celli] = 0.;
						}
					}
				}

				if (counter%(int(0.20*mesh.nCells())+1) == 0)
					Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

				counter++;

				// Output
				if (runTime.outputTime())
				{
					if (constPressureBatchReactor == true)
					{
						if (outputFormationRatesIndices.size() != 0)
						{
							#if OPENFOAM_VERSION >= 40
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].ref()[celli] = batchReactorHomogeneousConstantPressure.R()[outputFormationRatesIndices[i]+1] *
                                       	      		                                                   thermodynamicsMapXML->MW(outputFormationRatesIndices[i]);
							#else
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantPressure.R()[outputFormationRatesIndices[i]+1] *
                                       	      		                                                   thermodynamicsMapXML->MW(outputFormationRatesIndices[i]);
							#endif
						}
					}
					else
					{
						if (outputFormationRatesIndices.size() != 0)
						{
							#if OPENFOAM_VERSION >= 40
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].ref()[celli] = batchReactorHomogeneousConstantVolume.R()[outputFormationRatesIndices[i]+1] *
                                           	                                                           thermodynamicsMapXML->MW(outputFormationRatesIndices[i]);
							#else
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantVolume.R()[outputFormationRatesIndices[i]+1] *
                                           	                                                           thermodynamicsMapXML->MW(outputFormationRatesIndices[i]);
							#endif

						}
					}
				}
			}

			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
		}
	}
    */
    Info<< " * T gas min/max (after chemistry) = " << min(T).value() << ", " << max(T).value() << endl;
}
