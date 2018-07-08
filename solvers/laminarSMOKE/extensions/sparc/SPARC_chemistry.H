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
|   Authors: A. Cuoci, G. D'Alessio, A. Parente                           |
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
|   Copyright(C) 2018 A. Cuoci                                            |
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
	scalarField& cpuChemistryCells = cpuChemistry.ref();
	#else
	scalarField& TCells = T.internalField();
	scalarField& QCells = Q.internalField();
        scalarField& cpuChemistryCells = cpuChemistry.internalField();
	#endif

	const scalarField& clusterIndexCells = clusterIndex.ref();
	const scalarField& rhoCells = rho.internalField();
	const scalarField& cpCells = cp.internalField();
	const scalarField& mwCells = MWmix.internalField();
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
		unsigned int NCtot = thermodynamicsMapXML->NumberOfSpecies();
		Eigen::VectorXd yf(NCtot+1);

		Info <<" * Solving homogeneous chemistry with SPARC (OpenSMOKE++ solver)... "<<endl;
		{			
			unsigned int counter = 0;

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
				const unsigned int iCluster = clusterIndexCells[celli];
				const unsigned int NClocal  = vector_thermodynamicsMapXML[iCluster]->NumberOfSpecies();
				const unsigned int NRlocal  = vector_kineticsMapXML[iCluster]->NumberOfReactions(); 
				const unsigned int NEQ = NClocal+1;

				const double tStartLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

				//- Solving for celli:	
				if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry)
				{
					// Min and max values
					Eigen::VectorXd yMin(NEQ); for(unsigned int i=0;i<NEQ;i++) yMin(i) = 0.;  yMin(NClocal) = 200.;
					Eigen::VectorXd yMax(NEQ); for(unsigned int i=0;i<NEQ;i++) yMax(i) = 1.;  yMax(NClocal) = 6000.;
					Eigen::VectorXd y0local(NEQ);
					Eigen::VectorXd yflocal(NEQ);

					for(unsigned int i=0;i<NClocal;i++)
						y0local(i) = Y[classifier_indices[iCluster](i)].internalField()[celli];
					y0local(NClocal) = TCells[celli];

					for(unsigned int i=0;i<NCtot;i++)
						yf(i) = Y[i].internalField()[celli];
					yf(NCtot) = TCells[celli];

					if (constPressureBatchReactor == true)
					{
						// Set reactor
						vector_batchReactorHomogeneousConstantPressure[iCluster].SetReactor(thermodynamicPressure, cpCells[celli], mwCells[celli]);
						vector_batchReactorHomogeneousConstantPressure[iCluster].SetEnergyEquation(energyEquation);
					
						// Set initial conditions
						vector_odeSolverConstantPressure[iCluster].SetInitialConditions(t0, y0local);

						// Additional ODE solver options
						//if (celli == 0)
						{
							// Set linear algebra options
							vector_odeSolverConstantPressure[iCluster].SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
							vector_odeSolverConstantPressure[iCluster].SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

							// Set relative and absolute tolerances
							vector_odeSolverConstantPressure[iCluster].SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
							vector_odeSolverConstantPressure[iCluster].SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

							// Set minimum and maximum values
							vector_odeSolverConstantPressure[iCluster].SetMinimumValues(yMin);
							vector_odeSolverConstantPressure[iCluster].SetMaximumValues(yMax);
						}
					
						// Solve
						OdeSMOKE::OdeStatus status = vector_odeSolverConstantPressure[iCluster].Solve(t0+DeltaTCells[celli]);
						vector_odeSolverConstantPressure[iCluster].Solution(yflocal);

						if (status == -6)	// Time step too small
						{
							Info << "Constant pressure reactor: " << celli << endl;
							Info << " * T: " << TCells[celli] << endl;
							for(unsigned int i=0;i<NClocal;i++)
							 	Info << " * " << vector_thermodynamicsMapXML[iCluster]->NamesOfSpecies()[i] << ": " << y0local(i) << endl;
						}

						QCells[celli] = vector_batchReactorHomogeneousConstantPressure[iCluster].QR();
					}
					else
					{
						// Set reactor pressure
						vector_batchReactorHomogeneousConstantVolume[iCluster].SetReactor(vCells[celli], thermodynamicPressure, rhoCells[celli]);
						vector_batchReactorHomogeneousConstantVolume[iCluster].SetEnergyEquation(energyEquation);
					
						// Set initial conditions
						vector_odeSolverConstantVolume[iCluster].SetInitialConditions(t0, y0local);

						// Additional ODE solver options
						//if (celli == 0)
						{
							// Set linear algebra options
							vector_odeSolverConstantVolume[iCluster].SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
							vector_odeSolverConstantVolume[iCluster].SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

							// Set relative and absolute tolerances
							vector_odeSolverConstantVolume[iCluster].SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
							vector_odeSolverConstantVolume[iCluster].SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

							// Set minimum and maximum values
							vector_odeSolverConstantVolume[iCluster].SetMinimumValues(yMin);
							vector_odeSolverConstantVolume[iCluster].SetMaximumValues(yMax);
						}
					
						// Solve
						OdeSMOKE::OdeStatus status = vector_odeSolverConstantVolume[iCluster].Solve(t0+DeltaTCells[celli]);
						vector_odeSolverConstantVolume[iCluster].Solution(yflocal);

						if (status == -6)	// Time step too small
						{
							Info << "Constant volume reactor: " << celli << endl;
							Info << " * T: " << TCells[celli] << endl;
							for(unsigned int i=0;i<NClocal;i++)
							 	Info << " * " << vector_thermodynamicsMapXML[iCluster]->NamesOfSpecies()[i] << ": " << y0local(i) << endl;
						}

						QCells[celli] = vector_batchReactorHomogeneousConstantVolume[iCluster].QR();
					}
					
					for(unsigned int i=0;i<NClocal;i++)
						yf(classifier_indices[iCluster](i)) = yflocal(i);
					yf(NCtot) = yflocal(NClocal);
					
				}
				else
				{
					for(unsigned int i=0;i<NCtot;i++)
						yf(i) = Y[i].internalField()[celli];
					yf(NCtot) = TCells[celli];
				}

				// Check mass fractions
				normalizeMassFractions(yf, celli, massFractionsTol, vc_main_species);

				if (strangAlgorithm != STRANG_COMPACT)
				{
					// Assign mass fractions
					#if OPENFOAM_VERSION >= 40
					for(int i=0;i<NCtot;i++)
						Y[i].ref()[celli] = yf(i);
					#else
					for(int i=0;i<NCtot;i++)
						Y[i].internalField()[celli] = yf(i);
					#endif

					//- Allocating final values: temperature
					if (energyEquation == true)
						TCells[celli] = yf(NCtot);
				}
				/*
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

				const double tEndLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuChemistryCells[celli] = (tEndLocal-tStartLocal)*1000.;

				
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
	
    Info<< " * T gas min/max (after chemistry) = " << min(T).value() << ", " << max(T).value() << endl;
}
