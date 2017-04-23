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
	const scalarField& vCells = mesh.V();

	if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
		const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies()+2;
		
		// Min and max values
		Eigen::VectorXd yMin(NEQ); for(unsigned int i=0;i<NEQ;i++) yMin(i) = 0.;  yMin(NC) = 200.;	yMin(NC+1) = 0.;
		Eigen::VectorXd yMax(NEQ); for(unsigned int i=0;i<NEQ;i++) yMax(i) = 1.;  yMax(NC) = 6000.;	yMax(NC+1) = 1e16;
		Eigen::VectorXd y0(NEQ);
		Eigen::VectorXd yf(NEQ);

		Info <<" * Solving homogeneous chemistry (OpenSMOKE solver + ISAT)... "<<endl;
		{			
			unsigned int counter = 0;
			unsigned int nAddHOM = 0;
			unsigned int nGrowHOM = 0;
			unsigned int nRetHOM = 0;

			double cpuTimeRet      = 0.;
			double cpuTimeDI       = 0.;
			double cpuTimeGrowth   = 0.;
			double cpuTimeAddition = 0.;

			
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
				//- Solving for celli:	
				{
					for(unsigned int i=0;i<NC;i++)
						y0(i) = Y[i].internalField()[celli];
					y0(NC)   = TCells[celli];
					y0(NC+1) = tf-t0;
					
					if (constPressureBatchReactor == true)
					{
						// ISAT Algorithm
						{
							for(unsigned int i=0;i<NEQ;i++)
								phiISAT_HOM(i) = y0(i)*scalingFactors_ISAT(i);

							chemComp *phi0base = NULL;
							if(isat_HOM->retrieve(phiISAT_HOM, phi0base)) 
							{					
								double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
								
								// makes interpolation
								isat_HOM->interpol(phiISAT_HOM, RphiISAT_HOM, phi0base);
					
								//check negative value
								for(unsigned int i=0;i<NEQ;i++)
									yf(i) = std::max(RphiISAT_HOM(i), 0.)/scalingFactors_ISAT(i);

								nRetHOM++;	
								
								double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
								
								cpuTimeRet += (t2-t1);
							} 
							else 
							{			
								// Direct integration
								{
									double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
		
									// Set reactor
									batchReactorHomogeneousConstantPressure.SetReactor(thermodynamicPressure);
									batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);
									batchReactorHomogeneousConstantPressure.SetISAT(true);
						
									// Set initial conditions
									odeSolverConstantPressure().SetInitialConditions(t0, y0);

									// Additional ODE solver options
									if (nGrowHOM == 0 && nAddHOM == 0)
									{
										// Set linear algebra options
										odeSolverConstantPressure().SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
										odeSolverConstantPressure().SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

										// Set relative and absolute tolerances
										odeSolverConstantPressure().SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
										odeSolverConstantPressure().SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

										// Set minimum and maximum values
										odeSolverConstantPressure().SetMinimumValues(yMin);
										odeSolverConstantPressure().SetMaximumValues(yMax);
									}

									OdeSMOKE::OdeStatus status = odeSolverConstantPressure().Solve(tf);
									odeSolverConstantPressure().Solution(yf);

									// Move the solution from DI to ISAT
									for(unsigned int i=0;i<NEQ;i++)
										RphiISAT_HOM(i) = std::max(yf(i), 0.)*scalingFactors_ISAT(i);

									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

									cpuTimeDI += (t2-t1);
								}

								// Growth
								if(isat_HOM->grow(phiISAT_HOM, RphiISAT_HOM, phi0base)) 
								{
									double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();

									phi0base->growEOA(phiISAT_HOM);
									nGrowHOM++;

									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

									cpuTimeGrowth += (t2-t1);
								} 
								// Addition
								else
								{
									double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
									
									// compute mapping gradient
									calcMappingGradient(	phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, scalingFactors_ISAT, 
												luSolver_ISAT, (tf-t0), numberSubSteps_ISAT, &odeSolverConstantPressure());
							
									// add a new leaf 
									bool flag = isat_HOM->add(phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, phi0base); 
						
									if(flag == false)	
										Info << "ISAT Error - Addition process failed..." << endl;
									
									nAddHOM++;
									
									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
									
									cpuTimeAddition += (t2-t1);
								}
							}	
				
							// check if is necessary to balance
							isat_HOM->cleanAndBalance();
						}
					}
					else
					{
						Info << "ISAT can be used only with constant pressure reactors" << endl;
						abort();
					}
				}
			
				// Check mass fractions
				normalizeMassFractions(yf, NC, celli, massFractionsTol);

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

				if (counter%(int(0.20*mesh.nCells())+1) == 0)
					Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

				counter++;

				// Output
				if (runTime.outputTime())
				{
					if (constPressureBatchReactor == true)
					{
						QCells[celli] = batchReactorHomogeneousConstantPressure.QR();
						if (outputFormationRatesIndices.size() != 0)
						{
							#if OPENFOAM_VERSION >= 40
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].ref()[celli] = 	batchReactorHomogeneousConstantPressure.R()[outputFormationRatesIndices[i]+1] *
														thermodynamicsMapXML->MW(outputFormationRatesIndices[i]);
							#else
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = 	batchReactorHomogeneousConstantPressure.R()[outputFormationRatesIndices[i]+1] *
														thermodynamicsMapXML->MW(outputFormationRatesIndices[i]);
							#endif
						}
					}
					else
					{
						Info << "ISAT can be used only with constant pressure reactors" << endl;
						abort();
					}
				}

			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;

			if(isatCheck == true) 
			{
				Info << endl;
				Info << " ********* ISAT HOM stats **********" << endl;
				
				Info << "   Direct Integration : " << isat_HOM->nAdd()+isat_HOM->nGrow()  << " (" << nAddHOM+nGrowHOM << ")" << " (" << (nAddHOM+nGrowHOM)/double(mesh.nCells())*100. << "%)" << endl;
				Info << "      Add             : " << isat_HOM->nAdd()  << " (" << nAddHOM  << ")" << " (" << nAddHOM/double(mesh.nCells())*100. << "%)" << endl;
				Info << "      Grow            : " << isat_HOM->nGrow() << " (" << nGrowHOM << ")" << " (" << nGrowHOM/double(mesh.nCells())*100. << "%)" << endl;
				Info << "   Retrieve           : " << isat_HOM->nUse()  << " (" << nRetHOM  << ")" << " (" << nRetHOM/double(mesh.nCells())*100. << "%)" << endl;
				Info << endl;				
		
				const double cpuTimeIntegration = cpuTimeDI + cpuTimeGrowth + cpuTimeAddition;
				Info << "   CPU Integration  : " << cpuTimeIntegration  << " (" << cpuTimeIntegration/(tEnd-tStart)*100. << "%)" << endl;
				Info << "     CPU DI         : " << cpuTimeDI           << " (" << cpuTimeDI/(tEnd-tStart)*100.          << "%)" << endl;				
				Info << "     CPU Growth     : " << cpuTimeGrowth       << " (" << cpuTimeGrowth/(tEnd-tStart)*100.         << "%)" << endl;
				Info << "     CPU Addition   : " << cpuTimeAddition     << " (" << cpuTimeAddition/(tEnd-tStart)*100.    << "%)" << endl;
				Info << "   CPU Retrieve     : " << cpuTimeRet          << " (" << cpuTimeRet/(tEnd-tStart)*100.         << "%)" << endl;
				Info << endl;

				Info << "      BTS  : " << isat_HOM->nBTS()  << endl;
				Info << "      MRU  : " << isat_HOM->nMRU()  << endl;
				Info << "      MFU  : " << isat_HOM->nMFU()  << endl << endl;
				Info << endl;
			}
		}
	}
	else if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	{
		Info << "ISAT can be used only in conjuction with the OpenSMOKE++ ODE solver" << endl;
		abort();
	}
    
    Info<< " * T gas min/max (after chemistry) = " << min(T).value() << ", " << max(T).value() << endl;
}
