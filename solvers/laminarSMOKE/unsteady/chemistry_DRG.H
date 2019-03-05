/*-----------------------------------------------------------------------*\
|                                                                         |
|                                                   |
|                                                   | 
|                                   |
|                                   |
|                                  |
|                                  |
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
	#if OPENFOAM_VERSION >= 40
	scalarField& drgSpeciesCells = drg_important_species.ref();
	scalarField& drgReactionsCells = drg_important_reactions.ref();
	scalarField& cpuReactionsCells = cpu_reaction_times.ref();
	#else
	scalarField& drgSpeciesCells = drg_important_species.internalField();
	scalarField& drgReactionsCells = drg_important_reactions.internalField();
	scalarField& cpuReactionsCells = cpu_reaction_times.internalField();
	#endif


	if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
		
		// Min and max values
		Eigen::VectorXd yMin;
		Eigen::VectorXd yMax;
		
		// Initial and final values
		Eigen::VectorXd y0;
		Eigen::VectorXd yf;

		// Auxiliary vectors
		Eigen::VectorXd yff(NC+1);
		OpenSMOKE::OpenSMOKEVectorDouble omega_(NC);
		OpenSMOKE::OpenSMOKEVectorDouble x_(NC);
		OpenSMOKE::OpenSMOKEVectorDouble c_(NC);

		Info <<" * Solving homogeneous chemistry with DRG (OpenSMOKE solver)... "<<endl;
		{	
			unsigned int counter = 0;
			unsigned int counter_skipped = 0;
			
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			forAll(TCells, celli)
			{
				//- Solving for celli:	
				double tLocalStart = OpenSMOKE::OpenSMOKEGetCpuTime();
				if (TCells[celli] > drg_minimum_temperature_for_chemistry)
				{
					for(unsigned int i=0;i<NC;i++)
						omega_[i+1] = Y[i].internalField()[celli];
							
					double MW_;
					thermodynamicsMapXML->MoleFractions_From_MassFractions(x_, MW_, omega_);
					double cTot_ = thermodynamicPressure/PhysicalConstants::R_J_kmol/TCells[celli];
					Product(cTot_, x_, &c_);

					for(unsigned int i=0;i<drg_temperature.size();i++)
						if(TCells[celli] < drg_temperature[i])
						{
							drg->SetEpsilon(drg_epsilon[i]);
							break;
						}

					drg->Analysis(TCells[celli], thermodynamicPressure, c_);
						
					unsigned int NEQ = drg->number_important_species()+1;
					y0.resize(NEQ);
					yf.resize(NEQ);
					yMax.resize(NEQ);	yMax.setConstant(1.); yMax(NEQ-1) = 6000.; 
					yMin.resize(NEQ);	yMin.setConstant(0.); yMin(NEQ-1) = 200.;

					for (unsigned int i=0;i<drg->number_important_species();++i)	
					{
						const unsigned int j = drg->indices_important_species()[i]+1;
						y0(i) = omega_[j];
					}
					y0(NEQ-1) = TCells[celli];
					
					if (constPressureBatchReactor == true)
					{
						// Set reactor
						batchReactorHomogeneousConstantPressure.SetReactor(thermodynamicPressure);
						batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);

						// Set DRG
						batchReactorHomogeneousConstantPressure.SetDRG(drg);
						batchReactorHomogeneousConstantPressure.SetMassFractions(omega_);
						
						// Set initial conditions
						OdeSMOKE::MultiValueSolver<methodGearConstantPressure> odeSolverConstantPressureDRG;
						odeSolverConstantPressureDRG.SetReactor(&batchReactorHomogeneousConstantPressure);
						odeSolverConstantPressureDRG.SetInitialConditions(t0, y0);

						// Additional ODE solver options
						if (celli == 0)
						{
							// Set linear algebra options
							odeSolverConstantPressureDRG.SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
							odeSolverConstantPressureDRG.SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

							// Set relative and absolute tolerances
							odeSolverConstantPressureDRG.SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
							odeSolverConstantPressureDRG.SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

							// Set minimum and maximum values
							odeSolverConstantPressureDRG.SetMinimumValues(yMin);
							odeSolverConstantPressureDRG.SetMaximumValues(yMax);
						}
						
						// Solve
						OdeSMOKE::OdeStatus status = odeSolverConstantPressureDRG.Solve(tf);
						odeSolverConstantPressureDRG.Solution(yf);
					}
					else
					{
						Info << "DRG Analysis can be used only with constant pressure reactors" << endl;
						abort();
					}

					for(unsigned int i=0;i<NC;i++)
					yff(i) = omega_[i+1];
					yff(NC) = yf(NEQ-1);
				
					for (unsigned int i=0;i<drg->number_important_species();++i)	
					{
						const unsigned int j = drg->indices_important_species()[i];
						yff(j) = yf(i);
					}

					// Check mass fractions
					normalizeMassFractions(yff, celli, massFractionsTol, vc_main_species);

					// Assign mass fractions
					#if OPENFOAM_VERSION >= 40
					for(int i=0;i<NC;i++)
						Y[i].ref()[celli] = yff(i);
					#else
					for(int i=0;i<NC;i++)
						Y[i].internalField()[celli] = yff(i);
					#endif

					//- Allocating final values: temperature
					if (energyEquation == true)
						TCells[celli] = yff(NC);

					// Fill the relevant DRG fields
					drgSpeciesCells[celli]   = drg->number_important_species();
					drgReactionsCells[celli] = drg->number_important_reactions();
				}
				else
				{
					drgSpeciesCells[celli]   = 0;
					drgReactionsCells[celli] = 0;
					counter_skipped++;
				}	
				double tLocalEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

				cpuReactionsCells[celli] = (tLocalEnd - tLocalStart)*1.e6;		

				if (counter%(int(0.20*mesh.nCells())+1) == 0)
					Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << "(" << counter_skipped << ")" << endl;

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
								FormationRates[i].ref()[celli] = batchReactorHomogeneousConstantPressure.R()[outputFormationRatesIndices[i]+1] *
                                       	      		                                                   thermodynamicsMapXML->MW()[outputFormationRatesIndices[i]+1];
							#else
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantPressure.R()[outputFormationRatesIndices[i]+1] *
                                       	      		                                                   thermodynamicsMapXML->MW()[outputFormationRatesIndices[i]+1];
							#endif	
						}
					}
					else
					{
						QCells[celli] = batchReactorHomogeneousConstantVolume.QR();
						if (outputFormationRatesIndices.size() != 0)
						{
							#if OPENFOAM_VERSION >= 40
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].ref()[celli] = batchReactorHomogeneousConstantVolume.R()[outputFormationRatesIndices[i]+1] *
                                           	                                                           thermodynamicsMapXML->MW()[outputFormationRatesIndices[i]+1];
							#else
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantVolume.R()[outputFormationRatesIndices[i]+1] *
                                           	                                                           thermodynamicsMapXML->MW()[outputFormationRatesIndices[i]+1];
							#endif	
						}
					}
				}
			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
		}
	}
	else if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	{
		Info << "DRG can be used only in conjuction with the OpenSMOKE++ ODE solver" << endl;
		abort();
	}
    
    Info<< " * T gas min/max (after chemistry) = " << min(T).value() << ", " << max(T).value() << endl;
}
