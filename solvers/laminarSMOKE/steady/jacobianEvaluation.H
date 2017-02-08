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

if (homogeneousReactions == true && (speciesEquations == true || energyEquation == true))
{
	
	if (implicitSourceTerm == true && jacobianUpdate == 1)
	{
		Info<< "Jacobian evaluation (continous)... " << endl;

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			const unsigned int NE = NC+1+hmom_number_moments;

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField();
			const scalarField& HMOM_Norm_N0Cells  = HMOM_Norm_N0.internalField();
			const scalarField& HMOM_Norm_M00Cells = HMOM_Norm_M00.internalField();
			const scalarField& HMOM_Norm_M10Cells = HMOM_Norm_M10.internalField();
			const scalarField& HMOM_Norm_M01Cells = HMOM_Norm_M01.internalField();
	
			Eigen::VectorXd J(NE);
			OpenSMOKE::OpenSMOKEVectorDouble Source(NE);	
			OpenSMOKE::OpenSMOKEVectorDouble y(NE);

			forAll(TCells, celli)
			{
				for(int i=0;i<NC;i++)
					y[i+1] = Y[i].internalField()[celli];
				y[NC+1] = TCells[celli];
				y[NC+2] = HMOM_Norm_M00Cells[celli];
				y[NC+3] = HMOM_Norm_M10Cells[celli];
				y[NC+4] = HMOM_Norm_M01Cells[celli];
				y[NC+5] = HMOM_Norm_N0Cells[celli];

				linear_model.reactionSourceTerms(	*thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], Source);


				if (sparseJacobian == false)
					linear_model.reactionJacobian( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J );
				else
					linear_model.reactionJacobianSparse( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J, energyEquation, speciesEquations);

				for(int i=0;i<NE;i++)
					sourceImplicit[i].internalField()[celli] = J(i);
			
				for(int i=0;i<NE;i++)
					sourceExplicit[i].internalField()[celli] = Source[i+1] - J(i)*y[i+1];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;

	}
	else if (implicitSourceTerm == true && jacobianUpdate != 1)
	{
		Info<< "Jacobian evaluation (discrete)... " << endl;

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			const unsigned int NE = NC+1+hmom_number_moments;

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
			const scalarField& HMOM_Norm_N0Cells  = HMOM_Norm_N0.internalField();
			const scalarField& HMOM_Norm_M00Cells = HMOM_Norm_M00.internalField();
			const scalarField& HMOM_Norm_M10Cells = HMOM_Norm_M10.internalField();
			const scalarField& HMOM_Norm_M01Cells = HMOM_Norm_M01.internalField();
	
			Eigen::VectorXd J(NE);
			OpenSMOKE::OpenSMOKEVectorDouble Source(NE);	
			OpenSMOKE::OpenSMOKEVectorDouble y(NE);

			forAll(TCells, celli)
			{
				for(int i=0;i<NC;i++)
					y[i+1] = Y[i].internalField()[celli];
				y[NC+1] = TCells[celli];
				y[NC+2] = HMOM_Norm_M00Cells[celli];
				y[NC+3] = HMOM_Norm_M10Cells[celli];
				y[NC+4] = HMOM_Norm_M01Cells[celli];
				y[NC+5] = HMOM_Norm_N0Cells[celli];

				linear_model.reactionSourceTerms(	*thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], Source);

				if (jacobianCounter == jacobianUpdate)
				{
					if (sparseJacobian == false)
						linear_model.reactionJacobian( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J );
					else
						linear_model.reactionJacobianSparse( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J, energyEquation, speciesEquations);

					Jstore[celli] = J;
				}

				for(int i=0;i<NE;i++)
					sourceImplicit[i].internalField()[celli] = Jstore[celli](i);
			
				for(int i=0;i<NE;i++)
					sourceExplicit[i].internalField()[celli] = Source[i+1] - Jstore[celli](i)*y[i+1];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;

		if (jacobianCounter == jacobianUpdate)
			jacobianCounter = 0;

		jacobianCounter++;
	}
	else
	{
		Info<< "Source terms evaluation... " << endl;

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			const unsigned int NE = NC+1+hmom_number_moments;

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
			const scalarField& HMOM_Norm_N0Cells  = HMOM_Norm_N0.internalField();
			const scalarField& HMOM_Norm_M00Cells = HMOM_Norm_M00.internalField();
			const scalarField& HMOM_Norm_M10Cells = HMOM_Norm_M10.internalField();
			const scalarField& HMOM_Norm_M01Cells = HMOM_Norm_M01.internalField();
	
			OpenSMOKE::OpenSMOKEVectorDouble Source(NE);	
			OpenSMOKE::OpenSMOKEVectorDouble y(NE);

			forAll(TCells, celli)
			{
				for(int i=0;i<NC;i++)
					y[i+1] = Y[i].internalField()[celli];
				y[NC+1] = TCells[celli];
				y[NC+2] = HMOM_Norm_M00Cells[celli];
				y[NC+3] = HMOM_Norm_M10Cells[celli];
				y[NC+4] = HMOM_Norm_M01Cells[celli];
				y[NC+5] = HMOM_Norm_N0Cells[celli];

				linear_model.reactionSourceTerms( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], Source );

				for(int i=0;i<NE;i++)
					sourceImplicit[i].internalField()[celli] = 0.;
			
				for(int i=0;i<NE;i++)
					sourceExplicit[i].internalField()[celli] = Source[i+1];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;	
	}
}
