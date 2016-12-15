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
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "readGravitationalAcceleration.H"
	#include "createBasicFields.H"
	#include "readOptions.H"
	#include "createGasFields.H"
	#include "createChemicalFields.H"
	#include "createLaminarSMOKEThermoClass.H"
	#include "createDRGFields.H"
	#include "createFvOptions.H"
	#include "createRadiationModel.H"
	#include "memoryAllocation.H"
	#include "properties.H"
	#include "createAdditionalFields.H"
	#include "initContinuityErrs.H"
	#include "readTimeControls.H"
	#include "compressibleCourantNo.H"
	#include "setInitialDeltaT.H"
	pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"


		runTime++;
		Info<< "Time = " << runTime.timeName() << nl << endl;

		if (momentumEquations == false)
		{		
			if (strangAlgorithm == STRANG_MOMENTUM_TRANSPORT_REACTION)
			{
				#include "Policy_TransportReaction.H"
			}
			else if (strangAlgorithm == STRANG_MOMENTUM_REACTION_TRANSPORT)
			{
				#include "Policy_ReactionTransport.H"
			}
			else if (strangAlgorithm == STRANG_COMPACT)
			{
				#include "Policy_Compact.H"
			}
		}
		else
		{
			if (strangAlgorithm == STRANG_MOMENTUM_TRANSPORT_REACTION)
			{
				#include "Policy_MomentumTransportReaction.H"
			}
			else if (strangAlgorithm == STRANG_MOMENTUM_REACTION_TRANSPORT)
			{
				#include "Policy_MomentumReactionTransport.H"
			}
			else if (strangAlgorithm == STRANG_COMPACT)
			{
				#include "Policy_MomentumCompact.H"
			}		
		}

		// Passive scalars
            	#include "zMixEqn.H"
		#include "tauEqn.H"

		// Local post processing
		#include "localPostProcessing.H"
	
		runTime.write();

		
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
