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
	#if OPENFOAM_VERSION >= 60
		#include "rhoEqn.H"
	#elif OPENFOAM_VERSION >= 40
		if (pimple.nCorrPIMPLE() <= 1)
		{
		    #include "rhoEqn.H"
		}
	#else
		#include "rhoEqn.H"
	#endif

	// --- Pressure-velocity PIMPLE corrector loop
	while (pimple.loop())
	{
		#include "UEqn.H"

		// Chemical step
		{
			double t0 = runTime.value() - runTime.deltaT().value();
			double tf = runTime.value();
			#include "chemistry.H"
			#include "properties.H"
		}

		#include "YEqn.H"
		#include "TEqn.H"

		// --- Pressure corrector loop
		while (pimple.correct())
		{
			#if OPENFOAM_VERSION >= 40
			if (pimple.consistent())
			{
				#include "pcEqn.H"
			}
			else
			#endif
			{
				#include "pEqn.H"
			}
		}
	}

	#if OPENFOAM_VERSION <= 30
	rho = Rho(p,psi);
	#endif
}
