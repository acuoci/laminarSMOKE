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
	volScalarField rAU(1.0/UEqn.A());
	volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
	volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
	tUEqn.clear();

	bool closedVolume = false;

	surfaceScalarField phiHbyA("phiHbyA", fvc::flux(rho*HbyA));
	MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

	volScalarField rhorAtU("rhorAtU", rho*rAtU);

	// Update the pressure BCs to ensure flux consistency
	constrainPressure(p, rho, U, phiHbyA, rhorAtU, MRF);

	if (simple.transonic())
	{
	    surfaceScalarField phid
	    (
		"phid",
		(fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
	    );

	    phiHbyA +=
		fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf()
	      - fvc::interpolate(p)*phid;

	    HbyA -= (rAU - rAtU)*fvc::grad(p);

	    while (simple.correctNonOrthogonal())
	    {
		fvScalarMatrix pEqn
		(
		    fvc::div(phiHbyA)
		  + fvm::div(phid, p)
		  - fvm::laplacian(rhorAtU, p)
		 ==
		    fvOptions(psi, p, rho.name())
		);

		// Relax the pressure equation to maintain diagonal dominance
		pEqn.relax();

		#if DEVVERSION == 1
		pEqn.setReference
		(
		    pressureControl.refCell(),
		    pressureControl.refValue()
		);
		#else
		pEqn.setReference(pRefCell, pRefValue);
		#endif

		pEqn.solve();

		if (simple.finalNonOrthogonalIter())
		{
		    phi = phiHbyA + pEqn.flux();
		}
	    }
	}
	else
	{
	    closedVolume = adjustPhi(phiHbyA, U, p);

	    phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf();
	    HbyA -= (rAU - rAtU)*fvc::grad(p);

	    while (simple.correctNonOrthogonal())
	    {
		fvScalarMatrix pEqn
		(
		    fvc::div(phiHbyA)
		  - fvm::laplacian(rhorAtU, p)
		  ==
		    fvOptions(psi, p, rho.name())
		);

		#if DEVVERSION == 1
		pEqn.setReference
		(
		    pressureControl.refCell(),
		    pressureControl.refValue()
		);
		#else
		pEqn.setReference(pRefCell, pRefValue);
		#endif

		pEqn.solve();

		if (simple.finalNonOrthogonalIter())
		{
		    phi = phiHbyA + pEqn.flux();
		}
	    }
	}

	// The incompressibe form of the continuity error check is appropriate for
	// steady-state compressible also.
	#include "continuityErrs.H"

	// Explicitly relax pressure for momentum corrector
	p.relax();

	U = HbyA - rAtU*fvc::grad(p);
	U.correctBoundaryConditions();
	fvOptions.correct(U);

	#if DEVVERSION == 1
	pressureControl.limit(p);
	#endif

	// For closed-volume cases adjust the pressure and density levels
	// to obey overall mass continuity
	if (closedVolume)
	{
	    p += (initialMass - fvc::domainIntegrate(psi*p))
		/fvc::domainIntegrate(psi);
	}

	p.correctBoundaryConditions();

	// Recalculate density from the relaxed pressure
	//rho = thermo.rho();
    	rho = Rho(p,psi);
    	rho = max(rho, rhoMin);
    	rho = min(rho, rhoMax);

	if (!simple.transonic())
	{
	    rho.relax();
	}

	Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

	if (includeDpDt == true)
	{
		DpDt = fvc::div(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
	}
}
