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
	volScalarField rAU(1.0/UEqn().A());
	volScalarField rAtU(1.0/(1.0/rAU - UEqn().H1()));

	volVectorField HbyA("HbyA", U);
	HbyA = rAU*UEqn().H();

	UEqn.clear();

	bool closedVolume = false;

	if (simple.transonic())
	{
	    surfaceScalarField phid
	    (
		"phid",
		fvc::interpolate(psi)
	       *(fvc::interpolate(HbyA) & mesh.Sf())
	    );

	    MRF.makeRelative(fvc::interpolate(psi), phid);

	    surfaceScalarField phic
	    (
		"phic",
		fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf()
	    );

	    HbyA -= (rAU - rAtU)*fvc::grad(p);

	    volScalarField rhorAtU("rhorAtU", rho*rAtU);

	    while (simple.correctNonOrthogonal())
	    {
		fvScalarMatrix pEqn
		(
		    fvm::div(phid, p)
		  + fvc::div(phic)
		  - fvm::laplacian(rhorAtU, p)
		 ==
		    fvOptions(psi, p, rho.name())
		);

		// Relax the pressure equation to maintain diagonal dominance
		pEqn.relax();

		pEqn.setReference(pRefCell, pRefValue);

		pEqn.solve();

		if (simple.finalNonOrthogonalIter())
		{
		    phi == phic + pEqn.flux();
		}
	    }
	}
	else
	{
	    surfaceScalarField phiHbyA
	    (
		"phiHbyA",
		fvc::interpolate(rho*HbyA) & mesh.Sf()
	    );

	    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

	    closedVolume = adjustPhi(phiHbyA, U, p);

	    phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf();
	    HbyA -= (rAU - rAtU)*fvc::grad(p);

	    volScalarField rhorAtU("rhorAtU", rho*rAtU);

	    while (simple.correctNonOrthogonal())
	    {
		fvScalarMatrix pEqn
		(
		    fvc::div(phiHbyA)
		  - fvm::laplacian(rhorAtU, p)
		  ==
		    fvOptions(psi, p, rho.name())
		);

		pEqn.setReference(pRefCell, pRefValue);

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

	// For closed-volume cases adjust the pressure and density levels
	// to obey overall mass continuity
	if (closedVolume)
	{
	    p += (initialMass - fvc::domainIntegrate(psi*p))
		/fvc::domainIntegrate(psi);
	}

	// Recalculate density from the relaxed pressure
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

