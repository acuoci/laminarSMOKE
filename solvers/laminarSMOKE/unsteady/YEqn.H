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

// Mass Fluxes
#include "fluxes.H"

// Discretization schemes
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi)")
    )
);


{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
   	volScalarField& Yi = Y[i];
	 
	fvScalarMatrix YiEqn
	(
		fvm::ddt(rho, Yi)
	      + mvConvection->fvmDiv(phi, Yi)
	      - fvm::laplacian(rho*alpha/LewisNumbers(i), Yi)
		== 
              - fvm::div(Jc,Yi, "div(Jc,Yi)")
	      + RR[i]
	      + fvOptions(rho, Yi)
	);


	// Solve
	YiEqn.relax();
	fvOptions.constrain(YiEqn);
	#if OPENFOAM_VERSION >= 1000
	YiEqn.solve("Yi");
	#else
	YiEqn.solve(mesh.solver("Yi"));
	#endif
	fvOptions.correct(Yi);	
	
	// Sum of mass fractions
    	Yi.max(0.0);
	Yt += Yi;
    }

    one_minus_Yt = scalar(1.) - Yt;

    const scalarField& YtCells = Yt.internalField();
    forAll(YtCells, celli)
    {
	const double one_plus_alpha = 1./YtCells[celli];
	for (label i=0; i<Y.size(); i++)
		#if OPENFOAM_VERSION >= 40
		Y[i].ref()[celli] *= one_plus_alpha;
		#else
		Y[i].internalField()[celli] *= one_plus_alpha;
		#endif
    }

    // Boundaries
    forAll(Yt.boundaryField(), patchi)
    {
	#if OPENFOAM_VERSION >= 40
	fvPatchScalarField& pYt = Yt.boundaryFieldRef()[patchi];
	#else
	fvPatchScalarField& pYt = Yt.boundaryField()[patchi];
	#endif

	forAll(pYt, facei)
	{
		const double one_plus_alpha = 1./pYt[facei];
		for (label i=0; i<Y.size(); i++)
			#if OPENFOAM_VERSION >= 40
			Y[i].boundaryFieldRef()[patchi][facei] *= one_plus_alpha;
			#else
			Y[i].boundaryField()[patchi][facei] *= one_plus_alpha;
			#endif
			
	}
    }

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	
    Info << "Transport equations of species solved in " << tEnd - tStart << " s " << endl;
}

