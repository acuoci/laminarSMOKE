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

//- Creating field rho [kg/m3]
Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	Rho(p,psi)
);

#include "compressibleCreatePhi.H"

#if   OPENFOAM_VERSION == 30 
mesh.setFluxRequired(p.name());
#endif

#if   OPENFOAM_VERSION == 40 
#if STEADYSTATE != 1
pressureControl pressureControl(p, rho, pimple.dict(), false);
#else
pressureControl pressureControl(p, rho, simple.dict());
#endif
mesh.setFluxRequired(p.name());
#endif

#if STEADYSTATE != 1

Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt
(
	IOobject
	(
		"DpDt",
		runTime.timeName(),
		mesh
	),
	mesh,
	dimensionedScalar("DpDt", p.dimensions()/dimTime, 0)
);

Info<< "Creating DeltaT field" << endl;
volScalarField DeltaT
(
	IOobject
	(
		"DeltaT",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("one", dimTime, 0),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

#else

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);

#endif

// Kinetic energy
Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

// Interpolation scheme
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
	fields.add(Y[i]);
}

// Spark
spark.Setup(mesh);

	// Check the gravity
	if ( (mag(g)).value() != 0.)
	{
		Info << endl;
		Info << "------------------------------------------------------------------" << endl;
		Info << " ERROR:   gravity is turned on                                    " << endl;
		Info << "          Please consider that this solver is not conceived       " << endl;
 		Info << "          for bouyancy driven flows. If you think gravity         " << endl;
 		Info << "          plays a non-negligible role, consider the possibility   " << endl;
  		Info << "          to use the laminarBouyantPimpleSMOKE solver (for        " << endl;
  		Info << "          unsteady flows) or the laminarBouyantSimpleSMOKE solver " << endl;
  		Info << "          (for steady-state flows)                                " << endl;
		Info << "------------------------------------------------------------------" << endl;
		Info << endl;
		abort();
	}
