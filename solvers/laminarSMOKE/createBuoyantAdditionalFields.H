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
	#endif

	// Kinetic energy
	Info<< "Creating field kinetic energy K\n" << endl;
	volScalarField K("K", 0.5*magSqr(U));

	#if   OPENFOAM_VERSION >= 40

	#include "readhRef.H"
	#include "gh.H"
	#include "readpRef.H"

	volScalarField p_rgh
	(
	    IOobject
	    (
		"p_rgh",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	);

	mesh.setFluxRequired(p_rgh.name());

	// Force p_rgh to be consistent with p
	p_rgh = p - rho*gh;

	#if STEADYSTATE == 1
	#endif

	#elif   OPENFOAM_VERSION == 30

		#include "readhRef.H"
		#include "gh.H"

		volScalarField p_rgh
		(
		    IOobject
		    (
			"p_rgh",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		    ),
		    mesh
		);

		// Force p_rgh to be consistent with p
		p_rgh = p - rho*gh;

		mesh.setFluxRequired(p_rgh.name());

		#if STEADYSTATE == 1

			label pRefCell = 0;
			scalar pRefValue = 0.0;

			setRefCell
			(
				p,
				p_rgh,
				simple.dict(),
				pRefCell,
				pRefValue
			);

		#endif
	
	#else

		Info<< "Calculating field g.h\n" << endl;
		volScalarField gh("gh", g & mesh.C());

		surfaceScalarField ghf("ghf", g & mesh.Cf());

		volScalarField p_rgh
		(
			IOobject
			(
			    "p_rgh",
			    runTime.timeName(),
			    mesh,
			    IOobject::MUST_READ,
			    IOobject::AUTO_WRITE
			),
			mesh
		);

		// Force p_rgh to be consistent with p
		p_rgh = p - rho*gh;

		#if STEADYSTATE == 1

			label pRefCell = 0;
			scalar pRefValue = 0.0;

			setRefCell
			(
				p,
				p_rgh,
				mesh.solutionDict().subDict("SIMPLE"),
				pRefCell,
				pRefValue
			);

		#endif

	#endif

	// Interpolation scheme
	multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
	forAll(Y, i)
	{
		fields.add(Y[i]);
	}

	// Spark
	spark.Setup(mesh);

	// Check the gravity
	if ( (mag(g)).value() == 0.)
	{
		Info << endl;
		Info << "------------------------------------------------------------------" << endl;
		Info << " ERROR:   gravity is turned off                                   " << endl;
		Info << "          Please consider that this solver is conceived           " << endl;
 		Info << "          for flows in which gravity play a major role.           " << endl;
 		Info << "          If you think gravity plays a negligible role, consider  " << endl;
  		Info << "          the possibility to use the laminarPimpleSMOKE solver    " << endl;
		Info << "------------------------------------------------------------------" << endl;
		Info << endl;
		abort();
	}
