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
	volScalarField Dhmom
	(
	    IOobject
	    (
		"Dhmom",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("Dhmom", dimensionSet(0., 2., -1., 0., 0.), 1.e-8)
	);

	{
		fvScalarMatrix HMOM_Norm_N0Eqn
		(
			fvm::ddt(rho,HMOM_Norm_N0)
		      + fvm::div(phi, HMOM_Norm_N0) 
		      - fvm::laplacian(rho*Dhmom, HMOM_Norm_N0) 
		        == 
		        // source terms
		     	fvOptions(rho, HMOM_Norm_N0)
		);

		// Solve
		HMOM_Norm_N0Eqn.relax();
		fvOptions.constrain(HMOM_Norm_N0Eqn);
		HMOM_Norm_N0Eqn.solve(mesh.solver("HMOM_Norm_N0"));

		HMOM_Norm_N0.max(1e-32);
	}

	{
		fvScalarMatrix HMOM_Norm_M00Eqn
		(
			fvm::ddt(rho,HMOM_Norm_M00)
		      + fvm::div(phi, HMOM_Norm_M00) 
		      - fvm::laplacian(rho*Dhmom, HMOM_Norm_M00) 
		        == 
		        // source terms
		     	fvOptions(rho, HMOM_Norm_M00)
		);

		// Solve
		HMOM_Norm_M00Eqn.relax();
		fvOptions.constrain(HMOM_Norm_M00Eqn);
		HMOM_Norm_M00Eqn.solve(mesh.solver("HMOM_Norm_M00"));

		HMOM_Norm_M00.max(1e-32);
	}

	{
		fvScalarMatrix HMOM_Norm_M10Eqn
		(
			fvm::ddt(rho,HMOM_Norm_M10)
		      + fvm::div(phi, HMOM_Norm_M10) 
		      - fvm::laplacian(rho*Dhmom, HMOM_Norm_M10) 
		        == 
		        // source terms
		     	fvOptions(rho, HMOM_Norm_M10)
		);

		// Solve
		HMOM_Norm_M10Eqn.relax();
		fvOptions.constrain(HMOM_Norm_M10Eqn);
		HMOM_Norm_M10Eqn.solve(mesh.solver("HMOM_Norm_M10"));

		HMOM_Norm_M10.max(1e-32);
	}

	{
		fvScalarMatrix HMOM_Norm_M01Eqn
		(
			fvm::ddt(rho,HMOM_Norm_M01)
		      + fvm::div(phi, HMOM_Norm_M01) 
		      - fvm::laplacian(rho*Dhmom, HMOM_Norm_M01) 
		        == 
		        // source terms
		     	fvOptions(rho, HMOM_Norm_M01)
		);

		// Solve
		HMOM_Norm_M01Eqn.relax();
		fvOptions.constrain(HMOM_Norm_M01Eqn);
		HMOM_Norm_M01Eqn.solve(mesh.solver("HMOM_Norm_M01"));

		HMOM_Norm_M01.max(1e-32);
	}

}
