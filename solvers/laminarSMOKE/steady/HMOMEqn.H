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
		unsigned int index = Y.size()+4;

		fvScalarMatrix HMOM_Norm_N0Eqn
		(
		        fvm::div(phi, HMOM_Norm_N0) 
		      - fvm::laplacian(rho*Dhmom, HMOM_Norm_N0) 
		        == 
		        // source terms
			rho*sourceExplicit[index] + fvm::Sp(rho*sourceImplicit[index],HMOM_Norm_N0) + 
		     	fvOptions(rho, HMOM_Norm_N0)
		);

		// Add thermophoretic effect
		if (thermophoreticEffect == true)
			HMOM_Norm_N0Eqn -= fvc::laplacian(0.55*mu/T*HMOM_Norm_N0, T, "laplacian(teta,hmom)");

		// Solve
		HMOM_Norm_N0Eqn.relax();
		fvOptions.constrain(HMOM_Norm_N0Eqn);
		HMOM_Norm_N0Eqn.solve(mesh.solver("HMOM_Norm_N0"));

		HMOM_Norm_N0.max(hmom_minimum_normalized_value);
		HMOM_Norm_N0.min(1e0);
	}

	{
		unsigned int index = Y.size()+1;

		fvScalarMatrix HMOM_Norm_M00Eqn
		(
			fvm::div(phi, HMOM_Norm_M00) 
		      - fvm::laplacian(rho*Dhmom, HMOM_Norm_M00) 
		        == 
		        // source terms
			rho*sourceExplicit[index] + fvm::Sp(rho*sourceImplicit[index],HMOM_Norm_M00) + 
		     	fvOptions(rho, HMOM_Norm_M00)
		);

		// Add thermophoretic effect
		if (thermophoreticEffect == true)
			HMOM_Norm_M00Eqn -= fvc::laplacian(0.55*mu/T*HMOM_Norm_M00, T, "laplacian(teta,hmom)");

		// Solve
		HMOM_Norm_M00Eqn.relax();
		fvOptions.constrain(HMOM_Norm_M00Eqn);
		HMOM_Norm_M00Eqn.solve(mesh.solver("HMOM_Norm_M00"));

		HMOM_Norm_M00.max(hmom_minimum_normalized_value);
		HMOM_Norm_M00.min(1e0);
	}

	{
		unsigned int index = Y.size()+2;

		fvScalarMatrix HMOM_Norm_M10Eqn
		(
			fvm::div(phi, HMOM_Norm_M10) 
		      - fvm::laplacian(rho*Dhmom, HMOM_Norm_M10) 
		        == 
		        // source terms
			rho*sourceExplicit[index] + fvm::Sp(rho*sourceImplicit[index],HMOM_Norm_M10) + 
		     	fvOptions(rho, HMOM_Norm_M10)
		);

		// Add thermophoretic effect
		if (thermophoreticEffect == true)
			HMOM_Norm_M10Eqn -= fvc::laplacian(0.55*mu/T*HMOM_Norm_M10, T, "laplacian(teta,hmom)");

		// Solve
		HMOM_Norm_M10Eqn.relax();
		fvOptions.constrain(HMOM_Norm_M10Eqn);
		HMOM_Norm_M10Eqn.solve(mesh.solver("HMOM_Norm_M10"));

		HMOM_Norm_M10.max(hmom_minimum_normalized_value);
		HMOM_Norm_M10.min(1e0);
	}

	{
		unsigned int index = Y.size()+3;

		fvScalarMatrix HMOM_Norm_M01Eqn
		(
			fvm::div(phi, HMOM_Norm_M01) 
		      - fvm::laplacian(rho*Dhmom, HMOM_Norm_M01) 
		        == 
		        // source terms
			rho*sourceExplicit[index] + fvm::Sp(rho*sourceImplicit[index],HMOM_Norm_M01) + 
		     	fvOptions(rho, HMOM_Norm_M01)
		);

		// Add thermophoretic effect
		if (thermophoreticEffect == true)
			HMOM_Norm_M01Eqn -= fvc::laplacian(0.55*mu/T*HMOM_Norm_M01, T, "laplacian(teta,hmom)");

		// Solve
		HMOM_Norm_M01Eqn.relax();
		fvOptions.constrain(HMOM_Norm_M01Eqn);
		HMOM_Norm_M01Eqn.solve(mesh.solver("HMOM_Norm_M01"));

		HMOM_Norm_M01.max(hmom_minimum_normalized_value);
		HMOM_Norm_M01.min(1e0);
	}

}
