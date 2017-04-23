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

if (calculateLocalStrainRate == true)
{
	Info<< "Reading zMix field\n" << endl;
	volScalarField zMix
	(
	    IOobject
	    (
		"zMix",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	);

	volScalarField strain_rate
	(
	    IOobject
	    (
		"strain_rate",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionSet(0, 0, -1, 0, 0)
	);

	dimensionedScalar epsilon_mag_grad_z
	(
		"epsilon_mag_grad_z",
		dimensionSet(0,-1,0,0,0,0,0),
		scalar(1.e-9)
	); 

	// Calculate
	volVectorField grad_z = fvc::grad(zMix);
	volScalarField mag_grad_z_uncorrected = mag(grad_z); 
	volScalarField mag_grad_z = mag(grad_z) + epsilon_mag_grad_z; 
	strain_rate = -( sqr(grad_z) && fvc::grad(U) )/sqr(mag_grad_z)+ fvc::div(U); 

	// Clean
	{
		const scalarField& zMixCells = zMix.internalField();
		const scalarField& mag_grad_z_uncorrectedCells = mag_grad_z_uncorrected.internalField();
		#if OPENFOAM_VERSION >= 40
		scalarField& strain_rateCells = strain_rate.ref();
		#else
		scalarField& strain_rateCells = strain_rate.internalField();
		#endif      

		forAll(zMixCells, celli)
		{	
			if (zMixCells[celli] <= 1e-3 || zMixCells[celli] >= 0.999)
				 strain_rateCells[celli] = 0.;

			if (mag_grad_z_uncorrectedCells[celli] < 1.e-9)
				strain_rateCells[celli] = 0.;
		}
	}
	
	// Write
	strain_rate.write();
}
