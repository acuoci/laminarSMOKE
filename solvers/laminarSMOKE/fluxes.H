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
	// Mass diffusion
	forAll (Y,i)
	{
		if (mwCorrectionInDiffusionFluxes == true)
		{
			volScalarField& Xi = X[i];
			volScalarField& Dmixi = Dmix[i];
			dimensionedScalar MWi("MWi", dimensionSet(1,0,0,0,-1,0,0),scalar(thermodynamicsMapXML->MW(i)) ); 

			J[i] = ( fvc::interpolate ( -rho*Dmixi*MWi/MWmix ) ) * ( fvc::interpolate ( fvc::grad(Xi) ) & mesh.Sf() );
		}
		else
		{
			volScalarField& Yi = Y[i];
			volScalarField& Dmixi = Dmix[i];

			J[i] = ( fvc::interpolate (-rho*Dmixi) ) * ( fvc::interpolate ( fvc::grad(Yi) ) & mesh.Sf() );
		}
	}

	// Soret effect
	forAll (Y,i)
	{
		if (soretEffect == true)
		{
			if (soretEffectList[i] == true)
			{
				J[i] += ( fvc::interpolate ( -rho*Dsoret[indexSoret[i]]/T ) ) *
					( fvc::interpolate (  fvc::grad(T) ) & mesh.Sf () ) ;
			}
		}
	}

	// Thermophoretic effect
	forAll (Y,i)
	{
		if (thermophoreticEffect == true)
		{
			if (thermophoreticEffectList[i] == true)
			{
				volScalarField& Yi = Y[i];

				J[i] += ( fvc::interpolate ( -0.55*mu/T*Yi ) ) *
					( fvc::interpolate ( fvc::grad(T) ) & mesh.Sf () ) ;
			}
		}
	}

	// Correction flux
	Jc *= 0.;
	forAll (Y,i)
	{
		Jc -= J[i];
	}

	// Corrected fluxes
	forAll (Y,i)
	{
		J[i] += Jc*fvc::interpolate(Y[i]);
	}
}
