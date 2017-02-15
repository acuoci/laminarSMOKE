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


if(speciesEquations == true)
{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    
    volScalarField Yt = 0.0*Y[0];

    for (label j=0; j<Y.size(); j++)
    {
	label i = species_order[j];

        if (i != inertIndex)
        {
           	volScalarField& Yi = Y[i];
	    	volScalarField& Dmixi = Dmix[i];
		dimensionedScalar Mi( "Mi", dimensionSet(1,0,0,0,-1,0,0),thermodynamicsMapXML->MW()[i+1] ); 

		// Corrections
		{
			sumDiffusionCorrections *= 0.;
			for (label k=0; k<Y.size(); k++)
			{
				volScalarField& Yk = Y[k];	
				dimensionedScalar Mk( "Mk", dimensionSet(1,0,0,0,-1,0,0), thermodynamicsMapXML->MW()[k+1] );

				if (i!=k)
				sumDiffusionCorrections += fvc::laplacian(rho*MWmix/Mk*Dmixi*Yi, Yk);
			}
		}

		fvScalarMatrix YiEqn
		(
		    	fvm::div(phi, Yi, "div(phi,Yi)")
		      - fvm::laplacian(rho*Dmixi, Yi) 
                        == 
	              - fvm::laplacian(rho*Dmixi*MWmix*Yi/Mi, Yi)
		      - sumDiffusionCorrections + 
                      - fvm::div(Jc,Yi, "div(Jc,Yi)") + 
		        sourceExplicit[i] + fvm::Sp(sourceImplicit[i],Yi) +						
		     	fvOptions(rho, Yi)
		);

		// Add Soret effect
		if (soretEffect == true)
		{ 
			if (soretEffectList[i] == true)
				YiEqn -= fvc::laplacian(rho*Dsoret[indexSoret[i]]/T, T, "laplacian(teta,Yi)");
		}

		// Add thermophoretic effect
		if (thermophoreticEffect == true)
		{
			if (thermophoreticEffectList[i] == true)
				YiEqn -= fvc::laplacian(0.55*mu/T*Yi, T, "laplacian(teta,Yi)");
		}

		// Solve
		YiEqn.relax();
		fvOptions.constrain(YiEqn);
		YiEqn.solve(mesh.solver(exceptional_species[i]));
		fvOptions.correct(Yi);	
			
		// Sum of mass fractions
            	Yi.max(0.0);
           	Yt += Yi;
        }
    }

    Info << "Inert species is " << Y[inertIndex].name() << " with local index equal to " << inertIndex << endl;
    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

    if (species_order_policy == SPECIES_ORDER_POLICY_SWEEP)
    	std::reverse(species_order.begin(),species_order.end());
    else if (species_order_policy == SPECIES_ORDER_POLICY_RANDOM_SHUFFLE)
    	std::random_shuffle ( species_order.begin(), species_order.end() );
    else if (species_order_policy == SPECIES_ORDER_POLICY_ROTATE)
	std::rotate(species_order.begin(),species_order.begin()+1,species_order.end());

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	
    Info << "Transport equations of species solved in " << tEnd - tStart << " s " << endl;
}
