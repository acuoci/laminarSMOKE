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

volScalarField rho
(
	IOobject
	(
		"rho",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(1, -3, 0, 0, 0)
);

{
	unsigned int ns = Y.size();
	OpenSMOKE::OpenSMOKEVectorDouble y(ns);
	OpenSMOKE::OpenSMOKEVectorDouble x(ns);
	Eigen::VectorXd y_eigen(ns);
	Eigen::VectorXd x_eigen(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	      scalarField& rhoCells = rho.internalField();

	forAll(TCells, celli)
	{	
		double mw;

		// Extract the mean mass fractions
		for(unsigned int i=0;i<ns;i++)
			y[i+1] = Y[i].internalField()[celli];
		const double sum = y.SumElements();
		for(unsigned int i=0;i<ns;i++)
			y[i+1] /= sum;
		

		// Molecular weight
		thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

		// Eigen vectors
		for(unsigned int i=0;i<ns;i++)
		{
			y_eigen(i) = y[i+1];
			x_eigen(i) = x[i+1];
		}

		// Density
		rhoCells[celli] = pCells[celli]*mw/PhysicalConstants::R_J_kmol/TCells[celli];
	}


	forAll(T.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT  = T.boundaryField()[patchi];
		const fvPatchScalarField& pp  = p.boundaryField()[patchi];
		      fvPatchScalarField& prho  = rho.boundaryField()[patchi];


		forAll(pT, facei)
		{
			double mw;

			// Extract the mass fractions
			for(unsigned int i=0;i<ns;i++)
				y[i+1] = Y[i].boundaryField()[patchi][facei];
			const double sum = y.SumElements();
			for(unsigned int i=0;i<ns;i++)
				y[i+1] /= sum;

			// Molecular weight
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

			// Eigen vectors
			for(unsigned int i=0;i<ns;i++)
			{
				y_eigen(i) = y[i+1];
				x_eigen(i) = x[i+1];
			}

			// Density
			prho[facei] = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];
		}
	}
}
