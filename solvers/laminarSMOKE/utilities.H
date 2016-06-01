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

/*
#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "energyJumpFvPatchScalarField.H"
#include "energyJumpAMIFvPatchScalarField.H"
*/
tmp<volScalarField> Rho(volScalarField& p, volScalarField& psi)
{
	return p*psi;
}

void normalizeMassFractions(double* omega, const label celli, const double massFractionsTol, const unsigned int n)
{
	const double eps = 1.e-3;

	double sumFractions = 0.;
    	for(int i=0; i < n; i++)
		sumFractions += omega[i];

	if (sumFractions > 1.+eps || sumFractions < 1.-eps)
	{
		Info << "WARNING: sum of mass-fractions = " << sumFractions << " in cell " << celli <<endl;

		if (time == 0)
		{
		    Info << "\nFATAL ERROR: sum of Yi is not 1" <<endl;
		    Info << "Check internal field on cell: " << celli <<endl;
		    abort();
		}    
	}
	
	for(int i=0; i < n; i++)
		omega[i] /= sumFractions;
}

void normalizeMassFractions(Eigen::VectorXd& omega_plus_additional_variables, const int N, const label celli, const double massFractionsTol)
{
    double sumFractions = 0.;
    for(int i=0; i < N; i++)	
	sumFractions += omega_plus_additional_variables(i);

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of mass-fractions = " << sumFractions << " in cell " << celli <<endl;
        
        if (time == 0)
        {
            Info << "\nFATAL ERROR: sum of Yi is not 1" <<endl;
            Info << "Check internal field on cell: " << celli <<endl;
            abort();
        }    
    }
	
    for(int i=0; i < N; i++)
	omega_plus_additional_variables(i) /= sumFractions;
}


void normalizeMassFractions(Eigen::VectorXd& omega_plus_temperature, const label celli, const double massFractionsTol)
{
	double sumFractions = 0.;
    for(int i=0; i < omega_plus_temperature.size()-1; i++)
            sumFractions += omega_plus_temperature(i);

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of mass-fractions = " << sumFractions << " in cell " << celli <<endl;

        //for(int i=0; i < omega_plus_temperature.size()-1; i++)
        //    cout << i << ")\t" << omega_plus_temperature(i) << endl;
        
        if (time == 0)
        {
            Info << "\nFATAL ERROR: sum of Yi is not 1" <<endl;
            Info << "Check internal field on cell: " << celli <<endl;
            abort();
        }    
    }
	
	 for(int i=0; i < omega_plus_temperature.size()-1; i++)
            omega_plus_temperature(i) /= sumFractions;
}

void normalizeMassFractions(OpenSMOKE::OpenSMOKEVectorDouble& omega_plus_temperature, const label celli, const double massFractionsTol)
{
    double sumFractions = 0.;
    for(int i=1; i <= omega_plus_temperature.Size()-1; i++)
            sumFractions += omega_plus_temperature[i];

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of mass-fractions = " << sumFractions << " in cell " << celli <<endl;

        //for(int i=1; i <= omega_plus_temperature.Size()-1; i++)
        //    cout << i << ")\t" << omega_plus_temperature[i] << endl;
        
        if (time == 0)
        {
            Info << "\nFATAL ERROR: sum of Yi is not 1" <<endl;
            Info << "Check internal field on cell: " << celli <<endl;
            abort();
        }    
    }
	
	 for(int i=1; i <= omega_plus_temperature.Size()-1; i++)
            omega_plus_temperature[i] /= sumFractions;
}
/*
Foam::wordList heBoundaryBaseTypes(volScalarField& T_)
{
    const volScalarField::GeometricBoundaryField& tbf = T_.boundaryField();

    wordList hbt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
        if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpFvPatchScalarField&>(tbf[patchi]);

            hbt[patchi] = pf.interfaceFieldType();
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpAMIFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpAMIFvPatchScalarField&>
                (
                    tbf[patchi]
                );

            hbt[patchi] = pf.interfaceFieldType();
        }
    }

    return hbt;
}


Foam::wordList heBoundaryTypes(volScalarField& T_)
{
    const volScalarField::GeometricBoundaryField& tbf = T_.boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientEnergyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedEnergyFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpAMIFvPatchScalarField::typeName;
        }
        else if (tbf[patchi].type() == "energyRegionCoupledFvPatchScalarField")
        {
            hbt[patchi] = "energyRegionCoupledFvPatchScalarField";
        }
    }

    return hbt;
}

void heBoundaryCorrection(volScalarField& h)
{
    volScalarField::GeometricBoundaryField& hbf = h.boundaryField();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}
*/
