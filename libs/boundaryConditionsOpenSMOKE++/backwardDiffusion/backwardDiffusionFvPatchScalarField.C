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

#include "backwardDiffusionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedUserDefinedFvPatchScalarField(p, iF)
{
	alfa()    = 0.0;
	beta()    = 0.0;
        eta()     = 0.0;
	omega0()  = 0.0;
	rho0()    = 0.0;
        epsilon() = 0.0;
	#if OPENFOAM_VERSION >= 40
	nameInternal_ = internalField().name();
	#else
        nameInternal_ = dimensionedInternalField().name();
	#endif
        
}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const backwardDiffusionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedUserDefinedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedUserDefinedFvPatchScalarField(p, iF)
{
    // Name of field
    #if OPENFOAM_VERSION >= 40
    nameInternal_ = internalField().name();
    #else
    nameInternal_ = dimensionedInternalField().name();
    #endif
    
    // Set the nominal value
    omega0() = scalarField("omega0", dict, p.size());
    rho0()   = scalarField("rho0", dict, p.size());

    // Fixed value condition is forced
    alfa() = 1000.;
    eta()  = 1.;

    // Calculating epsilon
    epsilon() = 0;

    // Calculating beta
    const double Dmix = 1e-10;
    beta() = Dmix*this->patch().deltaCoeffs();

    // Read value if available
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const backwardDiffusionFvPatchScalarField& tppsf
)
:
    mixedUserDefinedFvPatchScalarField(tppsf)
{}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const backwardDiffusionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedUserDefinedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::backwardDiffusionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void Foam::backwardDiffusionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedUserDefinedFvPatchScalarField::rmap(ptf, addr);

//    const backwardDiffusionFvPatchScalarField& tiptf =
//        refCast<const backwardDiffusionFvPatchScalarField>(ptf);
}


void Foam::backwardDiffusionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Calculating alfa
    const volVectorField& U = db().lookupObject<volVectorField>("U");
    tmp<vectorField> n = patch().nf();
    alfa() = -(n & U.boundaryField()[patchi]);
    
    // Calculating eta
    const volScalarField& rho = db().lookupObject<volScalarField>("rho");
    eta() = rho0() / rho.boundaryField()[patchi];

    #if OPENFOAM_VERSION >= 40
    nameInternal_ = internalField().name();
    #else
    nameInternal_ = dimensionedInternalField().name();
    #endif

    bool soretEffect  = db().foundObject<volScalarField>("gas_Dsoret_" + nameInternal_);

    // Calculating epsilon
    if (soretEffect == true)
    {
        const volScalarField& Dsoret = db().lookupObject<volScalarField>("gas_Dsoret_" + nameInternal_);
    	const volScalarField& T = db().lookupObject<volScalarField>("T");
    	epsilon() = -T.boundaryField()[patchi].snGrad() / T.boundaryField()[patchi] * Dsoret.boundaryField()[patchi];
    }
    else
    { 
	epsilon() = 0.;
    }

    // Calculating beta
    #if OPENFOAM_VERSION >= 40
    nameInternal_ = internalField().name();
    #else
    nameInternal_ = dimensionedInternalField().name();
    #endif
    const volScalarField& Dmix = db().lookupObject<volScalarField>("gas_Dmix_" + nameInternal_);
    beta() = Dmix.boundaryField()[patchi]*this->patch().deltaCoeffs();

    if (debug)
    {
    }

    mixedUserDefinedFvPatchScalarField::updateCoeffs();
}

void Foam::backwardDiffusionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    omega0().writeEntry("omega0", os);
    rho0().writeEntry("rho0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        backwardDiffusionFvPatchScalarField
    );
}

// ************************************************************************* //
