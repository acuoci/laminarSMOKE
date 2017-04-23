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

#include "speciesWallFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesWallFvPatchScalarField::speciesWallFvPatchScalarField
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
	epsilon() = 0.0;
	omega0()  = 0.0;
	rho0()    = 0.0;
	#if OPENFOAM_VERSION >= 40
	nameInternal_ = internalField().name();
	#else
	nameInternal_ = dimensionedInternalField().name();
	#endif
}


Foam::speciesWallFvPatchScalarField::speciesWallFvPatchScalarField
(
    const speciesWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedUserDefinedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::speciesWallFvPatchScalarField::speciesWallFvPatchScalarField
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
    omega0() = 0;
    rho0()   = 0;
    alfa()   = 0.;
    eta()    = 0.;

    // Calculating beta
    const double Dmix = 1e-10;
    beta() = Dmix*this->patch().deltaCoeffs();
 
    // Calculating epsilon
    epsilon() = 0;
 
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


Foam::speciesWallFvPatchScalarField::speciesWallFvPatchScalarField
(
    const speciesWallFvPatchScalarField& tppsf
)
:
    mixedUserDefinedFvPatchScalarField(tppsf)
{}


Foam::speciesWallFvPatchScalarField::speciesWallFvPatchScalarField
(
    const speciesWallFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedUserDefinedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speciesWallFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void Foam::speciesWallFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedUserDefinedFvPatchScalarField::rmap(ptf, addr);

//    const speciesWallFvPatchScalarField& tiptf =
//        refCast<const speciesWallFvPatchScalarField>(ptf);
}


void Foam::speciesWallFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Name of field
    #if OPENFOAM_VERSION >= 40
    nameInternal_ = internalField().name();
    #else
    nameInternal_ = dimensionedInternalField().name();
    #endif

    // Theve variables are kept equal to 0
    omega0() = 0.;
    rho0()   = 0.;
    eta()    = 0.;
    alfa()   = 0.;	
    
    // Index of patch
    const label patchi = patch().index();

    // Calculating beta    
    const volScalarField& Dmix = db().lookupObject<volScalarField>("gas::Dmix_" + nameInternal_);
    beta() = Dmix.boundaryField()[patchi]*this->patch().deltaCoeffs();
    
    // Calculating epsilon
    bool soretEffect = db().foundObject<volScalarField>("gas::Dsoret_" + nameInternal_);
    if (soretEffect == true)
    {
    	const volScalarField& Dsoret = db().lookupObject<volScalarField>("gas::Dsoret_" + nameInternal_);
    	const volScalarField& T = db().lookupObject<volScalarField>("T");
    	epsilon() = -T.boundaryField()[patchi].snGrad() / T.boundaryField()[patchi] * Dsoret.boundaryField()[patchi];
    }

    if (debug)
    {
    }

    mixedUserDefinedFvPatchScalarField::updateCoeffs();
}

void Foam::speciesWallFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        speciesWallFvPatchScalarField
    );
}

// ************************************************************************* //
