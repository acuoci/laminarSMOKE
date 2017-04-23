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

#include "depositionWallFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::depositionWallFvPatchScalarField::depositionWallFvPatchScalarField
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


Foam::depositionWallFvPatchScalarField::depositionWallFvPatchScalarField
(
    const depositionWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedUserDefinedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::depositionWallFvPatchScalarField::depositionWallFvPatchScalarField
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


Foam::depositionWallFvPatchScalarField::depositionWallFvPatchScalarField
(
    const depositionWallFvPatchScalarField& tppsf
)
:
    mixedUserDefinedFvPatchScalarField(tppsf)
{}


Foam::depositionWallFvPatchScalarField::depositionWallFvPatchScalarField
(
    const depositionWallFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedUserDefinedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::depositionWallFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void Foam::depositionWallFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedUserDefinedFvPatchScalarField::rmap(ptf, addr);

//    const depositionWallFvPatchScalarField& tiptf =
//        refCast<const depositionWallFvPatchScalarField>(ptf);
}


void Foam::depositionWallFvPatchScalarField::updateCoeffs()
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
    omega0()  = 0.;
    rho0()    = 0.;
    eta()     = 0.;
    alfa()    = 0.;
    epsilon() = 0.;	
    
    // Index of patch
    const label patchi = patch().index();

    // Calculating beta    
    const volScalarField& Dmix = db().lookupObject<volScalarField>("gas::Dmix_" + nameInternal_);
    beta() = Dmix.boundaryField()[patchi]*this->patch().deltaCoeffs();

    // Calculating alfa (only if thermophoretic effect is on)  
    if (nameInternal_.substr(0,3) == "BIN"  || nameInternal_.substr(0,3) == "bin")
    {  
	const volScalarField& mu = db().lookupObject<volScalarField>("gas::mu");
        const volScalarField& rho = db().lookupObject<volScalarField>("rho");
	const volScalarField& T = db().lookupObject<volScalarField>("T");
    	alfa() = 0.55*mu.boundaryField()[patchi]/rho.boundaryField()[patchi]/T.boundaryField()[patchi]*T.boundaryField()[patchi].snGrad();
    }
    
    // Calculating epsilon
    bool soretEffect  = db().foundObject<volScalarField>("gas::Dsoret_" + nameInternal_);
    if (soretEffect == true)
    {
	Info << nameInternal_ << " Soret " << endl;
    	const volScalarField& Dsoret = db().lookupObject<volScalarField>("gas::Dsoret_" + nameInternal_);
    	const volScalarField& T = db().lookupObject<volScalarField>("T");
    	epsilon() = -T.boundaryField()[patchi].snGrad() / T.boundaryField()[patchi] * Dsoret.boundaryField()[patchi];
    }

    if (debug)
    {
    }

    mixedUserDefinedFvPatchScalarField::updateCoeffs();
}

void Foam::depositionWallFvPatchScalarField::write(Ostream& os) const
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
        depositionWallFvPatchScalarField
    );
}

// ************************************************************************* //
