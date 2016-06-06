/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::radiation::OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a grey-diffuse condition for radiative
    heat flux, \c Qr, for use with the view factor model

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        Qro          | external radiative heat flux | yes    |
        emissivityMode | emissivity mode: solidRadiation or lookup | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            greyDiffusiveRadiationViewFactor;
        Qro             uniform 0;
        emissivityMode  solidRadiation;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::radiationCoupledBase
    Foam::radiation::radiationModel
    Foam::radiation::viewFactor
    Foam::fixedValueFvPatchField

SourceFiles
    OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField_H
#define OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField_H

#include "radiationCoupledBase.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
Class OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField declaration
\*---------------------------------------------------------------------------*/

class OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField,
    public radiationCoupledBase
{
    // Private data

        //- External radiative heat flux
        scalarField Qro_;


public:

    //- Runtime type information
    TypeName("OpenSMOKEgreyDiffusiveRadiationViewFactor");


    // Constructors

        //- Construct from patch and internal field
        OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        // OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField onto a new patch
        OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new OpenSMOKEgreyDiffusiveViewFactorFixedValueFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return external radiative heat flux
            const scalarList& Qro()
            {
                return Qro_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
