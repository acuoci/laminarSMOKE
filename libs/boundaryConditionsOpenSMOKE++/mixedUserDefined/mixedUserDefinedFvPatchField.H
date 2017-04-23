/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::mixedUserDefinedFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a base class for 'mixedUserDefined' type boundary
    conditions, i.e. conditions that mix fixed value and patch-normal gradient
    conditions.

    The respective contributions from each is determined by a weight field:

        \f[
            x_p = w x_p + (1-w) \left(x_c + \frac{\nabla_\perp x}{\Delta}\right)
        \f]

    where
    \vartable
        x_p   | patch values
        x_c   | patch internal cell values
        w     | weight field
        \Delta| inverse distance from face centre to internal cell centre
        w     | weighting (0-1)
    \endvartable


    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        valueFraction | weight field           | yes         |
        refValue     | fixed value             | yes         |
        refGrad      | patch normal gradient   | yes         |
    \endtable

Note
    This condition is not usually applied directly; instead, use a derived
    mixedUserDefined condition such as \c inletOutlet

SeeAlso
    Foam::inletOutletFvPatchField

SourceFiles
    mixedUserDefinedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mixedUserDefinedFvPatchField_H
#define mixedUserDefinedFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mixedUserDefinedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixedUserDefinedFvPatchField
:
    public fvPatchField<Type>
{
    // Private data
    scalarField alfa_;
    scalarField beta_;
    scalarField eta_;
    scalarField epsilon_;
    scalarField omega0_;
    scalarField rho0_;

public:

    //- Runtime type information
    TypeName("mixedUserDefined");


    // Constructors

        //- Construct from patch and internal field
        mixedUserDefinedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixedUserDefinedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given mixedUserDefinedFvPatchField onto a new patch
        mixedUserDefinedFvPatchField
        (
            const mixedUserDefinedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixedUserDefinedFvPatchField
        (
            const mixedUserDefinedFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixedUserDefinedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mixedUserDefinedFvPatchField
        (
            const mixedUserDefinedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixedUserDefinedFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }


        // Return defining fields

            virtual scalarField& alfa()
            {
                return alfa_;
            }

            virtual const scalarField& alfa() const
            {
                return alfa_;
            }

            virtual scalarField& beta()
            {
                return beta_;
            }

            virtual const scalarField& beta() const
            {
                return beta_;
            }

            virtual scalarField& eta()
            {
                return eta_;
            }

            virtual const scalarField& eta() const
            {
                return eta_;
            }

	    virtual scalarField& epsilon()
            {
                return epsilon_;
            }

            virtual const scalarField& epsilon() const
            {
                return epsilon_;
            }

            virtual scalarField& omega0()
            {
                return omega0_;
            }

            virtual const scalarField& omega0() const
            {
                return omega0_;
            }

            virtual scalarField& rho0()
            {
                return rho0_;
            }

            virtual const scalarField& rho0() const
            {
                return rho0_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
		#if OPENFOAM_VERSION >= 40
		const Pstream::commsTypes commsType =Pstream::commsTypes::blocking
		#else
                const Pstream::commsTypes commsType = Pstream::blocking
		#endif
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixedUserDefinedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
