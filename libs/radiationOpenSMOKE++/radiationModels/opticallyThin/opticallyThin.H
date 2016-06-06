/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::radiation::opticallyThin

SourceFiles
    opticallyThin.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelopticallyThin_H
#define radiationModelopticallyThin_H

#include "OpenSMOKEradiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class opticallyThin Declaration
\*---------------------------------------------------------------------------*/

class opticallyThin
:
    public OpenSMOKEradiationModel
{
    // Private data

        //- Total radiative heat flux [W/m2]
        volScalarField Qr_;

        //- Absorption coefficient
        volScalarField a_;

        //- Emission coefficient
        volScalarField e_;

        //- Emission contribution
        volScalarField E_;

	//- Ambient temperature
	scalar ambientTemperature_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        opticallyThin(const opticallyThin&);

        //- Disallow default bitwise assignment
        void operator=(const opticallyThin&);


public:

    //- Runtime type information
    TypeName("opticallyThin");


    // Constructors

        //- Construct from components
        opticallyThin(const volScalarField& T);

        //- Construct from components
        opticallyThin(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~opticallyThin();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh> > Ru() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
