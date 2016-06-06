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

Namespace
    Foam::radiation

Description
    Namespace for radiation modelling

Class
    Foam::radiation::OpenSMOKEradiationModel

Description
    Top level model for radiation modelling

SourceFiles
    OpenSMOKEradiationModel.C
    OpenSMOKEradiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef OpenSMOKEradiationModel_H
#define OpenSMOKEradiationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "fvMatricesFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace radiation
{

// Forward declaration of classes
class OpenSMOKEabsorptionEmissionModel;
class OpenSMOKEscatterModel;

/*---------------------------------------------------------------------------*\
                   Class OpenSMOKEradiationModel Declaration
\*---------------------------------------------------------------------------*/

class OpenSMOKEradiationModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the temperature field
        const volScalarField& T_;

        //- Radiation model on/off flag
        Switch radiation_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  radiation solver iteration
        label solverFreq_;

        //- Flag to enable radiation model to be evaluated on first iteration
        bool firstIter_;


        // References to the radiation sub-models

            //- Absorption/emission model
            autoPtr<OpenSMOKEabsorptionEmissionModel> absorptionEmission_;

            //- Scatter model
            autoPtr<OpenSMOKEscatterModel> scatter_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        OpenSMOKEradiationModel(const OpenSMOKEradiationModel&);

        //- Disallow default bitwise assignment
        void operator=(const OpenSMOKEradiationModel&);


public:

    //- Runtime type information
    TypeName("OpenSMOKEradiationModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            OpenSMOKEradiationModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            OpenSMOKEradiationModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& T
            ),
            (dict, T)
        );


    // Constructors

        //- Null constructor
        OpenSMOKEradiationModel(const volScalarField& T);

        //- Construct from components
        OpenSMOKEradiationModel(const word& type, const volScalarField& T);

        //- Construct from components
        OpenSMOKEradiationModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& T
        );


    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<OpenSMOKEradiationModel> New(const volScalarField& T);

        //- Return a reference to the selected radiation model
        static autoPtr<OpenSMOKEradiationModel> New
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~OpenSMOKEradiationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve radiation equation(s)
            virtual void calculate() = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Radiation model on/off flag
            const Switch radiation() const
            {
                return radiation_;
            }

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh> > Ru() const = 0;

	    //- Temperature source term
            virtual tmp<fvScalarMatrix> divq(volScalarField& T) const;

	    //- Qloss [W/m3]
	    virtual void Qloss( volScalarField& T, volScalarField& Qloss);

	    /*
            //- Energy source term
            virtual tmp<fvScalarMatrix> Sh(fluidThermo& thermo) const;
            */

            //- Access to absorptionEmission model
            const OpenSMOKEabsorptionEmissionModel& absorptionEmission() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToRadiationRunTimeSelectionTables(model)                           \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        OpenSMOKEradiationModel,                                              \
        model,                                                                \
        dictionary                                                            \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        OpenSMOKEradiationModel,                                              \
        model,                                                                \
        T                                                                     \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
