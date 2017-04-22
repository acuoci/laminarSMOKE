/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "opticallyThin.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "OpenSMOKEabsorptionEmissionModel.H"
#include "OpenSMOKEscatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(opticallyThin, 0);
        addToRadiationRunTimeSelectionTables(opticallyThin);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::opticallyThin::opticallyThin(const volScalarField& T)
:
    OpenSMOKEradiationModel(typeName, T),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
    ambientTemperature_(readScalar(coeffs_.lookup("ambientTemperature")))
{}


Foam::radiation::opticallyThin::opticallyThin(const dictionary& dict, const volScalarField& T)
:
    OpenSMOKEradiationModel(typeName, dict, T),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
    ambientTemperature_(readScalar(coeffs_.lookup("ambientTemperature")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::opticallyThin::~opticallyThin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::opticallyThin::read()
{
    if (OpenSMOKEradiationModel::read())
    {
        // nothing to read

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::opticallyThin::calculate()
{
    a_ = absorptionEmission_->a();
    e_ = absorptionEmission_->e();
    E_ = absorptionEmission_->E();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::opticallyThin::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4.0*absorptionEmission_->aCont()*physicoChemical::sigma
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> > Foam::radiation::opticallyThin::Ru() const
{
    #if OPENFOAM_VERSION >= 40
    tmp<volScalarField> Tenv4
    (
        new volScalarField
        (
            IOobject
            (
                "Tenv4",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("Tenv4", dimensionSet(0, 0, 0, 4, 0), pow(ambientTemperature_,4.)),
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );
    #else
    tmp<volScalarField> Tenv4
    (
        new volScalarField
        (
            IOobject
            (
                "Tenv4",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("Tenv4", dimensionSet(0, 0, 0, 4, 0), pow(ambientTemperature_,4.)),
            zeroGradientFvPatchVectorField::typeName
        )
    );
    #endif

    #if OPENFOAM_VERSION >= 40
    	const DimensionedField<scalar, volMesh> E  = absorptionEmission_->ECont()()();
    	const DimensionedField<scalar, volMesh> a  = absorptionEmission_->aCont()()();
    	const DimensionedField<scalar, volMesh> T4 = Tenv4()();
    #else
   	const DimensionedField<scalar, volMesh> E  = absorptionEmission_->ECont()().dimensionedInternalField();
    	const DimensionedField<scalar, volMesh> a  = absorptionEmission_->aCont()().dimensionedInternalField();
    	const DimensionedField<scalar, volMesh> T4 = Tenv4().dimensionedInternalField();
    #endif

    return 4.0*a*physicoChemical::sigma*T4 -4.0*E;
}


// ************************************************************************* //
