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

\*---------------------------------------------------------------------------*/

#include "OpenSMOKEabsorptionEmissionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(OpenSMOKEabsorptionEmissionModel, 0);
        defineRunTimeSelectionTable(OpenSMOKEabsorptionEmissionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::OpenSMOKEabsorptionEmissionModel::OpenSMOKEabsorptionEmissionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::OpenSMOKEabsorptionEmissionModel::~OpenSMOKEabsorptionEmissionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::a(const label bandI) const
{
    return aDisp(bandI) + aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::aCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "aCont",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::aDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "aDisp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::e(const label bandI) const
{
    return eDisp(bandI) + eCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::eCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "eCont",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::eDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "eDisp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::E(const label bandI) const
{
    return EDisp(bandI) + ECont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::ECont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ECont",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::OpenSMOKEabsorptionEmissionModel::EDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "EDisp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}


Foam::label Foam::radiation::OpenSMOKEabsorptionEmissionModel::nBands() const
{
    return pTraits<label>::one;
}


const Foam::Vector2D<Foam::scalar>&
Foam::radiation::OpenSMOKEabsorptionEmissionModel::bands(const label n) const
{
    return Vector2D<scalar>::one;
}


bool Foam::radiation::OpenSMOKEabsorptionEmissionModel::isGrey() const
{
    return false;
}


void Foam::radiation::OpenSMOKEabsorptionEmissionModel::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aj
) const
{
    a = this->a();
    aj[0] =  a;
}


// ************************************************************************* //
