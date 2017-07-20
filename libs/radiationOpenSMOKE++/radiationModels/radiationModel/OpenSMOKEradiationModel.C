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

#include "OpenSMOKEradiationModel.H"
#include "OpenSMOKEabsorptionEmissionModel.H"
#include "OpenSMOKEscatterModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(OpenSMOKEradiationModel, 0);
        defineRunTimeSelectionTable(OpenSMOKEradiationModel, T);
        defineRunTimeSelectionTable(OpenSMOKEradiationModel, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiation::OpenSMOKEradiationModel::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    #if DEVVERSION == 1
    if (io.typeHeaderOk<IOdictionary>(true))
    #else
    if (io.headerOk())
    #endif
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::radiation::OpenSMOKEradiationModel::initialise()
{
    if (radiation_)
    {
        solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

        absorptionEmission_.reset
        (
            OpenSMOKEabsorptionEmissionModel::New(*this, mesh_).ptr()
        );

        scatter_.reset(OpenSMOKEscatterModel::New(*this, mesh_).ptr());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::OpenSMOKEradiationModel::OpenSMOKEradiationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(NULL),
    scatter_(NULL)
{}


Foam::radiation::OpenSMOKEradiationModel::OpenSMOKEradiationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(lookupOrDefault("radiation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(NULL),
    scatter_(NULL)
{
    if (readOpt() == IOobject::NO_READ)
    {
        radiation_ = false;
    }

    initialise();
}


Foam::radiation::OpenSMOKEradiationModel::OpenSMOKEradiationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(lookupOrDefault("radiation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(NULL),
    scatter_(NULL)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::OpenSMOKEradiationModel::~OpenSMOKEradiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::OpenSMOKEradiationModel::read()
{
    if (regIOobject::read())
    {
        lookup("radiation") >> radiation_;
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::OpenSMOKEradiationModel::correct()
{
    if (!radiation_)
    {
        return;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }
}

Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::OpenSMOKEradiationModel::divq( volScalarField& T ) const
{
    return
    (
        Ru() - fvm::Sp(Rp()*pow3(T), T)
    );
}

void Foam::radiation::OpenSMOKEradiationModel::Qloss( volScalarField& T, volScalarField& Qloss)
{
	// TODO: Ru is an external source
	//       Thus it can be neglected in the expression below, unless the user
	//	 explicitly defined it
        Qloss = Rp()*pow4(T);	// - Ru();
}

/* 
// To be used for the enthalpy equation
Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::OpenSMOKEradiationModel::Sh
(
    fluidThermo& thermo
) const
{
    volScalarField& he = thermo.he();
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_ - 4.0*he/Cpv)
    );
}

*/

const Foam::radiation::OpenSMOKEabsorptionEmissionModel&
Foam::radiation::OpenSMOKEradiationModel::absorptionEmission() const
{
    if (!absorptionEmission_.valid())
    {
        FatalErrorIn
        (
            "const Foam::radiation::OpenSMOKEabsorptionEmissionModel&"
            "Foam::radiation::OpenSMOKEradiationModel::absorptionEmission() const"
        )
            << "Requested radiation absorptionEmission model, but model is "
            << "not activate" << abort(FatalError);
    }

    return absorptionEmission_();
}


// ************************************************************************* //
