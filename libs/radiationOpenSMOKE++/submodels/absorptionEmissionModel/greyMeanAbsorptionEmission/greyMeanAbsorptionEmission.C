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

#include "greyMeanAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyMeanAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            OpenSMOKEabsorptionEmissionModel,
            greyMeanAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyMeanAbsorptionEmission::greyMeanAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    OpenSMOKEabsorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    speciesNames_(0),
    specieIndex_(label(0)),
    lookUpTablePtr_(),
    thermo_(mesh.lookupObject<laminarSMOKEthermoClass>("laminarSMOKEthermoClass")),
    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff"))),
    Yj_(nSpecies_)
{
    label nFunc = 0;
    const dictionary& functionDicts = dict.subDict(typeName + "Coeffs");

    forAllConstIter(dictionary, functionDicts, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();
        speciesNames_.insert(key, nFunc);
        const dictionary& dict = iter().dict();
        coeffs_[nFunc].initialise(dict);
        nFunc++;
    }

    if (coeffsDict_.found("lookUpTableFileName"))
    {
        const word name = coeffsDict_.lookup("lookUpTableFileName");
        if (name != "none")
        {
            lookUpTablePtr_.set
            (
                new interpolationLookUpTable<scalar>
                (
                    fileName(coeffsDict_.lookup("lookUpTableFileName")),
                    mesh.time().constant(),
                    mesh
                )
            );

            if (!mesh.foundObject<volScalarField>("ft"))
            {
                FatalErrorIn
                (
                    "Foam::radiation::greyMeanAbsorptionEmission(const"
                    "dictionary& dict, const fvMesh& mesh)"
                )   << "specie ft is not present to use with "
                    << "lookUpTableFileName " << nl
                    << exit(FatalError);
            }
        }
    }

    // Check that all the species on the dictionary are present in the
    // look-up table and save the corresponding indices of the look-up table

    label j = 0;
    forAllConstIter(HashTable<label>, speciesNames_, iter)
    {
        if (!lookUpTablePtr_.empty())
        {
            if (lookUpTablePtr_().found(iter.key()))
            {
                label index = lookUpTablePtr_().findFieldIndex(iter.key());

                Info<< "specie: " << iter.key() << " found on look-up table "
                    << " with index: " << index << endl;

                specieIndex_[iter()] = index;
            }
            else if (mesh.foundObject<volScalarField>(iter.key()))
            {
                volScalarField& Y =
                    const_cast<volScalarField&>
                    (
                        mesh.lookupObject<volScalarField>(iter.key())
                    );
                Yj_.set(j, &Y);
                specieIndex_[iter()] = 0;
                j++;
                Info<< "specie: " << iter.key() << " is being solved" << endl;
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::radiation::greyMeanAbsorptionEmission(const"
                    "dictionary& dict, const fvMesh& mesh)"
                )   << "specie: " << iter.key()
                    << " is neither in look-up table: "
                    << lookUpTablePtr_().tableName()
                    << " nor is being solved" << nl
                    << exit(FatalError);
            }
        }
        else if (mesh.foundObject<volScalarField>(iter.key()))
        {
            volScalarField& Y =
                const_cast<volScalarField&>
                (
                    mesh.lookupObject<volScalarField>(iter.key())
                );

            Yj_.set(j, &Y);
            specieIndex_[iter()] = 0;
            j++;
        }
        else
        {
            FatalErrorIn
            (
                "Foam::radiation::greyMeanAbsorptionEmission(const"
                "dictionary& dict, const fvMesh& mesh)"
            )   << " there is not lookup table and the specie" << nl
                << iter.key() << nl
                << " is not found " << nl
                << exit(FatalError);

        }
    }

	gas_correction_coefficient_  = coeffsDict_.lookupOrDefault<scalar>(word("gasCorrectionCoefficient"),  scalar(1.));
	Info << "Gas correction coefficient: " << gas_correction_coefficient_ << endl;

    	word soot_radiation(coeffsDict_.lookup("sootRadiation"));

	if (soot_radiation == "none")
		soot_planck_coefficient_ = SOOT_RADIATION_PLANCK_COEFFICIENT_NONE;
	else if (soot_radiation == "Smooke")
		soot_planck_coefficient_ = SOOT_RADIATION_PLANCK_COEFFICIENT_SMOOKE;
	else if (soot_radiation == "Kent")
		soot_planck_coefficient_ = SOOT_RADIATION_PLANCK_COEFFICIENT_KENT;
	else if (soot_radiation == "Sazhin")
		soot_planck_coefficient_ = SOOT_RADIATION_PLANCK_COEFFICIENT_SAZHIN;
	else
	{
		FatalErrorIn( "Foam::radiation::greyMeanAbsorptionEmission::greyMeanAbsorptionEmission")
	    		<< "Wrong sootRadiation model. Available models: none | Smooke | Kent | Sazhin " << abort(FatalError);
	}

	soot_correction_coefficient_ = coeffsDict_.lookupOrDefault<scalar>(word("sootCorrectionCoefficient"), scalar(1.));
	Info << "Soot radiation model: " << soot_radiation << endl;
	Info << " - Correction coeff.: " << soot_correction_coefficient_ << endl;

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyMeanAbsorptionEmission::~greyMeanAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::aCont(const label bandI) const
{
    const laminarSMOKEthermoClass& mixture = dynamic_cast<const laminarSMOKEthermoClass&>(thermo_);

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();


    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    #if OPENFOAM_VERSION >= 40
    	scalarField& a = ta.ref().primitiveFieldRef();
    #else
	scalarField& a = ta().internalField();
    #endif
    
    forAll(a, cellI)
    {
        forAllConstIter(HashTable<label>, speciesNames_, iter)
        {
            label n = iter();
            scalar Xipi = 0.0;
            if (specieIndex_[n] != 0)
            {
                //Specie found in the lookUpTable.
                const volScalarField& ft =
                    mesh_.lookupObject<volScalarField>("ft");

                const List<scalar>& Ynft = lookUpTablePtr_().lookUp(ft[cellI]);
                //moles x pressure [atm]
                Xipi = Ynft[specieIndex_[n]]*paToAtm(p[cellI]);
            }
            else
            {
                scalar invWt = 0.0;
                forAll (mixture.Y(), s)
                {
                    invWt += mixture.Y(s)[cellI]/mixture.W(s);
                }

		label index = mixture.species_index(iter.key());
                scalar Xk = mixture.Y(index)[cellI]/(mixture.W(index)*invWt);

                Xipi = Xk*paToAtm(p[cellI]);
            }

            const absorptionCoeffs::coeffArray& b = coeffs_[n].coeffs(T[cellI]);

            scalar Ti = T[cellI];
            // negative temperature exponents
            if (coeffs_[n].invTemp())
            {
                Ti = 1.0/T[cellI];
            }
            a[cellI] +=
                Xipi
               *(
                    ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                  + b[0]
                );
        }

	a[cellI] *= gas_correction_coefficient_;
    }

    // Soot contribution
    if (soot_planck_coefficient_ != SOOT_RADIATION_PLANCK_COEFFICIENT_NONE)
    {
	const scalarField& fvsoot = T.db().objectRegistry::lookupObject<volScalarField>("soot_fv_large").internalField();

	forAll(a, cellI)
    	{
		const scalar fvi = fvsoot[cellI];
		const scalar Ti = T[cellI];

		double aSoot = 0.;
															
		if (soot_planck_coefficient_ == SOOT_RADIATION_PLANCK_COEFFICIENT_SMOOKE)
			aSoot = 1307.*fvi*Ti;						// [1/m]	(Smooke et al. Combustion and Flame 2009)
		else if (soot_planck_coefficient_ == SOOT_RADIATION_PLANCK_COEFFICIENT_KENT)
			aSoot = 2262.*fvi*Ti;						// [1/m]	(Kent al. Combustion and Flame 1990)
		else if (soot_planck_coefficient_ == SOOT_RADIATION_PLANCK_COEFFICIENT_SAZHIN)
			aSoot= 1232.*fvi*(1. + 4.8e-4*(Ti - 2000.));			// [1/m]	(Sazhin, Fluent 1994)	

		a[cellI] += aSoot*soot_correction_coefficient_;
        }
    }

    #if OPENFOAM_VERSION >= 40
    	ta.ref().correctBoundaryConditions();
    #else
	ta().correctBoundaryConditions();
    #endif

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::eCont(const label bandI) const
{
   return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    if (mesh_.foundObject<volScalarField>("dQ"))
    {
        const volScalarField& dQ = mesh_.lookupObject<volScalarField>("dQ");

    #if OPENFOAM_VERSION >= 40
    	if (dQ.dimensions() == dimEnergy/dimTime)
        {
            E.ref().primitiveFieldRef() = EhrrCoeff_*dQ/mesh_.V();
        }
        else if (dQ.dimensions() == dimEnergy/dimTime/dimVolume)
        {
            E.ref().primitiveFieldRef() = EhrrCoeff_*dQ;
        }
        else
        {
            if (debug)
            {
                WarningIn
                (
                    "tmp<volScalarField>"
                    "radiation::greyMeanAbsorptionEmission::ECont"
                    "("
                        "const label"
                    ") const"
                )
                    << "Incompatible dimensions for dQ field" << endl;
            }
        }
    #else
	if (dQ.dimensions() == dimEnergy/dimTime)
        {
            E().internalField() = EhrrCoeff_*dQ/mesh_.V();
        }
        else if (dQ.dimensions() == dimEnergy/dimTime/dimVolume)
        {
            E().internalField() = EhrrCoeff_*dQ;
        }
        else
        {
            if (debug)
            {
                WarningIn
                (
                    "tmp<volScalarField>"
                    "radiation::greyMeanAbsorptionEmission::ECont"
                    "("
                        "const label"
                    ") const"
                )
                    << "Incompatible dimensions for dQ field" << endl;
            }
        }
    #endif
        
    }
    else
    {
        WarningIn
        (
            "tmp<volScalarField>"
            "radiation::greyMeanAbsorptionEmission::ECont"
            "("
                "const label"
            ") const"
        ) << "dQ field not found in mesh" << endl;
    }

    return E;
}


// ************************************************************************* //
