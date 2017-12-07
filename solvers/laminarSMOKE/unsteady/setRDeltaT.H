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

{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    // Maximum flow Courant number
    scalar maxCo(readScalar(pimpleDict.lookup("maxCo")));

    // Maximum time scale
    scalar maxDeltaT(pimpleDict.lookupOrDefault<scalar>("maxDeltaT", GREAT));

    // Smoothing parameter (0-1) when smoothing iterations > 0
    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    // Damping coefficient (1-0)
    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    // Maximum change in cell temperature per iteration
    // (relative to previous value)
    scalar alphaTemp(pimpleDict.lookupOrDefault("alphaTemp", 0.05));

    // Maximum change in cell concentration per iteration
    // (relative to reference value)
    scalar alphaY(pimpleDict.lookupOrDefault("alphaY", 1.0));

    Info<< "Time scales min/max:" << endl;

    // Cache old reciprocal time scale field
    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Flow time scale
    {
        rDeltaT.ref() =
        (
            fvc::surfaceSum(mag(phi))()()
           /((2*maxCo)*mesh.V()*rho())
        );

        // Limit the largest time scale
        rDeltaT.max(1/maxDeltaT);

        Info<< "    Flow        = "
            << 1/gMax(rDeltaT.primitiveField()) << ", "
            << 1/gMin(rDeltaT.primitiveField()) << endl;
    }

    // Heat release rate time scale
    if (alphaTemp < 1)
    {
        volScalarField::Internal rDeltaTT
        (
            mag(Q)/(alphaTemp*rho*cp*T)
        );

        Info << "    Temperature = "
             << 1/(gMax(rDeltaTT.field()) + VSMALL) << ", "
             << 1/(gMin(rDeltaTT.field()) + VSMALL) << endl;

        rDeltaT.ref() = max(rDeltaT(), rDeltaTT);
    }

    // Reaction rate time scale
    if (alphaY < 1)
    {
        dictionary Yref(pimpleDict.subDict("Yref"));

        volScalarField::Internal rDeltaTY
        (
            IOobject
            (
                "rDeltaTY",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("rDeltaTY", rDeltaT.dimensions(), 0)
        );

        bool foundY = false;

// TODO
/*
        forAll(Y, i)
        {
            if (i != inertIndex)
            {
                volScalarField& Yi = Y[i];

                if (Yref.found(Yi.name()))
                {
                    foundY = true;
                    scalar Yrefi = readScalar(Yref.lookup(Yi.name()));

                    rDeltaTY.field() = max
                    (
                        mag
                        (
                            reaction->R(Yi)().source()
                           /((Yrefi*alphaY)*(rho*mesh.V()))
                        ),
                        rDeltaTY
                    );
                }
            }
        }

        if (foundY)
        {
            Info << "    Composition = "
                 << 1/(gMax(rDeltaTY.field()) + VSMALL) << ", "
                 << 1/(gMin(rDeltaTY.field()) + VSMALL) << endl;

            rDeltaT.ref() = max(rDeltaT(), rDeltaTY);
        }
        else
        {
            IOWarningIn(args.executable().c_str(), Yref)
                << "Cannot find any active species in Yref " << Yref
                << endl;
        }
*/
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    // Spatially smooth the time scale field
    if (rDeltaTSmoothingCoeff < 1)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (scalar(1.0) - rDeltaTDampingCoeff)*rDeltaT0
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "    Overall     = "
        << 1/gMax(rDeltaT.primitiveField())
        << ", " << 1/gMin(rDeltaT.primitiveField()) << endl;
}

// ************************************************************************* //
