Info<< "Creating field number of active species\n" << endl;
volScalarField drg_important_species
(
    IOobject
    (
        "drg_important_species",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputDRG
    ),
    mesh,
    dimensionedScalar("drg_important_species", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField drg_important_reactions
(
    IOobject
    (
        "drg_important_reactions",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputDRG
    ),
    mesh,
    dimensionedScalar("drg_important_reactions", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField cpu_reaction_times
(
    IOobject
    (
        "cpu_reaction_times",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        outputDRG
    ),
    mesh,
    dimensionedScalar("cpu_reaction_times", dimensionSet(0, 0, 1, 0, 0), 0.0)
);
