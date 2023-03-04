scalar timeFlag(readScalar(runTime.controlDict().lookup("timeFlag")));
Switch adjustFluxOut(runTime.controlDict().lookup("adjustFluxOut"));

Info << "\nReading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p_init
(
    IOobject
    (
        "p_init",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

#include "compressibleCreatePhi.H"

// Used for limiting debilitating spikes in rho
IOdictionary filterProperties
(   
    IOobject
    (
        "filterProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )   
);      
 
//autoPtr<LESfilter> filterTestPtr(LESfilter::New(U.mesh(), filterProperties));
//LESfilter& filterTest(filterTestPtr());

// Dummy fields declared just to use the turbulence model
volScalarField Zmix
(
    IOobject
    (   
        "Zmix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Zmix", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

volScalarField& alpha = mu;

// Create turbulence model
autoPtr<dynCompSansThermoSF::turbulenceModel> turbulence
(
    dynCompSansThermoSF::turbulenceModel::New(rho, mu, alpha, U, Zmix, phi)
);

// Used in adjustPhi()
volScalarField dmdt
(
    IOobject
    (
        "dmdt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dmdt", dimensionSet(1,0,-1,0,0), 0.0)
);

volScalarField sqrt_rhoTimeInterp
(
    IOobject
    (
        "sqrt_rhoTimeInterp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    sqrt(rho)

);

volScalarField sqrt_rhoTimeInterp_timeInterp
(
    IOobject
    (
        "sqrt_rhoTimeInterp_timeInterp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    sqrt(rho)
);

volScalarField rhoTimeInterp
(
    IOobject
    (
        "rhoTimeInterp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rho
);

// Blending factors
surfaceScalarField UBlendingFactor
(
    IOobject
    (
        "UBlendingFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::mag(mesh.Sf())/Foam::mag(mesh.Sf())*0.0
);

surfaceScalarField U_0BlendingFactor
(
    IOobject
    (
        "U_0BlendingFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::mag(mesh.Sf())/Foam::mag(mesh.Sf())*0.0
);

