Info << "Creating combustion model" << endl;

autoPtr<combustionModels::psiCombustionModel> combustion
(
    combustionModels::psiCombustionModel::New(mesh)
);

psiReactionThermo& thermo = combustion->thermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

Info << "Creating Mixture" << endl;

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));

if (!composition.contains(inertSpecie)) {
    FatalErrorIn(args.executable())
        << "Specified inert specie '" << inertSpecie << "' not found in "
        << "species list. Available species:" << composition.species()
        << exit(FatalError);
}

// CUSTOM PISO SELECTION
bool FULL_RHO_TRANSP(runTime.controlDict().lookupOrDefault("FULL_RHO_TRANSP", true));
bool NEGLECT_DtK(runTime.controlDict().lookupOrDefault("NEGLECT_DtK", false));
bool COMPPISO(runTime.controlDict().lookupOrDefault("CompressiblePISO", false));
bool ENTHCORR(runTime.controlDict().lookupOrDefault("EnthalpyCorrector", false));
bool BUOYANCY(runTime.controlDict().lookupOrDefault("BUOYANCY", false));


if (BUOYANCY) {
    Info << "BUOYANCY MODE ON" << endl;
}

autoPtr<volScalarField> p;
if(COMPPISO) {
    Info << "Compressible PISO MODE ON" << endl;
    // set it to a const. pressure for thermo calculations
    autoPtr<volScalarField> p_file(new volScalarField (
        IOobject (
            "p", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    ));
    p = p_file;
} else {
    Info << "Incompressible PISO MODE ON" << endl;
    // set it to a dyn. pressure for piso calculations
    autoPtr<volScalarField> p_file(new volScalarField(
        IOobject (
           "pd", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::AUTO_WRITE
        ),
        mesh
    ));
    p = p_file;
    // create separate constant pressure for
    // the equation of state
    volScalarField const_p (
        IOobject (
            "p", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("p", dimPressure, 1e5)
    );
    thermo.p() = const_p;
}

// volScalarField& p = thermo.p();
const auto& T   = thermo.T();
const auto& psi = thermo.psi();

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i) {fields.add(Y[i]);}
fields.add(thermo.he());

Info<< "Reading field rho" << endl;
volScalarField rho (
    IOobject (
        "rho", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    thermo.rho()
);


volScalarField thermoRhoOld (
    IOobject (
        "rho", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    thermo.rho()
);


Info<< "Reading field U" << endl;
volVectorField U (
    IOobject (
        "U", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence (
    compressible::turbulenceModel::New (rho, U, phi, thermo)
);

// new in OF24
// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());

volScalarField mu (
    IOobject (
        "mu", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    turbulence->mu()
);


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p().dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

bool DEBUG = runTime.controlDict().lookupOrDefault("Debug",false);
Info << "Debug Mode: " << DEBUG << endl;

clock_t t_end_U, t_begin_U;
double elapsed_secs_U;

IOdictionary LesCoeffs (
    IOobject  (
        "constant/LESProperties",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const scalar Sc(LesCoeffs.lookupOrDefault<scalar>("Sc", 0.7));
const scalar Sct(LesCoeffs.lookupOrDefault<scalar>("Sct", 0.4));
const bool setPrandtl = runTime.controlDict().lookupOrDefault("setPrandtl", false);
const scalar Pr(LesCoeffs.lookupOrDefault<scalar>("Pr", 0.7));
const scalar Prt(LesCoeffs.lookupOrDefault<scalar>("Prt", 0.4));
