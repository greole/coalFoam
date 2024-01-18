// update rho from thermodynamics, moved out of piso loop as suggested by Jarczy 2013
// since thermo is not updated during PISO loop ?
rho = thermo.rho();
clock_t t_begin_P = std::clock();
Info <<  "=== Start Pressure Corrector ===" << endl;

volScalarField rAU(1.0/UEqn.A());
U = rAU*UEqn.H();

surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        (fvc::interpolate(rho*HbyA) & mesh.Sf())
      + rhorAUf*fvc::ddtCorr(rho, U, phi)
    )
);

fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p())
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, p())
     ==
        coalParcels.Srho()
      + fvOptions(psi, p(), rho.name())
    );

    fvOptions.constrain(pEqn);

    pEqn.solve(mesh.solver(p().select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter()) {phi = phiHbyA + pEqn.flux();}
}

#include "rhoEqn_Sp.H"
#include "compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(p());
U.correctBoundaryConditions();
fvOptions.correct(U);

K = 0.5*magSqr(U);

if (thermo.dpdt()) {dpdt = fvc::ddt(p());}

clock_t t_end_P = std::clock();
double elapsed_secs_P = double(t_end_P - t_begin_P) / CLOCKS_PER_SEC;
Info <<  "=== Pressure Corrector [Done]: " << elapsed_secs_P << " s ===\n" << endl;
