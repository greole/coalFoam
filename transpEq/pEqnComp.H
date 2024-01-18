// rho = thermo.rho();1        // update rho from thermodynamics, moved out of piso loop as suggested byJarczy 2013
clock_t t_begin_P = std::clock();
Info <<  "=== Start Pressure Corrector ===" << endl;

volScalarField rAU(1.0/UEqn.A());
U = rAU*UEqn.H();

phi = fvc::interpolate(rho)*(           // Kg/s  
    (fvc::interpolate(U) & mesh.Sf())
  + fvc::ddtPhiCorr(rAU, rho, U, phi)   // ddtPhiCor for backwardScheme 
);                                      // defined by backwardDdtScheme.C
                                        // L 670
for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++) {
    fvScalarMatrix pdEqn (       // assume a constant pressure (p) since psi*p=rho neglect effect of p on rho
            fvc::ddt(rho)        // ContiErrs: sum local = ~e-15, global = ~e-17, cumulative = 7.59582e-15 (1e-05,100 TS)
          //  fvc::ddt(psi,p)    // ContiErrs: Temp convection test: sum local ~e-14, global = ~e-15 Final  3.69453e-13 (1e-05,100 TS) 
          //  fvm::ddt(psi,pd)   // unstable
          // - psi*fvm::ddt(pd)  // unstable
          + fvc::div(phi)
          - fvm::laplacian(rho*rAU, p)
          == coalParcels.Srho()
    );

    pdEqn.solve (
        mesh.solver(p().select(pimple.finalInnerIter(corr, nonOrth)))
    );
    if (nonOrth == pimple.nNonOrthCorr()) {  // update flux in last iteration
            phi += pdEqn.flux();             // fvMatrix.C L.830 
    }
}

#include "rhoEqn_Sp.H"
#include "compressibleContinuityErrs.H"

// if (pdLimit) {
//     if (DEBUG) {
//         Info << "max pd: "
//              << max(p()).value()
//              << "min pd: "
//              << min(p()).value()
//              << endl;
//     }
//     p().max(pdLowerLimit);
//     p().min(pdUpperLimit);
// }

U -= rAU*fvc::grad(p);     // correct u with dynamic pressure
U.correctBoundaryConditions();
//DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), pd);

clock_t t_end_P = std::clock();
double elapsed_secs_P = double(t_end_P - t_begin_P) / CLOCKS_PER_SEC;
Info <<  "=== Pressure Corrector [Done]: " << elapsed_secs_P << " s ===\n" << endl;
