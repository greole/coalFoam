    t_begin_U = std::clock();
    Info <<  "\n=== Start Momentum predictor ===" << endl;

    scalar BT = 1.0;
    if (!ENTHCORR) {BT = 0.0;}
    dimensionedVector g2(g);
    g2[2] = dimensionedScalar("g2", g.dimensions(), 9.81);

    fvVectorMatrix UEqn (
        // Programmers Guide P 37 - ddt(rho,U) =
        // \frac{\partial \overbar{\rho}\tilde{U}}{\partial t}
        fvm::ddt(rho, U)
        // phi units [1 0 -1 0 0 0 0] = kg/s
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
      // + BT*g2*rho.dimensionedInternalField()*(1.0-rho.dimensionedInternalField()/dimensionedScalar("rho0",dimDensity,1.205))
      + BT*g2*(rho.dimensionedInternalField()-dimensionedScalar("rho0",dimDensity,1.205))
      + coalParcels.SU(U)
      + fvOptions(rho, U)
    );

    fvOptions.constrain(UEqn);

    // predictor step u^*(p^n), using the dynamic pressure
    // fvc::grad(pd) uses conservative gauss intergration Programmers guide P36
    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p()));

        fvOptions.correct(U);

        K = 0.5*magSqr(U);
    }

    // update mag strain rate field
    // volTensorField gradU = fvc::grad(U);
    // S = mag(symm(gradU))*Foam::sqrt(2.0);

    t_end_U = std::clock();
    elapsed_secs_U = double(t_end_U - t_begin_U) / CLOCKS_PER_SEC;

    Info <<  "=== Momentum predictor [Done]: " << elapsed_secs_U << " s ===\n" << endl;
