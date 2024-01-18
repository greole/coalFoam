{
    clock_t t_begin_h = std::clock();
    Info <<  "=== Start Sensible Enthalpy Transport ====" << endl;
    auto dt = runTime.deltaT();

    volScalarField& he = thermo.he();

    scalar DtK = 1.0;
    if (NEGLECT_DtK) {DtK = 0.0;}

    scalar DH = 1.0;
    if (!ENTHCORR) {DH = 0.0;}

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + (fvc::ddt(rho, K) + fvc::div(phi, K)) * DtK
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p(),
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
      + fvc::laplacian(turbulence->alphaEff()*T, thermo.Cp()) * DH
   // - fvc::laplacian(turbulence->alphaEff()*thermo.Cp(), T)
     ==
        combustion->Sh()
      + coalParcels.Sh(he)
      + radiation->Sh(thermo)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

    Info << "T gas min/max  = "
         << min(T).value() << ", "
         << max(T).value() << endl;

    clock_t t_end_h = std::clock();
    double elapsed_secs_h = double(t_end_h - t_begin_h) / CLOCKS_PER_SEC;

    Info <<  "=== Sensible Enthalpy Transport [Done]: "
         << " exp diff corrector " << ENTHCORR
         << " [Done]: "
         << elapsed_secs_h
         << " s ====\n"
         << endl;
}

