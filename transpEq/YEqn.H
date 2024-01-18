tmp<fv::convectionScheme<scalar> > mvConvection (// ???
    fv::convectionScheme<scalar>::New (
        mesh, fields, phi, mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    clock_t t_begin_Y = std::clock();
    Info <<  "=== Start Species Transport ===" << endl;

    combustion->correct();
    dQ = combustion->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->mut()/Sct + turbulence->mu()/Sc, Yi)
              ==
                coalParcels.SYi(i, Yi)
              + combustion->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);

    clock_t t_end_Y = std::clock();
    double elapsed_secs_Y = double(t_end_Y-t_begin_Y)/CLOCKS_PER_SEC;
    Info << "=== Species Transport "
         << "Sc " << Sc << " Sct " << Sct
         << " [Done]: "
         << elapsed_secs_Y << " s ===\n\n";
}

