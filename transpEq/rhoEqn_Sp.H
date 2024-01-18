/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    rhoEqn

Description
    Solve the continuity for density.

\*---------------------------------------------------------------------------*/

{
    clock_t t_begin_rho = std::clock();

    if (FULL_RHO_TRANSP) {
    Info << "\n=== Start FULL Conti Equation ===" << endl;
        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho)
          + fvc::div(phi)
          ==
            coalParcels.Srho(rho)
          + fvOptions(rho)
        );

        fvOptions.constrain(rhoEqn);

        rhoEqn.solve();

        fvOptions.correct(rho);
    } else {
    Info << "\n=== Start PART Conti Equation ===" << endl;
        rho.internalField() +=
               coalParcels.Srho(rho)().source()
               * runTime.deltaTValue();
    }

    clock_t t_end_rho = std::clock();
    double elapsed_secs_rho = double(t_end_rho - t_begin_rho) / CLOCKS_PER_SEC;
    Info << "=== Conti Equation [Done]:" << elapsed_secs_rho << " s ===\n" << endl;
}
// ************************************************************************* //
