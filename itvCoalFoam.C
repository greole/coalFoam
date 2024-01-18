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

Application
    coalChemistryFoam

Description
    Transient solver for:
    - compressible,
    - turbulent flow,
    with
    - coal and limestone parcel injections,
    - energy source, and
    - combustion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"

// LEGACY INCLUDES
//#include "hCombustionThermo.H"
//#include "LESModel.H"
//#include "psiChemistryModel.H"
//#include "chemistrySolver.H"
//#include "reactingMixture.H"

// OF240 INCLUDES
#include "basicThermoCloud.H"
#include "coalCloud.H"
#include "psiCombustionModel.H"
#include "fvIOoptionList.H"

#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
// // CUSTOM INCLUDES
#include "ListListOps.H"
#include <ctime>
// // CUSTOM CLASSES
#include "particle.H"
#include "coalParcel.H"
#include "probesContainer.H"
#include "particleProbe.H"
//#include "EBU.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[]) {
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    //#include "readChemistryProperties.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    // Next is new in OF24
    #include "createFvOptions.H"
    #include "createClouds.H"
    #include "createRadiationModel.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info << "\n---- Starting time loop ----\n" << endl;

    while (runTime.run()) {
        auto cpu0   = runTime.elapsedCpuTime();
        auto clock0 = runTime.elapsedClockTime();

        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info << "Time = " << runTime.timeName()
             << "\nTime index: " << runTime.timeIndex()
             <<  nl << endl;

        // Removed in OF24?
        // #include "chemistry.H"     // store reaction rates,
        #include "evolveParticles.H"
        #include "sampleParticles.H"

        // --- Pressure-velocity PIMPLE corrector loop
        #include "rhoEqn_Sp.H" // transport \tilda{rho^n} due to
                               // \tilda{phi^n}, particle mass sources

        while (pimple.loop()) {
            #include "UEqn.H"      // mom. pred. u* = f(rho^n,)

            // Note: in std piso formulation Y and hsEqn would be
            // solved after pressure corretor, however since Poisson
            // Eqn depends on ddt(rho) T,rho should be updated first.

            #include "YEqn.H"     // y*
            if (!setPrandtl) {    // h*
                #include "EEqn.H"
            }
            else {
                #include "EEqnExpPrandtl.H"
            }


            // rho = thermo.rho();   // account for temperature changes

            if (COMPPISO) {
                Info << "Compressible Piso Loop" << endl; // PISO loop
                while (pimple.correct()) {
                    #include "pCompEqn.H"
                }
            } else {
                Info << "Incompressible Piso Loop" << endl; // PISO loop
                while (pimple.correct()) {
                    #include "pIncompEqn.H"
                }
            }

            if (pimple.turbCorr()) turbulence->correct();

        }

        rho = thermo.rho();
        mu = turbulence->mu();


        #include "writeData.H"
        auto cpu   = runTime.elapsedCpuTime();
        auto clock = runTime.elapsedClockTime();

        Info << "---- ExecutionTime = " << cpu << " s" << " delta " << cpu-cpu0
             << "  ClockTime = " << clock << " s ----" << " delta " << clock - clock0
             << nl << endl;
    }

    Info << "End\n" << endl;
    return(0);
}
// ************************************************************************* //
