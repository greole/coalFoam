/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/



#include "coalCloud.H"
#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelPatchInteractionModels.H"
#include "makeReactingMultiphaseParcelStochasticCollisionModels.H" // MP variant
#include "makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"

// Reacting multiphase
#include "makeMyReactingMultiphaseParcelDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

#include "makeCoalParcelSurfaceReactionModels.H"
#include "makeCoalParcelThermolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(coalCloud, 0);

    makeParcelCloudFunctionObjects(coalCloud);

    // Kinematic sub-models
    makeThermoParcelForces(coalCloud);
    makeParcelDispersionModels(coalCloud);
    makeReactingMultiphaseParcelInjectionModels(coalCloud);
    makeParcelPatchInteractionModels(coalCloud);
    makeReactingMultiphaseParcelStochasticCollisionModels
    (
        coalCloud
    );
    makeReactingParcelSurfaceFilmModels(coalCloud);

    // Thermo sub-models
    makeParcelHeatTransferModels(coalCloud);

    // Reacting sub-models
    makeReactingMultiphaseParcelCompositionModels
    (
        coalCloud
    );
    makeReactingParcelPhaseChangeModels(coalCloud);

    // Reacting multiphase sub-models
    makeReactingMultiphaseParcelDevolatilisationModels
    (
        coalCloud
    );
    makeReactingMultiphaseParcelSurfaceReactionModels
    (
        coalCloud
    );

    makeCoalParcelSurfaceReactionModels(coalCloud);

    makeThermolysisModel(coalCloud);
    makeThermolysisModels(coalCloud);

    defineTemplateTypeNameAndDebug(coalParcel, 0);
    defineTemplateTypeNameAndDebug(Cloud<coalParcel>, 0);
}


// ************************************************************************* //
