// // #include "basicReactingMultiphaseCloud.H"
#include "coalCloud.H"

#include "BiniDispersionForce.H"

#include "ParticleForce.H" // thermo variant

namespace Foam {

    // makeParticleForceModel(basicReactingMultiphaseCloud);
    // makeParticleForceModel(coalCloud);
    // makeParticleForceModelType(BiniDispersionForce, basicReactingMultiphaseCloud);
    makeParticleForceModelType(BiniDispersionForce, coalCloud);
}
