// // #include "basicReactingMultiphaseCloud.H"
#include "coalCloud.H"

#include "myGravity.H"

#include "ParticleForce.H" // thermo variant

namespace Foam {

    // makeParticleForceModel(basicReactingMultiphaseCloud);
    // makeParticleForceModel(coalCloud);
    // makeParticleForceModelType(myGravity, basicReactingMultiphaseCloud);
    makeParticleForceModelType(myGravity, coalCloud);
}
