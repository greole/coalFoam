#ifndef makeReactingMultiphaseParcelInjectionModels_H
#define makeReactingMultiphaseParcelInjectionModels_H

#include "relativeVelocity.H"
// #include "basicReactingMultiphaseCloud.H"
#include "coalCloud.H"

// #define makeReactingMultiphaseParcelInjectionModels(CloudType)                \
//                                                                               \
//     makeInjectionModel(CloudType);                                            \
//     makeInjectionModelType(relativeVelocity, CloudType);
//
//
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
//
namespace Foam {
    makeInjectionModelType(relativeVelocity, coalCloud);
}

#endif
