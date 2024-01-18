/*---------------------------------------------------------------------------*\
Class
    probesContainer

Description
    Stores particle data


SourceFiles
    probesContainer.H

     Author: Gregor
     Todo:    make it a template class (particle type)
               introduce a second class mean value !

\*---------------------------------------------------------------------------*/

#ifndef probesContainer_H
#define probesContainer_H
#include "fvCFD.H"
#include "particleProbe.H"
#include <map>
#include <algorithm>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

class probesContainer {
public:

    probesContainer(Time& runTime_, fvMesh& mesh_, const coalCloud& cC);

    // bool register_(dictionary& probe);

    bool register_(
        const coalCloud& cC,
        word master,
        word type,
        vector pos,
        vector normal,
        scalar probe_size,
        scalar sampleFreq,
        scalar radius
        );

    bool sample(coalCloud& coalParcels);

    bool write();


private:

    Foam::Time& runTime;

    fvMesh& mesh;

    const IOdictionary particleSampleDict;

    const scalar timeStart;

    std::map<word, std::vector<particleProbe*> > probes;
};

}

#endif
