/*---------------------------------------------------------------------------*\
Class
    particleProbe

Description
    Stores particle data


SourceFiles
    particleProbe.H

     Author: Gregor
     Todo:    make it a template class (particle type)
               introduce a second class mean value !

\*---------------------------------------------------------------------------*/

#ifndef particleProbe_H
#define particleProbe_H
#include "fvCFD.H"
#include "coalParcel.H"
#include "coalCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

class particleProbe {

public:
    // collect all values in a struct for simple
    // return
    struct avg_values {
        scalar position = 0.0;
        vector avg_vel = Foam::vector::zero;
        vector vel_var = Foam::vector::zero;
        scalar TMean = 0.0;
        scalar ddtTemp = 0.0;
        vector gas_phase_vel = Foam::vector::zero;
        scalar nParticles = 0.0;
        scalar nParcel = 0.0;
        scalarField ddtSolid = Foam::scalarField(5,0);
        scalarField ddtGas = Foam::scalarField(5,0);
        wordList gasNames;
        wordList solidNames;
    };

    particleProbe (
        word probeName_,
        word probeType_,
        vector position_,
        scalar deltax_,
        Foam::Time& runTime_,
        const coalCloud& cC,
        scalar sampleFreq_
    );
    // radial probe
    particleProbe (
        word probeName_,
        word probeType_,
        vector position_,
        vector probeNormal_,
        scalar probeRadius_,
        scalar deltax_,
        Foam::Time& runTime_,
        const coalCloud& cC,
        scalar sampleFreq_
    );

    bool insertParticleData(const coalParcel& pP);

    void writeOutData();

    word getProbeName()     {return name;}

    word getProbeType()     {return type;}

    vector getProbePos()    {return position;}

    scalar getProbeRadius() {return probeRadius;}

    vector getVar()         {return vel_var;}

    vector getVel()         {return avg_vel;}

    scalar getT()           {return TMean;}

    scalar getNParticle()   {return n;}

    scalar getNParcel()     {return nParcel;}

    vector component_sqr(vector vec);

    // test if particle is inside probe volume
    bool isInside(vector pos_);

    // average
    vector average(vector valuevector, vector oldvaluevector);

    void printProbeName(){Info << "ProbeName: " <<  name << endl;}

    bool sample();

    avg_values get_averages();

    private:
        // TODO store on vector
        scalar deltax;               // probe size
        scalar probeRadius;          // probe radius for radial probes
        word name;                   // Probe name
        word type;                   // Probe type
        scalar TMean;                // mean temperature
        scalar ddtTemp ;             // change in temperature
        vector position;             // Probe position
        vector avg_vel;              // averaged particle velocity
        vector vel_var;              // particle velocity variances
        vector Uc_;                  // A vector for the gas-phase velocity
        scalar n;                    // used for averaging
        scalar nParcel;
        scalar alpha;                // factor for averaging
        scalar beta;                 // factor for averaging
        scalarField ddtSolid;
        scalarField ddtGas;
        Foam::Time& runTime;
        scalar sampleFreq;           // sample every n seconds
        scalar remSampleTime;        // time left to next sample event

        wordList gasNames;

        wordList solidNames;

        bool adjustSize;
        avg_values probe_averages;
};

}

#endif
