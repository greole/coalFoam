/*---------------------------------------------------------------------------*\

Class
    particleProbe

Description
    Stores particle data


SourceFiles
    particleProbe.C

    Author: Gregor

\*---------------------------------------------------------------------------*/

#include "particleProbe.H"
#include "Time.H"
#include "OFstream.H"


Foam::particleProbe::particleProbe (
    word probeName_,
    word probeType_,
    vector position_,
    scalar deltax_,
    Foam::Time& runTime_,
    const coalCloud& cC,
    scalar sampleFreq_=1
) :
    runTime(runTime_),
    sampleFreq(sampleFreq_)
{
    name          = probeName_;
    type          = probeType_;
    position      = position_;
    avg_vel       = vector::zero;
    vel_var       = vector::zero;
    deltax        = deltax_;
    remSampleTime = 0.0;
    nParcel       = 0.0;
    n             = 0.0;
    TMean         = 0.0;
    ddtTemp       = 0.0;
    alpha         = 0.0;
    beta          = 0.0;

    probe_averages.position = position[2];
    probe_averages.avg_vel = avg_vel;
    probe_averages.vel_var = vel_var;
    probe_averages.TMean = 0.0;
    probe_averages.ddtTemp = 0.0;
    probe_averages.gas_phase_vel = Foam::vector::zero;
    probe_averages.nParticles = n;
    probe_averages.nParcel = nParcel;

    const label idGas = cC.composition().idGas();
    gasNames = cC.composition().componentNames(idGas);
    const label idSolid = cC.composition().idSolid();
    solidNames = cC.composition().componentNames(idSolid);

    ddtSolid = Foam::scalarField(solidNames.size(), 0);
    ddtGas = Foam::scalarField(gasNames.size(), 0);
}


Foam::particleProbe::particleProbe (
    word probeName_,
    word probeType_,
    vector position_,
    vector probeNormal_,
    scalar probeRadius_,
    scalar deltax_,
    Foam::Time& runTime_,
    const coalCloud& cC,
    scalar sampleFreq_=1
) :
    runTime(runTime_),
    sampleFreq(sampleFreq_)
{
    name          = probeName_;
    type          = probeType_;
    position      = position_;
    avg_vel       = vector::zero;
    vel_var       = vector::zero;
    deltax        = deltax_;
    probeRadius   = probeRadius_;
    remSampleTime = 0.0;
    nParcel       = 0.0;
    n             = 0.0;
    TMean         = 0.0;
    ddtTemp       = 0.0;
    alpha         = 0.0;
    beta          = 0.0;
    // FIXME make it DRY
    probe_averages.position = position[2];
    probe_averages.avg_vel = avg_vel;
    probe_averages.vel_var = vel_var;
    probe_averages.TMean = 0.0;
    probe_averages.ddtTemp = 0.0;
    probe_averages.gas_phase_vel = Foam::vector::zero;
    probe_averages.nParticles = n;
    probe_averages.nParcel = nParcel;
    probe_averages.ddtSolid = Foam::scalarField(5, 0);
    probe_averages.ddtGas = Foam::scalarField(5, 0);

    const label idGas = cC.composition().idGas();
    gasNames = cC.composition().componentNames(idGas);
    const label idSolid = cC.composition().idSolid();
    solidNames = cC.composition().componentNames(idSolid);

    ddtSolid = Foam::scalarField(solidNames.size(), 0);
    ddtGas = Foam::scalarField(gasNames.size(), 0);
}


bool Foam::particleProbe::sample(){
    remSampleTime--;
    if (remSampleTime <= 0){
        remSampleTime = sampleFreq;
        return true;}
    else {return false;}
}


bool Foam::particleProbe::insertParticleData (
        const coalParcel& pP
    ) {
    if (!isInside(pP.position())) {return false;}
    nParcel++;
    scalar pM = pP.nParticle(); // new particle Mass
    n         = n + pM;         // new total Mass
    scalar rn = 1.0/n;
    alpha     = (n-pM)*rn;
    beta      = pM*rn;

    // Velocities
    // adjusting variance for old mean value
    vel_var += component_sqr(avg_vel);
    // calc average particle velocity

    if (type == "pointProbe") {
        avg_vel = alpha*avg_vel + beta*pP.U();
    }
    else {
        auto pi2 = 1.5707963267948966;
        const auto x = pP.position()[0]+SMALL;
        const auto y = pP.position()[1];
        auto theta = (Foam::atan(y/x));

        if (x < 0.0 && y > 0) theta = pi2*2 + theta;
        if (x < 0.0 && y < 0) theta = pi2*2 + theta;
        if (x > 0.0 && y < 0) theta = pi2*4 + theta;

        const auto uswirl = Foam::cos(theta)*pP.U()[1] - Foam::sin(theta)*pP.U()[0];
        const auto urad   = Foam::sin(theta)*pP.U()[1] - Foam::cos(theta)*pP.U()[0];
        avg_vel[0] = alpha*avg_vel[0] + beta*urad;
        avg_vel[1] = alpha*avg_vel[1] + beta*uswirl;
        avg_vel[2] = alpha*avg_vel[2] + beta*pP.U()[2];
    }

    // calc average particle velocity variance
    vel_var = alpha*vel_var
            + beta*(component_sqr(pP.U()))
            - component_sqr(avg_vel);
    // Temperatures
    // calc average particle velocity
    TMean = alpha*TMean + beta*pP.T();
    ddtTemp = alpha*ddtTemp + beta*pP.getDdtTemp();

    auto ddtSolid_ = pP.getDdtSolid();
    auto ddtGas_ = pP.getDdtGas();


    ddtSolid = alpha*ddtSolid + beta*ddtSolid_;
    ddtGas = alpha*ddtGas + beta*ddtGas_;
    return true;
}

inline vector Foam::particleProbe::component_sqr(vector vec) {
    return vector {vec[0]*vec[0], vec[1]*vec[1], vec[2]*vec[2]};} // RVO

bool Foam::particleProbe::isInside (vector particle_pos) {
    if (type=="pointProbe" && deltax >= mag(particle_pos - position))
        {return true;}
    if (type=="radialProbe"
            && deltax >= std::abs(probeRadius - pow(
                    pow(particle_pos[0],2)
                  + pow(particle_pos[1],2), 0.5))
            && deltax >= std::abs(particle_pos[2] - position[2]))
       {return true;}
    return false;
}

particleProbe::avg_values particleProbe::get_averages() {
    scalar pos;
    if (type == "radialProbe") {pos = probeRadius;}
    else {pos = position[2];}
        probe_averages.position = pos;
        probe_averages.avg_vel = avg_vel;
        probe_averages.vel_var = vel_var;
        probe_averages.TMean = TMean;
        probe_averages.ddtTemp = ddtTemp;
        probe_averages.gas_phase_vel= vector::zero;
        probe_averages.nParticles = n ;
        probe_averages.nParcel = nParcel;
        probe_averages.ddtSolid = ddtSolid;
        probe_averages.ddtGas = ddtGas;
        probe_averages.gasNames = gasNames;
        probe_averages.solidNames = solidNames;
    return probe_averages;
}
