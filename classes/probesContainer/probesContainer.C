/*---------------------------------------------------------------------------*\

Class
    particleProbe

Description
    Stores particle data


SourceFiles
    particleProbe.C

    Author: Gregor

\*---------------------------------------------------------------------------*/

#include "probesContainer.H"

probesContainer::probesContainer(
        Foam::Time& runTime_, fvMesh& mesh_, const coalCloud& cC):
    runTime(runTime_),
    mesh(mesh_),
    particleSampleDict(IOobject (
       "system/particleSampleDict",
       runTime,
       IOobject::MUST_READ,
       IOobject::NO_WRITE
    )),
    timeStart(particleSampleDict.lookupOrDefault("timeStart", 0.0))
{
    Info << "Constructing Probes " << endl;
    for(auto probe: particleSampleDict.toc()) {
        // NOTE we need to check if the particleSampleDict entry is
        // a key or an actual probe sub dict
        if (probe == "timeStart") {}
        else {
            auto probeSubDict = particleSampleDict.subDict(probe);
            if (probeSubDict.found("type")) {
                word probeType = probeSubDict.lookup("type");
                if (probeType.find("multipoint") != std::string::npos) {
                   word stype    = probeType.erase(0, 5);

                   word master   = probeSubDict["master"];
                   label nProbes = readLabel(probeSubDict.lookup("nProbes"));
                   vector start  = probeSubDict["startPosition"];
                   vector end    = probeSubDict["endPosition"];
                   vector delta  = 1.0/float(nProbes-1) * (end - start);

                   for (int i=0; i < nProbes; i++) {
                        register_(
                            cC, master, stype,
                            delta * float(i),
                            vector(probeSubDict["normal"]),
                            readScalar(probeSubDict["probeSize"]),
                            readScalar(probeSubDict["sampleFreq"]),
                            0.0
                        );
                   }
                }

                else if (probeType.find("multiradial") != std::string::npos) {

                   word stype    = probeType.erase(0, 5);

                   word master   = probeSubDict["master"];
                   label nProbes = readLabel(probeSubDict.lookup("nProbes"));
                   vector pos    = probeSubDict["position"];
                   scalar start  = readScalar(probeSubDict["startRadius"]);
                   scalar end    = readScalar(probeSubDict["endRadius"]);
                   scalar delta  = 1.0/float(nProbes - 1)*(end - start);

                   for (int i=0; i < nProbes; i++) {
                        register_(
                            cC, master, stype, pos,
                            vector(probeSubDict["normal"]),
                            readScalar(probeSubDict["probeSize"]),
                            readScalar(probeSubDict["sampleFreq"]),
                            delta*float(i)
                        );
                   }

                }

                // else {register_(probeSubDict);}
            }
        }
    }
}

bool probesContainer::register_ (
        const coalCloud& cC,
        Foam::word master,
        Foam::word type,
        Foam::vector pos,
        Foam::vector normal,
        Foam::scalar probe_size,
        Foam::scalar sampleFreq,
        Foam::scalar radius
    ) {

    if (mesh.findCell(pos) == -1) {return false;}

    if (type != "radialProbe"){
        probes[master].push_back(
            new particleProbe(
                "probe", type, pos, probe_size, runTime, cC, sampleFreq
            )
        );
    } else {
        probes[master].push_back(
            new particleProbe(
               "probe", type, pos, normal, radius, probe_size, runTime, cC, sampleFreq
            )
        );
    }

    return true;
}

// bool probesContainer::register_ (
//         Foam::dictionary& probe
//     ) {
//     vector pos(probe["position"]);
//
//     if (mesh.findCell(pos) == -1) {return false;}
//
//     word master(probe.lookup("master"));
//     word type(probe.lookup("type"));
//     vector normal(probe["normal"]);
//
//     scalar freq = readScalar(probe.lookup("sampleFreq"));
//     scalar size = readScalar(probe.lookup("probeSize"));
//
//     scalar radius = 0.0;
//
//     return register_(master, type, pos, normal, size, freq, );
// }

bool probesContainer::sample(coalCloud& coalParcels) {
    if (runTime.timeOutputValue() >= timeStart &&
        coalParcels.solution().active()) {

        Info << "Sampling particles" << endl;
        for (auto& probes_per_type: probes) {
            for (auto probe: probes_per_type.second) {
                if (probe->sample()) {
                    forAllConstIter(coalCloud, coalParcels, iter) {
                        probe->insertParticleData(iter());
                    }
                }
            }
        }

    }
    return true;
}

bool probesContainer::write() {
    for (auto& probes_per_type: probes) {
        fileName samplePath(
            runTime.path()/runTime.timeName()/"lagrangian/probes"
        );
        Foam::mkDir(samplePath);

        word fn = "_UMean_uuMean_vvMean_wwMean_TMean_ddt(T)_nParticles_nParcels";

        word ddt = "";
        word ddt_ = "";
        auto avgs0 = probes_per_type.second[0]->get_averages();
        for (auto &ddtSolid_: avgs0.solidNames) {
            ddt += "_ddt(" + ddtSolid_ + ")";
            ddt_ += " ddt(" + ddtSolid_ + ") ";
        }
        for (auto &ddtGas_: avgs0.gasNames) {
            ddt += "_ddt(" + ddtGas_ + ")";
            ddt_ += " ddt(" + ddtGas_ + ") ";
        }
        fn += ddt;
        fn += ".dat";

        fileName sampleName(samplePath/(probes_per_type.first + fn));

        OFstream samples(sampleName);
        samples << "# position "
                << " uMean"
                << " vMean"
                << " wMean"
                << " uVar "
                << " vVar "
                << " wVar "
                << " TMean"
                << " ddtTemp "
                << " nParticles"
                << " nParcels "
                << ddt_ << nl;

        for (auto probe: probes_per_type.second) {
            auto avgs = probe->get_averages();
            samples << avgs.position << " "
                    << avgs.avg_vel[0] << " "
                    << avgs.avg_vel[1] << " "
                    << avgs.avg_vel[2] << " "
                    << avgs.vel_var[0] << " "
                    << avgs.vel_var[1] << " "
                    << avgs.vel_var[2] << " "
                    << avgs.TMean << " "
                    << avgs.ddtTemp << " "
                    << avgs.nParticles << " "
                    << avgs.nParcel << " ";

                    for (auto &ddtSolid_: avgs.ddtSolid) {samples << ddtSolid_ << " ";}
                    for (auto &ddtGas_: avgs.ddtGas) {samples << ddtGas_ << " ";}

                    samples << nl;
        }
    }
    return true;
}


