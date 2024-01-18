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

\*---------------------------------------------------------------------------*/

#include "RFRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RFRate<CloudType>::RFRate
(
    const dictionary& dict,
    CloudType& owner
)
:
    ThermolysisModel<CloudType>(dict, owner, typeName),
    Y_(scalarList(2, 0.0)),
    // emmitToGasphase(this->coeffDict().lookupOrDefault<bool>("emmitToGasphase", false)),
    A_(scalarList(2, 0.0)),
    Ta_(scalarList(2, 0.0)),
    a0_(scalarList(2, 0.0)),
    a1_(scalarList(2, 0.0)),
    a2_(scalarList(2, 0.0)),
    a3_(scalarList(2, 0.0)),
    Tlow_(scalarList(2, 0.0)),
    hs_(scalarList(2, 0.0)),
    SINGLE_SPECIES_MODE(this->coeffDict().lookupOrDefault("SINGLE_SPECIES_MODE", false))
{
    // Fill Reaction Lists
    Info << "SINGLE_SPECIES_MODE: "  << SINGLE_SPECIES_MODE << endl;
    int i = 0;
    auto dict_ = this->coeffDict();
    for(auto spec: dict_.toc()) {
        if (spec != "Tlow" && spec != "SINGLE_SPECIES_MODE") {
            auto specSubDict = dict_.subDict(spec);
            Y_[i] = readScalar(specSubDict.lookup("Y"));
            Ta_[i] = readScalar(specSubDict.lookup("Ta"));
            A_[i] = readScalar(specSubDict.lookup("A"));
            a0_[i] = readScalar(specSubDict.lookup("a0"));
            a1_[i] = readScalar(specSubDict.lookup("a1"));
            a2_[i] = readScalar(specSubDict.lookup("a2"));
            a3_[i] = readScalar(specSubDict.lookup("a3"));
            Tlow_[i] = readScalar(specSubDict.lookup("Tlow"));
            hs_[i] = readScalar(specSubDict.lookup("hs"));
            i++;
        }
    }
}


template<class CloudType>
Foam::RFRate<CloudType>::RFRate
(
    const RFRate<CloudType>& srm
)
:
    ThermolysisModel<CloudType>(srm.owner_),
    Y_(srm.Y_),
    // emmitToGasphase(srm.emmitToGasphase),
    A_(srm.A_),
    Ta_(srm.Ta_),
    a0_(srm.a0_),
    a1_(srm.a1_),
    a2_(srm.a2_),
    a3_(srm.a3_),
    Tlow_(srm.Tlow_),
    hs_(srm.hs_),
    SINGLE_SPECIES_MODE(srm.SINGLE_SPECIES_MODE)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RFRate<CloudType>::~RFRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::RFRate<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::RFRate<CloudType>::calculate
(
    const scalar dt,
    const label cellI,
    const scalar d,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar mass,
    scalarField& YGas,
    scalarField& YLiquid,
    scalarField& YSolid,
    scalarField& YMixture,
    const scalar N,
    scalarField& dMassGas,
    scalarField& dMassLiquid,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier,
    scalarField& dHs
) const
{
    scalarList dmG_ = scalarList(2,0.0);
    scalar dHsTot = 0;

    const label idSolid = CloudType::parcelType::SLD;
    const label idGas = CloudType::parcelType::GAS;

    for (int i=0; i<2; i++) {
        if (T >= Tlow_[i]) {

        scalar Ysolid_ = YMixture[idSolid];
        scalar Ygas_   = YMixture[idGas];
        scalar mRaw    = mass*YSolid[2]*Ysolid_; // mass unconsumed raw
        scalar mAsh    = mass*YSolid[0]*Ysolid_; // mass unconsumed raw

        double m = (1-(mRaw+mAsh)/mass);
        double Z = a0_[i] + a1_[i]*m + a2_[i]*m*m + a3_[i]*m*m*m;

        scalar omega = A_[i]*exp(-(Ta_[i] + Z)/T);
        scalar dm = omega*dt*mRaw;

        scalar dmGas  = Y_[i]*dm;
        scalar dmChar = (1.0 - Y_[i])*dm;

        // Reset Index for Single Species Mode
        const int j = (SINGLE_SPECIES_MODE) ? 0 : i;

        dMassGas[j]   -= dmGas;
        dMassSolid[1] -= dmChar;
        dMassSolid[2] += (dmChar + dmGas);

        // Info << "dmGas: " << dmGas << " dmChar " << dmChar
        //      << " Z " << Z  << " T " << T
        //      << " A[i] " << A_[i] << " E[i] " << E_[i]
        //      << " j " << j
        //      << " dt " << dt << " omega " << omega
        //      << " m0 " << m
        //      << " mRaw1 " <<  mRaw << endl;
        //
        dHs[0] += dmGas * hs_[i];
        }

    }
    return 0.0;
}


// ************************************************************************* //
