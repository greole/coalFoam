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

#include "competingRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::competingRate<CloudType>::competingRate
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
            Ta_[i] = readScalar(specSubDict.lookup("E"));
            A_[i] = readScalar(specSubDict.lookup("A"));
            Tlow_[i] = readScalar(specSubDict.lookup("Tlow"));
            hs_[i] = readScalar(specSubDict.lookup("hs"));
            i++;
        }
    }
}


template<class CloudType>
Foam::competingRate<CloudType>::competingRate
(
    const competingRate<CloudType>& srm
)
:
    ThermolysisModel<CloudType>(srm.owner_),
    Y_(srm.Y_),
    // emmitToGasphase(srm.emmitToGasphase),
    A_(srm.A_),
    Ta_(srm.Ta_),
    Tlow_(srm.Tlow_),
    hs_(srm.hs_),
    SINGLE_SPECIES_MODE(srm.SINGLE_SPECIES_MODE)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::competingRate<CloudType>::~competingRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::competingRate<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::competingRate<CloudType>::calculate
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

    const label idSolid = CloudType::parcelType::SLD;
    const label idGas = CloudType::parcelType::GAS;

    for (int i=0; i<2; i++) {
        if (T >= Tlow_[i]) {

        scalar Ysolid_ = YMixture[idSolid];
        scalar Ygas_   = YMixture[idGas];
        scalar mRaw    = mass*YSolid[2]*Ysolid_;

        scalar omega = A_[i]*exp(-Ta_[i]/T);
        scalar dm = omega*dt*mRaw;

        scalar dmGas  = Y_[i]*dm;
        scalar dmChar = (1.0 - Y_[i])*dm;

        // Reset Index for Single Species Mode
        const int j = (SINGLE_SPECIES_MODE) ? 0 : i;

        dMassGas[j]   -= dmGas;
        dMassSolid[1] -= dmChar;
        dMassSolid[2] += (dmChar + dmGas);
        }

    }
    return 0.0;
}


// ************************************************************************* //
