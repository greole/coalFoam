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

#include "singleKineticRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::singleKineticRate<CloudType>::singleKineticRate
(
    const dictionary& dict,
    CloudType& owner
)
:
    ThermolysisModel<CloudType>(dict, owner, typeName),
    Y_(readScalar(this->coeffDict().lookup("Y"))),
    // emmitToGasphase(this->coeffDict().lookupOrDefault<bool>("emmitToGasphase", false)),
    A_(readScalar(this->coeffDict().lookup("A"))),
    Ta_(readScalar(this->coeffDict().lookup("Ta"))),
    Tlow_(readScalar(this->coeffDict().lookup("Tlow"))),
    hs_(readScalar(this->coeffDict().lookup("h")))
{}


template<class CloudType>
Foam::singleKineticRate<CloudType>::singleKineticRate
(
    const singleKineticRate<CloudType>& srm
)
:
    ThermolysisModel<CloudType>(srm.owner_),
    Y_(srm.Y_),
    // emmitToGasphase(srm.emmitToGasphase),
    A_(srm.A_),
    Ta_(srm.Ta_),
    Tlow_(srm.Tlow_),
    hs_(srm.hs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::singleKineticRate<CloudType>::~singleKineticRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::singleKineticRate<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::singleKineticRate<CloudType>::calculate
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
    // do nothing if below
    // lower temperature limit
    if (T < Tlow_) return 0.0;
    const label idSolid = CloudType::parcelType::SLD;
    const label idGas = CloudType::parcelType::GAS;
    const scalar Ysolid_ = YMixture[idSolid];
    const scalar Ygas_ = YMixture[idGas];

    // Info << "specie::RR " <<  specie::RR
    //      << " T  " <<  T << endl;
    scalar omega = A_*exp(-Ta_/T);
    scalar dm = omega*dt*YSolid[2]*Ysolid_*mass;

    scalar dmGas  = Y_*dm;
    scalar dmChar = (1.0 - Y_)*dm;

    dMassGas[0]   -= dmGas;
    dMassSolid[1] -= dmChar;
    dMassSolid[2] += (dmChar + dmGas);

    dHs[0] = dmGas * hs_;

    return 0.0;
}


// ************************************************************************* //
