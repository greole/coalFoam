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

#include "myGravity.H"

//#define INC() std::getenv("$FOAM_SRC")  "/turbulenceModels/compressible/turbulenceModel/turbulenceModel.H"
#include TURBMODEL
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myGravity<CloudType>::myGravity
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    mesh(this->owner().mesh())
{}


template<class CloudType>
Foam::myGravity<CloudType>::myGravity(const myGravity& gf)
:
    ParticleForce<CloudType>(gf),
    mesh(this->owner().mesh())
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myGravity<CloudType>::~myGravity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::myGravity<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const vector g = {0.0, 0.0, -9.81};

    forceSuSp value(vector::zero, 0.0);

    // Explicit Force
    value.Su() = mass*g;

    return value;
}


// ************************************************************************* //
