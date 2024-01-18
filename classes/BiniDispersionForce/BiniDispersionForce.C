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

#include "BiniDispersionForce.H"

//#define INC() std::getenv("$FOAM_SRC")  "/turbulenceModels/compressible/turbulenceModel/turbulenceModel.H"
#include TURBMODEL
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BiniDispersionForce<CloudType>::BiniDispersionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    mesh(this->owner().mesh()),
    turbModel(mesh.lookupObject<compressible::turbulenceModel>("turbulenceModel")),
    g_(owner.g().value()),
    positionOffset(readScalar(this->coeffs().lookup("positionOffset"))),
    alpha_(readScalar(this->coeffs().lookup("alpha"))),
    C0_(readScalar(this->coeffs().lookup("C0"))),
    taut_(this->coeffs().lookupOrDefault("Method", word("default"))),
    delta_(&mesh.lookupObject<compressible::LESModel>("turbulenceModel").delta())
{}


template<class CloudType>
Foam::BiniDispersionForce<CloudType>::BiniDispersionForce(const BiniDispersionForce& gf)
:
    ParticleForce<CloudType>(gf),
    mesh(this->owner().mesh()),
    turbModel(mesh.lookupObject<Foam::compressible::turbulenceModel>("turbulenceModel")),
    g_(gf.g_),
    positionOffset(gf.positionOffset),
    alpha_(gf.alpha_),
    C0_(gf.C0_),
    taut_(gf.taut_),
    delta_(gf.delta_)
    // k_(gf.k_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BiniDispersionForce<CloudType>::~BiniDispersionForce()
{}


template<class CloudType>
Foam::vector Foam::BiniDispersionForce<CloudType>::rnd()
{
    cachedRandom& rnd = this->owner().rndGen();
    return rnd.sample01<vector>()*2.0-vector::one;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::BiniDispersionForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{

    if (p.position()[2] < positionOffset) {return forceSuSp(vector::zero, 0.0);}



    // const scalar Cd = 24.0 + 4.0*cbrt(Re*Re); // Re^(2/3)
    // Stokes Flow for Re -> 0
    // Loth p.30 Eqn 1.46b
    const scalar Cd = 24.0/(Re+SMALL);

    const scalar UrelMag = mag(p.U() - p.Uc());

    // Bini & Jones 2007  eqn 3.3
    const scalar tau_p = 4*p.d()/(3*Cd*UrelMag);

    const scalar Delta = delta_()[p.cell()];
    const label cellI = p.cell();

    // FIXME this is only valid for basic Smag model
    const scalar cellM = mesh.V()[cellI]*p.rhoc();
    const scalar k_SGS = cellM*Foam::sqr((mesh.lookupObject<const volScalarField>("muSgs")[cellI]/p.rhoc())/(Delta*0.094));

    const scalar k_SGS_SQRT = Foam::sqrt(k_SGS) + SMALL;


    scalar tau_t;

    if (taut_ == "default" ) {
           tau_t = Foam::pow(tau_p, 2.0*alpha_)
               / Foam::pow(Delta/k_SGS_SQRT, 2.0*alpha_-1.0);
    }
    if (taut_ == "particle" ) {
           tau_t = Foam::pow(Delta, 1.0/3.0) / mag(p.U());
    }
    if (taut_ == "ksgs" ) {
           tau_t = Foam::pow(Delta, 1.0/3.0) / k_SGS_SQRT;
    }

    const vector rd = const_cast<BiniDispersionForce<CloudType> *>(this)->rnd();

    forceSuSp value(vector::zero, 0.0);

    // Explicit Force
    value.Su() = Foam::sqrt(C0_*k_SGS/tau_t)*mass/dt*rd;

    return value;
}


// ************************************************************************* //
