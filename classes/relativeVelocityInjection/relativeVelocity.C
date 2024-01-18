/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "relativeVelocity.H"
#include "TimeDataEntry.H"
#include "distributionModel.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::relativeVelocity<CloudType>::relativeVelocity
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    relativeVelocityBase(owner.mesh(), this->coeffDict().lookup("patchName")),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    U0_(this->coeffDict().lookup("Urel")),
    flowRateProfile_
    (
        TimeDataEntry<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()
        )
    ),
    sizeDistribution_
    (
        distributionModels::distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    )
{
    duration_ = owner.db().time().userTimeToTime(duration_);

    relativeVelocityBase::updateMesh(owner.mesh());

    // Set total volume/mass to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0.0, duration_);
}


template<class CloudType>
Foam::relativeVelocity<CloudType>::relativeVelocity
(
    const relativeVelocity<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    relativeVelocityBase(im),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    U0_(im.U0_),
    flowRateProfile_(im.flowRateProfile_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::relativeVelocity<CloudType>::~relativeVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::relativeVelocity<CloudType>::updateMesh()
{
    relativeVelocityBase::updateMesh(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::relativeVelocity<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::relativeVelocity<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar nParcels = (time1 - time0)*parcelsPerSecond_;

        cachedRandom& rnd = this->owner().rndGen();

        label nParcelsToInject = floor(nParcels);

        // Inject an additional parcel with a probability based on the
        // remainder after the floor function
        if
        (
            nParcelsToInject > 0
         && (
               nParcels - scalar(nParcelsToInject)
             > rnd.globalPosition(scalar(0), scalar(1))
            )
        )
        {
            ++nParcelsToInject;
        }

        return nParcelsToInject;
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::relativeVelocity<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return flowRateProfile_.integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::relativeVelocity<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
)
{
    relativeVelocityBase::setPositionAndCell
    (
        this->owner().mesh(),
        this->owner().rndGen(),
        position,
        cellOwner,
        tetFaceI,
        tetPtI
    );
}


template<class CloudType>
void Foam::relativeVelocity<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    const volVectorField& U =  this->owner().U();
    // construct UInterp
    autoPtr<interpolation<vector> > UInterp = interpolation<vector>::New
    (
            "cell",
            U
    );

    // set particle velocity
    // parcel.U() = U0_  + this->owner().U()[parcel.cell()];
    parcel.U() = U0_ + UInterp().interpolate(parcel.position(), parcel.cell());

    // set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::relativeVelocity<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::relativeVelocity<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
