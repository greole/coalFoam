/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::relativeVelocity

Description
    Patch injection
    - User specifies
      - Total mass to inject
      - Name of patch
      - Injection duration
      - Initial parcel velocity
      - Injection volume flow rate
    - Parcel diameters obtained by distribution model
    - Parcels injected randomly across the patch

SourceFiles
    relativeVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef relativeVelocity_H
#define relativeVelocity_H

#include "InjectionModel.H"
#include "relativeVelocityBase.H"
#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class TimeDataEntry;

class distributionModel;

/*---------------------------------------------------------------------------*\
                       Class relativeVelocity Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class relativeVelocity
:
    public InjectionModel<CloudType>,
    public relativeVelocityBase
{
    // Private data

        //- Injection duration [s]
        scalar duration_;

        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Initial parcel velocity [m/s]
        const vector U0_;

        //- Flow rate profile relative to SOI []
        const TimeDataEntry<scalar> flowRateProfile_;

        //- Parcel size distribution model
        const autoPtr<distributionModels::distributionModel> sizeDistribution_;


public:

    //- Runtime type information
    TypeName("relativeVelocity");


    // Constructors

        //- Construct from dictionary
        relativeVelocity
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        relativeVelocity(const relativeVelocity<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType> > clone() const
        {
            return autoPtr<InjectionModel<CloudType> >
            (
                new relativeVelocity<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~relativeVelocity();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFaceI,
                label& tetPtI
            );

            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "relativeVelocity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
