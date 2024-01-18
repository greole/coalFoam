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

Class
    coalCloud

Description
    Cloud class to introduce coal parcels

\*---------------------------------------------------------------------------*/

#ifndef coalCloud_H
#define coalCloud_H

#include "thermolysisModel.H"

#include "Cloud.H"
#include "KinematicCloud.H"
#include "ThermoCloud.H"
#include "ReactingCloud.H"
#include "ReactingMultiphaseCloud.H"
#include "coalParcel.H"

#include "typeInfo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes

template<class CloudType>
class thermolysisModel;


// class myReactingMultiphaseCloudNT
// {
//     // Private Member Functions
//
//         //- Disallow default bitwise copy construct
//         myReactingMultiphaseCloudNT(const myReactingMultiphaseCloudNT&);
//
//         //- Disallow default bitwise assignment
//         void operator=(const myReactingMultiphaseCloudNT&);
//
//
// public:
//
//     //- Runtime type information
//     TypeName("myReactingMultiphaseCloudNT");
//
//     // Constructors
//
//         //- Null constructor
//         myReactingMultiphaseCloudNT(){};
//
//
//     //- Destructor
//     virtual ~myReactingMultiphaseCloudNT(){};
// };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
class myReactingMultiphaseCloud
:
    public CloudType
{
    public:

        // Public typedefs

            //- Type of cloud this cloud was instantiated for
            typedef CloudType cloudType;

            //- Type of parcel the cloud was instantiated for
            typedef typename CloudType::particleType parcelType;

            typedef myReactingMultiphaseCloud<CloudType> myreactingMultiphaseCloudType;
private:

    autoPtr<myReactingMultiphaseCloud<CloudType> > cloudCopyPtr_;

protected:

    typename parcelType::constantProperties constProps_;

    autoPtr
    <
        ThermolysisModel<myReactingMultiphaseCloud<CloudType> >
    >
    thermolysisModel_;

    void setModels();

public:

        //- Return const access to devolatilisation model
        inline const ThermolysisModel
        <
            myReactingMultiphaseCloud<CloudType>
        >&
        thermolysis() const {return thermolysisModel_;};

        //- Return reference to devolatilisation model
        inline ThermolysisModel
        <
           myReactingMultiphaseCloud<CloudType>
        >&
        thermolysis() {return thermolysisModel_();};


        myReactingMultiphaseCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        myReactingMultiphaseCloud
        (
            myReactingMultiphaseCloud<CloudType>& c,
            const word& name
        );

        void cloudReset(myReactingMultiphaseCloud<CloudType>& c);

        //- Copy constructor with new name - creates bare cloud
        myReactingMultiphaseCloud
        (
            const fvMesh& mesh,
            const word& name,
            const myReactingMultiphaseCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new myReactingMultiphaseCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new myReactingMultiphaseCloud(this->mesh(), name, *this)
            );
        }


        void evolve();

        //- Destructor
        virtual ~myReactingMultiphaseCloud();

        TypeName("CoalCloud");

};

typedef myReactingMultiphaseCloud
<
    ReactingMultiphaseCloud
    <
        ReactingCloud
        <
            ThermoCloud
            <
                KinematicCloud
                <
                    Cloud
                    <
                        coalParcel
                    >
                >
            >
        >
    >
> coalCloud;

}


#include "coalCloud.C"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

 // ************************************************************************* //
