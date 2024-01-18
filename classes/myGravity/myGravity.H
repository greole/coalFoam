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
    Foam::myGravity

Description
    Calculates particle gravity force

SourceFiles
    myGravityI.H
    myGravity.C

\*---------------------------------------------------------------------------*/

#ifndef myGravity_H
#define myGravity_H

#include "fvCFD.H"
#include <math.h>
#include "ParticleForce.H"
#include "basicReactingMultiphaseCloud.H"
#include TURBMODEL
#include LESMODEL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

    // forwad declaration in namespace to avoid
    // unneeded warnings
    namespace compressible {
        class turbulenceModel;
    }


    class fvMesh;

// class ParticleForce<basicReactingMultiphaseCloud>;

/*---------------------------------------------------------------------------*\
                       Class myGravity Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class myGravity
:
    public ParticleForce<CloudType>
{
    // Private data

        const fvMesh& mesh;

        // start position for dispersion z-coordinate

public:

    //- Runtime type information
    TypeName("myGravity");


    // Constructors

        //- Construct from mesh
        myGravity
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        myGravity(const myGravity& gf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType> > clone() const
        {
            return autoPtr<ParticleForce<CloudType> >
            (
                new myGravity<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~myGravity();


    // Member Functions

        // Access

            //- Return the acceleration due to gravity
            inline const vector& g() const;

            inline vector rnd();


        // Evaluation

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myGravityI.H"

#ifdef NoRepository
    #include "myGravity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //::
