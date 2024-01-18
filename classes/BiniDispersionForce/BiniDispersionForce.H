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
    Foam::BiniDispersionForce

Description
    Calculates particle gravity force

SourceFiles
    BiniDispersionForceI.H
    BiniDispersionForce.C

\*---------------------------------------------------------------------------*/

#ifndef BiniDispersionForce_H
#define BiniDispersionForce_H

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
                       Class BiniDispersionForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BiniDispersionForce
:
    public ParticleForce<CloudType>
{
    // Private data

        const fvMesh& mesh;

        const compressible::turbulenceModel& turbModel;
        //- Reference to the acceleration due to gravity
        const vector& g_;

        // start position for dispersion z-coordinate
        const scalar positionOffset;

        const scalar alpha_;

        const scalar C0_;

        const word taut_;

        const autoPtr<const volScalarField> delta_;

        // const autoPtr<const volScalarField> k_;

public:

    //- Runtime type information
    TypeName("BiniDispersion");


    // Constructors

        //- Construct from mesh
        BiniDispersionForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        BiniDispersionForce(const BiniDispersionForce& gf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType> > clone() const
        {
            return autoPtr<ParticleForce<CloudType> >
            (
                new BiniDispersionForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~BiniDispersionForce();


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

#include "BiniDispersionForceI.H"

#ifdef NoRepository
    #include "BiniDispersionForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //::
