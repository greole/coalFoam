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
    Foam::RFRate

Description
    Dummy surface reaction model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef RFRate_H
#define RFRate_H

#include "thermolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class RFRate Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RFRate
:
    public ThermolysisModel<CloudType>
{

    // ratio between formed Yield and formed Char
    scalarList Y_;

    // emmit formed Yield to Gasphase directly
    // const bool emmitToGasphase;

    // reaction rate coefficitient [-]
    scalarList A_;

    // reaction rate coefficitient [J/kmol]
    scalarList Ta_;

    //
    scalarList a0_;
    scalarList a1_;
    scalarList a2_;
    scalarList a3_;

    scalarList Tlow_;

    scalarList hs_;

    bool SINGLE_SPECIES_MODE;

public:

    //- Runtime type information
    TypeName("RFRate");


    // Constructors

        //- Construct from dictionary
        RFRate(const dictionary& dict, CloudType& owner);

        //- Construct copy
        RFRate(const RFRate<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<ThermolysisModel<CloudType> > clone() const
        {
            return autoPtr<ThermolysisModel<CloudType> >
            (
                new RFRate<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RFRate();


    // Member Functions

        //- Flag to indicate whether model activates devolatisation model
        virtual bool active() const;

        //- Update surface reactions
        virtual scalar calculate
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
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
#   include "RFRate.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
