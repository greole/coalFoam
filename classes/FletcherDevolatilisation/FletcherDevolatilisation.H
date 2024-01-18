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
    Foam::FletcherDevolatilisation

Description
    Single kinetic rate devolatisation model.
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = A1.exp(- E/R.T)

    Where:
        kappa = rate constant
        A1    = activation energy (user input)
        E     = pre-exponential factor (user input)
        R     = universal gas constant
        T     = temperature

    Usage:

        FletcherDevolatilisationCoeffs
        {
            volatileData
            (
                (CxHyOz0 12 0.5 0.1 5(c1 c2 c3 c4 c5))   // (name A1 E Y coeffs)
                (CxHyOz1 12 0.5 0.2 5(c1 c2 c3 c4 c5))   // (name A1 E Y coeffs)
            );

            volatileResidualCoeff 1e-6;
        }

\*---------------------------------------------------------------------------*/

#ifndef FletcherDevolatilisation_H
#define FletcherDevolatilisation_H

// #include "fvCFD.H"
#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class FletcherDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FletcherDevolatilisation
:
    public DevolatilisationModel<CloudType>
{
    // Helper class to store specie-local volatile data
    class volatileData
    {
        // Private data

            //- Specie name
            word name_;

            //- Activation energy
            scalar A1_;

            //- Pre-exponential factor
            scalar E_;

            //- Char conversion ratio
            scalar Y_;

            //- Correction Coeffs
            // scalarField c_;
            scalar c0_;
            scalar c1_;
            scalar c2_;
            scalar c3_;
            scalar c4_;

            scalarField c_;

    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
                name_(word::null),
                A1_(0.0),
                E_(0.0),
                Y_(0.0),
                c0_(0.0),
                c1_(0.0),
                c2_(0.0),
                c3_(0.0),
                c4_(0.0),
                c_(scalarField(5, 0.0))
            {}

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is),
                A1_(readScalar(is)),
                E_(readScalar(is)),
                Y_(readScalar(is)),
                c0_(readScalar(is)),
                c1_(readScalar(is)),
                c2_(readScalar(is)),
                c3_(readScalar(is)),
                c4_(readScalar(is)),
                c_(scalarField(5,0))
            {
                c_[0] = c0_;
                c_[1] = c1_;
                c_[2] = c2_;
                c_[3] = c3_;
                c_[4] = c4_;
            }

            //- Construct as copy
            volatileData(const volatileData& vd)
            :
                name_(vd.name_),
                A1_(vd.A1_),
                E_(vd.E_),
                Y_(vd.Y_),
                c0_(vd.c0_),
                c1_(vd.c1_),
                c2_(vd.c2_),
                c3_(vd.c3_),
                c4_(vd.c4_),
                c_(vd.c_)
            {}


        //- Destructor
        ~volatileData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the activation energy
                scalar A1() const
                {
                    return A1_;
                }

                //- Return const access to the pre-exponential factor
                scalar E() const
                {
                    return E_;
                }

                scalar Y() const
                {
                    return Y_;
                }

                scalarField c() const
                {
                    return c_;
                }



        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_  >> vd.A1_  >> vd.E_ >> vd.Y_
                    >> vd.c0_ >> vd.c1_ >> vd.c2_ >> vd.c3_ >> vd.c4_ ;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.A1_ << token::SPACE
                    << vd.E_ << token::SPACE
                    << vd.Y_ << token::SPACE
                    << vd.c0_ << vd.c1_ << vd.c2_ << vd.c3_ << vd.c4_
                    << token::END_LIST;

                return os;
            }
    };


    // Private data

        // Model constants

            //- List of volatile data - (name A1 E)
            List<volatileData> volatileData_;

            //- List of initial volatile mass fractions
            List<scalar> YVolatile0_;

            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;


public:

    //- Runtime type information
    TypeName("FletcherDevolatilisation");


    // Constructors

        //- Construct from dictionary
        FletcherDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        FletcherDevolatilisation
        (
            const FletcherDevolatilisation<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DevolatilisationModel<CloudType> > clone() const
        {
            return autoPtr<DevolatilisationModel<CloudType> >
            (
                new FletcherDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~FletcherDevolatilisation();


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const scalar age,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalarField& YGasEff,
            const scalarField& YLiquidEff,
            const scalarField& YSolidEff,
            label& canCombust,
            scalarField& dMassDV
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
#   include "FletcherDevolatilisation.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
