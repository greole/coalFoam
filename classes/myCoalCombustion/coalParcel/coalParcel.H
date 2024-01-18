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
    coalParcel

Description
    Definition of coal parcel

SourceFiles
    CoalParcel.C

\*---------------------------------------------------------------------------*/

#ifndef coalParcel_H
#define coalParcel_H

#include "contiguous.H"
#include "particle.H"
#include "KinematicParcel.H"
#include "ThermoParcel.H"
#include "ReactingParcel.H"
// #include "ReactingMultiphaseParcel.H"

#include "Cloud.H"
#include "SLGThermo.H"

#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// NOTE fwd declaration needed
// for operator fwd declaration
template<class Type>
class myCoalParcel;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const myCoalParcel<Type>&
);

template<class Type>
class myCoalParcel :
    public Type
{

public:

    // IDs of phases in ReacingParcel phase list (Y)

        static const label GAS;
        static const label LIQ;
        static const label SLD;


    //- Class to hold reacting multiphase particle constant properties
    class constantProperties
    :
        public Type::constantProperties
    {
        // Private data

            demandDrivenEntry<scalar> fixedCp_;

            //- Devolatilisation activation temperature [K]
            demandDrivenEntry<scalar> TDevol_;

            //- Latent heat of devolatilisation [J/kg]
            demandDrivenEntry<scalar> LDevol_;

            //- Fraction of enthalpy retained by parcel due to surface
            //  reactions
            demandDrivenEntry<scalar> hRetentionCoeff_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access

            inline scalar fixedCp() const;

            //- Return const access to the devolatilisation temperature
            inline scalar TDevol() const;

            //- Return const access to the latent heat of devolatilisation
            inline scalar LDevol() const;

            //- Return const access to the fraction of enthalpy retained by
            //  parcel due to surface reactions
            inline scalar hRetentionCoeff() const;
    };

    private:

        scalar Y_ash; // kg_ash/kg_dry

        scalar Y_char; // kg_char/kg_dry

        scalarField ddtSolid;

        scalarField ddtLiquid;

        scalarField ddtGas;

        scalar ddtTemp;

    // Private Member Functions

        //- Return the mixture effective specific heat capacity
        template<class TrackData>
        scalar CpEff
        (
            TrackData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective sensible enthalpy
        template<class TrackData>
        scalar HsEff
        (
            TrackData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective latent heat
        template<class TrackData>
        scalar LEff
        (
            TrackData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Update the mass fractions (Y, YGas, YLiquid, YSolid)
        scalar updateMassFractions
        (
            const scalar mass0,
            const scalarField& dMassGas,
            const scalarField& dMassLiquid,
            const scalarField& dMassSolid
        );


protected:

    // Protected data

        // Parcel properties

            //- Mass fractions of gases []
            scalarField YGas_;

            //- Mass fractions of liquids []
            scalarField YLiquid_;

            //- Mass fractions of solids []
            scalarField YSolid_;

            //- Flag to identify if the particle can devolatilise and combust
            //  Combustion possible only after volatile content falls below
            //  threshold value.  States include:
            //  0 = can devolatilise, cannot combust but can change
            //  1 = can devolatilise, can combust
            // -1 = cannot devolatilise or combust, and cannot change
            label canCombust_;


    // Protected Member Functions

        //- Calculate Devolatilisation
        template<class TrackData>
        void calcDevolatilisation
        (
            TrackData& td,
            const scalar dt,           // timestep
            const scalar age,          // age
            const scalar Ts,           // surface temperature
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const scalar mass0,        // mass (initial on injection)
            const scalarField& YGasEff,// gas component mass fractions
            const scalarField& YLiquidEff,// liquid component mass fractions
            const scalarField& YSolidEff,// solid component mass fractions
            label& canCombust,          // 'can combust' flag
            scalarField& dMassDV,      // mass transfer - local to particle
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        ) const;

        //- Calculate surface reactions
        template<class TrackData>
        void calcSurfaceReactions
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label canCombust,     // 'can combust' flag
            const scalar N,            // flux of species emitted from particle
            const scalarField& YMix,   // mixture mass fractions
            const scalarField& YGas,   // gas-phase mass fractions
            const scalarField& YLiquid,// liquid-phase mass fractions
            const scalarField& YSolid, // solid-phase mass fractions
            scalarField& dMassSRGas,   // gas-phase mass transfer - local
            scalarField& dMassSRLiquid,// liquid-phase mass transfer - local
            scalarField& dMassSRSolid, // solid-phase mass transfer - local
            scalarField& dMassSRCarrier, // carrier phase mass transfer
            scalar& Sh,                // explicit particle enthalpy source
            scalar& dhsTrans           // sensible enthalpy transfer to carrier
        ) const;


    public:

    // Constructors

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const myCoalParcel<Type>&
        );

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline myCoalParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );


        //- Construct from components
        inline myCoalParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const scalarField& Y0,
            const scalarField& YGas0,
            const scalarField& YLiquid0,
            const scalarField& YSolid0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        myCoalParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        myCoalParcel(const myCoalParcel& p);

        //- Construct as a copy
        myCoalParcel
        (
            const myCoalParcel& p,
            const polyMesh& mesh
        );

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new myCoalParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new myCoalParcel(*this, mesh));
        }

        // Inheritance of constructor
        // using Type::Type;

        TypeName("CoalParcel");

        //- String representation of properties
        AddToPropertyList
        (
            Type,
            " nGas(Y1..YN)"
          + " nLiquid(Y1..YN)"
          + " nSolid(Y1..YN)"
        );

    // Member Functions

        // Access

            //- Return const access to mass fractions of gases
            inline const scalarField& YGas() const;

            //- Return const access to mass fractions of liquids
            inline const scalarField& YLiquid() const;

            //- Return const access to mass fractions of solids
            inline const scalarField& YSolid() const;

            //- Return const access to the canCombust flag
            inline label canCombust() const;

        // Edit

            //- Return access to mass fractions of gases
            inline scalarField& YGas();

            //- Return access to mass fractions of liquids
            inline scalarField& YLiquid();

            //- Return access to mass fractions of solids
            inline scalarField& YSolid();

            //- Return access to the canCombust flag
            inline label& canCombust();


        // Main calculation loop

            //- Set cell values
            template<class TrackData>
            void setCellValues
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct cell values using latest transfer information
            template<class TrackData>
            void cellValueSourceCorrection
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            template<class TrackData>
            void calcThermolysis
            (
                TrackData& td,
                const scalar dt,
                const label cellI,
                const scalar d,
                const scalar T,
                const scalar mass,
                const label canCombust,
                const scalar N,
                scalarField& YMix,
                scalarField& YGas,
                scalarField& YLiquid,
                scalarField& YSolid,
                scalarField& dMassSRGas,
                scalarField& dMassSRLiquid,
                scalarField& dMassSRSolid,
                scalarField& dMassSRCarrier,
                scalar& Sh,
                scalar& dhsTrans
             );

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Read - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);

        // template<class TrackData>
        // void calc
        // (
        //     TrackData& td,
        //     const scalar dt,
        //     const label cellI
        // );

        const scalarField getDdtSolid() const {return ddtSolid;}

        const scalarField getDdtLiquid() const {return ddtLiquid;}

        const scalarField getDdtGas() const {return ddtGas;}

        scalar getDdtTemp() const {return ddtTemp;}

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<myCoalParcel<Type> > operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<myCoalParcel<Type> >
                (
                    new myCoalParcel<Type>(mesh_, is, true)
                );
            }
        };
    };

    // typedef myCoalParcel
    // <
    //     ReactingMultiphaseParcel
    //     <
    //         ReactingParcel
    //         <
    //             ThermoParcel
    //             <
    //                 KinematicParcel
    //                 <
    //                     particle
    //                 >
    //             >
    //         >
    //     >
    // > coalParcel;

    typedef myCoalParcel
    <
        ReactingParcel
        <
            ThermoParcel
            <
                KinematicParcel
                <
                    particle
                >
            >
        >
    > coalParcel;

    template<>
    inline bool contiguous<coalParcel>()
    {
        return false;
    }

}

#include "coalParcelI.H"
#include "coalParcel.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
