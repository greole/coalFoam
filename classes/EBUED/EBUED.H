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
    Foam::combustionModels::EBUED

Description
    Partially stirred reactor combustion model.  The model calculates a finite
    rate, based on both turbulence and chemistry time scales.  Depending on
    mesh resolution, the Cmix parameter can be used to scale the turbulence
    mixing time scale.

SourceFiles
    EBUED.C

Based on
    Large-eddy simulation of turbulent combustion using
    different combustion models
    Zhou, LX and Hu, LY and Wang, F
    Fuel 2008

\*---------------------------------------------------------------------------*/

#ifndef EBUED_H
#define EBUED_H

#include LAMINAR
#include "fvCFD.H" // for fvc::grad()
#include "Reaction.H"
#include "basicChemistryModel.H"
#include "reactingMixture.H"
#include "thermoPhysicsTypes.H"
#include "psiThermoCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class EBUED Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class EBUED
:
    public laminar<Type>
{

    // Private data
        // EBUED Coeffs
        const scalar C1_;

        const scalar C2_;

        const scalar C3_;

        const reactingMixture<gasHThermoPhysics>& reactions_;

        // stoiciometric mass based Ox/Fuel ratios
        const scalarField betas;

        //const bool DEBUEDG;
        const PtrList<volScalarField>& Y_;

        volScalarField YProd;
        //- List of fuel reaction rate per Reaction [kg/m3/s]
        PtrList<volScalarField> omegas;

    // Protected Member Functions


    // Private Member Functions

        //- Disallow copy construct
        EBUED(const EBUED&);

        //- Disallow default bitwise assignment
        void operator=(const EBUED&);


public:

    //- Runtime type information
    TypeName("EBUED");


    // Constructors

        //- Construct from components
        EBUED(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~EBUED();


    // Member Functions

        // Return fuel consumption vector per reaction
        PtrList<volScalarField>& getOmegas() {return omegas;};

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix.
            virtual tmp<fvScalarMatrix> R(volScalarField& Y);

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //-  Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

            // Compute the stoichiometric coefficient
            scalar beta(const Reaction<gasHThermoPhysics>& R);
            // construct the list of betas
            scalarField constrBetas(
                const reactingMixture<gasHThermoPhysics>&);

            // wrapper for testing
            void updateR(
                const volScalarField& rho,
                const volScalarField& S);

            //- Write access to chemical source terms
            //  (e.g. for multi-chemistry model)
            void calcR_(const label speciesI);

            // Compute fuel consumption rate of reaction R
            inline volScalarField omega(
                const volScalarField& YFuel,
                const volScalarField& YOx,
                const volScalarField& YProd,
                const volScalarField& rho,
                const volScalarField& S,
                const scalar C1,
                const scalar C2,
                const scalar C3,
                const scalar beta
            );

            PtrList<volScalarField>& constrOmegas(
                const reactingMixture<gasHThermoPhysics>& Rs,
                const PtrList<volScalarField>& Y,
                const volScalarField& rho,
                const volScalarField& S
            );

            const reactingMixture<gasHThermoPhysics>& Reactions(){return reactions_;}

            volScalarField RR(label speciesI) const;

            volScalarField& computeProdsMassFraction(
                const Reaction<gasHThermoPhysics>& R,
                const PtrList<volScalarField>& Y
            );
    // I-O

            //- Update properties from given dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "EBUED.C"

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
