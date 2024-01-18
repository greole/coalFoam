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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::myCoalParcel<Type>::constantProperties::
constantProperties()
:
    Type::constantProperties(),
    fixedCp_(this->dict_, 0.0),
    TDevol_(this->dict_, 0.0),
    LDevol_(this->dict_, 0.0),
    hRetentionCoeff_(this->dict_, 0.0)
{}


template<class Type>
inline Foam::myCoalParcel<Type>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    Type::constantProperties(cp),
    fixedCp_(cp.fixedCp_),
    TDevol_(cp.TDevol_),
    LDevol_(cp.LDevol_),
    hRetentionCoeff_(cp.hRetentionCoeff_)
{}


template<class Type>
inline Foam::myCoalParcel<Type>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    Type::constantProperties(parentDict),
    fixedCp_(this->dict_, "fixedCp"),
    TDevol_(this->dict_, "TDevol"),
    LDevol_(this->dict_, "LDevol"),
    hRetentionCoeff_(this->dict_, "hRetentionCoeff")
{}


template<class Type>
inline Foam::myCoalParcel<Type>::myCoalParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    Type(mesh, position, cellI, tetFaceI, tetPtI),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{
}


template<class Type>
inline Foam::myCoalParcel<Type>::myCoalParcel
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
)
:
    Type
    (
        mesh,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
        constProps
    ),
    YGas_(YGas0),
    YLiquid_(YLiquid0),
    YSolid_(YSolid0),
    canCombust_(0)
{
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class Type>
inline Foam::scalar
Foam::myCoalParcel<Type>::constantProperties::fixedCp() const
{
    return fixedCp_.value();
}

template<class Type>
inline Foam::scalar
Foam::myCoalParcel<Type>::constantProperties::TDevol() const
{
    return TDevol_.value();
}


template<class Type>
inline Foam::scalar
Foam::myCoalParcel<Type>::constantProperties::LDevol() const
{
    return LDevol_.value();
}


template<class Type>
inline Foam::scalar
Foam::myCoalParcel<Type>::constantProperties::
hRetentionCoeff() const
{
    scalar value = hRetentionCoeff_.value();

    if ((value < 0) || (value > 1))
    {
        FatalErrorIn
        (
            "myCoalParcel<Type>::constantProperties::"
            "constantProperties"
        )   << "hRetentionCoeff must be in the range 0 to 1" << nl
            << exit(FatalError) << endl;
    }

    return value;
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class Type>
inline const Foam::scalarField& Foam::myCoalParcel<Type>::
YGas() const
{
    return YGas_;
}


template<class Type>
inline const Foam::scalarField& Foam::myCoalParcel<Type>::
YLiquid() const
{
    return YLiquid_;
}


template<class Type>
inline const Foam::scalarField& Foam::myCoalParcel<Type>::
YSolid() const
{
    return YSolid_;
}


template<class Type>
inline Foam::label
Foam::myCoalParcel<Type>::canCombust() const
{
    return canCombust_;
}


template<class Type>
inline Foam::scalarField& Foam::myCoalParcel<Type>::YGas()
{
    return YGas_;
}


template<class Type>
inline Foam::scalarField& Foam::myCoalParcel<Type>::YLiquid()
{
    return YLiquid_;
}


template<class Type>
inline Foam::scalarField& Foam::myCoalParcel<Type>::YSolid()
{
    return YSolid_;
}


template<class Type>
inline Foam::label& Foam::myCoalParcel<Type>::canCombust()
{
    return canCombust_;
}


// ************************************************************************* //
