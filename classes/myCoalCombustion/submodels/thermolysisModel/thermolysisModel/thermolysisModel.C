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

#include "thermolysisModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermolysisModel<CloudType>::ThermolysisModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner),
    dMass_(0.0)
{}


template<class CloudType>
Foam::ThermolysisModel<CloudType>::ThermolysisModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    dMass_(0.0)
{}


template<class CloudType>
Foam::ThermolysisModel<CloudType>::ThermolysisModel
(
    const ThermolysisModel<CloudType>& srm
)
:
    CloudSubModelBase<CloudType>(srm),
    dMass_(srm.dMass_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermolysisModel<CloudType>::~ThermolysisModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ThermolysisModel<CloudType>::calculate
(
    const scalar,
    const label,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    scalarField&,
    scalarField&,
    scalarField&,
    scalarField&,
    const scalar,
    scalarField&,
    scalarField&,
    scalarField&,
    scalarField&,
    scalarField&
) const
{
    notImplemented
    (
        "Foam::scalar Foam::thermolysisModel<CloudType>::calculate"
        "("
            "const scalar, "
            "const label, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalarField&, "
            "const scalarField&, "
            "const scalarField&, "
            "const scalarField&, "
            "const scalar, "
            "scalarField&, "
            "scalarField&, "
            "scalarField&, "
            "scalarField&"
        ") const"
    );

    return 0.0;
}


// template<class CloudType>
// void Foam::thermolysisModel<CloudType>::addTothermolysisMass
// (
//     const scalar dMass
// )
// {
//     dMass_ += dMass;
// }


template<class CloudType>
void Foam::ThermolysisModel<CloudType>::info(Ostream& os)
{
    // const scalar mass0 = this->template getBaseProperty<scalar>("mass");
    // const scalar massTotal = mass0 + returnReduce(dMass_, sumOp<scalar>());
    //
    // Info<< "    Mass transfer surface reaction  = " << massTotal << nl;
    //
    // if (this->outputTime())
    // {
    //     this->setBaseProperty("mass", massTotal);
    //     dMass_ = 0.0;
    // }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermolysisModelNew.C"

// ************************************************************************* //

