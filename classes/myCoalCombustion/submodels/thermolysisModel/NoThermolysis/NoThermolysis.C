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

\*---------------------------------------------------------------------------*/

#include "NoThermolysis.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoThermolysis<CloudType>::NoThermolysis
(
    const dictionary&,
    CloudType& owner
)
:
    ThermolysisModel<CloudType>(owner)
{}


template<class CloudType>
Foam::NoThermolysis<CloudType>::NoThermolysis
(
    const NoThermolysis<CloudType>& srm
)
:
    ThermolysisModel<CloudType>(srm.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoThermolysis<CloudType>::~NoThermolysis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoThermolysis<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::NoThermolysis<CloudType>::calculate
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
    // do nothing
    return 0.0;
}


// ************************************************************************* //
