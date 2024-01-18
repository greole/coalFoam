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

#include "coalParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::string Foam::myCoalParcel<Type>::propertyList_ =
    Foam::myCoalParcel<Type>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::myCoalParcel<Type>::myCoalParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    Type(mesh, is, readFields),
    Y_ash(0),
    Y_char(0),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{
    // NOTE used for parallel transfer
    if (readFields)
    {
        // Pout << "Ptransfer " << endl;
        DynamicList<scalar> Yg;
        DynamicList<scalar> Yl;
        DynamicList<scalar> Ys;

        is >> Yg >> Yl >> Ys;

        YGas_.transfer(Yg);
        YLiquid_.transfer(Yl);
        YSolid_.transfer(Ys);

        // scale the mass fractions
        const scalarField& YMix = this->Y_;
        YGas_ /= YMix[GAS] + ROOTVSMALL;
        YLiquid_ /= YMix[LIQ] + ROOTVSMALL;
        YSolid_ /= YMix[SLD] + ROOTVSMALL;
    }

    // Check state of Istream
    is.check
    (
        "myCoalParcel<Type>::myCoalParcel"
        "("
            "const polyMesh&, "
            "Istream&, "
            "bool"
        ")"
    );
        // Pout << "[done] " << endl;
}


template<class Type>
template<class CloudType>
void Foam::myCoalParcel<Type>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    Type::readFields(c);
}


template<class Type>
template<class CloudType, class CompositionType>
void Foam::myCoalParcel<Type>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    if (!c.size())
    {
        return;
    }

    Type::readFields(c, compModel);

    // Get names and sizes for each Y...
    const label idGas = compModel.idGas();
    const wordList& gasNames = compModel.componentNames(idGas);
    const label idLiquid = compModel.idLiquid();
    const wordList& liquidNames = compModel.componentNames(idLiquid);
    const label idSolid = compModel.idSolid();
    const wordList& solidNames = compModel.componentNames(idSolid);
    const wordList& stateLabels = compModel.stateLabels();

    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<myCoalParcel<Type> >, c, iter)
    {
        myCoalParcel<Type>& p = iter();
        p.YGas_.setSize(gasNames.size(), 0.0);
        p.YLiquid_.setSize(liquidNames.size(), 0.0);
        p.YSolid_.setSize(solidNames.size(), 0.0);
    }

    // Populate YGas for each parcel
    forAll(gasNames, j)
    {
        IOField<scalar> YGas
        (
            c.fieldIOobject
            (
                "Y" + gasNames[j] + stateLabels[idGas],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<myCoalParcel<Type> >,
            c,
            iter
        )
        {
            myCoalParcel<Type>& p = iter();
            p.YGas_[j] = YGas[i++]/(p.Y()[GAS] + ROOTVSMALL);
        }
    }
    // Populate YLiquid for each parcel
    forAll(liquidNames, j)
    {
        IOField<scalar> YLiquid
        (
            c.fieldIOobject
            (
                "Y" + liquidNames[j] + stateLabels[idLiquid],
                 IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<myCoalParcel<Type> >,
            c,
            iter
        )
        {
            myCoalParcel<Type>& p = iter();
            p.YLiquid_[j] = YLiquid[i++]/(p.Y()[LIQ] + ROOTVSMALL);
        }
    }
    // Populate YSolid for each parcel
    forAll(solidNames, j)
    {
        IOField<scalar> YSolid
        (
            c.fieldIOobject
            (
                "Y" + solidNames[j] + stateLabels[idSolid],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<myCoalParcel<Type> >,
            c,
            iter
        )
        {
            myCoalParcel<Type>& p = iter();
            p.YSolid_[j] = YSolid[i++]/(p.Y()[SLD] + ROOTVSMALL);
        }
    }
}


template<class Type>
template<class CloudType>
void Foam::myCoalParcel<Type>::writeFields(const CloudType& c)
{
    Type::writeFields(c);
}


template<class Type>
template<class CloudType, class CompositionType>
void Foam::myCoalParcel<Type>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    Type::writeFields(c, compModel);

    label np = c.size();

    // Write the composition fractions
    if (np > 0)
    {
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            IOField<scalar> YGas
            (
                c.fieldIOobject
                (
                    "Y" + gasNames[j] + stateLabels[idGas],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<myCoalParcel<Type> >,
                c,
                iter
            )
            {
                const myCoalParcel<Type>& p0 = iter();
                YGas[i++] = p0.YGas()[j]*p0.Y()[GAS];
            }

            YGas.write();
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            IOField<scalar> YLiquid
            (
                c.fieldIOobject
                (
                    "Y" + liquidNames[j] + stateLabels[idLiquid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<myCoalParcel<Type> >,
                c,
                iter
            )
            {
                const myCoalParcel<Type>& p0 = iter();
                YLiquid[i++] = p0.YLiquid()[j]*p0.Y()[LIQ];
            }

            YLiquid.write();
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            IOField<scalar> YSolid
            (
                c.fieldIOobject
                (
                    "Y" + solidNames[j] + stateLabels[idSolid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<myCoalParcel<Type> >,
                c,
                iter
            )
            {
                const myCoalParcel<Type>& p0 = iter();
                YSolid[i++] = p0.YSolid()[j]*p0.Y()[SLD];
            }

            YSolid.write();
        }

        // write ddtYsolids
        IOField<scalarField> ddtYSolid
        (
            c.fieldIOobject
            (
                "ddtYsolid",
                IOobject::NO_READ
            ),
            np
        );

        label is = 0;
        forAllConstIter
        (
            typename Cloud<myCoalParcel<Type> >,
            c,
            iter
        )
        {
            const myCoalParcel<Type>& p0 = iter();
            ddtYSolid[is++] = p0.getDdtSolid();
        }

        ddtYSolid.write();

        // write ddtYGas
        IOField<scalarField> ddtYGas
        (
            c.fieldIOobject
            (
                "ddtYGas",
                IOobject::NO_READ
            ),
            np
        );

        label ig = 0;
        forAllConstIter
        (
            typename Cloud<myCoalParcel<Type> >,
            c,
            iter
        )
        {
            const myCoalParcel<Type>& p0 = iter();
            ddtYGas[ig++] = p0.getDdtGas();
        }

        ddtYGas.write();
        //
        // write ddtYLiquid
        IOField<scalarField> ddtYLiquid
        (
            c.fieldIOobject
            (
                "ddtYLiquid",
                IOobject::NO_READ
            ),
            np
        );

        label il = 0;
        forAllConstIter
        (
            typename Cloud<myCoalParcel<Type> >,
            c,
            iter
        )
        {
            const myCoalParcel<Type>& p0 = iter();
            ddtYGas[il++] = p0.getDdtLiquid();
        }

        ddtYLiquid.write();
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const myCoalParcel<Type>& p
)
{
    // NOTE needed for parallel transfer
    //
    // TODO  add new quantities
    // Pout << "Transfer "
    //     << " p.Y()[0] " << p.Y()[0]
    //     << " p.Y()[1] " << p.Y()[1]
    //     << " p.Y()[2] " << p.Y()[2]
    //     << endl;
    scalarField YGasLoc(p.YGas()*p.Y()[0]);
    scalarField YLiquidLoc(p.YLiquid()*p.Y()[1]);
    scalarField YSolidLoc(p.YSolid()*p.Y()[2]);
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const Type&>(p)
            << token::SPACE << YGasLoc
            << token::SPACE << YLiquidLoc
            << token::SPACE << YSolidLoc;
    }
    else
    {
        os  << static_cast<const Type&>(p);
        os  << YGasLoc << YLiquidLoc << YSolidLoc;
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const myCoalParcel<Type>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
