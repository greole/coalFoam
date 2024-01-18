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
    Foam::thermolysisModel

Description
    Templated surface reaction model class

SourceFiles
    thermolysisModel.C
    thermolysisModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef thermolysisModel_H
#define thermolysisModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class thermolysisModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermolysisModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //-Runtime type information
    TypeName("thermolysisModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ThermolysisModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        ThermolysisModel(CloudType& owner);

        //- Construct from dictionary
        ThermolysisModel
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& type
        );

        //- Construct copy
        ThermolysisModel(const ThermolysisModel<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<ThermolysisModel<CloudType> > clone() const
        {
            return autoPtr<ThermolysisModel<CloudType> >
            (
                new ThermolysisModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ThermolysisModel();


    //- Selector
    static autoPtr<ThermolysisModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& cloud
    );


    // Member Functions

        //- Update surface reactions
        //  Returns the heat of reaction
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

        //- Add to devolatilisation mass
        // void addTothermolysisMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermolysisModel(CloudType)                      \
                                                             \
    typedef CloudType::myreactingMultiphaseCloudType         \
        myreactingMultiphaseCloudType;                       \
    defineNamedTemplateTypeNameAndDebug                      \
    (                                                        \
        ThermolysisModel<CloudType>,                         \
        0                                                    \
    );                                                       \
    defineTemplateRunTimeSelectionTable                      \
    (                                                        \
        ThermolysisModel<CloudType>,                         \
        dictionary                                           \
    );


#define makeThermolysisModelType(SS, CloudType)              \
                                                             \
    typedef CloudType::myreactingMultiphaseCloudType         \
        myreactingMultiphaseCloudType;                       \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);   \
                                                             \
    ThermolysisModel<CloudType>::                            \
        adddictionaryConstructorToTable<SS<CloudType> >      \
        add##SS##CloudType##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
#   include "thermolysisModel.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
