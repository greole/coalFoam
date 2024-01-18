#include "coalCloud.H"
#include TURBMODEL
#include LESMODEL

namespace Foam
{


template<class CloudType>
void Foam::myReactingMultiphaseCloud<CloudType>::setModels()
{
    thermolysisModel_.reset
    (
        ThermolysisModel<myReactingMultiphaseCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}

template<class CloudType>
void Foam::myReactingMultiphaseCloud<CloudType>::cloudReset
(
    myReactingMultiphaseCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    thermolysisModel_.reset(c.thermolysisModel_.ptr());
}

template<class CloudType>
void Foam::myReactingMultiphaseCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<myReactingMultiphaseCloud<CloudType> > td(*this);

        this->solve(td);
    }
}


template<class CloudType>
Foam::myReactingMultiphaseCloud<CloudType>::myReactingMultiphaseCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, readFields),
    cloudCopyPtr_(NULL),
    constProps_(this->particleProperties()),
    thermolysisModel_(NULL)
{
    Info << "myReactingMultiphaseCloud1" << endl;
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        this->resetSourceTerms();
    }
}

template<class CloudType>
Foam::myReactingMultiphaseCloud<CloudType>::myReactingMultiphaseCloud
(
    myReactingMultiphaseCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    cloudCopyPtr_(NULL),
    constProps_(c.constProps_),
    thermolysisModel_(c.thermolysisModel_->clone())
{
    Info << "myReactingMultiphaseCloud2" << endl;
}


template<class CloudType>
Foam::myReactingMultiphaseCloud<CloudType>::myReactingMultiphaseCloud
(
    const fvMesh& mesh,
    const word& name,
    const myReactingMultiphaseCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    cloudCopyPtr_(NULL),
    constProps_(),
    thermolysisModel_(NULL)
{
    Info << "myReactingMultiphaseCloud2" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myReactingMultiphaseCloud<CloudType>::~myReactingMultiphaseCloud()
{}

}

