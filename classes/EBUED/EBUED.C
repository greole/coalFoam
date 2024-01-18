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

#include "EBUED.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
EBUED<Type>::EBUED
(
    const word& modelType,
    const fvMesh& mesh
)
:
    laminar<Type>(modelType, mesh),
    C1_(readScalar(this->coeffs().lookup("C1"))),
    C2_(readScalar(this->coeffs().lookup("C2"))),
    C3_(readScalar(this->coeffs().lookup("C3"))),
    reactions_(dynamic_cast<
        const reactingMixture<gasHThermoPhysics>&>(
            this->chemistryPtr_->thermo())),
    betas(constrBetas(reactions_)),
    Y_(this->thermo().composition().Y()),
    YProd(Y_[0])
{
    omegas.setSize(reactions_.size());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
EBUED<Type>::~EBUED()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
volScalarField& EBUED<Type>::computeProdsMassFraction(
    const Reaction<gasHThermoPhysics>& R,
    const PtrList<volScalarField>& Y
) {
    // Product concentrations are determining the overall rates
    // instead of temperature, if we have multiple products and
    // reactions forming these products one reference product
    // concentration has to be found. Thus the maximum is taken
    // fuel consumption vector
    //volScalarField YProds(Y[0]);
    YProd = YProd*0.0;

    const label nP = R.rhs().size(); // Number of product species
    for (label s = 0; s < nP; s++) {
        const label speciesIndex = R.rhs()[s].index;
        YProd += Y[speciesIndex];
        // if (speciesIndex != inertIndex) { // skip if inert
        // }
    }
    return YProd;
}

template<class Type>
Foam::volScalarField EBUED<Type>::omega(
    const volScalarField& YFuel,
    const volScalarField& YOx,
    const volScalarField& YProd,
    const volScalarField& rho,
    const volScalarField& S,
    const scalar C1,
    const scalar C2,
    const scalar C3,
    const scalar beta
)
{
    scalar rbeta = 1/beta;
    volScalarField min1 = min(YFuel, YOx*rbeta);
    volScalarField min2 = min(min1, C2*YProd/(1+beta));

    // compute maximum amount of products that can be formed
    volScalarField maxYProd = min((1+beta)*YFuel, (1+rbeta)*YOx) + YProd;

    // xi is the ratio of current to maximum formable products
    volScalarField xi = -YProd*C3/(maxYProd+SMALL);

    return C1*rho*S*min2*exp(xi);
}

template<class Type>
Foam::scalar EBUED<Type>::beta(
    const Reaction<gasHThermoPhysics>& R
)
{
    const scalar nFuel = R.lhs()[0].stoichCoeff;
    const scalar nOx = R.lhs()[1].stoichCoeff;
    // first species on lhs side is assumed to be the fuel
    const scalar fuelIndex = R.lhs()[0].index;
    // second species on lhs side is assumed to be the ox
    const scalar oxIndex = R.lhs()[1].index;
    const scalar Wfuel = this->thermo().composition().W(fuelIndex);
    const scalar Wox = this->thermo().composition().W(oxIndex);
    // mass stoich. coefficient fuel
    scalar sFuel = nFuel * Wfuel;
    // mass stoich. coefficient Ox
    scalar sOx   = nOx * Wox;
    return sOx/sFuel;
}

template<class Type>
Foam::scalarField EBUED<Type>::constrBetas(
    const reactingMixture<gasHThermoPhysics>& Rs)
{
    scalarField betas_ (Rs.size(), 0);

    forAll(Rs, i) {betas_[i] = beta(Rs[i]);}

    return betas_;
}

template<class Type>
inline PtrList<volScalarField>& EBUED<Type>::constrOmegas(
    const reactingMixture<gasHThermoPhysics>& Rs,
    const PtrList<volScalarField>& Y,
    const volScalarField& rho,
    const volScalarField& S
) {
    forAll(Rs, i) {
        auto& R = Rs[i];
        const scalar fuelIndex = R.lhs()[0].index;
        const scalar oxIndex = R.lhs()[1].index;
        const volScalarField& YProds = computeProdsMassFraction(R, Y);
        auto omeg = omega(Y[fuelIndex], Y[oxIndex], YProds, rho, S, C1_, C2_, C3_, betas[i]);
        // omegas.set(i, &omeg);
        omegas.set(i, new volScalarField(omeg));
    }
    return omegas;
}

template<class Type>
void EBUED<Type>::correct()
{
    if (this->active()) {
        auto& rho =  this->turbulence().rho();
        auto& U =  this->turbulence().U();
        const auto S = mag(symm(fvc::grad(U)))*sqrt(2.0);
        updateR(rho, S);
    }
}


// for simpler testing an extra method has been created
// where an specific S can be passed as argument
// omegas and the RR matrix are updated as side effect,
template<class Type>
void EBUED<Type>::updateR(
    const volScalarField& rho,
    const volScalarField& S)
{
    // NOTE constrOmegas changes omegas as side effect
    constrOmegas(reactions_, Y_, rho, S);
    // Info << "correct()" << endl;
    forAll(Y_,i) {
        const label specieI = this->thermo().composition().species()[Y_[i].name()];
        calcR_(specieI);
    }
}

// // FIXME update RR_
template<class Type>
void EBUED<Type>::calcR_(const label speciesI)
{
    auto W = this->thermo().composition().W(speciesI);
    auto& RR_ = this->chemistryPtr_->RR(speciesI);
    RR_ *= 0.0;
    forAll(reactions_, Ri) {
        scalar nFwd = 0;
        scalar nRev = 0;

        auto& R = reactions_[Ri];
        const scalar fuelIndex = R.lhs()[0].index;
        auto Wfuel = this->thermo().composition().W(fuelIndex);

        forAll(R.lhs(), i) {
            if (R.lhs()[i].index == speciesI) {nFwd += R.lhs()[i].stoichCoeff;}
        }

        forAll(R.rhs(), i) {
            if (R.rhs()[i].index == speciesI) {nRev += R.rhs()[i].stoichCoeff;}
        }
        // kg/m3/s = mol_i/mol_fuel * (kg/mol)/(kg_f/mol_f) * omega
        RR_ += (-nFwd + nRev)*W/Wfuel *  omegas[Ri];
    }
}

// TODO never called ?
template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
EBUED<Type>::R(volScalarField& Y)
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    fvScalarMatrix& Su = tSu();

    if (this->active())
    {
        // Info << "R()" << endl;
        const label specieI = this->thermo().composition().species()[Y.name()];
        Su += this->chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
EBUED<Type>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
EBUED<Type>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );


    if (this->active())
    {
        scalarField& Sh = tSh();

        forAll(Y_, i)
        {
            forAll(Sh, cellI)
            {
                const scalar hi = this->thermo().composition().Hc(i);
                Sh[cellI] -= hi*this->chemistryPtr_->RR(i)[cellI];
            }
        }
    }

    return tSh;
}


template<class Type>
bool EBUED<Type>::read()
{
    // if (laminar<Type>::read())
    // {
    //     this->coeffs().lookup("C1") >> C1_;
    //     this->coeffs().lookup("C2") >> C2_;
    //     this->coeffs().lookup("turbulentReaction") >> turbulentReaction_;
    //     return true;
    // }
    // else
    // {
    //     return false;
    // }
    return true;
}


// ************************************************************************* //
