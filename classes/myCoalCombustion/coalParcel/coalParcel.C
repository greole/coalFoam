#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// template<class Type>
// template<class TrackData>
// void Foam::myCoalParcel<Type>::calc
// (
//     TrackData& td,
//     const scalar dt,
//     const label cellI
// )
// {
//     const scalar rdt = 1/dt;
//     const scalar mass0   = this->mass();
//     scalarField dSolid0  = this->YSolid_*mass0;
//     scalarField dLiquid0 = this->YLiquid_*mass0;
//     scalarField dGas0    = this->YGas_*mass0;
//     scalar Temp0         = this->T_;
//
//     Type::calc(td, dt, cellI);
//
//     const scalar massi = this->mass();
//     ddtSolid  = rdt*(this->YSolid_*massi  - dSolid0);
//     ddtLiquid = rdt*(this->YLiquid_*massi - dLiquid0);
//     ddtGas    = rdt*(this->YGas_*massi    - dGas0);
//     ddtTemp   = rdt*(this->T_ - Temp0);
// }

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::label Foam::myCoalParcel<Type>::GAS(0);

template<class Type>
const Foam::label Foam::myCoalParcel<Type>::LIQ(1);

template<class Type>
const Foam::label Foam::myCoalParcel<Type>::SLD(2);


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class Type>
template<class TrackData>
Foam::scalar Foam::myCoalParcel<Type>::CpEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().Cp(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().Cp(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().Cp(idS, YSolid_, p, T);
}

// template<class Type>
// template<class TrackData>
// Foam::scalar Foam::myCoalParcel<Type>::CpCorr
// (
//     const scalar& T
// ) const
// {
//     // TODO dont hardcode coeffs
//     const scalar A = 1.08;
//     const scalar B = 3.81e-3;
//     const scalar C = 5.48e-7;
//     const scalar D = -7.83e-9;
//     const scalar E = 3.99e-12;
//     const scalar T2 = T*T;
//
//     return A + B*T + C*T2 + D*T2*T + E*T2*T2;
// }
//
//

template<class Type>
template<class TrackData>
Foam::scalar Foam::myCoalParcel<Type>::HsEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().Hs(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().Hs(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().Hs(idS, YSolid_, p, T);
}


template<class Type>
template<class TrackData>
Foam::scalar Foam::myCoalParcel<Type>::LEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().L(idS, YSolid_, p, T);
}


template<class Type>
Foam::scalar Foam::myCoalParcel<Type>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, ROOTVSMALL);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = 1.0 - YMix[GAS] - YMix[LIQ];

    return massNew;
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class Type>
template<class TrackData>
void Foam::myCoalParcel<Type>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    Type::setCellValues(td, dt, cellI);
}


template<class Type>
template<class TrackData>
void Foam::myCoalParcel<Type>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // Re-use correction from reacting parcel
    Type::cellValueSourceCorrection(td, dt, cellI);
}

template<class Type>
template<class TrackData>
void Foam::myCoalParcel<Type>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();

    // Info << "calc" << endl;

    const scalar rdt = 1/dt;
    const scalar mass0_  = this->mass();
    scalarField& YMix    = this->Y_;

    scalarField dSolid0  = this->YSolid_*mass0_*YMix[SLD];
    scalarField dLiquid0 = this->YLiquid_*mass0_*YMix[LIQ];
    scalarField dGas0    = this->YGas_*mass0_*YMix[GAS];
    scalar Temp0         = this->T_;

    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();

    const scalar pc = this->pc_;

    const label idG = composition.idGas();
    const label idL = composition.idLiquid();
    const label idS = composition.idSolid();


    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(td, cellI, T0, Ts, rhos, mus, Prs, kappas);
    scalar Res = this->Re(U0, d0, rhos, mus);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = vector::zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);

    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(composition.carrier().species().size(), 0.0);

    // Calc mass and enthalpy transfer due to phase change
    this->calcPhaseChange
    (
        td,
        dt,
        cellI,
        Res,
        Prs,
        Ts,
        mus/rhos,
        d0,
        T0,
        mass0,
        idL,
        YMix[LIQ],
        YLiquid_,
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs
    );

    // thermolysis
    // ~~~~~~~~~~~~~~~~~
    // Change in carrier phase composition due to surface reactions
    scalarField dMassTHGas(YGas_.size(), 0.0);
    scalarField dMassTHLiquid(YLiquid_.size(), 0.0);
    scalarField dMassTHSolid(YSolid_.size(), 0.0);
    scalarField dMassTHCarrier(composition.carrier().species().size(), 0.0);

    calcThermolysis(
        td,
        dt,
        cellI,
        d0,
        T0,
        mass0,
        canCombust_,
        Ne,
        YMix,
        YGas_,
        YLiquid_,
        YSolid_,
        dMassTHGas,
        dMassTHLiquid,
        dMassTHSolid,
        dMassTHCarrier,
        Sh,
        dhsTrans
    );

    if (YSolid_[2]*YMix[SLD] < 1.0e-15) {
        canCombust_ = true;
    }

    scalarField dMassGas0(dMassTHGas);
    scalarField dMassLiquid0(dMassTHLiquid);
    scalarField dMassSolid0(dMassTHSolid);
    scalar mass01 =
        updateMassFractions(mass0, dMassGas0, dMassLiquid0, dMassSolid0);

    // Devolatilisation
    // ~~~~~~~~~~~~~~~~

    // Mass transfer due to devolatilisation
    scalarField dMassDV(YGas_.size(), 0.0);

    // Calc mass and enthalpy transfer due to devolatilisation
    calcDevolatilisation
    (
        td,
        dt,
        this->age_,
        Ts,
        d0,
        T0,
        mass0,
        this->mass0_,
        YMix[GAS]*YGas_,
        YMix[LIQ]*YLiquid_,
        YMix[SLD]*YSolid_,
        canCombust_,
        dMassDV,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // Change in carrier phase composition due to surface reactions
    scalarField dMassSRGas(YGas_.size(), 0.0);
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField dMassSRCarrier(composition.carrier().species().size(), 0.0);


    // Calc mass and enthalpy transfer due to surface reactions
    calcSurfaceReactions
    (
        td,
        dt,
        cellI,
        d0,
        T0,
        mass0,
        canCombust_,
        Ne,
        YMix,
        YGas_,
        YLiquid_,
        YSolid_,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier,
        Sh,
        dhsTrans
    );


    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField dMassGas(dMassDV + dMassSRGas);
    scalarField dMassLiquid(dMassPC + dMassSRLiquid);
    scalarField dMassSolid(dMassSRSolid);
    scalar mass1 =
        updateMassFractions(mass0, dMassGas, dMassLiquid, dMassSolid);

    scalar fixedCp = td.cloud().constProps().fixedCp();

    if (fixedCp > 0) {
        this->Cp_ = fixedCp;
    }
    if (fixedCp == 0) {
        this->Cp_ = CpEff(td, pc, T0, idG, idL, idS);
    }
    if (fixedCp < 0) {
        const scalar A = 1.08;
        const scalar B = 3.81e-3;
        const scalar C = 5.48e-7;
        const scalar D = -7.83e-9;
        const scalar E = 3.99e-12;

        // Correllation for Deg Celsius
        const scalar Tmod = T0-273.15;

        if (Tmod > 1700.0) {this->Cp_ = 1100.0;}
        else {
            const scalar T2 = Tmod*Tmod;
            this->Cp_ = (A + B*Tmod + C*T2 + D*T2*Tmod + E*T2*T2)*1000.0;
        }
    }

    // Update particle density or diameter
    if (td.cloud().constProps().constantVolume())
    {
        this->rho_ = mass1/this->volume();
    }
    else
    {
        this->d_ = cbrt(mass1/this->rho_*6.0/pi);
    }

    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < td.cloud().constProps().minParcelMass())
    {
        td.keepParticle = false;

        if (td.cloud().solution().coupled())
        {
            scalar dm = np0*mass0;

            // Absorb parcel into carrier phase
            forAll(YGas_, i)
            {
                label gid = composition.localToGlobalCarrierId(GAS, i);
                td.cloud().rhoTrans(gid)[cellI] += dm*YMix[GAS]*YGas_[i];
            }
            forAll(YLiquid_, i)
            {
                label gid = composition.localToGlobalCarrierId(LIQ, i);
                td.cloud().rhoTrans(gid)[cellI] += dm*YMix[LIQ]*YLiquid_[i];
            }
/*
            // No mapping between solid components and carrier phase
            forAll(YSolid_, i)
            {
                label gid = composition.localToGlobalCarrierId(SLD, i);
                td.cloud().rhoTrans(gid)[cellI] += dm*YMix[SLD]*YSolid_[i];
            }
*/
            td.cloud().UTrans()[cellI] += dm*U0;

            td.cloud().hsTrans()[cellI] += dm*HsEff(td, pc, T0, idG, idL, idS);

            td.cloud().phaseChange().addToPhaseChangeMass(np0*mass1);
        }

        return;
    }

    // Correct surface values due to emitted species
    this->correctSurfaceValues(td, cellI, Ts, Cs, rhos, mus, Prs, kappas);
    Res = this->Re(U0, this->d_, rhos, mus);


    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            td,
            dt,
            cellI,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );


    // this->Cp_ = CpEff(td, pc, this->T_, idG, idL, idS);


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(td, dt, cellI, Res, mus, mass1, Su, dUTrans, Spu);


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (td.cloud().solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(YGas_, i)
        {
            scalar dm = np0*dMassGas[i];
            label gid = composition.localToGlobalCarrierId(GAS, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            td.cloud().rhoTrans(gid)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*U0;
            td.cloud().hsTrans()[cellI] += dm*hs;
        }
        forAll(YLiquid_, i)
        {
            scalar dm = np0*dMassLiquid[i];
            label gid = composition.localToGlobalCarrierId(LIQ, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            td.cloud().rhoTrans(gid)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*U0;
            td.cloud().hsTrans()[cellI] += dm*hs;
        }
/*
        // No mapping between solid components and carrier phase
        forAll(YSolid_, i)
        {
            scalar dm = np0*dMassSolid[i];
            label gid = composition.localToGlobalCarrierId(SLD, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            td.cloud().rhoTrans(gid)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*U0;
            td.cloud().hsTrans()[cellI] += dm*hs;
        }
*/
        forAll(dMassSRCarrier, i)
        {
            scalar dm = np0*dMassSRCarrier[i];
            scalar hs = composition.carrier().Hs(i, pc, T0);
            td.cloud().rhoTrans(i)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*U0;
            td.cloud().hsTrans()[cellI] += dm*hs;
        }

        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;
        td.cloud().UCoeff()[cellI] += np0*Spu;

        // Update sensible enthalpy transfer
        td.cloud().hsTrans()[cellI] += np0*dhsTrans;
        td.cloud().hsCoeff()[cellI] += np0*Sph;

        // Update radiation fields
        if (td.cloud().radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(this->T_);
            td.cloud().radAreaP()[cellI] += dt*np0*ap;
            td.cloud().radT4()[cellI] += dt*np0*T4;
            td.cloud().radAreaPT4()[cellI] += dt*np0*ap*T4;
        }
    }

    const scalar massi = this->mass();
    ddtSolid  = rdt*(this->YSolid_*massi*YMix[SLD]  - dSolid0);
    ddtLiquid = rdt*(this->YLiquid_*massi*YMix[LIQ] - dLiquid0);
    ddtGas    = rdt*(this->YGas_*massi*YMix[GAS]    - dGas0);
    ddtTemp   = rdt*(this->T_ - Temp0);
}

template<class Type>
template<class TrackData>
void Foam::myCoalParcel<Type>::calcThermolysis
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
)
{
    scalarField dHs = scalarField(1,0.0);
    // Total mass of volatiles evolved
    td.cloud().thermolysis().calculate
    (
        dt,
        cellI,
        d,
        T,
        this->Tc_,
        this->pc_,
        this->rhoc_,
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier,
        dHs
    );

    scalar dMassTot = sum(dMassSRGas);
    // Pout << "dHs: " << dHs << endl;

    Sh -= sum(dHs)/dt;
}


template<class Type>
template<class TrackData>
void Foam::myCoalParcel<Type>::calcDevolatilisation
(
    TrackData& td,
    const scalar dt,
    const scalar age,
    const scalar Ts,
    const scalar d,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
) const
{
    // Check that model is active
    if (!td.cloud().devolatilisation().active())
    {
        return;
    }

    // Initialise demand-driven constants
    (void)td.cloud().constProps().fixedCp();
    (void)td.cloud().constProps().TDevol();
    (void)td.cloud().constProps().LDevol();

    // Check that the parcel temperature is within necessary limits for
    // devolatilisation to occur
    if (T < td.cloud().constProps().TDevol() || canCombust == -1)
    {
        return;
    }

    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();


    // Total mass of volatiles evolved
    td.cloud().devolatilisation().calculate
    (
        dt,
        age,
        mass0,
        mass,
        T,
        YGasEff,
        YLiquidEff,
        YSolidEff,
        canCombust,
        dMassDV
    );

    scalar dMassTot = sum(dMassDV);

    td.cloud().devolatilisation().addToDevolatilisationMass
    (
        this->nParticle_*dMassTot
    );

    Sh -= dMassTot*td.cloud().constProps().LDevol()/dt;

    // Update molar emissions
    if (td.cloud().heatTransfer().BirdCorrection())
    {
        // Molar average molecular weight of carrier mix
        const scalar Wc =
            max(SMALL, this->rhoc_*specie::RR*this->Tc_/this->pc_);

        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));

        forAll(dMassDV, i)
        {
            const label id = composition.localToGlobalCarrierId(GAS, i);
            const scalar Cp = composition.carrier().Cp(id, this->pc_, Ts);
            const scalar W = composition.carrier().W(id);
            const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);

            // Dab calc'd using API vapour mass diffusivity function
            const scalar Dab =
                3.6059e-3*(pow(1.8*Ts, 1.75))
               *sqrt(1.0/W + 1.0/Wc)
               /(this->pc_*beta);

            N += Ni;
            NCpW += Ni*Cp*W;
            Cs[id] += Ni*d/(2.0*Dab);
        }
    }
}


template<class Type>
template<class TrackData>
void Foam::myCoalParcel<Type>::calcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans
) const
{

    // Check that model is active
    if (!td.cloud().surfaceReaction().active())
    {
        return;
    }

    // Initialise demand-driven constants
    (void)td.cloud().constProps().hRetentionCoeff();
    (void)td.cloud().constProps().TMax();

    // Check that model is active
    if (canCombust != 1)
    {
        return;
    }


    // Update surface reactions
    const scalar hReaction = td.cloud().surfaceReaction().calculate
    (
        dt,
        cellI,
        d,
        T,
        this->Tc_,
        this->pc_,
        this->rhoc_,
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier
    );

    td.cloud().surfaceReaction().addToSurfaceReactionMass
    (
        this->nParticle_
       *(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    );

    const scalar xsi = min(T/td.cloud().constProps().TMax(), 1.0);
    const scalar coeff =
        (1.0 - xsi*xsi)*td.cloud().constProps().hRetentionCoeff();

    Sh += coeff*hReaction/dt;

    dhsTrans += (1.0 - coeff)*hReaction;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::myCoalParcel<Type>::myCoalParcel
(
    const myCoalParcel<Type>& p
)
:
    Type(p),
    Y_ash(p.Y_ash),
    Y_char(p.Y_char),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{
}


template<class Type>
Foam::myCoalParcel<Type>::myCoalParcel
(
    const myCoalParcel<Type>& p,
    const polyMesh& mesh
)
:
    Type(p, mesh),
    Y_ash(p.Y_ash),
    Y_char(p.Y_char),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "coalParcelIO.C"

// ************************************************************************* //
