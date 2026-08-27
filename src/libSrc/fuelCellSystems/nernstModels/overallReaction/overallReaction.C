/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
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

#include "overallReaction.H"
#include "phaseModel.H"

#include "constants.H"

// Rgas and F are defined at file scope in standard.C, which is included into
// the same translation unit (NernstModels.C) ahead of this file (matching how
// fixedValue.C reuses them).

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::nernstModels::overallReaction<Thermo, OtherThermo>::overallReaction
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const dictionary& dict
)
:
    NernstModel<Thermo, OtherThermo>(phase1, phase2, dict),
    phase1_(phase1),
    phase2_(phase2),
    deltaH0_(dict.lookupOrDefault<scalar>("deltaH0", 285830.0)),
    deltaS0_(dict.lookupOrDefault<scalar>("deltaS0", 163.3)),
    nElectrons_(dict.lookupOrDefault<scalar>("nElectrons", 2.0)),
    split_(dict.lookupOrDefault<scalar>("split", 1.0))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::nernstModels::overallReaction<Thermo, OtherThermo>::~overallReaction()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::nernstModels::overallReaction<Thermo, OtherThermo>::correct()
{
    this->deltaH() *= 0.0;
    this->deltaS() *= 0.0;

    const word water = phaseModel::water;

    scalarField& nernst = this->operator()();
    scalarField Qrxn(nernst.size(), 1.0);

    const scalarField& p = this->thermo_.p();
    const scalarField& T = this->thermo_.T();

    scalarField pRef = p/this->pRef().value();

    //- Concentration/pressure correction on the local gas species only
    //- (identical form to the standard/fixedValue models). The overall-
    //- reaction EMF below replaces the half-cell enthalpy/entropy sums, so
    //- no per-species formation enthalpy is accumulated here.
    forAllConstIter(HashTable<scalar>, this->rxnList(), iter)
    {
        const word& nameI = iter.key();

        if (nameI != "e" && nameI != water)
        {
            scalar stoiCoeffI = this->rxnList()[nameI];

            const scalarField& X = phase1_.X(nameI);

            Qrxn *= Foam::pow(Foam::max(X, this->residualY())*pRef, stoiCoeffI);
        }
    }

    //- Overall water-splitting EMF on a single common basis:
    //-   E_cell(T) = (deltaH0 - T*deltaS0)/(n*F)
    //- split across electrodes by the dictionary fraction split_, so that
    //- nernst_air - nernst_fuel = E_cell + concentration corrections.
    const scalar n = this->rxnList()["e"];

    forAll(nernst, cellI)
    {
        const scalar Ecell =
            (deltaH0_ - T[cellI]*deltaS0_)/(nElectrons_*F);

        nernst[cellI] =
            split_*Ecell
          + Rgas*T[cellI]*Foam::log(Qrxn[cellI])/n/F;
    }

    Info<< "Nernst " << this->operator()().mesh().name()
        << ": min = " << Foam::min(this->operator()().primitiveField())
        << ", mean = " << Foam::average(this->operator()().primitiveField())
        << ", max = " << Foam::max(this->operator()().primitiveField())
        << endl;
}

// ************************************************************************* //
