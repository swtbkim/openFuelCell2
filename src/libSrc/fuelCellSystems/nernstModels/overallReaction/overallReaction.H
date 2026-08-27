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

Class
    Foam::overallReaction.C

Description
    Cell EMF from ONE overall water-splitting reaction, split across the two
    electrodes, instead of composing two independently-referenced half-cells.

    The standard nernstModel builds a separate half-cell potential per
    electrode (nernst_air, nernst_fuel) from each species' own janaf absolute
    enthalpy zero, and the cell EMF emerges as nernst_air - nernst_fuel
    (ButlerVolmer.C: eta = phiE - phiI - nernst, applied per electrode). For
    alkaline water electrolysis the gas-phase O2 (air) and H2 (fuel) janaf
    references do not share a common elemental zero, so the two half-cells sum
    to ~2.8 V instead of the ~1.23 V overall water-splitting Gibbs (see
    docs/benchmarks/B1_results.md "EMF INVESTIGATION").

    This model fixes that by computing the overall-reaction EMF once from a
    single, internally-consistent thermodynamic basis:

    \f[
        E_{cell}(T) = \frac{\Delta H^0 - T\,\Delta S^0}{n F}
    \f]

    for the overall reaction H2O(l) -> H2(g) + 1/2 O2(g), with dictionary-
    supplied deltaH0/deltaS0 (defaults for liquid-water splitting:
    deltaH0 = 285830 J/mol, deltaS0 = 163.3 J/mol/K, n = 2), giving
    E_cell ~ 1.18-1.23 V over 298-353 K.

    Partition across electrodes (the consuming code uses nernst per electrode
    and the cell EMF = nernst_air - nernst_fuel):

    \f[
        nernst = split\,E_{cell}
               + \frac{R T}{n F}\ln\prod_i X_i^{w_i}
    \f]

    where the concentration/pressure term is identical to the standard model
    (gas species in the local RxnList, n = signed RxnList["e"]). The chosen
    partition is split = 1 on the air (O2/anode) electrode and split = 0 on the
    fuel (H2/cathode) electrode, so that

        nernst_air - nernst_fuel = E_cell + (concentration corrections)

    matches the overall EMF exactly while each electrode keeps its own
    concentration correction on its own gas species. split is a dictionary
    parameter so other conventions are possible, but the air=1 / fuel=0 split
    is the default and the validated choice for AWE; in electrolysis it drives
    each electrode's eta = phiE - phiI - nernst with the correct sign.

SourceFiles
    overallReaction.C

\*---------------------------------------------------------------------------*/

#ifndef overallReaction_H
#define overallReaction_H

#include "NernstModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseModel;

namespace nernstModels
{

/*---------------------------------------------------------------------------*\
                        Class overallReaction Declaration
\*---------------------------------------------------------------------------*/
template<class Thermo, class OtherThermo>
class overallReaction
:
    public NernstModel<Thermo, OtherThermo>
{
    //- Phase 1
    const phaseModel& phase1_;

    //- Phase 2
    const phaseModel& phase2_;

    //- Overall-reaction enthalpy change [J/mol] (liquid-water splitting)
    scalar deltaH0_;

    //- Overall-reaction entropy change [J/mol/K]
    scalar deltaS0_;

    //- Number of electrons in the overall reaction
    scalar nElectrons_;

    //- Fraction of the cell EMF assigned to this electrode
    scalar split_;

public:
    //- Runtime type information
    TypeName("overallReaction");

    // Constructors

    overallReaction
    (
        const phaseModel& phase1,
        const phaseModel& phase2,
        const dictionary& dict
    );

    //- Destructor
        virtual ~overallReaction();


    // Member Functions

        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace nernstModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "overallReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
