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

#include "kohSigma.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{
    defineTypeNameAndDebug(kohSigma, 0);

    addToRunTimeSelectionTable
    (
        sigmaModel,
        kohSigma,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmaModels::kohSigma::kohSigma
(
    const fvMesh& mesh,
    const dictionary& sigmaDictionary
)
:
    sigmaModel(mesh, sigmaDictionary),
    w_(sigmaDictionary_.get<scalar>("w")),
    TName_(sigmaDictionary_.getOrDefault<word>("T", "T")),
    epsilon_(sigmaDictionary_.getOrDefault<scalar>("epsilon", 1.0)),
    tau_(sigmaDictionary_.getOrDefault<scalar>("tau", 1.0))
{
    // Sanity guard: warn if epsilon/tau are both at defaults but user may
    // have intended a diaphragm correction (tau=0 would divide by zero).
    if (tau_ <= 0)
    {
        FatalErrorInFunction
            << "kohSigma: tau (tortuosity) must be > 0; got tau = " << tau_
            << exit(FatalError);
    }

    // Sanity check: verify correlation gives physically reasonable result
    // at the reference point w=6.9 mol/L (~30 wt% KOH), T=353.15 K (80 C).
    // Expected: ~100-140 S/m.  A value far outside this window indicates a
    // unit-convention error (the S/cm vs S/m 100x trap).
    const scalar wRef = 6.9;
    const scalar Tref = 353.15;
    const scalar sigRef = 100.0 *
    (
        -2.041*wRef
        - 0.0028*wRef*wRef
        + 0.005332*wRef*Tref
        + 207.2*wRef/Tref
        + 0.001043*wRef*wRef*wRef
        - 0.0000003*wRef*wRef*Tref*Tref
    );

    if (sigRef < 50.0 || sigRef > 300.0)
    {
        FatalErrorInFunction
            << "kohSigma sanity check FAILED: Gilliam correlation at "
            << "w=6.9 mol/L, T=353.15 K gives " << sigRef
            << " S/m; expected ~100-140 S/m.  Check unit convention."
            << exit(FatalError);
    }
    else
    {
        Info << "kohSigma: sanity check OK  "
             << "(w=6.9 mol/L, T=353.15 K -> sigma = "
             << sigRef << " S/m, expected ~100-140 S/m)" << nl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigmaModels::kohSigma::~kohSigma()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigmaModels::kohSigma::correct
(
    volScalarField& sigmaField
) const
{
    const scalarField& T =
        mesh_.lookupObject<volScalarField>(TName_).internalField();

    // Porosity/tortuosity correction factor (1.0 if no diaphragm params set)
    const scalar porFactor = epsilon_ / tau_;

    // Floor value: keep sigma positive outside the validity window
    const scalar sigFloor = 1.0e-6;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const scalar Ti = T[cellI];

            // Gilliam et al. (2007): correlation coefficients are in S/cm;
            // multiply by 100 to obtain S/m.
            const scalar sigSm = 100.0 *
            (
                -2.041*w_
                - 0.0028*w_*w_
                + 0.005332*w_*Ti
                + 207.2*w_/Ti
                + 0.001043*w_*w_*w_
                - 0.0000003*w_*w_*Ti*Ti
            );

            sigmaField[cellI] = Foam::max(sigSm * porFactor, sigFloor);
        }
    }
}

// ************************************************************************* //
