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
    Foam::sigmaModels::kohSigma

Description
    Ionic conductivity of aqueous KOH electrolyte using the Gilliam et al.
    (2007) correlation as a function of molarity w [mol/L] and temperature
    T [K]:

        sigma [S/m] = 100 * (
            -2.041*w - 0.0028*w^2 + 0.005332*w*T + 207.2*w/T
            + 0.001043*w^3 - 0.0000003*w^2*T^2 )

    The factor of 100 converts from S/cm (the native unit of the correlation
    coefficients as published) to S/m (OpenFOAM SI units).

    An optional diaphragm porosity/tortuosity correction is applied when
    epsilon and tau are present in the dictionary:

        sigma_eff = sigma * (epsilon / tau)

    Valid range: w = 0-12 mol/L, T = 273-373 K. Values outside this window
    are clipped to a small positive floor (1e-6 S/m).

    Dictionary example (free electrolyte zone):
    \verbatim
        sigmaModel  kohSigma;
        kohSigmaCoeffs
        {
            cellZone  ccl;
            w         6.9;     // mol/L (~30 wt% KOH)
            T         T;       // temperature field name (default "T")
        }
    \endverbatim

    Dictionary example (diaphragm zone with Bruggeman-type correction):
    \verbatim
        sigmaModel  kohSigma;
        kohSigmaCoeffs
        {
            cellZone  membrane;
            w         6.9;
            T         T;
            epsilon   0.55;    // diaphragm porosity [-]
            tau       1.35;    // tortuosity [-]; sigma_eff = sigma*epsilon/tau
        }
    \endverbatim

SourceFiles
    kohSigma.C

References
    Gilliam et al. (2007), Int. J. Hydrogen Energy, 32, 359-364.

\*---------------------------------------------------------------------------*/

#ifndef kohSigma_H
#define kohSigma_H

#include "sigmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{

/*---------------------------------------------------------------------------*\
                          Class kohSigma Declaration
\*---------------------------------------------------------------------------*/

class kohSigma
:
    public sigmaModel
{
    // Private data

        //- KOH molarity [mol/L] (constant in v1)
        scalar w_;

        //- Temperature field name
        word TName_;

        //- Diaphragm porosity [-] (1.0 if not set -> no correction)
        scalar epsilon_;

        //- Diaphragm tortuosity [-] (1.0 if not set -> no correction)
        scalar tau_;


public:

    //- Runtime type information
    TypeName("kohSigma");


    // Constructors

        //- Construct from components
        kohSigma
        (
            const fvMesh& mesh,
            const dictionary& sigmaDictionary
        );


    //- Destructor
    virtual ~kohSigma();


    // Member Functions

        //- Correct the sigma field over the cellZone
        virtual void correct(volScalarField& sigmaField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sigmaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
