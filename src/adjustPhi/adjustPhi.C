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

#include "adjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorFvsPatchFields.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::adjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    volScalarField& p
)
{
    if (p.needReference())
    {
        // p coefficients should not be updated here
        // p.boundaryField().updateCoeffs();

        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;

        surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

//            if (!isA<processorFvsPatchScalarField>(phip))
	    if (!phip.coupled())
            {
                if (Up.fixesValue() && !isA<inletOutletFvPatchVectorField>(Up))
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }
                }
                else
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            adjustableMassOut += phip[i];
                        }
                    }
                }
            }
        }

        // Calculate the total flux in the domain, used for normalisation
        scalar totalFlux = vSmall + sum(mag(phi)).value();

        reduce(massIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());

        scalar massCorr = 1.0;
        scalar magAdjustableMassOut = mag(adjustableMassOut);

        if
        (
            magAdjustableMassOut > vSmall
         && magAdjustableMassOut/totalFlux > small
        )
        {
            massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }
        else if (mag(fixedMassOut - massIn)/totalFlux > 1e-8)
        {
            FatalErrorInFunction
                << "Continuity error cannot be removed by adjusting the"
                   " outflow.\nPlease check the velocity boundary conditions"
                   " and/or run potentialFoam to initialise the outflow." << nl
                << "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << massIn << nl
                << "Specified mass outflow  : " << fixedMassOut << nl
                << "Adjustable mass outflow : " << adjustableMassOut << nl
                << exit(FatalError);
        }

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = bphi[patchi];

//            if (!isA<processorFvsPatchScalarField>(phip))
	    if (!phip.coupled())
            {
                if
                (
                    !Up.fixesValue()
                 || isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                    }
                }
            }
        }

        return mag(massIn)/totalFlux < small
            && mag(fixedMassOut)/totalFlux < small
            && mag(adjustableMassOut)/totalFlux < small;
    }
    else
    {
        return false;
    }
}


bool Foam::adjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    volScalarField& p,
    volScalarField dmdt
)
{
    if (p.needReference())
    {
        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;
        scalar dmdtGlobal = 0.0;

        surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();

        forAll(dmdt.mesh().cells(), celli)
        {
            dmdtGlobal += dmdt[celli];
        }

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

//            if (!isA<processorFvsPatchScalarField>(phip))
	    if (!phip.coupled())
            {
                if
                (
                    Up.fixesValue()
                 && !isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }
                }
                else
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            adjustableMassOut += phip[i];
                        }
                    }
                }
            }
        }

        // Calculate the total flux in the domain, used for normalisation
        scalar totalFlux = vSmall + sum(mag(phi)).value();

        reduce(massIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());
        reduce(dmdtGlobal, sumOp<scalar>());

        scalar massCorr = 1.0;
        scalar magAdjustableMassOut = mag(adjustableMassOut);
 
        if
        (
            magAdjustableMassOut > vSmall
         && magAdjustableMassOut/totalFlux > small
        )
        {
            massCorr = (massIn - dmdtGlobal - fixedMassOut)/adjustableMassOut;
        }
        else if (mag(fixedMassOut - massIn)/totalFlux > 1e-8)
        {
            FatalErrorInFunction
                << "Continuity error cannot be removed by adjusting the"
                   " outflow.\nPlease check the velocity boundary conditions"
                   " and/or run potentialFoam to initialise the outflow." << nl
                << "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << massIn << nl
                << "Specified mass outflow  : " << fixedMassOut << nl
                << "Adjustable mass outflow : " << adjustableMassOut << nl
                << exit(FatalError);
        }

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = bphi[patchi];

            if (!isA<processorFvsPatchScalarField>(phip))
            {
                if
                (
                    !Up.fixesValue()
                 || isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                    }
                }
            }
        }

        return mag(massIn)/totalFlux < small
            && mag(fixedMassOut)/totalFlux < small
            && mag(adjustableMassOut)/totalFlux < small;
    }
    else
    {
        return false;
    }
}


bool Foam::adjustPhi2
(
    surfaceScalarField& phi,
    const volVectorField& U,
    volScalarField& p,
    volScalarField dmdt
)
{
    if (p.needReference())
    {
        scalar fixedMassIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassIn = 0.0;
        scalar adjustableMassOut = 0.0;
        scalar dmdtGlobal = 0.0;

        surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();

        forAll(dmdt.mesh().cells(), celli)
        {
            dmdtGlobal += dmdt[celli];
        }

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

//            if (!isA<processorFvsPatchScalarField>(phip))
	    if (!phip.coupled())
            {
                if
                (
                    Up.fixesValue()
                 && !isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            fixedMassIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }
                }
                else
                {
                    forAll(phip, i)
                    {   
                        if (phip[i] < 0.0)
                        {
                            adjustableMassIn -= phip[i];
                        }
                        else
                        {
                            adjustableMassOut += phip[i];
                        }
                    }
                }
            }
        }
        
        // Calculate the total flux in the domain, used for normalisation
        scalar totalFlux = vSmall + sum(mag(phi)).value();

        reduce(fixedMassIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassIn, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());
        reduce(dmdtGlobal, sumOp<scalar>());

        scalar massCorr1 = 1.0;
        scalar massCorr2 = 1.0;
        scalar magAdjustableMassOut = mag(adjustableMassOut);

        if
        (
            magAdjustableMassOut > vSmall
         && magAdjustableMassOut/totalFlux > small
        )
        {
            scalar b = fixedMassIn - fixedMassOut - dmdtGlobal;
            massCorr1 = mag(-b + sqrt(b*b - 4.0*adjustableMassOut*adjustableMassIn))/2.0/adjustableMassOut;
            massCorr2 = mag(-b - sqrt(b*b - 4.0*adjustableMassOut*adjustableMassIn))/2.0/adjustableMassOut;
        }
        else if (mag(fixedMassOut - fixedMassIn)/totalFlux > 1e-8)
        {
            FatalErrorInFunction
                << "Continuity error cannot be removed by adjusting the"
                   " outflow.\nPlease check the velocity boundary conditions"
                   " and/or run potentialFoam to initialise the outflow." << nl
		<< "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << fixedMassIn << nl
                << "Adjustable mass inflow : " << adjustableMassIn << nl
                << "Specified mass outflow  : " << fixedMassOut << nl
                << "Adjustable mass outflow : " << adjustableMassOut << nl
                << exit(FatalError);
        }

        scalar massCorr = max(massCorr1, massCorr2);

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = bphi[patchi];

            if (!isA<processorFvsPatchScalarField>(phip))
            {
                if
                (
                    !Up.fixesValue()
                 || isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    forAll(phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                        else
                        {
                            phip[i] /= massCorr;
                        }
                    }
                }
            }
        }

        return mag(fixedMassIn)/totalFlux < small
            && mag(adjustableMassIn)/totalFlux < small
            && mag(fixedMassOut)/totalFlux < small
            && mag(adjustableMassOut)/totalFlux < small;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
