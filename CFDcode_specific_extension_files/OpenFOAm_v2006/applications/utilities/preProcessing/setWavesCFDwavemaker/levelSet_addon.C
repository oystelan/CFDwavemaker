/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

#include "levelSet.H"
#include "levelSet_addon.H"
#include "cut.H"
#include "polyMeshTetDecomposition.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::levelSetFraction
(
    const volScalarField& levelC,
    const pointScalarField& levelP,
    const bool above
)
{
    const fvMesh& mesh = levelC.mesh();

    tmp<volScalarField> tResult
    (
        volScalarField::New
        (
            "levelSetFraction",
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    volScalarField& result = tResult.ref();

    result.primitiveFieldRef() =
        levelSetFraction
        (
            mesh,
            levelC.primitiveField(),
            levelP.primitiveField(),
            above
        );

    forAll(mesh.boundary(), patchi)
    {
        result.boundaryFieldRef()[patchi] =
            levelSetFraction
            (
                mesh.boundary()[patchi],
                levelC.boundaryField()[patchi],
                levelP.boundaryField()[patchi].patchInternalField()(),
                above
            );
    }

    return tResult;
}


// ************************************************************************* //
