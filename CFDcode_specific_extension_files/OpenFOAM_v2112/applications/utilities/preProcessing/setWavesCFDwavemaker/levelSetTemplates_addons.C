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

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::levelSetAverage
(
    const fvMesh& mesh,
    const scalarField& levelC,
    const scalarField& levelP,
    const Field<Type>& positiveC,
    const Field<Type>& positiveP,
    const Field<Type>& negativeC,
    const Field<Type>& negativeP
)
{
    tmp<Field<Type>> tResult(new Field<Type>(mesh.nCells(), Zero));
    Field<Type>& result = tResult.ref();

    forAll(result, cI)
    {
        const List<tetIndices> cellTetIs =
            polyMeshTetDecomposition::cellTetIndices(mesh, cI);

        scalar v = 0;
        Type r = Zero;

        forAll(cellTetIs, cellTetI)
        {
            const triFace triIs = cellTetIs[cellTetI].faceTriIs(mesh);

            const FixedList<point, 4>
                tet =
            {
                mesh.cellCentres()[cI],
                mesh.points()[triIs[0]],
                mesh.points()[triIs[1]],
                mesh.points()[triIs[2]]
            };
            const FixedList<scalar, 4>
                level =
            {
                levelC[cI],
                levelP[triIs[0]],
                levelP[triIs[1]],
                levelP[triIs[2]]
            };
            const cut::volumeIntegrateOp<Type>
                positive = FixedList<Type, 4>
                ({
                    positiveC[cI],
                    positiveP[triIs[0]],
                    positiveP[triIs[1]],
                    positiveP[triIs[2]]
                    });
            const cut::volumeIntegrateOp<Type>
                negative = FixedList<Type, 4>
                ({
                    negativeC[cI],
                    negativeP[triIs[0]],
                    negativeP[triIs[1]],
                    negativeP[triIs[2]]
                    });

            v += cut::volumeOp()(tet);

            r += tetCut(tet, level, positive, negative);
        }

        result[cI] = r / v;
    }

    return tResult;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::levelSetAverage
(
    const volScalarField& levelC,
    const pointScalarField& levelP,
    const GeometricField<Type, fvPatchField, volMesh>& positiveC,
    const GeometricField<Type, pointPatchField, pointMesh>& positiveP,
    const GeometricField<Type, fvPatchField, volMesh>& negativeC,
    const GeometricField<Type, pointPatchField, pointMesh>& negativeP
)
{
    const fvMesh& mesh = levelC.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tResult
    (
        GeometricField<Type, fvPatchField, volMesh>::New
        (
            positiveC.name() + ":levelSetAverage",
            mesh,
            dimensioned<Type>("0", positiveC.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& result = tResult.ref();

    result.primitiveFieldRef() =
        levelSetAverage
        (
            mesh,
            levelC.primitiveField(),
            levelP.primitiveField(),
            positiveC.primitiveField(),
            positiveP.primitiveField(),
            negativeC.primitiveField(),
            negativeP.primitiveField()
        );

    forAll(mesh.boundary(), patchi)
    {
        result.boundaryFieldRef()[patchi] =
            levelSetAverage
            (
                mesh.boundary()[patchi],
                levelC.boundaryField()[patchi],
                levelP.boundaryField()[patchi].patchInternalField()(),
                positiveC.boundaryField()[patchi],
                positiveP.boundaryField()[patchi].patchInternalField()(),
                negativeC.boundaryField()[patchi],
                negativeP.boundaryField()[patchi].patchInternalField()()
            );
    }

    return tResult;
}




// ************************************************************************* //
