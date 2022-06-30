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

Description:
    Functions to create an averaged field from a discontinuous field defined by
    a level-set.

SourceFiles:
    levelSet.C
    levelSetTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef levelSet_addons_H
#define levelSet_addons_H

#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Calculate the average value of two fields, one on each side of a level set
//  The level set and the fields are given both on the points and cell-centres.
    template<class Type>
    tmp<Field<Type>> levelSetAverage
    (
        const fvMesh& mesh,
        const scalarField& levelC,
        const scalarField& levelP,
        const Field<Type>& positiveC,
        const Field<Type>& positiveP,
        const Field<Type>& negativeC,
        const Field<Type>& negativeP
    );

    //- As the above overload, but both in cells and on patch faces
    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh>> levelSetAverage
    (
        const volScalarField& levelC,
        const pointScalarField& levelP,
        const GeometricField<Type, fvPatchField, volMesh>& positiveC,
        const GeometricField<Type, pointPatchField, pointMesh>& positiveP,
        const GeometricField<Type, fvPatchField, volMesh>& negativeC,
        const GeometricField<Type, pointPatchField, pointMesh>& negativeP
    );


    //- As the above overload, but both in cells and on patch faces
    tmp<volScalarField> levelSetFraction
    (
        const volScalarField& levelC,
        const pointScalarField& levelP,
        const bool above
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "levelSetTemplates_addons.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
