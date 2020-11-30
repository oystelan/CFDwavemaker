/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

Application
    setWaves

Description
    Applies wave models to the entire domain for case initialisation using
    level sets for second-order accuracy.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "levelSet.H"
#include "levelSet_addon.H"
#include "pointFields.H"
#include "timeSelector.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "waveAlphaFvPatchScalarField.H"
#include "waveVelocityFvPatchVectorField.H"
#include "CFDwavemaker.h"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, false);

    #include "addDictOption.H"
    #include "addRegionOption.H"

    argList::addOption
    (
        "alpha",
        "name",
        "name of the volume fraction field, default is \"alpha\""
    );

    argList::addOption
    (
        "U",
        "name",
        "name of the velocity field, default is \"U\""
    );

    argList::addBoolOption
    (
        "gas",
        "the volume fraction field is that that of the gas above the wave"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createNamedMesh.H"

    const word dictName("setWavesCFDwavemakerDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictName << "\n" << endl;
    IOdictionary setWavesCFDwavemakerDict(dictIO);

    //#include "readGravitationalAcceleration.H"

    const pointMesh& pMesh = pointMesh::New(mesh);

    // Parse options
    const word alphaName = setWavesCFDwavemakerDict.lookupOrDefault<word>("alpha", "alpha");
    const word UName = setWavesCFDwavemakerDict.lookupOrDefault<word>("U", "U");
    const bool liquid = setWavesCFDwavemakerDict.lookupOrDefault<bool>("liquid", true);

    // Get the wave models
    //const waveSuperposition& waves = waveSuperposition::New(mesh);

    wave_Initialize(); // initialize CFDwavemaker and read input file located in constants folder

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        const scalar t = runTime.value();

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.readUpdate();

        // Read the fields which are to be set
        volScalarField alpha
        (
            IOobject
            (
                alphaName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );
        volVectorField U
        (
            IOobject
            (
                UName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        // Create modelled fields on both cells and points
        volScalarField h
        (
            IOobject("h", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(dimLength, 0)
        );
        pointScalarField hp
        (
            IOobject("hp", runTime.timeName(), mesh),
            pMesh,
            dimensionedScalar(dimLength, 0)
        );
        volVectorField uGas
        (
            IOobject("uGas", runTime.timeName(), mesh),
            mesh,
            dimensionedVector(dimVelocity, vector::zero)
        );
        pointVectorField uGasp
        (
            IOobject("uGasp", runTime.timeName(), mesh),
            pMesh,
            dimensionedVector(dimVelocity, vector::zero)
        );
        volVectorField uLiq
        (
            IOobject("uLiq", runTime.timeName(), mesh),
            mesh,
            dimensionedVector(dimVelocity, vector::zero)
        );
        pointVectorField uLiqp
        (
            IOobject("uLiqp", runTime.timeName(), mesh),
            pMesh,
            dimensionedVector(dimVelocity, vector::zero)
        );

        // Cell centres and points
        const pointField& ccs = mesh.cellCentres();
        const pointField& pts = mesh.points();


        // loop over primitive fields and set height and velocity at all ccs
        scalarField x_ccs(ccs.component(0));
        scalarField y_ccs(ccs.component(1));
        scalarField z_ccs(ccs.component(2));
        scalarField height_above_wave_ccs(x_ccs.size());
        scalarField xvel_liq_ccs(x_ccs.size());
        scalarField yvel_liq_ccs(x_ccs.size());
        scalarField zvel_liq_ccs(x_ccs.size());
        scalarField xvel_gas_ccs(x_ccs.size());
        scalarField yvel_gas_ccs(x_ccs.size());
        scalarField zvel_gas_ccs(x_ccs.size());

        forAll(x_ccs, i)
        {
            height_above_wave_ccs[i] = z_ccs[i] - wave_SurfElev(x_ccs[i], y_ccs[i], t);
            xvel_liq_ccs[i] = wave_VeloX(x_ccs[i], y_ccs[i], z_ccs[i], t);
            yvel_liq_ccs[i] = wave_VeloY(x_ccs[i], y_ccs[i], z_ccs[i], t);
            zvel_liq_ccs[i] = wave_VeloZ(x_ccs[i], y_ccs[i], z_ccs[i], t);
            xvel_gas_ccs[i] = 0.0;
            yvel_gas_ccs[i] = 0.0;
            zvel_gas_ccs[i] = 0.0;
        }

        // set primitive fields
        h.primitiveFieldRef() = height_above_wave_ccs;
        uLiq.primitiveFieldRef() = zip(xvel_liq_ccs, yvel_liq_ccs, zvel_liq_ccs);
        uGas.primitiveFieldRef() = zip(xvel_gas_ccs, yvel_gas_ccs, zvel_gas_ccs);


        // loop over primitive fields and set height and velocity at all pts
        const scalarField x_pts(pts.component(0));
        const scalarField y_pts(pts.component(1));
        const scalarField z_pts(pts.component(2));
        scalarField height_above_wave_pts(x_pts.size());
        scalarField xvel_liq_pts(x_pts.size());
        scalarField yvel_liq_pts(x_pts.size());
        scalarField zvel_liq_pts(x_pts.size());
        scalarField xvel_gas_pts(x_pts.size());
        scalarField yvel_gas_pts(x_pts.size());
        scalarField zvel_gas_pts(x_pts.size());

        forAll(x_pts, i)
        {
            height_above_wave_pts[i] = z_pts[i] - wave_SurfElev(x_pts[i], y_pts[i], t);
            xvel_liq_pts[i] = wave_VeloX(x_pts[i], y_pts[i], z_pts[i], t);
            yvel_liq_pts[i] = wave_VeloY(x_pts[i], y_pts[i], z_pts[i], t);
            zvel_liq_pts[i] = wave_VeloZ(x_pts[i], y_pts[i], z_pts[i], t);
            xvel_gas_pts[i] = 0.0;
            yvel_gas_pts[i] = 0.0;
            zvel_gas_pts[i] = 0.0;
        }

        // set primitive fields
        hp.primitiveFieldRef() = height_above_wave_pts;
        uLiqp.primitiveFieldRef() = zip(xvel_liq_pts, yvel_liq_pts, zvel_liq_pts);
        uGasp.primitiveFieldRef() = zip(xvel_gas_pts, yvel_gas_pts, zvel_gas_pts);



        // Internal field
        //h.primitiveFieldRef() = waves.height(t, ccs);
        //hp.primitiveFieldRef() = waves.height(t, pts);
        //uGas.primitiveFieldRef() = waves.UGas(t, ccs);
        //uGasp.primitiveFieldRef() = waves.UGas(t, pts);
        //uLiq.primitiveFieldRef() = waves.ULiquid(t, ccs);
        //uLiqp.primitiveFieldRef() = waves.ULiquid(t, pts);

        // Boundary fields
        forAll(mesh.boundary(), patchj)
        {
            const pointField& fcs = mesh.boundary()[patchj].Cf();

            
            const scalarField x_patch(fcs.component(0));
            const scalarField y_patch(fcs.component(1));
            const scalarField z_patch(fcs.component(2));
            scalarField height_sf_patch(x_patch.size());
            scalarField xvel_liq_patch(x_patch.size());
            scalarField yvel_liq_patch(x_patch.size());
            scalarField zvel_liq_patch(x_patch.size());
            scalarField xvel_gas_patch(x_patch.size());
            scalarField yvel_gas_patch(x_patch.size());
            scalarField zvel_gas_patch(x_patch.size());

            forAll(x_patch, i)
            {
                height_sf_patch[i] = z_patch[i] - wave_SurfElev(x_patch[i], y_patch[i], t);
                xvel_liq_patch[i] = wave_VeloX(x_patch[i], y_patch[i], z_patch[i], t);
                yvel_liq_patch[i] = wave_VeloY(x_patch[i], y_patch[i], z_patch[i], t);
                zvel_liq_patch[i] = wave_VeloZ(x_patch[i], y_patch[i], z_patch[i], t);
                xvel_gas_patch[i] = 0.0;
                yvel_gas_patch[i] = 0.0;
                zvel_gas_patch[i] = 0.0;
            }


            h.boundaryFieldRef()[patchj] = height_sf_patch;
            uGas.boundaryFieldRef()[patchj] = zip(xvel_gas_patch, yvel_gas_patch, zvel_gas_patch);
            uLiq.boundaryFieldRef()[patchj] = zip(xvel_liq_patch, yvel_liq_patch, zvel_liq_patch);
            //h.boundaryFieldRef()[patchj] = waves.height(t, fcs);
            //uGas.boundaryFieldRef()[patchj] = waves.UGas(t, fcs);
            //uLiq.boundaryFieldRef()[patchj] = waves.ULiquid(t, fcs);
        }





        // Calculate the fields
        volScalarField alphaNoBCs(levelSetFraction(h, hp, !liquid));
        volVectorField UNoBCs(levelSetAverage(h, hp, uGas, uGasp, uLiq, uLiqp));

        // Set the wave and non-wall fixed-value patch fields
        forAll(mesh.boundary(), patchi)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];

            fvPatchScalarField& alphap = alpha.boundaryFieldRef()[patchi];
            fvPatchVectorField& Up = U.boundaryFieldRef()[patchi];
            if
            (
               !isA<wallPolyPatch>(patch)
             || isA<waveAlphaFvPatchScalarField>(alphap)
             || isA<waveVelocityFvPatchVectorField>(Up)
            )
            {
                alphap == alphaNoBCs.boundaryField()[patchi];
                Up == UNoBCs.boundaryField()[patchi];
            }
        }

        // Set the internal fields and all non-fixed value patch fields
        alpha = alphaNoBCs;
        U = UNoBCs;

        // Output
        Info<< "Writing " << alpha.name() << nl;
        alpha.write();
        Info<< "Writing " << U.name() << nl << endl;
        U.write();
    }
    wave_Cleanup();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
