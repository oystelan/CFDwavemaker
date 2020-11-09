/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "CFDwavemakerWaveModel.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "polyPatch.H"
#include "fvPatchFields.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace waveModels
    {
        defineTypeNameAndDebug(CFDwavemakerModel, 0);
        addToRunTimeSelectionTable
        (
            waveModel,
            CFDwavemakerModel,
            patch
        );
    }
}


void Foam::waveModels::CFDwavemakerModel::set_face_center_coordinates()
{
    // Local face centres
    const vectorField& Cf = patch_.faceCentres();
    const vectorField CfLocal(Rgl_ & Cf);
    x_ = CfLocal.component(0);
    y_ = CfLocal.component(1);

}


void Foam::waveModels::CFDwavemakerModel::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    std::cout << "setLevel i CFDwavemodel. Do nothing" << std::endl;
    /*forAll(level, paddlei)
    {
        scalar eta = wave_SurfElev(xPaddle_[paddlei], yPaddle_[paddlei], t);

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }*/
}




void Foam::waveModels::CFDwavemakerModel::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{



    forAll(U_, facei)
    {
        // Calculate surface elevation
        scalar surface_elev = wave_SurfElev(x_[facei], y_[facei], t);

        //std::cout << "Elev: " << surface_elev << " z: " << z_[facei] << " vz: " << wave_VeloZ(x_[facei], y_[facei], z_[facei], t) << std::endl;

        // Set volume fraction
        const scalar zMin0 = zMin_[facei] - zMin0_;
        const scalar zMax0 = zMax_[facei] - zMin0_;

        if (zMax0 < surface_elev)
        {
            alpha_[facei] = 1.0;
        }
        else if (zMin0 > surface_elev)
        {
            alpha_[facei] = 0.0;
        }
        else
        {
            alpha_[facei] = (surface_elev - zMin0) / (zMax0 - zMin0);
        }
       
        // Set velocity components
        const vector Uf = vector(wave_VeloX(x_[facei], y_[facei], z_[facei], t) ,
            wave_VeloY(x_[facei], y_[facei], z_[facei] , t),
            wave_VeloZ(x_[facei], y_[facei], z_[facei] , t));
        //std::cout << fraction << " z: " << z_[facei] << "elev: " << level[paddlei] << std::endl;
        U_[facei] = alpha_[facei]*Uf*tCoeff;



    }
}


void Foam::waveModels::CFDwavemakerModel::setAlpha(const scalarField& level)
{
    
    // override setAlpha
    std::cout << "setalpha i CFDwavemodel. Do nothing" << std::endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::CFDwavemakerModel::CFDwavemakerModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    irregularWaveModel(dict, mesh, patch, false)
{
    
    if (readFields)
    {
        readDict(dict);
    }
    // load waveinput.dat file
    int initcheck = wave_Initialize();
}

// 

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //
Foam::waveModels::CFDwavemakerModel::~CFDwavemakerModel() {
    wave_Cleanup();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::CFDwavemakerModel::readDict
(
    const dictionary& overrideDict
)
{
    set_face_center_coordinates();
    if (irregularWaveModel::readDict(overrideDict))
    {
        

        return true;
    }

    return false;
}


void Foam::waveModels::CFDwavemakerModel::info(Ostream& os) const
{
    irregularWaveModel::info(os);
    /*
    os  << "    Wave periods    : " << irregWavePeriods_.size() << nl
        << "    Wave heights    : " << irregWaveHeights_.size() << nl
        << "    Wave phases     : " << irregWavePhases_.size() << nl
        << "    Wave lengths    : " << irregWaveLengths_.size() << nl
        << "    Wave directions : " << irregWaveDirs_.size() << nl;
    */
}


// ************************************************************************* //
