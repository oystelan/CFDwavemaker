/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "CFDwavemaker_OF1912.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// standard callable functions for static linking
extern "C" int wave_Initialize(); //function needs to be called at startup (reading wave input file)
extern "C" int wave_Cleanup(); // Needs to be called before program end to clear variables stored in memory
extern "C" double wave_VeloX(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
extern "C" double wave_VeloY(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
extern "C" double wave_VeloZ(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
extern "C" double wave_DynPres(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
extern "C" double wave_SurfElev(double, double, double); // input variables are {xpoint,ypoint,time}
extern "C" double wave_VFrac(double, double, double, double, double);// input variables are {xpoint,ypoint,zpoint,time, delta_cellsize}
    
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(landeCFDwavemaker, 0);

addToRunTimeSelectionTable
(
    externalWaveForcing,
    landeCFDwavemaker,
    externalWaveForcing
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


CFDwavemaker::CFDwavemaker
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    externalWaveForcing(io, rT, mesh)
{
    wave_Initialize();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CFDwavemaker::step()
{
    // Nothing to be done
}

// close statement
void CFDwavemaker::close(){
    wave_Cleanup();
}

scalar CFDwavemaker::eta
(
    const point& x,
    const scalar& time
) const
{
    return wave_SurfElev(x[0], x[1], time); // input variables are {xpoint,ypoint,time};
}


scalar CFDwavemaker::ddxPd
(
    const point& x,
    const scalar& time,
    const vector& unitVector
) const
{
    return 0.0;
}


scalar CFDwavemaker::p
(
    const point& x,
    const scalar& time
) const
{
    return 0.0;
}


vector CFDwavemaker::U
(
    const point& x,
    const scalar& time
) const
{
    // Map the solution of OF-format but still in rotated form
    vector U(vector::zero);
    U.x() = wave_VeloX(x[0], x[1], x[2], time); // input variables are {xpoint,ypoint,zpoint,time}
    U.y() = wave_VeloY(x[0], x[1], x[2], time); // input variables are {xpoint,ypoint,zpoint,time}
    U.z() = wave_VeloZ(x[0], x[1], x[2], time); // input variables are {xpoint,ypoint,zpoint,time}
    return U;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
