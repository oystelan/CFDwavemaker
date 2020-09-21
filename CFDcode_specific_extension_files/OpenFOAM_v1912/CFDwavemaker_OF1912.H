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

Class
    Foam::waveTheories::landeCFDwavemaker

Description


SourceFiles
    landeCFDwavemaker.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.
    Bo Terp Paulsen, Deltares, All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef CFDwavemaker_H
#define CFDwavemaker_H

#include "externalWaveForcing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class landeCFDwavemaker Declaration
\*---------------------------------------------------------------------------*/

class CFDwavemaker
:
    public externalWaveForcing
{

protected:

    // Protected data

    // Protected member functions

public:

    //- Runtime type information
    TypeName("CFDwavemaker");

    // Constructors

        //- from components
        CFDwavemaker
        (
            IOobject,
            Time&,
            const fvMesh&
        );


    // Destructor

        ~CFDwavemaker()
        {}


    // Member Functions
        virtual void step();


        virtual inline void close();


        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const;


        virtual scalar ddxPd
        (
            const point&,
            const scalar&,
            const vector&
        ) const;


        virtual vector U
        (
            const point&,
            const scalar&
        ) const;


        virtual scalar p
        (
            const point&,
            const scalar&
        ) const;


        virtual bool writeData(Foam::Ostream&) const
        {
            return true;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
