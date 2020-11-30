/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::fv::WaveDampToIncident

Description
    
    Add a forcing term (source/reaction) to the momentum equation through
    the fvOptions plugin framework. The user specifies start position and
    direction of growth for the damping zone in addition to the distance
    to maximum damping and the maximum daping (forcing) coefficient,
    lambda. The momentum equation then looks like

        LHS(U) == RHS - lambda * ramp * rho * alpha * (U - U_incident_wave)

    Notice that the added term is partially implicit (the unknown U is
    implicit, the incident wave velocity is obviously not).

    The implementation is inspired by acousticDampingSource and verticalDamping

Usage
    
    See the accompanying README.txt for instructions

SourceFiles
    WaveDampToIncident.C

\*---------------------------------------------------------------------------*/

#ifndef WaveDampToIncident_H
#define WaveDampToIncident_H

#include "tmp.H"
#include "cellSetOption.H"
#include "CFDwavemaker.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    template <class Type>
    class Function1;

    namespace fv
    {

        /*---------------------------------------------------------------------------*\
                   Class WaveDampToIncident Declaration
\*---------------------------------------------------------------------------*/

        class WaveDampToIncident
            : public option
        {
        private:
            // Private data

            //- Forcing coefficient [1/s]
            scalar penalty_;

            //- The ramping function
            autoPtr<Function1<scalar>> ramp_;

            //- Origins of the ramping coordinate
            vectorField origins_;

            //- Directions of increasing ramping coordinate
            vectorField directions_;

            //- The zone where damping/momentum forcing is applied []
            //  The field is zero for no damping, 1 for max damping
            volScalarField forcingZone_;

            //- The velocity in the wave field
            volVectorField waveVelocity_;

            //- Name of alpha (default alpha.water)
            word alphaName_;

            // Private Member Functions

            //- The forcing zone ramp field (0 in the middle, 1 at max damping/boundaries)
            void updateForcingZone();

            //- The reference (incident wave) velocity
            void updateWaveVelocity();

            //- Get the wave velocity at a given time and local coordinates. Local
            //  x is aligned with the direction of propagation, and z with negative
            //  gravity.
            virtual tmp<vectorField> velocity(
                const scalar t,
                const vectorField& xyz);

            // - Add the source term to the equation
            void addSourceTerm(
                const volScalarField &alpha,
                const volScalarField &rho,
                fvMatrix<vector> &eqn,
                const label fieldi);

        public:
            //- Runtime type information
            TypeName("WaveDampToIncident");

            // Constructors

            //- Construct from components
            WaveDampToIncident(
                const word &name,
                const word &modelType,
                const dictionary &dict,
                const fvMesh &mesh);

            //- Destructor
            virtual ~WaveDampToIncident()
            {
                //wave_Cleanup();
            }

            // Member Functions

            //- Add implicit contribution to momentum equation
            virtual void addSup(
                fvMatrix<vector> &eqn,
                const label fieldI);

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup(
                const volScalarField &rho,
                fvMatrix<vector> &eqn,
                const label fieldI);

            //- Add implicit contribution to phase momentum equation
            virtual void addSup(
                const volScalarField &alpha,
                const volScalarField &rho,
                fvMatrix<vector> &eqn,
                const label fieldI);

            //- Read dictionary at startup
            virtual bool read(const dictionary &dict);
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
