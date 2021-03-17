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

\*---------------------------------------------------------------------------*/

#include "tmp.H"
#include "fvm.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "meshTools.H"
#include "Function1.H"
#include "uniformDimensionedFields.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "WaveDampToIncident.H"
#include "waveVelocityFvPatchVectorField.H"
//#include "waveSuperposition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(WaveDampToIncident, 0);
    addToRunTimeSelectionTable(option, WaveDampToIncident, dictionary); 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


fv::WaveDampToIncident::WaveDampToIncident
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    ramp_(),
    origins_(),
    directions_(),
    // initialize empty forcing zone, zero everywhere
    forcingZone_
    (
        IOobject
        (
            name + "_waveForcingZone",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("forcingZone", dimless, 0.0)
    ),
    // initialize to no wave velocity, zero everywhere
    waveVelocity_
    (
        IOobject
        (
            name + "_waveForcingVelocity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("forcingVelocity", dimVelocity, vector(0, 0, 0))
    )
{
    // Read input from fvOptions
    read(dict);

    // Compute the forcing zone position before the time loop
    updateForcingZone();

    // Write the initial forcing zone
    forcingZone_.write();

    // load waveinput.dat file
    //int initcheck = wave_Initialize();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Calculate the ramp field (zero in centre of domain, 1 at maximum damping)
void fv::WaveDampToIncident::updateForcingZone() {
    scalarField& zone = forcingZone_.primitiveFieldRef();
    const vectorField& c = mesh_.cellCentres();

    zone = 0;
    if (!ramp_.valid()) {
        Info << "WaveDampToIncident:" << name_ << " the ramp function is not valid!" << endl;
        return;
    }

    // Iterate over all origins and directions given for the forcing zones
    // This allocates memory, should be optimized if run every time step
    forAll(origins_, i)
    {
        scalarField ramp_coordinates = scalarField((c - origins_[i]) & directions_[i]);
        zone = max(zone, ramp_->value(ramp_coordinates));
    }
    Info << "WaveDampToIncident:" << name_ << " updated forcing zone" << endl;
}

// Return the reference (incident wave) velocity
void fv::WaveDampToIncident::updateWaveVelocity() {
    const objectRegistry& db = mesh_.thisDb();
    const vectorField& centroids = mesh_.cellCentres();
    const scalar t = db.time().timeOutputValue();
    //const waveSuperposition& waves = waveSuperposition::New(db);
    
    // Calculate wave velocity field
    vectorField& Uref = waveVelocity_.primitiveFieldRef();
    //Uref = waves.ULiquid(t, centroids);
    Uref = velocity(t, centroids);
}

//- Get the wave velocity at a given time and local coordinates. Local
//  x is aligned with the direction of propagation, and z with negative
//  gravity.
tmp<vectorField> fv::WaveDampToIncident::velocity
(
    const scalar t,
    const vectorField& xyz
)
{
    //SpectralWaveData* swd = const_cast<SpectralWaveData*>(&swd_());
    //swd->UpdateTime(t);

    const scalarField x(xyz.component(0));
    const scalarField y(xyz.component(1));
    const scalarField z(xyz.component(2));
    scalarField xvel(x.size());
    scalarField yvel(x.size());
    scalarField zvel(x.size());

    forAll(x, i)
    {
        //vector_swd U = swd->GradPhi(x[i], y[i], z[i]);
        xvel[i] = wave_VeloX(x[i], y[i], z[i], t);
        yvel[i] = wave_VeloY(x[i], y[i], z[i], t);
        zvel[i] = wave_VeloZ(x[i], y[i], z[i], t);
        /*double zeta = scalar(swd->Elev(x[i], 0.0));
        if(zeta > z[i])
        {
            xvel[i] = U.x;
            zvel[i] = U.z;
        }
        else
        {
            xvel[i] = 0.0;
            zvel[i] = 0.0;
        }*/
    }

    return zip(xvel, yvel, zvel);
}

// Add source terms to the momentum equation
// This is called in the time loop by the addSup methods that are
// used by the fvOptions() machinery in OpenFOAM for extra source terms
void fv::WaveDampToIncident::addSourceTerm
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (mesh_.changing()) {
        updateForcingZone();
    }
    updateWaveVelocity();
    const volVectorField& U = eqn.psi();
    const volVectorField& Uref = waveVelocity_;
    const dimensionedScalar penalty("penalty", dimless/dimTime, penalty_);

    // Calculate and add forcing term to momentum equations
    fvMatrix<vector> dampingEqn
    (
         penalty * alpha * rho * (fvm::Sp(forcingZone_, U) - forcingZone_ * Uref)
    );
    eqn -= dampingEqn;

    Info << "WaveDampToIncident:" << name_ << " added source term at "
         << mesh_.time().timeName()
         << " for field " << fieldNames_[fieldi]  << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  

void fv::WaveDampToIncident::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    Info << "WaveDampToIncident:" << name_ << ": Calling addSup with two args" << endl;
    Info << "WaveDampToIncident:" << name_ << ": NOT IMPLEMENTED! MISSING RHO" << endl;
}


void fv::WaveDampToIncident::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    Info << "WaveDampToIncident:" << name_ << ": Calling addSup with three args" << endl;
    const volScalarField& alpha(mesh_.lookupObject<volScalarField>(alphaName_));
    addSourceTerm(alpha, rho, eqn, fieldI);
}


void fv::WaveDampToIncident::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    Info << "WaveDampToIncident:" << name_ << ": Calling addSup with four args" << endl;
    addSourceTerm(alpha, rho, eqn, fieldI);
}

// Read the fvoptions dictionary at startup time
bool fv::WaveDampToIncident::read(const dictionary& dict)
{
    if (!option::read(dict)) {
        return false;
    } 

    const bool foundRamp = coeffs_.found("ramp");
    const bool foundOgns = coeffs_.found("origins");
    const bool foundDirs = coeffs_.found("directions");
    const bool foundAll =
        foundRamp
        && (foundOgns && foundDirs);
        
    if (!foundAll)
    {
        FatalErrorInFunction
            << "WaveDampToIncident:" << name_ << " ERROR: "
            << "The ramping specification is incomplete. \"ramp\", "
            << "\"origins\" and \"directions\", must all be "
            << "specified in order to ramp the damping. "
            << exit(FatalError);
    }

    
    ramp_ = Function1<scalar>::New("ramp", coeffs_);      
    coeffs_.lookup("origins") >> origins_;
    coeffs_.lookup("directions") >> directions_;

    // Read the penalty
    penalty_ = readScalar(coeffs_.lookup("penalty")); 

    if (origins_.size() == 0 || directions_.size() == 0 || origins_.size() != directions_.size())
    {
        FatalErrorInFunction
                << "WaveDampToIncident:" << name_ << " ERROR: "
                << "The same (non-zero) number of origins and "
                << "directions must be provided" << exit(FatalError);
    }
    
    forAll(directions_, i)
    {
        directions_[i] /= mag(directions_[i]);
    }

    // The name of the VOF alpha field (0 in air, 1.0 in the water)
    alphaName_ = dict.lookupOrDefault<word>("alphaName", "alpha.water");

    // From "option", document here what they are fore when you find out!
    fieldNames_ = wordList(1, coeffs_.lookupOrDefault<word>("U", "U"));
    applied_.setSize(1, false);

    return true;
}

} // end namespace Foam
// ************************************************************************* //
