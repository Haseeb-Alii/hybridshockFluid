/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::solvers::hybridshockFluid

Description
    Solver module for density-based solution of compressible flow

    Based on central-upwind schemes of Kurganov and Tadmor with support for
    mesh-motion and topology change.

    Reference:
    \verbatim
        Greenshields, C. J., Weller, H. G., Gasparini, L.,
        & Reese, J. M. (2010).
        Implementation of semi‐discrete, non‐staggered central schemes
        in a colocated, polyhedral, finite volume framework,
        for high‐speed viscous flows.
        International journal for numerical methods in fluids, 63(1), 1-21.
    \endverbatim

SourceFiles
    hybridshockFluid.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::incompressibleFluid

\*---------------------------------------------------------------------------*/

#ifndef hybridshockFluid_H
#define hybridshockFluid_H

#include "fluidSolver.H"
#include "fluidMulticomponentThermo.H"
#include "compressibleMomentumTransportModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "multivariateScheme.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class hybridshockFluid Declaration
\*---------------------------------------------------------------------------*/

class hybridshockFluid
:
    public fluidSolver
{

protected:

    // Thermophysical properties

     // Thermophysical properties
        autoPtr<fluidMulticomponentThermo> thermoPtr_;
        
        fluidMulticomponentThermo& thermo_;

    // Composition

        basicSpecieMixture& composition;
        PtrList<volScalarField>& Y_;  
        
        
    //- Reference to the pressure field
         volScalarField& p_;
        

        //- The continuity density field
        volScalarField rho_;
        
        
         
    // Kinematic properties

        //- Velocity field
        volVectorField U_;
        
        volScalarField delta;
        
        //- Mass-flux field
        surfaceScalarField phi_;
              

        //- Kinetic energy field
        //  Used in the energy equation
        volScalarField K;
           

    // Momentum transport

        bool inviscid;

        //- Pointer to the momentum transport model
        autoPtr<compressible::momentumTransportModel> momentumTransport;


    // Thermophysical transport

        autoPtr<fluidMulticomponentThermophysicalTransportModel>
            thermophysicalTransport;


    // Controls

        word fluxScheme;


    // Cached temporary fields

        tmp<surfaceScalarField> pos;
        tmp<surfaceScalarField> neg;

        tmp<surfaceScalarField> rho_pos;
        tmp<surfaceScalarField> rho_neg;

        tmp<surfaceVectorField> rhoU_pos;
        tmp<surfaceVectorField> rhoU_neg;

        tmp<surfaceVectorField> U_pos;
        tmp<surfaceVectorField> U_neg;
        
        PtrList<surfaceScalarField> rhoYposList;
        PtrList<surfaceScalarField> rhoYnegList;

        tmp<surfaceScalarField> p_pos;
        tmp<surfaceScalarField> p_neg;
        
        tmp<surfaceScalarField> a_pos;
        tmp<surfaceScalarField> a_neg;
              

        tmp<surfaceScalarField> aSf;
        
        tmp<surfaceScalarField> CbPos;
        
        tmp<surfaceScalarField> CbNeg;
        
        tmp<surfaceScalarField> phiL;
         
        
        tmp<surfaceScalarField> aphiv_pos;
        tmp<surfaceScalarField> aphiv_neg;
        

        tmp<surfaceVectorField> devTau;
        
     

        //- Optional LTS reciprocal time-step field
        tmp<volScalarField> trDeltaT;


private:

    // Private Member Functions
	
	
        //- Set rDeltaT for LTS
        void setRDeltaT(const surfaceScalarField& amaxSf);

        //- Correct the cached Courant numbers
        void correctCoNum(const surfaceScalarField& amaxSf);

        //- Construct the continuity equation and correct the density
        void correctDensity();

        //- Interpolate field vf according to direction dir
        template<class Type>
        inline tmp<SurfaceField<Type>> interpolate
        (
            const VolField<Type>& vf,
            const surfaceScalarField& dir,
            const word& reconFieldName = word::null
        )
        {
            tmp<SurfaceField<Type>> tsf
            (
                fvc::interpolate
                (
                    vf,
                    dir,
                    "reconstruct("
                  + (reconFieldName != word::null ? reconFieldName : vf.name())
                  + ')'
                )
            );

            SurfaceField<Type>& sf = tsf.ref();

            sf.rename(vf.name() + '_' + dir.name());

            return tsf;
        }

        void fluxPredictor();

        void clearTemporaryFields();



public:

    // Public Data

        //- Reference to the fluid thermophysical properties
        const fluidMulticomponentThermo& thermo;
        
        //- Reference to the pressure field
        const volScalarField& p;
        
        //- Reference to the continuity density field
        const volScalarField& rho;
       
        //- Reference to the velocity field
        const volVectorField& U;
        
        //- Reference to the mass-flux field
        const surfaceScalarField& phi;
        
        //- Reference to the composition
        const PtrList<volScalarField>& Y;
        
                       

    //- Runtime type information
    TypeName("hybridshockFluid");


    // Constructors

        //- Construct from region mesh
        hybridshockFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        hybridshockFluid(const hybridshockFluid&) = delete;


    //- Destructor
    virtual ~hybridshockFluid();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();
        
        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        // Solving the specie transport equation
        virtual void massfractionpredictor();
        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();
        
        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hybridshockFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
