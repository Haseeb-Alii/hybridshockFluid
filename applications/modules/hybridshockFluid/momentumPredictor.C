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

\*---------------------------------------------------------------------------*/

#include "hybridshockFluid.H"
#include "fvmDdt.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::hybridshockFluid::momentumPredictor()
{
    volVectorField& U(U_);
    
    // Printing the maximum and minimum values during runtime
    Info<< "Max CbPosM: " << max(CbPos()) << ", Min CbPosM: " << min(CbPos()) << endl;
    
    surfaceVectorField phiUk
    (
        CbPos()*aphiv_pos()*rhoU_pos()  + CbNeg()*aphiv_neg()*rhoU_neg() +   (a_pos()*p_pos() + a_neg()*p_neg())*mesh.Sf()
    );
    
      surfaceVectorField phiUl
    (
        (1-CbPos())*aphiv_pos()*rhoU_pos() + (1-CbPos())*aSf()*rhoU_pos() + (1-CbNeg())*aphiv_neg()*rhoU_neg() - (1-CbNeg())*aSf()*rhoU_neg()
    );

    // Construct the divDevTau matrix first
    // so that the maxwellSlipU BC can access the explicit part
    tmp<fvVectorMatrix> divDevTau;
    if (!inviscid)
    {
        divDevTau = momentumTransport->divDevTau(U);
    }

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvc::div(phiUk) + fvc::div(phiUl) 
      ==
        fvModels().source(rho, U)
    );

    if (!inviscid)
    {
        UEqn += divDevTau();
    }

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    solve(UEqn);

    fvConstraints().constrain(U);
    
    K = 0.5*magSqr(U);

    if (!inviscid)
    {
        devTau = divDevTau->flux();
    }
    
     Info<< "Max U: " << max(U) << ", Min U: " << min(U) << endl;
}


// ************************************************************************* //
