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
#include "fvcDdt.H"
#include "multivariateScheme.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::hybridshockFluid::massfractionpredictor()
{ 
    forAll(Y, i)
    {
        if (composition.solve(i))
        {
            volScalarField& Yi = Y_[i];

            surfaceScalarField Ypos = interpolate(Yi, pos()); 
            surfaceScalarField Yneg = interpolate(Yi, neg());
         
           surfaceScalarField phiYi
           (
              "phiYi",
               CbPos()*aphiv_pos()*rho_pos()*Ypos +   CbNeg()*aphiv_neg()*rho_neg()*Yneg + (1-CbPos())*phiv_pos*Ypos + (1-CbNeg())*phiv_neg*Yneg
           );
           	
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi) + fvc::div(phiYi)  == fvModels().source(rho, Yi)
            );
	    if (!inviscid)
	    {
	    YiEqn += thermophysicalTransport->divj(Yi);    
	    }
	    
            YiEqn.relax();
            fvConstraints().constrain(YiEqn);
            YiEqn.solve("Yi");
            fvConstraints().constrain(Yi);
            Yi.correctBoundaryConditions();
        }
    }

    composition.normalise();
}


// ************************************************************************* //
