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
#include "fvcCurl.H"
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::hybridshockFluid::fluxPredictor()
{
    if (!pos.valid())
    {
        pos = surfaceScalarField::New
        (
            "pos",
            mesh,
            dimensionedScalar(dimless, 1)
        );

        neg = surfaceScalarField::New
        (
            "neg",
            mesh,
            dimensionedScalar(dimless, -1.0)
        );
    }

    rho_pos = interpolate(rho, pos());
    rho_neg = interpolate(rho, neg());

    const volVectorField rhoU(rho*U);
    rhoU_pos = interpolate(rhoU, pos(), U.name());
    rhoU_neg = interpolate(rhoU, neg(), U.name());

    U_pos = surfaceVectorField::New("U_pos", rhoU_pos()/rho_pos());
    U_neg = surfaceVectorField::New("U_neg", rhoU_neg()/rho_neg());

    const volScalarField& T = thermo.T();

    const volScalarField rPsi("rPsi", 1.0/thermo.psi());
    const surfaceScalarField rPsi_pos(interpolate(rPsi, pos(), T.name()));
    const surfaceScalarField rPsi_neg(interpolate(rPsi, neg(), T.name()));

    p_pos = surfaceScalarField::New("p_pos", rho_pos()*rPsi_pos);
    p_neg = surfaceScalarField::New("p_neg", rho_neg()*rPsi_neg);

  
    phiv_pos = surfaceScalarField::New("phiv_pos", U_pos() & mesh.Sf());
    phiv_neg = surfaceScalarField::New("phiv_neg", U_neg() & mesh.Sf());

    // Make fluxes relative to mesh-motion
    if (mesh.moving())
    {
        phiv_pos -= mesh.phi();
        phiv_neg -= mesh.phi();
    }

    const volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
    const surfaceScalarField cSf_pos
    (
        "cSf_pos",
        interpolate(c, pos(), T.name())*mesh.magSf()
    );
    const surfaceScalarField cSf_neg
    (
        "cSf_neg",
        interpolate(c, neg(), T.name())*mesh.magSf()
    );

    const dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0);

    const surfaceScalarField ap
    (
        "ap",
        max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
    );
    const surfaceScalarField am
    (
        "am",
        min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
    );

    a_pos = surfaceScalarField::New
    (
        "a_pos",
        fluxScheme == "Tadmor"
          ? surfaceScalarField::New("a_pos", mesh, 0.5)
          : ap/(ap - am)
    );

    a_neg = surfaceScalarField::New("a_neg", 1.0 - a_pos());

    phiv_pos *= a_pos();
    phiv_neg *= a_neg();

    aSf = surfaceScalarField::New
    (
        "aSf",
        fluxScheme == "Tadmor"
          ? -0.5*max(mag(am), mag(ap))
          : am*a_pos()
    );

/***************************************Shock Sensor*************************************/

/* Gradients and other intermediate fields */
    
   
    volScalarField divU = fvc::div(U);       // Divergence of U
    volVectorField curlU = fvc::curl(U);     // Vorticity of V
    
    const scalar epsilon1 = 1e-06;
    
 // Interpolate divU to faces
    faceDivPos = interpolate(divU, pos());
    faceDivNeg = interpolate(divU, neg());    
    
// Interpolate curlU to face 
    faceCurlPos = interpolate(curlU, pos());
    faceCurlNeg = interpolate(curlU, neg());
    
// Calculating the shock sensor
    
   DurcosPos = ((faceDivPos*faceDivPos)/((faceDivPos*faceDivPos)-(mag(faceCurlPos))+epsilon1));  
   DurcosNeg = ((faceDivNeg*faceDivNeg)/((faceDivNeg*faceDivNeg)-(mag(faceCurlNeg))+epsilon1)); 

    surfaceScalarField CbPos
    (
    "CbPos", 0.5*(1-tanh(2.5 + 10*(mesh.deltaCoeffs()/cSf_pos)))*DurcosPos
    ); 
    
   surfaceScalarField  CbNeg
    (
    "CbNeg", 0.5*(1-tanh(2.5 + 10*(mesh.deltaCoeffs()/cSf_neg)))*DurcosNeg
    ); 
    
   /***************************************Shock Sensor*************************************/ 
    
    aphiv_pos = surfaceScalarField::New("aphiv_pos", phiv_pos - aSf());
    aphiv_neg = surfaceScalarField::New("aphiv_neg", phiv_neg + aSf());

    phi_ = CbPos*aphiv_pos()*rho_pos() + (1-CbPos)*phiv_pos + CbNeg*aphiv_neg()*rho_neg() + (1- CbNeg)*phiv_neg; 

}

// ************************************************************************* //
