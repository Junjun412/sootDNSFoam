//rho n+1 estimated using linear time interpolation
rhoTimeInterp = (rho + rho.oldTime())/2.0;
rhoTimeInterp.boundaryFieldRef() = (rho.boundaryField() + rho.oldTime().boundaryField())/2.0;
//rho n estimated using linear time interpolation
rhoTimeInterp.oldTime() = (rho.oldTime() + rho.oldTime().oldTime())/2.0;
rhoTimeInterp.oldTime().boundaryFieldRef() =  (rho.oldTime().boundaryField() + rho.oldTime().oldTime().boundaryField())/2.0;

//rhoDoubleTimeInterp = 0.25*(rho+2.0*rho.oldTime()+rho.oldTime().oldTime());
//rhoDoubleTimeInterp.boundaryField() = 0.25*(rho.boundaryField()+2.0*rho.oldTime().boundaryField()+rho.oldTime().oldTime().boundaryField());


sqrt_rhoTimeInterp = sqrt(rhoTimeInterp);
sqrt_rhoTimeInterp.oldTime() = sqrt(rhoTimeInterp.oldTime());
sqrt_rhoTimeInterp_timeInterp = 0.5*(sqrt_rhoTimeInterp+sqrt_rhoTimeInterp.oldTime());
phi = fvc::interpolate(rhoTimeInterp*U) & mesh.Sf();

//Creation of x Ueqn just to get the array rAU 
fvVectorMatrix UEqn
(
    //fvm::ddt(rho, U)
    fvm::ddt(rhoTimeInterp,U)
  + (sqrt_rhoTimeInterp/sqrt_rhoTimeInterp_timeInterp)*fvm::div((phi+phi.oldTime())/4.0, U)
  + (sqrt_rhoTimeInterp.oldTime()/sqrt_rhoTimeInterp_timeInterp)*fvc::div((phi+phi.oldTime())/4.0, U.oldTime())
  - fvm::laplacian(0.5*turbulence->muEff(),U)
  - fvc::laplacian(0.5*turbulence->muEff()*U.oldTime())
  - fvc::div(turbulence->muEff()*0.5*dev2(Foam::T(fvc::grad(U))))
  - fvc::div(turbulence->muEff()*0.5*dev2(Foam::T(fvc::grad(U.oldTime()))))
);
volScalarField rAU(1.0/UEqn.A());

//Adjust the last face.
dmdt.ref() = (rho.internalField() - rho.oldTime().internalField())/runTime.deltaT().value()*mesh.V();
if(adjustFluxOut) 
{
  adjustPhi(phi, U, p_init, dmdt);
}

//Solve a poisson equation
fvScalarMatrix pEqn_init
(
     fvm::laplacian(rhoTimeInterp*rAU,p_init)
   - fvc::ddt(rho) 
   - fvc::div(phi)
);

  pEqn_init.setReference(pRefCell, pRefValue);
  pEqn_init.solve(mesh.solver(p_init.select(pimple.finalInnerIter())));

  // Correct only the face mass fluxes
  phi -= pEqn_init.flux();


//U -= (rAU)*fvc::grad(p_init);
//U.correctBoundaryConditions();


