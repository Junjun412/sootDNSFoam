//volScalarField psi = rho/p;
//volScalarField psiold = rho.oldTime()/p.oldTime();
volScalarField rAU(1.0/UEqn.A());
U = rAU*UEqn.H();

if (poissonPEqn)
{
  //rhoTimeInterp defined at n+1 and n
  rhoTimeInterp = (rho + rho.oldTime())/2.0;
  rhoTimeInterp.oldTime() = (rho.oldTime() + rho.oldTime().oldTime())/2.0;
  rhoTimeInterp.boundaryFieldRef() = (rho.boundaryField() + rho.oldTime().boundaryField())/2.0;
  rhoTimeInterp.oldTime().boundaryFieldRef() =  (rho.oldTime().boundaryField() + rho.oldTime().oldTime().boundaryField())/2.0;

  currentTimerType=nonReacting;
  #include "cpuTiming.H"
  
  phi = fvc::interpolate(rhoTimeInterp*U) & mesh.Sf();
  dmdt.ref() = (rho.internalField() - rho.oldTime().internalField())/runTime.deltaT()*mesh.V();

  // Related to drhodt at n+1
//  dmdt.ref() = (rho.internalField() - rho.oldTime().internalField())/runTime.deltaT()*mesh.V();
  adjustPhi(phi, U, p, dmdt);

  while (pimple.correctNonOrthogonal())
  {
      fvScalarMatrix pEqn
      (
     	 fvm::laplacian(rhoTimeInterp*rAU,p)
    	 - fvc::ddt(rho) 
    	 - fvc::div(phi)
         - src_rho         // source term from soot
       );

    pEqn.setReference(pRefCell, pRefValue);
//    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
      phi -= pEqn.flux();
    }
  }
}

else // use old pEqn
{
    volScalarField psi = rho/p;
    
    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
//          + fvc::ddtPhiCorr(rAU, rho, U, phi)
          + fvc::ddtCorr(rho, U, phi)
        );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*rAU, p)
          - fvm::ddt(psi, p)
          - fvc::div(phi)
        );

//        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
	pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

}

U -= (rAU)*fvc::grad(p);
//#include "ConvBCU.H"

currentTimerType=matrix;
#include "cpuTiming.H"
