const labelUList& owner = mesh.owner();
const labelUList& neighbour = mesh.neighbour();
const volScalarField magU = mag(U);
scalar maxU=gMax(magU);

forAll(owner, facei)
{
    label own = owner[facei];
    label nei = neighbour[facei];

    UBlendingFactor[facei] = 0.5*(magU[own] > Ulim) + 0.5*(magU[nei] > Ulim);
    U_0BlendingFactor[facei] = 0.5*(magU[own] > Ulim) + 0.5*(magU[nei] > Ulim);
}

// Solve the Momentum equation
// Equation written at time n+1/2
// U defined at time n and n+1
//rho n+1 estimated using linear time interpolation

rhoTimeInterp = 0.5*(rho+rho.oldTime());
rhoTimeInterp.boundaryFieldRef() = 0.5*(rho.boundaryField()+rho.oldTime().boundaryField());
//rho n estimated using linear time interpolation
rhoTimeInterp.oldTime() = 0.5*(rho.oldTime()+rho.oldTime().oldTime());
rhoTimeInterp.oldTime().boundaryFieldRef() = 0.5*(rho.oldTime().boundaryField()+rho.oldTime().oldTime().boundaryField());

//For Kinetic energy conservation
sqrt_rhoTimeInterp = sqrt(rhoTimeInterp);
sqrt_rhoTimeInterp.oldTime() = sqrt(rhoTimeInterp.oldTime());
sqrt_rhoTimeInterp_timeInterp = 0.5*(sqrt_rhoTimeInterp+sqrt_rhoTimeInterp.oldTime());

currentTimerType=nonReacting;
#include "cpuTiming.H"

fvVectorMatrix UEqn
(
    fvm::ddt(rhoTimeInterp,U)
  + (sqrt_rhoTimeInterp/sqrt_rhoTimeInterp_timeInterp)*fvm::div((phi+phi.oldTime())/4.0, U)
  + (sqrt_rhoTimeInterp.oldTime()/sqrt_rhoTimeInterp_timeInterp)*fvc::div((phi+phi.oldTime())/4.0, U.oldTime())
  - fvm::laplacian(0.5*turbulence->muEff(),U)
  - fvc::laplacian(0.5*turbulence->muEff(),U.oldTime())
  - fvc::div(turbulence->muEff()*0.5*dev2(Foam::T(fvc::grad(U))))
  - fvc::div(turbulence->muEff()*0.5*dev2(Foam::T(fvc::grad(U.oldTime()))))
  - (rhoTimeInterp-freestreamDensity)*g
//  - rhoTimeInterp*g
);

//Best Estimate of gradP at time n+1/2
if(predictGradp)
{
  solve(UEqn==-fvc::grad(p));
}
else
{
  solve(UEqn);
}

currentTimerType=matrix;
#include "cpuTiming.H"


