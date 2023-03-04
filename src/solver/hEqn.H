if(solvehEqn)
{
//++++++++++=====boudnary conditions (needs to match the Dirich b.c. of temperature)=====++++++++++
forAll(h.boundaryField(),patchI)
{
	forAll(h.boundaryField()[patchI],faceI)
	{
		doublereal Y_tmp[nsp];
		for(int j=0; j<Y.size(); j++)	
		{
			Y_tmp[j]=Y[j].boundaryField()[patchI][faceI];
		}
		cgas.setState_TPY(T.boundaryField()[patchI][faceI],press,Y_tmp);
		h.boundaryFieldRef()[patchI][faceI]=cgas.enthalpy_mass();
	}
}
currentTimerType=nonReacting;
#include "cpuTiming.H"

h.correctBoundaryConditions();
//++++++++++=====Preferential diffusion=====++++++++++
laplacePreffDiff*=0.0;
if(preffDiff)
{
	for (int i=0; i<nsp; i++)
	{
		laplacePreffDiff+=fvc::laplacian((lambda/CP-rhoTimeInterp*DIFF[i])*hi[i],Y[i])/2.0; //lambda, CP, DIFF, hi should be at mid_step
		laplacePreffDiff+=fvc::laplacian((lambda/CP-rhoTimeInterp*DIFF[i])*hi[i].oldTime(),Y[i].oldTime())/2.0;
	}
}
laplacePreffDiff.correctBoundaryConditions();
//++++++++++=====enthalpy transport equation=====++++++++++
#include "radiationSourceTerm.H"
if(radiationOnOff)
{

	fvScalarMatrix hEqn
	(	
		fvm::ddt(rho, h)
	        + fvm::div(phi/2.0, h)
	        + fvc::div(phi/2.0, h.oldTime())
	        + fvm::div(phi_correct/2.0, h)
	        + fvc::div(phi_correct/2.0, h.oldTime())
	        - fvm::laplacian(lambda/CP/2.0, h)
	        - fvc::laplacian(lambda/CP/2.0, h.oldTime())
		+ laplacePreffDiff
	        - fvm::laplacian(turbulence->mut()/2.0/0.72, h)
	        - fvc::laplacian(turbulence->mut()/2.0/0.72, h.oldTime())
	//	- radiation->Sh(h, CP)
		- radSRC   // add by JGUO
	);
	hEqn.solve();
	//radSRC = radiation->Sh(h, CP) & dummy;
	//radiation->correct();
}
else
{
        fvScalarMatrix hEqn
        (
                fvm::ddt(rho, h)
                + fvm::div(phi/2.0, h)
                + fvc::div(phi/2.0, h.oldTime())
                + fvm::div(phi_correct/2.0, h)
                + fvc::div(phi_correct/2.0, h.oldTime())
                - fvm::laplacian(lambda/CP/2.0, h)
                - fvc::laplacian(lambda/CP/2.0, h.oldTime())
                + laplacePreffDiff
                - fvm::laplacian(turbulence->mut()/2.0/0.72, h)
                - fvc::laplacian(turbulence->mut()/2.0/0.72, h.oldTime())
        );
        hEqn.solve();	
}

currentTimerType=matrix;
#include "cpuTiming.H"
}
