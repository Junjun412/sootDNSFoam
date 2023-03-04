// Time staggering as follows
//time labeling:                    ...| n+1 |n+1/2|                 n                          |                n-1/2                     |n-1|...
//time differentiation for scalar:        ...|  Y  |  dY/dt= (Y-Y.oldTime())/runTime().deltaT() |              Y.oldTime()                 |...
//time differentiation for vector:              ...|                 U                          | dU/dt=(U-U.oldTime())/runTime().deltaT() | U.oldTime() | ...

//++++++++++=====reaction source=====++++++++++
if(solveYEqn && CsrcOnOff>0)
{
	currentTimerType=matrix; //dump everythig since last mark into communication&matrix operation time
	#include "cpuTiming.H"
	#include "chemIntegration.H"
	currentTimerType=reacting;
	#include "cpuTiming.H"
}

//++++++++++=====correction vleocity=====++++++++++
U_correct=0*U;
if(preffDiff)
{
	for (int j=0; j<nsp; j++)
	{	
		U_correct+=DIFF[j]*fvc::grad(Y[j])/2.0; //DIFF should be at mid_step
		U_correct+=DIFF[j]*fvc::grad((Y[j]).oldTime())/2.0;
	}
	U_correct.correctBoundaryConditions();
}
phi_correct = fvc::interpolate(rhoTimeInterp*U_correct) & mesh.Sf();

//++++++++++=====species transport equation=====++++++++++
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div((phi|2),Yi)")
    )
);
tmp<fv::convectionScheme<scalar> > mvConvection_correct
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi_correct,
        mesh.divScheme("div((phi_correct|2),Yi)")
    )
);

if(solveYEqn)
{

volScalarField Ysum=0.0*Y[0];
for (int i=0; i<nsp; i++)
{
	volScalarField& Yi=Y[i];
	if(Yi.name()!=("Y_"+inertSpecies))
	{
		if(preffDiff)
		{
			fvScalarMatrix YiEqn
			(
        		fvm::ddt(rho, Yi)
                        + mvConvection->fvmDiv(phi/2.0, Yi)
                        + mvConvection->fvcDiv(phi/2.0, Yi.oldTime())
			+ mvConvection_correct->fvmDiv(phi_correct/2.0, Yi)
			+ mvConvection_correct->fvcDiv(phi_correct/2.0, Yi.oldTime())
        		- fvm::laplacian(rhoTimeInterp*DIFF[i]/2.0, Yi) 
        		- fvc::laplacian(rhoTimeInterp*DIFF[i]/2.0, Yi.oldTime())
        		- fvm::laplacian(turbulence->mut()/2.0/0.72, Yi)
        		- fvc::laplacian(turbulence->mut()/2.0/0.72, Yi.oldTime())
			-  CsrcOnOff*SRC[i]
			-  src_Y[i]
//                        -  PAHsrcOnOff*src_Y[i],
        		);
			YiEqn.solve("Yi");
		}
		else
		{
			fvScalarMatrix YiEqn
       	 		(
        		fvm::ddt(rho, Yi)
                        + mvConvection->fvmDiv(phi/2.0, Yi)
                        + mvConvection->fvcDiv(phi/2.0, Yi.oldTime())
        		- fvm::laplacian(lambda/CP/2.0, Yi) 
        		- fvc::laplacian(lambda/CP/2.0, Yi.oldTime())
        		- fvm::laplacian(turbulence->mut()/2.0/0.72, Yi)
        		- fvc::laplacian(turbulence->mut()/2.0/0.72, Yi.oldTime())
			- CsrcOnOff*SRC[i]
			- src_Y[i]
//                        -  PAHsrcOnOff*src_Y[i],
        		);
			YiEqn.solve("Yi");
		}
		Yi.correctBoundaryConditions();
		Ysum+=Yi;
	}
	else
	{
		Info << "Skip transport equation for inert species " << Yi.name() << endl;
	}
}
for (int i=0; i<nsp; i++)
{
	if(Y[i].name()==("Y_"+inertSpecies))
	{
		Y[i]=1.0-Ysum;
		Y[i].correctBoundaryConditions();
	}
}

currentTimerType=matrix;
#include "cpuTiming.H"

}
