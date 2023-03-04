/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    flameletPimpleFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "radiationModel.H"
#include "pimpleControl.H"
#include "turbulenceModel.H"
//#include "LESfilter.H"

//using namespace std;
#include "sootExtern.H"
#include "OFstream.H"
#include "cantera/IdealGasMix.h"
using namespace Cantera;
#include "chemFuncs.H" //functions used for detailed chemistry

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    pimpleControl pimple(mesh);
    
    Switch solveYEqn(runTime.controlDict().lookup("solveYEqn"));
    Switch solvehEqn(runTime.controlDict().lookup("solvehEqn"));
    Switch solveMomEqn(runTime.controlDict().lookup("solveMomEqn"));
    Switch includeSoot(runTime.controlDict().lookup("includeSoot"));
    scalar Ulim(readScalar(runTime.controlDict().lookup("Ulim")));
    scalar CsrcOnOff(readScalar(runTime.controlDict().lookup("CsrcOnOff"))); // dictionary for on/off of source terms
    scalar YsrcOnOff(readScalar(runTime.controlDict().lookup("YsrcOnOff"))); // dictionary for on/off of source terms from soot
    scalar MomSrcOnOff(readScalar(runTime.controlDict().lookup("MomSrcOnOff"))); // dictionary for on/off of source terms
    Switch predictGradp(runTime.controlDict().lookup("predictGradp"));

    #include "createCpuTiming.H" //monitor parallel computing time
    #include "createDetailedChemistry.H"
    #include "createFields.H"
//    #include "createRadiationModel.H"
    #include "initializeFlowField.H"
//    #include "chemThermoUpdate.H"
    //Guess old density values at first
    //rho n-1/2 = rho n+1/2
    rho.oldTime() = rho;
    rho.oldTime().correctBoundaryConditions();
    //rho n-3/2 = rho n+1/2
    rho.oldTime().oldTime() = rho;
    rho.oldTime().oldTime().correctBoundaryConditions();
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    int subiter = 0;
    int iter = 0;
    while (runTime.run())
    {
        bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
        scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
        scalar maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
        #include "readTimeControls.H"
//        #include "compressibleCourantNo.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
	Switch solveYEqn;
        Switch solvehEqn;
        Switch solveMomEqn;
	Switch includeSoot;
        solveYEqn = readBool(runTime.controlDict().lookup("solveYEqn"));
        solvehEqn = readBool(runTime.controlDict().lookup("solvehEqn"));
        solveMomEqn = readBool(runTime.controlDict().lookup("solveMomEqn"));
        includeSoot = readBool(runTime.controlDict().lookup("includeSoot"));
        scalar CsrcOnOff;
        scalar YsrcOnOff;
	scalar PAHsrcOnOff;
	scalar rhoSrcOnOff;
	scalar MomSrcOnOff;
        CsrcOnOff = readScalar(runTime.controlDict().lookup("CsrcOnOff"));
        PAHsrcOnOff = readScalar(runTime.controlDict().lookup("PAHsrcOnOff"));
	YsrcOnOff = readScalar(runTime.controlDict().lookup("YsrcOnOff"));
	rhoSrcOnOff = readScalar(runTime.controlDict().lookup("rhoSrcOnOff"));
        MomSrcOnOff = readScalar(runTime.controlDict().lookup("MomSrcOnOff"));

	iter++;
        if(iter==1)
        {
          //Solve a pressure equation to guarantuee that at first, mass is exactly conserved
          //#include "pEqn_init.H"
          phi.oldTime()=phi;
          phi.oldTime().oldTime()=phi;
          if(initialize) 
    	  {
		#include "initializeFlowField.H"
    	  }
	  if(equilibrate)
	  {
		#include "equilibrateFlowField.H"
          }
	}
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;
        
//        Info << "BeforeSoot = " << timer.cpuTimeIncrement() << nl << endl;
	#include "sootSrc.H"
//        Info << "AfterSoot = " << timer.cpuTimeIncrement() << nl << endl;

	int subiter = 0;
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            subiter ++;
            if (subiter ==1)
            {
              Info << "Guessing new value rho " << endl;
              //Guessing new value of rho and phi using linear extrapolation
              rho = 2.0*rho.oldTime() - rho.oldTime().oldTime();
              phi = 2.0*phi.oldTime() - phi.oldTime().oldTime();
              rho.correctBoundaryConditions();
            }
 
            turbulence->correct();        
            // Solver scalar equation
            // Equation ------> n+1
            // Y_i -----------> n+3/2
            // Y_i old -------> n+1/2 
            #include "YEqn.H"

            // Solver scalar equation
            // Equation ------> n+1
            // h -------------> n+3/2
            // h old ---------> n+1/2 
            #include "hEqn.H"

	    #include "MomEqn.H"

            // Relation between scalars and looked up variable
            // set state variable -----------> n+3/2
            // set state variable old -------> n+1/2
            #include "chemThermoUpdate.H"

            // Momentum Equation
            // Equation ------> n+1/2
            // U -------------> n+1
            // U old ---------> n 
            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }

        //Info << "Cumulative error = " << cumulativeError << endl << endl;
        //dmdtError = fvc::ddt(rho) + fvc::div(phi);

	if(iter==10) // CpuTiming for the first 10 steps
	{
		#include "writeCpuTiming.H"
	}

        if (runTime.deltaT().value() < timeFlag)
            runTime.writeAndEnd();

        runTime.write();
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
