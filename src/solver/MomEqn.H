if (includeSoot && solveMomEqn)
{
//tmp<fv::convectionScheme<scalar> > mvConvection
//(
//    fv::convectionScheme<scalar>::New
//    (
//        mesh,
//        momFields,
//        phi,
//        mesh.divScheme("div(phi,Mom)")
//    )
//);
	volVectorField gradT = fvc::grad(T);
	forAll(rho,idx)
	{
		double gradmag = mag(gradT[idx]);
		if (gradmag>5000.0)
		{
			gradT[idx] = gradT[idx]*5000.0/gradmag;
		}
	}
	volVectorField U55 = U - 0.5*0.55*mu/T/rho*gradT;

	volVectorField gradTold = fvc::grad(T.oldTime());

	forAll(rho,idx)
	{
		double gradmag = mag(gradTold[idx]);
		if (gradmag>5000.0)
		{
			gradTold[idx] = gradTold[idx]*5000.0/gradmag;
		}
	}
	U55 = U55 - 0.5*0.55*mu/T.oldTime()/rho.oldTime()*gradTold;
	//volVectorField U55 = U - 0.55*mu/T/rho*fvc::grad(T);

	surfaceScalarField phi55 = fvc::interpolate(U55) & mesh.Sf();
	surfaceScalarField rhophi55 = fvc::interpolate(rhoTimeInterp*U55) & mesh.Sf();

	if (thermophoresis)
	{
                for (int i=0; i<nMom; i++)
                {
                        volScalarField& Momi=Mom[i];
                        fvScalarMatrix MomEqn
                        (
        //		    fvm::ddt(Momi)
         //		 + fvm::div(phi55, Momi)
         	//	 - fvm::laplacian((turbulence->mut())/rho, Momi)
         //		 - MomSrc[i] * MomSrcOnOff
                        fvm::ddt(rho, Momi)
                        + fvm::div(rhophi55/2.0, Momi)
                        + fvc::div(rhophi55/2.0, Momi.oldTime())
//                        - fvm::laplacian(sootDIFF*rho/2.0, Momi)
//                        - fvc::laplacian(sootDIFF*rho/2.0, Momi.oldTime())
               //         - fvm::laplacian((turbulence->mut())/2.0, Momi)
               //         - fvc::laplacian((turbulence->mut())/2.0, Momi.oldTime())
                        - MomSrcOnOff * MomSrc[i] //* rhoTimeInterp
                        );
                        MomEqn.solve();
                }
	}
	
	else
	{
		for (int i=0; i<nMom; i++)
		{
			volScalarField& Momi=Mom[i];
			fvScalarMatrix MomEqn
                	(
                	fvm::ddt(rho, Momi)
                	+ fvm::div(phi/2.0, Momi)
                	+ fvc::div(phi/2.0, Momi.oldTime())
               // 	- fvm::laplacian(lambda/CP/2.0, Momi)
              //  	- fvc::laplacian(lambda/CP/2.0, Momi.oldTime())
               // 	- fvm::laplacian(turbulence->mut()/2.0, Momi)
               // 	- fvc::laplacian(turbulence->mut()/2.0, Momi.oldTime())
                	- MomSrcOnOff * MomSrc[i] //* rhoTimeInterp
                	);
			MomEqn.solve();
		}
	}
}
