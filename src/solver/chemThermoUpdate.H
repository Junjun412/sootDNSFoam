// Time staggering as follows
//time labeling:                    ...| n+1 |n+1/2|                 n                          |                n-1/2                     |n-1|...
//time differentiation for scalar:        ...|  Y  |  dY/dt= (Y-Y.oldTime())/runTime().deltaT() |              Y.oldTime()                 |...
//time differentiation for vector:              ...|                 U                          | dU/dt=(U-U.oldTime())/runTime().deltaT() | U.oldTime() | ...
//Note: oldTime() is "last_step"
Info << "Updating thermochemical fields:" << endl;

//-----Special treatment for RHO : RHO = RHO(SCn+3/2)-----
volScalarField rhoPrior = rho;
updateThermoField(rho, h, press, Y, current_step, cgas);
forAll(rho,j)
{
        double density = rho.internalField()[j];
        rho.ref()[j] = density + DensityScaling.value() * (freestreamDensity.value() - density);
}
if(smoothRho)
{
        volScalarField drho = rho - rhoPrior;
        scalar ave_drho = average(drho).value();
        scalar stdDev_drho = Foam::sqrt(average(drho*drho).value() - ave_drho*ave_drho);
        forAll(Zmix,j)
  	{
      		if (drho[j] > ave_drho + 10.0*stdDev_drho) rho[j] = rhoPrior[j] + ave_drho + 10.0*stdDev_drho;
      		if (drho[j] < ave_drho - 10.0*stdDev_drho) rho[j] = rhoPrior[j] + ave_drho - 10.0*stdDev_drho;
        }
}

//-----Treatment for T : not explicitly used in any transport equation, use whatever timestep you want-----
updateThermoField(T, h, press, Y, current_step, cgas);

//-----Special treatment for Scalar diffusion terms : {lambda,CP,DIFF,hi} = {lambda,CP,DIFF,hi}(1/2(SCn+3/2 + SCn+1/2))-----
updateThermoField(lambda, h, press, Y, middle_step, cgas);
updateThermoField(CP, h, press, Y, middle_step, cgas);
if(preffDiff)
{
	updateThermoField(DIFF, h, press, Y, middle_step, cgas);
	updateThermoField(hi, h, press, Y, middle_step, cgas);
}

//-----Special treatment for Momentum diffusion term : mu = mu(SCn+1/2)-----
updateThermoField(mu, h, press, Y, last_step, cgas);
