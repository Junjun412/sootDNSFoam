Info << "Perform chemistry ODE integration" << endl;
if (-1>0)
{
forAll(t_chem.internalField(),cellI)
{
        //-----First do time scales estimation-----
        //setState
        doublereal Y_tmp[nsp];
        for(int j=0;j<nsp;j++)
        {
                Y_tmp[j]=(Y[j].internalField()[cellI]+Y[j].oldTime().internalField()[cellI])/2.0;
        }
        cgas.setMassFractions(Y_tmp);
        cgas.setState_HP((h.internalField()[cellI]+h.oldTime().internalField()[cellI])/2.0,press);

      if(t_chem.mesh().C()[cellI].x()>0.009 && t_chem.mesh().C()[cellI].x()<0.0092 && t_chem.mesh().C()[cellI].y()<0.0001 && t_chem.mesh().C()[cellI].y()>0.0)
      {
              Info << t_chem.mesh().C()[cellI].x() << t_chem.mesh().C()[cellI].y() << endl;
              for(int i=0; i<nsp; i++) Info << cgas.speciesName(i) << " " << cgas.moleFraction(i) << endl;
              Info << "Tem " << cgas.temperature() << "  P " << cgas.pressure() << endl;
      }
}
}
forAll(t_chem.internalField(),cellI)
{
    if(t_chem.mesh().C()[cellI].x()>-0.0015)
    {
	//-----First do time scales estimation-----
	//setState
	doublereal Y_tmp[nsp];
	for(int j=0;j<nsp;j++)	
	{
		Y_tmp[j]=(Y[j].internalField()[cellI]+Y[j].oldTime().internalField()[cellI])/2.0;
	}
	cgas.setMassFractions(Y_tmp);
	cgas.setState_HP((h.internalField()[cellI]+h.oldTime().internalField()[cellI])/2.0,press);

	//estimate temperature source
	cgas.getNetProductionRates(wdot);
	doublereal eth[nsp];
	cgas.getPartialMolarEnthalpies(eth);
	scalar SRC_T=0.0;
	for(int i=0; i<nsp; i++)
	{
		SRC_T+=eth[i]*wdot[i];
	}
	if(SRC_T<0) SRC_T=-SRC_T;
	SRC_T=SRC_T/cgas.density()/cgas.cp_mass();
	//estimate species mass fraction source
	scalar SRC_Y=0.0;
	for(int i=0; i<nsp; i++)
	{
		scalar SRC_Yi=wdot[i]*mwt[i]/cgas.density();
		if(SRC_Yi<0) SRC_Yi=-SRC_Yi;
		SRC_Y= SRC_Yi>SRC_Y ? SRC_Yi:SRC_Y;
	}
	//minimal chemical time scale
	t_chem.ref()[cellI] = T.internalField()[cellI]/(SRC_T+ROOTVSMALL) < 1.0/(SRC_Y+ROOTVSMALL) ? \
	T.internalField()[cellI]/(SRC_T+ROOTVSMALL) : 1.0/(SRC_Y+ROOTVSMALL);
	t_chem.ref()[cellI] = t_chem.internalField()[cellI] < 1000 ? t_chem.internalField()[cellI] : 1000; //clip very large time scale, just for visualization
	
	//-----Then perform integration-----
	//if(runTime.deltaTValue()<1e-4*t_chem.internalField()[cellI]) //change this conditional statement to always false if you want to use stiff chemistry solver for the entire domain
	if(-1>0)
	{
		//Directly use forward time stepping when chemical time scale is large
		for(int i=0;i<nsp;i++)
		{
			SRC[i].ref()[cellI]=wdot[i]*mwt[i];
		}
		stiff_ode.ref()[cellI]=0.0;
	}
	else
	{
		//Otherwise, use CVODE to solve stiff chemistry
                doublereal Y_tmp[nsp];
                for(int i=0;i<nsp;i++)
                {
                        NV_Ith_S(y,i)=Y[i].oldTime().internalField()[cellI];
                        Y_tmp[i]=Y[i].oldTime().internalField()[cellI];
                }
                cgas.setMassFractions(Y_tmp);
                cgas.setState_HP(h.oldTime().internalField()[cellI],press);
		ode_time = T0;
		CVodeReInit(cvode_mem, ode_time, y);
		flag = CVode(cvode_mem, runTime.deltaTValue(), y, &ode_time, CV_NORMAL);
		doublereal* sun_data=N_VGetArrayPointer(y);
		for(int i=0;i<nsp;i++)
		{
			SRC[i].ref()[cellI]=(*(sun_data+i)-Y[i].oldTime().internalField()[cellI])/runTime.deltaTValue()*rhoTimeInterp.internalField()[cellI];
		}
		stiff_ode.ref()[cellI]=1.0;
	}
    }
}
