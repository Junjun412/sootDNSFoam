if (includeSoot)
{
  Info << "Calculate soot source term" << endl;

  double maxMom = 0.0;
// Initialize arrays
  forAll(rho,idx)
  {
    for(int i=0; i<nMom; i++)
    {
        MomSrc[i][idx] = 0.0;
    }
    for(int i=0; i<nsp; i++)
    {
        src_Y[i][idx] = 0.0;
    }

    src_Tem[idx]   = 0.0;

    sootND[idx]    = 0.0;	//Number density
    sootFV[idx]    = 0.0;	//volume fraction
    sootDiam[idx]  = 0.0;	//Particle Diameter
    sootAggr[idx]  = 0.0;

    sootNDsrcNucl[idx] = 0.0;
    sootNDsrcCoag[idx] = 0.0;
    sootNDsrcCond[idx] = 0.0;
    sootNDsrcSurf[idx] = 0.0;
    sootNDsrcOxid[idx] = 0.0;

    sootFVsrcNucl[idx] = 0.0;
    sootFVsrcCoag[idx] = 0.0;
    sootFVsrcCond[idx] = 0.0;
    sootFVsrcSurf[idx] = 0.0;
    sootFVsrcOxid[idx] = 0.0;
  }

    forAll(rho,idx)
    {
	double timeStep = runTime.deltaT().value(); 
        int    nSc      = 4 + 7 + 1; //8 PAH + OH, H2O, H, H2, C2H2, O2, CO,  Temperature
        double tempRho  = rhoTimeInterp.internalField()[idx];
        int    tempnMom = nMom;		//nMom is defined in createDetailedChem.H
        double MomArr[nMom];
        double MomSrcArr[nSc+nMom];	//Source term
        double tableScalars[nSc];	//Input: species and temperature
        double Y_tmp[nsp];
        double sum;

        sum = 0.0;
        for(int i=0; i<nMom; i++)     MomArr[i]    = Mom[i][idx];
        for(int i=0; i<nSc+nMom; i++) MomSrcArr[i] = 0.0;

	// Clip moments of MomArr, not Mom itself
//	#include "clip_moments_indiv.H"
        #include <cmath>

        if(MomArr[2] > maxMom) maxMom = MomArr[2];

        for(int i=0;i<nsp;i++)
        {
                Y_tmp[i]=Y[i].internalField()[idx];
                sum     = sum + Y_tmp[i] / mwt[i];
        }
        cgas.setMassFractions(Y_tmp);
	double tempWmix = 1.0 / sum / 1000.0;  // Mean molecular weight of the mixture (kg mol-1)

        tableScalars[0] =Y_tmp[cgas.speciesIndex("A2")];
//      tableScalars[1] =Y_tmp[cgas.speciesIndex("naphthyne")];
//	tableScalars[2] =Y_tmp[cgas.speciesIndex("P2")];
        tableScalars[1] =Y_tmp[cgas.speciesIndex("A2R5")];
        tableScalars[2] =Y_tmp[cgas.speciesIndex("A3")];
        tableScalars[3] =Y_tmp[cgas.speciesIndex("A4")];
//        tableScalars[6] =Y_tmp[cgas.speciesIndex("BAPYR")];
//        tableScalars[5] =Y_tmp[cgas.speciesIndex("A4")];
//        tableScalars[8] =Y_tmp[cgas.speciesIndex("BEPYREN")];
//        tableScalars[9] =Y_tmp[cgas.speciesIndex("BGHIPER")];
//        tableScalars[10]= Y_tmp[cgas.speciesIndex("ANTHAN")];
//        tableScalars[11]= Y_tmp[cgas.speciesIndex("CORONEN")];

        tableScalars[nSc-8]=Y_tmp[cgas.speciesIndex("OH")];
        tableScalars[nSc-7]=Y_tmp[cgas.speciesIndex("H2O")];
        tableScalars[nSc-6]=Y_tmp[cgas.speciesIndex("H")];
        tableScalars[nSc-5]=Y_tmp[cgas.speciesIndex("H2")];
        tableScalars[nSc-4]=Y_tmp[cgas.speciesIndex("C2H2")];
        tableScalars[nSc-3]=Y_tmp[cgas.speciesIndex("O2")];
        tableScalars[nSc-2]=Y_tmp[cgas.speciesIndex("CO")];

        tableScalars[nSc-1] = T.internalField()[idx];
		
        int nStat = 4 + 10;
        double statArr[nStat];	//Outpou
	int other = 9;
	double otherArr[other];

//        if(MomSrcCalls)
//        {
//        Info << "Call semi_drive for calculating the source term" << endl;
//        Info << "tableScalars=" << tableScalars << endl;
//        Info << "MomArr=" << MomArr << endl;
//        Info << "MomSrcArr=" << MomSrcArr << endl;
//        Info << "statArr=" << statArr << endl;
	mom_driver_(&timeStep,&nSc,tableScalars,&tempRho,&tempWmix,&tempnMom,MomArr,MomSrcArr,&nStat,statArr,&other,otherArr);
//        Info << "In sootSrc.H: After call mom_driver" << nl << endl;
//	}
	src_Tem[idx]    = MomSrcArr[nSc-1];
	wCoeff[idx]     = otherArr[0];
	oxCoeff[idx]    = otherArr[1]; 
	Dimer_conc[idx] = otherArr[2];
	sootStar[idx]   = otherArr[3];
        betaN[idx]      = otherArr[4];
        betaC[idx]      = otherArr[5];
        prodRate[idx]   = otherArr[6];
        aromconc[idx]   = otherArr[7];
	Cdimer[idx]      = otherArr[8];

        for(int i=0; i<nMom; i++)
        {
		MomSrc[i][idx] = MomSrcArr[nSc+i];
               // Info << "MomSrc[i]=" << MomSrc[i][idx] << endl;
                //If moment is already below zero, cut source not to go negative
                if (Mom[i][idx]<0.0) 
                {
              		 MomSrc[i][idx]=max(MomSrc[i][idx],0.0);
                       //  Info << "Mom < 0, MomSrc[" << i << "]=" << MomSrc[i][idx] << endl;
                }
               //If moment is positive but source brings it negative, cut source
                else
                {
                	if (Mom[i][idx]+runTime.deltaT().value()*MomSrc[i][idx] < 0.0)
               		{
               			 MomSrc[i][idx] = -Mom[i][idx]/runTime.deltaT().value();
                         //        Info << "Mom+MomSrc < 0, MomSrc[" << i << "]=" << MomSrc[i][idx] << endl;
                        }
                }
	}
		
	double rate_C2H2, rate_O2, rate_OH;
        rate_C2H2 = -timeStep*MomSrcArr[4]/max(1.0E-60, tableScalars[nSc-4]);
        rate_O2   = -timeStep*MomSrcArr[5]/max(1.0E-60, tableScalars[nSc-3]);
        rate_OH   = -timeStep*MomSrcArr[6]/max(1.0E-60, tableScalars[nSc-8]);
        MomSrcArr[4] = MomSrcArr[4] / max(1.0, rate_C2H2);
        MomSrcArr[5] = MomSrcArr[5] / max(1.0, rate_O2);
        MomSrcArr[6] = MomSrcArr[6] / max(1.0, rate_OH);

        for(int i=0; i<nsp; i++)
        {
               if     (cgas.speciesName(i)=="A2"       ) {if (!std::isnan(MomSrcArr[0]))  src_Y[i][idx] = MomSrcArr[0]*PAHsrcOnOff;}
//               else if(cgas.speciesName(i)=="naphthyne") {if (!std::isnan(MomSrcArr[1]))  src_Y[i][idx] = MomSrcArr[1]*PAHsrcOnOff;}
//               else if(cgas.speciesName(i)=="P2"       ) {if (!std::isnan(MomSrcArr[2]))  src_Y[i][idx] = MomSrcArr[2]*PAHsrcOnOff;}
               else if(cgas.speciesName(i)=="A2R5"     ) {if (!std::isnan(MomSrcArr[1]))  src_Y[i][idx] = MomSrcArr[1]*PAHsrcOnOff;}
               else if(cgas.speciesName(i)=="A3"       ) {if (!std::isnan(MomSrcArr[2]))  src_Y[i][idx] = MomSrcArr[2]*PAHsrcOnOff;}
               else if(cgas.speciesName(i)=="A4"       ) {if (!std::isnan(MomSrcArr[3]))  src_Y[i][idx] = MomSrcArr[3]*PAHsrcOnOff;}
//               else if(cgas.speciesName(i)=="BAPYR") {if (!std::isnan(MomSrcArr[6]))  src_Y[i][idx] = MomSrcArr[6]*PAHsrcOnOff;}
//               else if(cgas.speciesName(i)=="A4"       ) {if (!std::isnan(MomSrcArr[5]))  src_Y[i][idx] = MomSrcArr[5]*PAHsrcOnOff;}
//               else if(cgas.speciesName(i)=="BEPYREN") {if (!std::isnan(MomSrcArr[8]))  src_Y[i][idx] = MomSrcArr[8]*PAHsrcOnOff;}
//               else if(cgas.speciesName(i)=="BGHIPER") {if (!std::isnan(MomSrcArr[9]))  src_Y[i][idx] = MomSrcArr[9]*PAHsrcOnOff;}
//               else if(cgas.speciesName(i)=="ANTHAN" ) {if (!std::isnan(MomSrcArr[10])) src_Y[i][idx] = MomSrcArr[10]*PAHsrcOnOff;}
//               else if(cgas.speciesName(i)=="CORONEN") {if (!std::isnan(MomSrcArr[11])) src_Y[i][idx] = MomSrcArr[11]*PAHsrcOnOff;}
               else if(cgas.speciesName(i)=="C2H2"     ) {if (!std::isnan(MomSrcArr[4])) src_Y[i][idx] = MomSrcArr[4]*YsrcOnOff;}
//               else if(cgas.speciesName(i)=="O2"     ) {if (!std::isnan(MomSrcArr[13])) src_Y[i][idx] = MomSrcArr[13]*YsrcOnOff;}
//               else if(cgas.speciesName(i)=="OH"     ) {if (!std::isnan(MomSrcArr[14])) src_Y[i][idx] = MomSrcArr[14]*YsrcOnOff;}
//               else if(cgas.speciesName(i)=="CO"     ) {if (!std::isnan(MomSrcArr[15])) src_Y[i][idx] = MomSrcArr[15]*YsrcOnOff;}
//               else if(cgas.speciesName(i)=="H"      ) {if (!std::isnan(MomSrcArr[16])) src_Y[i][idx] = MomSrcArr[16]*YsrcOnOff;}
	       else src_Y[i][idx] = 0.0;
        }
        if (!std::isnan(MomSrcArr[9])) src_rho[idx]  = MomSrcArr[9]*PAHsrcOnOff*rhoSrcOnOff;


	sootND[idx]    = statArr[0];
        sootFV[idx]    = statArr[1];
        sootDiam[idx]  = statArr[2];
        sootAggr[idx]  = statArr[3];

        sootNDsrcNucl[idx] = statArr[4];
        sootNDsrcCoag[idx] = statArr[5];
        sootNDsrcCond[idx] = statArr[6];
        sootNDsrcSurf[idx] = statArr[7];
        sootNDsrcOxid[idx] = statArr[8];

        sootFVsrcNucl[idx] = statArr[9];
        sootFVsrcCoag[idx] = statArr[10];
        sootFVsrcCond[idx] = statArr[11];
        sootFVsrcSurf[idx] = statArr[12];
        sootFVsrcOxid[idx] = statArr[13];
    }
//    forAll(Zmix.internalField(),cellI)
//    {
//        for(int i=0; i<nMom; i++)
//        {
//              if(Zmix.mesh().C()[cellI].y()>0.003) MomSrc[i][cellI]=1e-20;
//        }
//	if(Zmix.mesh().C()[cellI].y()>0.003) sootFV[cellI]=1e-20;
//    }
}

