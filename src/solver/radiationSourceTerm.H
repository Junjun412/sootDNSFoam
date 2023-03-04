//-----Update moleFractions (for post processing or other usage)-----
forAll(XCO2.internalField(),cellI)
{
        //Update local gas massFractions
        double Y_tmp[nsp];
        for(int i=0;i<nsp;i++)
        {
                Y_tmp[i]=Y[i].internalField()[cellI];
        }
        cgas.setMassFractions(Y_tmp);
        //Extrate local gas moleFractions
        double X_tmp[nsp];
        cgas.getMoleFractions(X_tmp);
        for(int i=0;i<nsp;i++)
        {
                X[i].ref()[cellI]=X_tmp[i];
        }
        XCO2.ref()[cellI]=X_tmp[cgas.speciesIndex("CO2")];
        XH2O.ref()[cellI]=X_tmp[cgas.speciesIndex("H2O")];
        XCH4.ref()[cellI]=X_tmp[cgas.speciesIndex("CH4")];
        XCO.ref()[cellI] =X_tmp[cgas.speciesIndex("CO")];
//        C2H4.internalField()[cellI] =X_tmp[cgas.speciesIndex("C2H4")];
        element_C.ref()[cellI] = cgas.elementalMassFraction(0);
        element_H.ref()[cellI] = cgas.elementalMassFraction(1);
        element_O.ref()[cellI] = cgas.elementalMassFraction(2);
//        Zmix_bilger.internalField()[cellI] = (24.0/7.0*(cgas.elementalMassFraction('C')+cgas.elementalMassFraction('H'))-cgas.elementalMassFraction('O')+0.233)/(24.0/7.0*0.17609+0.233);
        Zmix_bilger.ref()[cellI] = (4.0*(cgas.elementalMassFraction(0)+cgas.elementalMassFraction(1))-cgas.elementalMassFraction(2)+0.488133)/(4.0*0.4149+0.488133);
//        X_PAH.ref()[cellI] =X_tmp[cgas.speciesIndex("A2")]+X_tmp[cgas.speciesIndex("A2R5")]+X_tmp[cgas.speciesIndex("A4")]+X_tmp[cgas.speciesIndex("P2")]+X_tmp[cgas.speciesIndex("A3")]+X_tmp[cgas.speciesIndex("BGHIF")]+X_tmp[cgas.speciesIndex("BAPYR")];
//	        X_CPAH.ref()[cellI] =10.0*X_tmp[cgas.speciesIndex("A2")]+12.0*X_tmp[cgas.speciesIndex("A2R5")]+16.0*X_tmp[cgas.speciesIndex("A4")]+12.0*X_tmp[cgas.speciesIndex("P2")]+14.0*X_tmp[cgas.speciesIndex("A3")]+18.0*X_tmp[cgas.speciesIndex("BGHIF")]+20.0*X_tmp[cgas.speciesIndex("BAPYR")];
}
//-----Calculating the Planck mean absorption coefficient for OTA-----
// The coefficients come from https://www.sandia.gov/TNF/radiation.html
forAll(a.internalField(),cellI)
{
	doublereal kCO2;
	doublereal kH2O;
	doublereal kCH4;
	doublereal kCO ;
        doublereal T_tmp;
	doublereal C;
	C = 1150.0;
        T_tmp = T.oldTime().ref()[cellI];
	T_tmp = min(T_tmp, 3000.0);
	T_tmp = max(T_tmp, 300.0 );

        kCO2=18.74100-1.2131e05/pow(T_tmp,1)+2.735e08/pow(T_tmp,2)-1.9405e11/pow(T_tmp,3)+5.631e13/pow(T_tmp,4)-5.8168e15/pow(T_tmp,5);
        kH2O=-0.23093-1.1239e03/pow(T_tmp,1)+9.4153e06/pow(T_tmp,2)-2.9988e09/pow(T_tmp,3)+5.1382e11/pow(T_tmp,4)-1.8684e10/pow(T_tmp,5); 
        kCH4= 6.63340-3.5686e-03*pow(T_tmp,1)+1.6682e-08*pow(T_tmp,2)+2.5611e-10*pow(T_tmp,3)-2.6558e-14*pow(T_tmp,4);
	if (T_tmp < 750.0)
        {
       	   kCO= 4.7869-6.9530e-02*pow(T_tmp,1)+2.95775e-04*pow(T_tmp,2)-4.25732e-07*pow(T_tmp,3)+2.02894e-10*pow(T_tmp,4);
        }
	else
	{
	   kCO= 10.09-1.183e-02*pow(T_tmp,1)+4.7753e-06*pow(T_tmp,2)-5.87209e-10*pow(T_tmp,3)-2.5334e-14*pow(T_tmp,4);
	}

        a.ref()[cellI]=XCO2.internalField()[cellI]*kCO2+XH2O.internalField()[cellI]*kH2O+XCH4.internalField()[cellI]*kCH4+XCO.internalField()[cellI]*kCO+C*sootFV.internalField()[cellI]*T_tmp;
}
// Calculating the radiation source term using OTA
forAll(radSRC.internalField(),cellI)
{
	doublereal Tsur = 300.0;
        doublereal sigma_SB = 5.669e-08;
        if(radiationOnOff)
	{
 	     radSRC.ref()[cellI]=-4.0*sigma_SB*a.internalField()[cellI]*(pow(T.oldTime().internalField()[cellI],4)-pow(Tsur,4));
	}
	else
	{
	      radSRC.ref()[cellI]=0.0;
	}
}

