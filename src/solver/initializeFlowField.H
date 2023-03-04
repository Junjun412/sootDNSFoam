Info << "Initialize enthalpy field" << endl;
forAll(rho.internalField(),cellI)
{
	//Enable this block if you want to initialze the whole domain with air, otherwise you need to make sure time folders "0" contains mass fraction that make physical sense.
	/*
	for(int i=0; i<nsp; i++)
	{
		if(Y[i].name()=="Y_O2")
		{
			Y[i].internalField()[cellI]=0.23;
		}
		else if(Y[i].name()=="Y_N2")
		{
			Y[i].internalField()[cellI]=0.77;
		}
	}
	*/
	//Initialize enthalpy field based on temperature field
	doublereal Y_tmp[Y.size()];
	for(int j=0; j<Y.size(); j++)	
	{
		Y_tmp[j]=Y[j].internalField()[cellI];
	}
	cgas.setState_TPY(T.internalField()[cellI],press,Y_tmp);
	h.ref()[cellI]=cgas.enthalpy_mass();
}
#include "chemThermoUpdate.H"
