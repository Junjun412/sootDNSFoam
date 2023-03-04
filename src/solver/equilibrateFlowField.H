Info << "Calculate equilibrium state" << endl;

// Compute equilibrium species composition
forAll(h.internalField(),cellI)
{
	doublereal Y_tmp[nsp];
	for(int i=0;i<nsp;i++)
	{
		Y_tmp[i]=Y[i].internalField()[cellI];
	}
	cgas.setMassFractions(Y_tmp);
	cgas.setState_HP(h.internalField()[cellI],press);
	cgas.equilibrate("TP","vcs");
	//cgas.equilibrate("HP","vcs",1e-10,999,999); //If simulation get stuck (probably at chemThermoUpdate), try use max_iter naumber for equilibrium calculation
	cgas.getMassFractions(Y_tmp);
	for(int i=0;i<nsp;i++)
	{
		//if(Zmix.mesh().C()[cellI].x()>0.002) Y[i].internalField()[cellI]=Y_tmp[i];
		Y[i].ref()[cellI]=Y_tmp[i];
	}
}
//Update thermo properties
for(int i=0;i<nsp;i++)
{
	Y[i].oldTime()=Y[i];
	Y[i].oldTime().oldTime()=Y[i];
}
#include "chemThermoUpdate.H"

