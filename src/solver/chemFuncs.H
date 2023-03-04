//+++++============================Thermochemical update=============================+++++
#include "cantera/transport.h"
enum time_step {current_step, middle_step, last_step};

void updateThermoField (
volScalarField& thermoField,
volScalarField enthalpy, scalar p, const PtrList<volScalarField>& massFrac,
time_step input_step, IdealGasMix ct_gas
)
{
	int nsp=massFrac.size();
	std::auto_ptr<Transport> ct_gas_transport(newTransportMgr("Mix", &ct_gas));
	scalar H;
	doublereal Y[nsp];
	int countErr=0;
	forAll(thermoField.internalField(),cellI)
	{
		//What time step to use to set the thermo state
		switch(input_step)
		{
			case current_step:
				H=enthalpy.internalField()[cellI];
				for(int j=0;j<nsp;j++)	
				{
					Y[j]=massFrac[j].internalField()[cellI];
				}
				break;
			case middle_step:
				H=(enthalpy.internalField()[cellI]+enthalpy.oldTime().internalField()[cellI])/2.0;
				for(int j=0;j<nsp;j++)	
				{
					Y[j]=(massFrac[j].internalField()[cellI]+massFrac[j].oldTime().internalField()[cellI])/2.0;
				}
				break;
			case last_step:
				H=enthalpy.oldTime().internalField()[cellI];
				for(int j=0;j<nsp;j++)	
				{
					Y[j]=massFrac[j].oldTime().internalField()[cellI];
				}
				break;
			default:
				Info << "updateThermo: Error! Cannot update with unknown timestep " << endl;
				Info << "Use \"current_step\", \"middle_step\", or \"last_step\" only!" << endl;
				break;
		}
		//Update thermochemical field
		ct_gas.setMassFractions(Y);
		bool skipUpdate=false;
		try
		{
			ct_gas.setState_HP(H,p);
		} 
		catch(CanteraError ctErr)
		{
			countErr=countErr+1;
			skipUpdate=true;
		}
		if(ct_gas.temperature()<295) skipUpdate=true;
		if(!skipUpdate)
		{	
		if (thermoField.name()=="rho")
		{
			if(cellI==0) Info << "Updating density field" << endl;
			thermoField.ref()[cellI]=ct_gas.density();
		}
		else if (thermoField.name()=="mu")
		{
			if(cellI==0) Info << "Updating viscosity field" << endl;
			thermoField.ref()[cellI]=ct_gas_transport->viscosity();
		}
		else if (thermoField.name()=="CP")
		{
			if(cellI==0) Info << "Updating heat capacity field" << endl;
			thermoField.ref()[cellI]=ct_gas.cp_mass();
		}
		else if (thermoField.name()=="lambda")
		{
			if(cellI==0) Info << "Updating heat conductivity field" << endl;
			thermoField.ref()[cellI]=ct_gas_transport->thermalConductivity();
		}
		else if (thermoField.name()=="T")
		{
			if(cellI==0) Info << "Updating temperature field" << endl;
			thermoField.ref()[cellI]=ct_gas.temperature();
		}
		else
		{
			Info << "updateThermo: Error! Cannot update unknown field " << thermoField.name() << endl;
		}
		}
	}
	if(countErr>0) Pout << "Warning: A certain number of cells cannot be updated by updateThermoField: " << countErr << endl;
	//Update boundary
	thermoField.correctBoundaryConditions();
}

//Function overload
void updateThermoField (
PtrList<volScalarField>& thermoFields,
volScalarField enthalpy, scalar p, const PtrList<volScalarField>& massFrac,
time_step input_step, IdealGasMix ct_gas
)
{
	int nsp=massFrac.size();
	std::auto_ptr<Transport> ct_gas_transport(newTransportMgr("Mix", &ct_gas));
	doublereal mwt[nsp];
	ct_gas.getMolecularWeights(mwt);
	scalar H;
	doublereal Y[nsp];
	int countErr=0;
	forAll(thermoFields[0].internalField(),cellI)
	{
		//What time step to use to set the thermo state
		switch(input_step)
		{
			case current_step:
				H=enthalpy.internalField()[cellI];
				for(int j=0;j<nsp;j++)	
				{
					Y[j]=massFrac[j].internalField()[cellI];
				}
				break;
			case middle_step:
				H=(enthalpy.internalField()[cellI]+enthalpy.oldTime().internalField()[cellI])/2.0;
				for(int j=0;j<nsp;j++)	
				{
					Y[j]=(massFrac[j].internalField()[cellI]+massFrac[j].oldTime().internalField()[cellI])/2.0;
				}
				break;
			case last_step:
				H=enthalpy.oldTime().internalField()[cellI];
				for(int j=0;j<nsp;j++)	
				{
					Y[j]=massFrac[j].oldTime().internalField()[cellI];
				}
				break;
			default:
				Info << "updateThermo: Error! Cannot update with unknown timestep " << endl;
				Info << "Use \"current_step\", \"middle_step\", or \"last_step\" only!" << endl;
				break;
		}
		//Update thermochemical field
		ct_gas.setMassFractions(Y);
		bool skipUpdate=false;
		try
		{
			ct_gas.setState_HP(H,p);
		} 
		catch(CanteraError ctErr)
		{
			countErr=countErr+1;
			skipUpdate=true;
		}
		if(ct_gas.temperature()<295) skipUpdate=true;
		if(!skipUpdate)
		{	
		if (thermoFields[0].name()=="h_"+ct_gas.speciesName(0))
		{
			if(cellI==0) Info << "Updating species enthalpy fields" << endl;
			doublereal tmp[nsp];
			ct_gas.getPartialMolarEnthalpies(tmp);
			for (int n=0;n<nsp;n++)
			{
				thermoFields[n].ref()[cellI]=tmp[n]/mwt[n];
			}
		}
		else if (thermoFields[0].name()=="DIFF_"+ct_gas.speciesName(0))
		{
			if(cellI==0) Info << "Updating species diffusivity fields" << endl;
			doublereal tmp[nsp];
			ct_gas_transport->getMixDiffCoeffs(tmp);
			for (int n=0;n<nsp;n++)
			{
				thermoFields[n].ref()[cellI]=tmp[n];
			}
		}
		else
		{
			Info << "updateThermo: Error! Cannot update unknown fields " << thermoFields[0].name() << endl;
		}
		}
	}
	//Update boundary
	for(int n=0; n<nsp; n++)
	{
		thermoFields[n].correctBoundaryConditions();
	}
	if(countErr>0) Pout << "Warning: A certain number of cells cannot be updated by updateThermoField: " << countErr << endl;
}

//+++++===========================Stiff chemistry solver=============================+++++
#include "sundials/sundials_types.h"     /* defs. of realtype, sunindextype      */
#include "cvodes/cvodes.h"               /* prototypes for CVODE fcts., consts.  */
#include "nvector/nvector_serial.h"      /* access to serial N_Vector            */
#include "sunmatrix/sunmatrix_dense.h"   /* access to dense SUNMatrix            */
#include "sunlinsol/sunlinsol_dense.h"   /* access to dense SUNLinearSolver      */
#include "cvodes/cvodes_direct.h"        /* access to CVDls interface            */
using namespace Cantera;

// Specify the right hand side function in y'=f(t,y)
// For details see declaration of "CVRhsFn" in "include/cantera/ext/cvodes/cvodes.h"
static int SRC_func(realtype time, N_Vector y, N_Vector ydot, void *user_data) 
{
   IdealGasMix* ct_gas = static_cast<IdealGasMix*>(user_data);
   // store the h,p before massfraction change   // Add in V1.2
   doublereal enthalpy=ct_gas->enthalpy_mass();
   doublereal pressure=ct_gas->pressure();
   // extract data from sundials n_vector and put into cantera friendly types
   int nsp = ct_gas->nSpecies();
   doublereal* sun_data = N_VGetArrayPointer(y);
   doublereal Y_tmp[nsp];
   for(int n=0; n<nsp; n++)
   {
	Y_tmp[n]=*(sun_data+n);
   }
//   doublereal pressure=*(sun_data+nsp);
//   doublereal enthalpy=*(sun_data+nsp+1);
   // setState 
   ct_gas->setMassFractions(Y_tmp);
   ct_gas->setState_HP(enthalpy,pressure);
   // compute source term
   doublereal wdot[nsp];
   ct_gas->getNetProductionRates(wdot);
   doublereal mwt[nsp];
   ct_gas->getMolecularWeights(mwt);
   doublereal rho = ct_gas->density();
   for(int n=0; n<nsp; n++)
   {
	NV_Ith_S(ydot,n)=wdot[n]*mwt[n]/rho;
   }
//   NV_Ith_S(ydot,nsp)=0.0;
//   NV_Ith_S(ydot,nsp+1)=0.0;
   return 0;
}

