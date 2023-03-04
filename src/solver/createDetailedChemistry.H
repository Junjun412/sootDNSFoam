IOdictionary inputParameters
(
    IOobject
    (
        "inputParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
bool smoothRho (readBool(inputParameters.lookup("smoothRho")));
dimensionedScalar DensityScaling (inputParameters.lookup("DensityScaling"));
dimensionedScalar freestreamDensity (inputParameters.lookup("freestreamDensity"));
bool poissonPEqn (readBool(inputParameters.lookup("poissonPEqn")));
bool thermophoresis (readBool(inputParameters.lookup("thermophoresis")));
bool radiationOnOff (readBool(inputParameters.lookup("radiationOnOff")));

if ( Pstream::master() ) printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
if ( Pstream::master() ) printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );

//+++++==============Read in related controls and chemical mechanism=================+++++
scalar press(readScalar(inputParameters.lookup("press"))); // this thermo pressure is different from the pressure solved in Poisson eqn.
bool initialize(readBool(inputParameters.lookup("initialize"))); // Initialize the enthalpy field according to current temperature and massfraction field
bool equilibrate(readBool(inputParameters.lookup("equilibrate"))); // Initialize the flow field with chemical equilibrium state
bool preffDiff(readBool(inputParameters.lookup("preffDiff")));   // Use non-unity Lewis number in species and enthalpy equation
Info << "createDetailedChemistry.H: Initializing Cantera"<< endl;
string canteraFile (word(inputParameters.lookup("canteraFile")));
string inertSpecies (word(inputParameters.lookup("inertSpecies")));

IdealGasMix cgas(canteraFile,"gas");
//number of species
int nsp = static_cast<int>(cgas.nSpecies());
// molecular weights (kg/kmol)
std::vector<doublereal> mwt;
mwt.reserve(nsp);
cgas.getMolecularWeights(mwt);
// heat of formation (J/kg)
std::vector<doublereal> hform;
hform.reserve(nsp);
for( int i=0; i<nsp; i++ ) hform[i] = cgas.Hf298SS(i) / mwt[i];
// net production rate (kmol/m^3/s)
doublereal wdot[nsp];

if ( Pstream::master() )
{
	cout << "\n";
	cout << "  canteraThermoFile : " << canteraFile << "\n";
	printf( "  number of species : %d\n", nsp );
	printf( "\n  | %12s---%20s---%19s | \n", "------------", "--------------------", "-------------------" );
	printf( "  | %12s | %20s | %19s | \n", "species name", 	"molec mass (kg/kmol)", "heat of form (J/kg)" );
	printf( "  | %12s | %20s | %19s | \n", "------------", "--------------------", "-------------------" );
	for( int i=0; i<nsp; i++ ) printf("  | %12s | %20.5f | %19.6e | \n", cgas.speciesName(i).c_str(), mwt[i], hform[i] );
	printf( "  | %12s---%20s---%19s | \n", "------------", "--------------------", "-------------------" );
}

//+++++=============variables related to stiff chemistry ode solver==================+++++
//general problem parameters for sundials
const realtype T0 = RCONST(0.0);                                            // time zero
const realtype reltol = RCONST(1.0e-8);                                     // tolerances
const realtype abstol = RCONST(1.0e-15);
//general problem variables for sundials
int flag;                                                                   // just to hold the error-checking flag, I don't intend to use it.
realtype ode_time = T0;	                                                    // time variable tracked by CVode
sunindextype NEQ = nsp;                                                   // number of dependent vars ({Y,p,h}: n_species+2).
N_Vector y = NULL;                                                          // empty vector for storing solution
SUNMatrix A = NULL;                                                         // empty dense matrix for solver
SUNLinearSolver LS = NULL;                                                  // empty dense linear solver

// Create serial vector of length NEQ for I.C.
y = N_VNew_Serial(NEQ);
//Create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration
void *cvode_mem = NULL;
cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
//Initialize the integrator memory and specify the user's rhs function, the inital time, and the initial dependent variable vector y.
flag = CVodeInit(cvode_mem, SRC_func, ode_time, y);
//Specify the scalar relative tolerance and vector absolute tolerances
flag = CVodeSStolerances(cvode_mem, reltol, abstol);
//Specify the maxinum numbers of ODE integration, default is 500. But don't hesitate to try much larger number if you get CVODE error message "mxstep steps taken before reaching tout."
flag = CVodeSetMaxNumSteps(cvode_mem, 3000);
//Attach the gas object to CVode
flag = CVodeSetUserData(cvode_mem, &cgas);
//Create dense SUNMatrix for use in linear solves
A = SUNDenseMatrix(NEQ, NEQ);
//Create dense SUNLinearSolver object for use by CVode
LS = SUNDenseLinearSolver(y, A);
//Attach the matrix and linear solver to CVode
flag = CVDlsSetLinearSolver(cvode_mem, LS, A);

//+++++================field variables related to chemistry==========================+++++
PtrList<volScalarField> Y(nsp);
PtrList<volScalarField> X(nsp);
PtrList<volScalarField> DIFF(nsp);
PtrList<volScalarField> SRC(nsp);
PtrList<volScalarField> hi(nsp);
PtrList<volScalarField> src_Y(nsp);  //Species source terms from soot formation
for ( int i=0; i<nsp; i++ ) {
	if(cgas.speciesName(i)=="O2"||cgas.speciesName(i)=="N2"||cgas.speciesName(i)=="CH4" \
	 ||cgas.speciesName(i)=="H2O"||cgas.speciesName(i)=="CO2"||cgas.speciesName(i)=="C2H4") // species that you want to initialize from time folder
	{
		Y.set (
			i,
			new volScalarField (
				IOobject (
					"Y_"+cgas.speciesName(i),
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);
	}
	else //species that are initialize with default settings, i.e., uniform 0
	{
		Y.set (
			i,
			new volScalarField (
				IOobject (
					"Y_"+cgas.speciesName(i),
					runTime.timeName(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE
				),
				mesh,
    				dimensionedScalar("Y_"+cgas.speciesName(i), dimensionSet(0,0,0,0,0,0,0), 0)
			)
		);
		//boundary conditions for default species, IMPORTANT: modify this part based on your geometry
		//boundary conditions for default species, IMPORTANT: modify this part based on your geometry
		//boundary conditions for default species, IMPORTANT: modify this part based on your geometry
            	const polyBoundaryMesh& patches = mesh.boundaryMesh();
            	forAll(patches,patchID)
		{
			if(patches[patchID].name()=="fuel"||patches[patchID].name()=="air"||patches[patchID].name()=="oxidant")
			{
				Y[i].boundaryFieldRef().set(patchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[patchID], Y[i]));
				fixedValueFvPatchScalarField& YPatch=refCast<fixedValueFvPatchScalarField>(Y[i].boundaryFieldRef()[patchID]);
				scalarField& YField = YPatch;
				YField=0.0;
			}
			else if((patches[patchID].name()=="axis"))
			{
				Y[i].boundaryFieldRef().set(patchID,fvPatchField<scalar>::New("empty", mesh.boundary()[patchID], Y[i]));
			}
			else if(patches[patchID].name()=="front"||patches[patchID].name()=="back")
			{
				Y[i].boundaryFieldRef().set(patchID,fvPatchField<scalar>::New("wedge", mesh.boundary()[patchID], Y[i]));
			}
			else
			{
				Y[i].boundaryFieldRef().set(patchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchID], Y[i]));
			}
		}
		Y[i].correctBoundaryConditions(); //This line is needed, otehrwise solution singularity
	}

        if(cgas.speciesName(i)=="A3"||cgas.speciesName(i)=="C2H2"||cgas.speciesName(i)=="A2R5" \
         ||cgas.speciesName(i)=="A2"||cgas.speciesName(i)=="A4") // species that you want to save
        {
                X.set (
                        i,
                        new volScalarField (
                                IOobject (
                                        "X_"+cgas.speciesName(i),
                                        runTime.timeName(),
                                        mesh,
                                        IOobject::NO_READ,
                                        IOobject::AUTO_WRITE
                                ),
                                mesh,
                                dimensionedScalar("X_"+cgas.speciesName(i), dimensionSet(0,0,0,0,0,0,0), 0)
                        )
                );
        }
        else //species you don't want to save
        {
                X.set (
                        i,
                        new volScalarField (
                                IOobject (
                                        "X_"+cgas.speciesName(i),
                                        runTime.timeName(),
                                        mesh,
                                        IOobject::NO_READ,
                                        IOobject::AUTO_WRITE
                                ),
                                mesh,
                                dimensionedScalar("X_"+cgas.speciesName(i), dimensionSet(0,0,0,0,0,0,0), 0)
                        )
                );
         }

        if(cgas.speciesName(i)=="C2H2"||cgas.speciesName(i)=="CO"||cgas.speciesName(i)=="A2")
        {
                src_Y.set (
                        i,
                        new volScalarField (
                                IOobject (
                                        "src_Y_"+cgas.speciesName(i),
                                        runTime.timeName(),
                                        mesh,
                                        IOobject::NO_READ,
                                        IOobject::AUTO_WRITE
                                ),
                                mesh,
                                dimensionedScalar("src_Y_"+cgas.speciesName(i), dimensionSet(1,-3,-1,0,0,0,0), 0)
                        )
                );
        }
        else //species you don't want to save
        {
                src_Y.set (
                        i,
                        new volScalarField (
                                IOobject (
                                        "src_Y_"+cgas.speciesName(i),
                                        runTime.timeName(),
                                        mesh,
                                        IOobject::NO_READ,
                                        IOobject::NO_WRITE
                                ),
                                mesh,
                                dimensionedScalar("src_Y_"+cgas.speciesName(i), dimensionSet(1,-3,-1,0,0,0,0), 0)
                        )
                );
         }

	DIFF.set (
		i,
		new volScalarField (
			IOobject (
				"DIFF_"+cgas.speciesName(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
    			dimensionedScalar("DIFF_"+cgas.speciesName(i), dimensionSet(0,2,-1,0,0,0,0), 1e-5)
		)
	);
	SRC.set (
		i,
		new volScalarField (
			IOobject (
				"SRC_"+cgas.speciesName(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
    			dimensionedScalar("SRC_"+cgas.speciesName(i), dimensionSet(1,-3,-1,0,0,0,0), 0.0) //Notice this is not wdot, but the actual source in YEqn
		)
	);
	hi.set (
		i,
		new volScalarField (
			IOobject (
				"h_"+cgas.speciesName(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
    			dimensionedScalar("h_"+cgas.speciesName(i), dimensionSet(0,2,-2,0,0,0,0), 0.0)
		)
	);
}
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
for (int i=0;i<nsp;i++) {
    	fields.add(Y[i]);
	if(Y[i].name() == ("Y_"+inertSpecies))
	{
		Info << inertSpecies << " is used as inert species" << endl;
	}
}



volScalarField rho
(
    IOobject
    (   
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh//,
    //dimensionedScalar("rho", dimensionSet(1,-3,0,0,0,0,0), 0.0)
);

volScalarField T
(
    IOobject
    (   
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField h
(
    IOobject
    (   
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh//,
    //dimensionedScalar("h", dimensionSet(0,2,-2,0,0,0,0), 0.0)
);


volScalarField mu
(
    IOobject
    (   
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh//,
    //dimensionedScalar("mu", dimensionSet(1,-1,-1,0,0,0,0), 0.0)
);

volScalarField lambda
(
    IOobject
    (   
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh//,
    //dimensionedScalar("lambda", dimensionSet(1,1,-3,-1,0,0,0), 0.0)
);

volScalarField CP
(
    IOobject
    (   
        "CP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh//,
    //dimensionedScalar("CP", dimensionSet(0,2,-2,-1,0,0,0), 0.0)
);

//Chemical timescale
volScalarField t_chem
(
    IOobject
    (   
        "t_chem",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("t_chem", dimensionSet(0,0,1,0,0,0,0), 0.0)
);
volScalarField stiff_ode
(
    IOobject
    (   
        "stiff_ode",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("stiff_ode", dimensionSet(0,0,1,0,0,0,0), 0.0)
);



//Preferential diffusion
volVectorField U_correct
(
    IOobject
    (
        "U_correct",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U_correct", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0))
    //0.0*linearInterpolate(rho*U) & mesh.Sf()
);

surfaceScalarField phi_correct
(
    IOobject
    (
        "phi_correct",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("phi_correct", dimensionSet(1,0,-1,0,0,0,0), 0.0)
    //0.0*linearInterpolate(rho*U) & mesh.Sf()
);

volScalarField laplacePreffDiff
(
    IOobject
    (
        "laplacePreffDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::laplacian(lambda/CP, h)
);

//Blending factor
surfaceScalarField hBlendingFactor
(
    IOobject
    (
        "hBlendingFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::mag(mesh.Sf())/Foam::mag(mesh.Sf())*1.0
);
surfaceScalarField h_0BlendingFactor
(   
    IOobject
    (   
        "h_0BlendingFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Foam::mag(mesh.Sf())/Foam::mag(mesh.Sf())*1.0
);

//absorption coefficient

volScalarField a
(
    IOobject (
            "a",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("a", dimensionSet(0,-1,0,0,0,0,0), 0)
);

//moleFractions output
volScalarField XCO2
(
    IOobject (
            "XCO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("XCO2", dimensionSet(0,0,0,0,0,0,0), 1E-10)
);

volScalarField XH2O
(
    IOobject (
            "XH2O",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("XH2O", dimensionSet(0,0,0,0,0,0,0), 1E-10)
);

volScalarField XCH4
(
    IOobject (
            "XCH4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("XCH4", dimensionSet(0,0,0,0,0,0,0), 0)
);

volScalarField XCO
(
    IOobject (
            "XCO",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("XCO", dimensionSet(0,0,0,0,0,0,0), 0)
);

// Radaition source term
volScalarField radSRC
(
    IOobject
    (
        "radSRC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("radSRC", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
);

int nMom = 4;
// Soot moment term array
PtrList<volScalarField> Mom(nMom);
// Soot moment source term array
PtrList<volScalarField> MomSrc(nMom);

for(int idx=0; idx<nMom; idx++)
{
    Mom.set
    (
         idx,
         new volScalarField
         (
          IOobject
          (
           "Mom_" + std::to_string(idx),
           runTime.timeName(),
           mesh,
           IOobject::MUST_READ,
           IOobject::AUTO_WRITE
           ),
          mesh
          )
         );
}

for(int idx=0; idx<nMom; idx++)
{
    MomSrc.set
    (
         idx,
         new volScalarField
         (
          IOobject
          (
           "MomSrc_"+ std::to_string(idx),
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
           ),
          mesh,
          dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
          )
         );
}

volScalarField src_Tem
(
    IOobject
    (
        "src_Tem",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("src_Tem", dimensionSet(0,0,0,0,0,0,0), 0.0) //unit?
);

volScalarField src_rho
(
    IOobject
    (
        "src_rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("src_rho", dimensionSet(1,-3,-1,0,0,0,0), 0.0) //unit?
);

// Soot postprocess
volScalarField sootND
(
    IOobject
    (
        "sootND",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootND", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootFV
(
    IOobject
    (
        "sootFV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootFV", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootDiam
(
    IOobject
    (
        "sootDiam",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootDiam", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootAggr
(
    IOobject
    (
        "sootAggr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootAggr", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootNDsrcNucl
(
    IOobject
    (
        "sootNDsrcNucl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sootNDsrcNucl", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootNDsrcCoag
(
    IOobject
    (
        "sootNDsrcCoag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sootNDsrcCoag", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootNDsrcCond
(
    IOobject
    (
        "sootNDsrcCond",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sootNDsrcCond", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootNDsrcSurf
(
    IOobject
    (
        "sootNDsrcSurf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sootNDsrcSurf", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootNDsrcOxid
(
    IOobject
    (
        "sootNDsrcOxid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sootNDsrcOxid", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootFVsrcNucl
(
    IOobject
    (
        "sootFVsrcNucl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootFVsrcNucl", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootFVsrcCoag
(
    IOobject
    (
        "sootFVsrcCoag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootFVsrcCoag", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootFVsrcCond
(
    IOobject
    (
        "sootFVsrcCond",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootFVsrcCond", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootFVsrcSurf
(
    IOobject
    (
        "sootFVsrcSurf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootFVsrcSurf", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootFVsrcOxid
(
    IOobject
    (
        "sootFVsrcOxid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sootFVsrcOxid", dimensionSet(0,0,0,0,0,0,0),0.0)
);


volScalarField sootDIFF
(
    IOobject
    (
        "sootDIFF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sootDIFF", dimensionSet(0,2,-1,0,0,0,0),10e-8)
);

volScalarField element_C
(
    IOobject
    (
        "element_C",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("element_C", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField element_H
(
    IOobject
    (
        "element_H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("element_H", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField element_O
(
    IOobject
    (
        "element_O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("element_O", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField Zmix_bilger
(
    IOobject
    (
        "Zmix_bilger",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Zmix_bilger", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField aromconc
(
      IOobject
      (
         "aromconc",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("aromconc", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField betaC
(
     IOobject
     (
       "betaC",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("betaC", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField prodRate
(
     IOobject
     (
       "prodRate",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("prodRate", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField sootStar
(
     IOobject
     (
       "sootStar",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("sootStar", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField betaN
(
     IOobject
     (
        "betaN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("betaN", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField oxCoeff
(
     IOobject
         (
	   "oxCoeff",
	   runTime.timeName(),
	   mesh,
	   IOobject::NO_READ,
	   IOobject::AUTO_WRITE
	 ),
	 mesh,
	 dimensionedScalar("oxCoeff", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField Dimer_conc
(
     IOobject
     (
        "Dimer_conc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("Dimer_conc", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField wCoeff
(
     IOobject
     (
       "wCoeff",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("wCoeff", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField X_PAH
(
       IOobject
       (
         "X_PAH",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
       ),
       mesh,
       dimensionedScalar("X_PAH", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField X_CPAH
(
        IOobject
        (
          "X_CPAH",
          runTime.timeName(),
	  mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("X_CPAH", dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField dummy
(
    IOobject
    (
        "dummy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dummy", dimensionSet(0,0,0,0,0), 1.0)
);

volScalarField Cdimer
(
    IOobject
    (
        "Cdimer",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Cdimer", dimensionSet(0,0,0,0,0,0,0),0.0)
);

if ( Pstream::master() ) printf( "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
