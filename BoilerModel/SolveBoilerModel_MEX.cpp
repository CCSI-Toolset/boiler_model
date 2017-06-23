//SolvBoilerModel_MEX.cpp
//MEX function to solve the boiler model
#include "mex.h"
#include "RadiantFurnace.h"

int SolveBoilerModel(int* piconf, double* pdinput, double* pdoutput)
{
	//piconf:		integer array for configuration inputs, typically fixed inputs
	//pdinput:		double array for inputs
	//pdoutput:		double array for outputs
	//return the number of product variables or -1 if failed
	CRadiantFurnace furn;				//an instance of boiler model
	furn.ifurn_type = 0;				//currently rear wall exit furnace only
	furn.npsd = 1;						//currently use single particle size distribution only
	furn.nsh = 1;						//currently use single superheater only
	furn.nck = 1;						//currently use single coal kinetic data set
	int i, j;
	int icount;							//count for input
	int istream;						//count for inlet stream
	int iparaview = piconf[0];			//flag to generate ParaView files
	int ioutput_desc = piconf[1];		//flag to generate output description file
	int imass_frac = piconf[2];			//flag to indication mass fraction for gas phase composition
	int nlevel_burner = piconf[3];		//number of burner levels
	int nlevel_ofa = piconf[4];			//number of OFA or boundary air levels, maximum is 10
	int nzone_burner = piconf[5];		//number of zones in burner region from hopper knuckle to noze bottom, usually > number of burner levels
	int nzone_exit = piconf[6];			//number of zones above tip of the nose, typically 1 or 2
	furn.nx = piconf[7];				//number of cells in x (depth) direction
	furn.nz = piconf[8];				//number of cells in z (width) direction
	furn.nzone_burn = nzone_burner;		//number of zones in burner region
	furn.nzone_nose = 1;				//always use 1 zone in nose region
	furn.nzone_uppr = nzone_exit;		//same as the exit zone number
	furn.nzone = 1 + furn.nzone_burn + furn.nzone_nose + furn.nzone_uppr;		//total number of zones
	icount = 9;
	for (i=0; i<furn.nzone; i++)
		furn.pncell_zone[i] = piconf[icount++];	//number of cells in each zone
	for (i=0; i<furn.nsh; i++)			//currently nsh is set to 1
	{
		furn.psh[i].npanel = piconf[icount++];
		furn.psh[i].npoints = piconf[icount++];
	}
	furn.ppsd[0].nbin = piconf[icount++];	//number of size bins for particle size distribution
	furn.iprop_zone = piconf[icount++];		//flag to indicate if enclosure wall properties depend on zone
	//end of integer configuration read

	//read double variables related to geometry
	icount = 0;
	furn.zwidth = pdinput[icount++];		//furnace width in z direction
	furn.xdepth = pdinput[icount++];		//furnace depth in x direction
	furn.xhop_bot0 = pdinput[icount++];		//x coordinate of the point of hopper bottom at front wall
	furn.xhop_bot1 = pdinput[icount++];		//x coordinate of the point of hopper bottom at rear wall
	furn.xnose_tip = pdinput[icount++];		//x coordinate of the point of nose tip
	furn.yheight = pdinput[icount++];		//furnace height in y direction
	furn.yhop_knuckle = pdinput[icount++];			//y coordinate of hopper knuckle
	furn.ynose_bot = pdinput[icount++];		//y coordinate of the bottom of the nose slope
	furn.ynose_tip = pdinput[icount++];		//y coordinate of the nose tip
	//y coordinate at the top of each zone
	//the first zone is always at the hopper knuckle
	for (i=1; i<furn.nzone_burn; i++)
		furn.pyzone[i] = pdinput[icount++];
	//the top of the last burner zone is always at the nose tip
	for (i=0; i<furn.nzone_uppr-1; i++)
		furn.pyzone[i+furn.nzone-furn.nzone_uppr] = pdinput[icount++];
	//the top of the last exit zone is always at the furnace roof

	//superheater geometry and indices
	for (i=0; i<furn.nsh; i++)
	{
		for (j=0; j<furn.psh[i].npoints; j++)
		{
			furn.psh[i].xp[j] = pdinput[icount++];
			furn.psh[i].yp[j] = pdinput[icount++];
		}
		furn.psh[i].UpdateExtraPoint();
		furn.psh[i].iwall = i+1;
		furn.psh[i].iwater = i+1;
	}
	//location of primary and secondary air stream in each burner level
	for (i=0; i<nlevel_burner; i++)
	{
		furn.pyinlet[i] = pdinput[icount++];
		furn.pyinlet[nlevel_burner+i] = furn.pyinlet[i];
	}
	//location of OFA ports
	for (i=0; i<nlevel_ofa; i++)
		furn.pyinlet[2*nlevel_burner+i] = pdinput[icount++];

	//particle size distribution
	for (i=0; i<furn.ppsd[0].nbin; i++)
	{
		furn.ppsd[0].adp[i] = pdinput[icount++];
		furn.ppsd[0].mf[i] = pdinput[icount++];
	}

	//coal particle density
	furn.ppsd[0].denp = pdinput[icount++];

	//coal reactivity
	furn.pck[0].fswell = pdinput[icount++];
	furn.pck[0].alpha = pdinput[icount++];
	furn.pck[0].ecoco2 = pdinput[icount++];
	furn.pck[0].acoco2 = pdinput[icount++];
	furn.pck[0].echar_o2 = pdinput[icount++];
	furn.pck[0].achar_o2 = pdinput[icount++];
	furn.pck[0].echar_h2o = pdinput[icount++];
	furn.pck[0].achar_h2o = pdinput[icount++];
	furn.pck[0].echar_co2 = pdinput[icount++];
	furn.pck[0].achar_co2 = pdinput[icount++];

	//set stream data
	CMaterialStream	pa;
	CMaterialStream sa;
	CMaterialStream ofa[10];
	CMaterialStream flue_gas;
	//set primary air with coal stream
	pa.name = "Primary";
	pa.bgas = true;
	pa.bsolid = true;
	pa.bliquid = false;
	pa.temp = pdinput[icount++];
	pa.pres = pdinput[icount++];
	pa.gas.mflow = pdinput[icount++];
	pa.gas.ns = furn.nspecies;
	for (i=0; i<furn.nspecies; i++)
		pa.gas.pgas[i] = CIdealGas::GetGasByIndex(furn.ispecies[i]);
	for (i=0; i<furn.nspecies; i++)
	{
		if (imass_frac)		//based on mass fraction
		{
			pa.gas.imol = 1;
			pa.gas.fsp_mass[i] = pdinput[icount++];
		}
		else	//default imol = 0
			pa.gas.fsp_mole[i] = pdinput[icount++];
	}
	pa.gas.InitArraysFromFsp();
	pa.gas.CalcAllProperties();
	pa.coal.mflow = pdinput[icount++];
	pa.coal.fcmass[0][0] = pdinput[icount++];
	pa.coal.fcmass[0][1] = pdinput[icount++];
	pa.coal.fcmass[0][2] = pdinput[icount++];
	pa.coal.fcmass[0][3] = pdinput[icount++];
	pa.coal.fcmass[0][4] = pdinput[icount++];
	pa.coal.fcmass[0][5] = pdinput[icount++];
	pa.coal.fcmass[0][6] = pdinput[icount++];
	pa.coal.fcmass[0][7] = pdinput[icount++];
	pa.coal.hhv[0] = pdinput[icount++];
	pa.coal.CalcAllProperties();
	pa.UpdateProperties();
	//set secondary air stream
	sa.name = "Secondary";
	sa.bgas = true;
	sa.bsolid = false;
	sa.bliquid = false;
	sa.temp = pdinput[icount++];
	sa.pres = pdinput[icount++];
	sa.gas.mflow = pdinput[icount++];
	sa.gas.ns = furn.nspecies;
	for (i=0; i<furn.nspecies; i++)
		sa.gas.pgas[i] = CIdealGas::GetGasByIndex(furn.ispecies[i]);
	for (i=0; i<furn.nspecies; i++)
	{
		if (imass_frac)
		{
			sa.gas.imol = 1;
			sa.gas.fsp_mass[i] = pdinput[icount++];
		}
		else
			sa.gas.fsp_mole[i] = pdinput[icount++];
	}
	sa.gas.InitArraysFromFsp();
	sa.gas.CalcAllProperties();
	sa.UpdateProperties();
	//set ofa streams
	for (i=0; i<nlevel_ofa; i++)
	{
		ofa[i] = sa;
		ofa[i].gas.mflow = pdinput[icount++];
		ofa[i].UpdateProperties();
	}
	//inlet streams for combustion
	istream = 0;
	for (i=0; i<nlevel_burner; i++)
		furn.pstm_comb[istream++] = &pa;
	for (i=0; i<nlevel_burner; i++)
		furn.pstm_comb[istream++] = &sa;
	for (i=0; i<nlevel_ofa; i++)
		furn.pstm_comb[istream++] = &ofa[i];
	furn.ns_comb = istream;
	//assume ns_watr always 2 (1 for feed water and 1 for platen superheater)
	furn.ns_watr = 2;
	//outlet streams
	flue_gas.name = "Flue_Gas";
	furn.pstm_flue = &flue_gas;
	//wall properties
	//uniform enclosure wall properties
	if (furn.iprop_zone)		//specify properties in each zone
	{
		for (i=0; i<furn.nzone; i++)
			furn.pemis_zone[i] = pdinput[icount++];
		for (i=0; i<furn.nzone; i++)
			furn.prwall_zone[i] = pdinput[icount++];
		for (i=0; i<furn.nzone; i++)
			furn.ptback_zone[i] = pdinput[icount++];

	}
	else
	{
		furn.pemis_zone[0] = pdinput[icount++];
		furn.prwall_zone[0] = pdinput[icount++];
		furn.ptback_zone[0] = pdinput[icount++];
	}
	//superheater wall properties
	furn.pemis_sh[0] = pdinput[icount++];
	furn.prwall_sh[0] = pdinput[icount++];
	furn.ptback_sh[0] = pdinput[icount++];
	//exit wall properties
	furn.pemis_sh[1] = pdinput[icount++];
	furn.ptback_sh[1] = pdinput[icount++];
	//effectiveness factors
	furn.ef_kag = pdinput[icount++];
	furn.ef_kap = pdinput[icount++];
	furn.ef_rxn = pdinput[icount++];

	//finally solve the boiler model
	if (furn.Solve())		//failed to solve boiler model
		return -1;

	if (iparaview)
		furn.WriteVtkFiles("boiler_model.vtk");
	if (ioutput_desc)
		furn.WriteResultDescriptionFile("boiler_model_output_description.txt");
	std::vector<double> results;
	furn.PrepareResults(results);
	j = results.size();
	for (i=0; i<j; i++)
		pdoutput[i] = results[i];
	return j;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int *inIntMatrix;				/* 1xN integer input matrix */
    double *inDoubleMatrix;         /* 1xN double input matrix */
    double *outDoubleMatrix;        /* 1xN double output matrix */
    size_t ncols = 100;             /* size of double output matrix */

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is integer array */
    if( !mxIsInt16(prhs[0]) ||
         mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])<2 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notInteger","First input must be an integer array.");
    }

    /* check that number of rows in first input argument is 1 */
    if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","First input must be a row vector.");
    }

    /* check that number of rows in second input argument is 1 */
	    if(mxGetM(prhs[1])!=1) {
	        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Second input must be a row vector.");
    }

    /* get the value of the integer array input  */
    inIntMatrix = mxGetIr(prhs[0]);

    /* create a pointer to the real data in the input matrix  */
    inDoubleMatrix = mxGetPr(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outDoubleMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    SolveBoilerModel(inIntMatrix,inDoubleMatrix,outDoubleMatrix);
}
