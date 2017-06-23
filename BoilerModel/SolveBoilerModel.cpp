//SolvBoilerModel.cpp
//MEX function to solve the boiler model

#include "RadiantFurnace.h"
#include "SolveBoilerModel.h"

int SolveBoilerModel(int niinput, int ndinput, int* piconf, double* pdinput, double* pdoutput)
{
	//niinput:		number of integer inputs
	//ndinput:		number of double inputs
	//piconf:		integer array for configuration inputs, typically fixed inputs
	//pdinput:		double array for inputs
	//pdoutput:		double array for outputs
	//return the number of product variables or -1 if failed
	CRadiantFurnace furn;				//an instance of boiler model
	furn.npsd = 1;						//currently use single particle size distribution only
	furn.nck = 1;						//currently use single coal kinetic data set
	int i, j;
	int icount;							//count for input
	int istream;						//count for inlet stream
	int iparaview = piconf[0];			//flag to generate ParaView files
	int ioutput_desc = piconf[1];		//flag to generate output description file
	int imass_frac = piconf[2];			//flag to indication mass fraction for gas phase composition
	furn.ifurn_type = piconf[3];		//furnace type, 0=rear exit, 1=top exit
	int nlevel_burner = piconf[4];		//number of burner levels
	int nlevel_ofa = piconf[5];			//number of OFA or boundary air levels, maximum is 10
	int nzone_burner = piconf[6];		//number of zones in burner region from hopper knuckle to noze bottom, usually > number of burner levels
	int nzone_nose = piconf[7];			//number of zones in nose region from nose bottom to nose tip
	int nzone_exit = piconf[8];			//number of zones above tip of the nose, typically 1 or 2
	int nfuels;							//number of solid fuels in priamry stream
	furn.nx = piconf[9];				//number of cells in x (depth) direction
	furn.nz = piconf[10];				//number of cells in z (width) direction
	furn.nzone_burn = nzone_burner;		//number of zones in burner region
	furn.nzone_nose = nzone_nose;		//number of zones in nose region
	furn.nzone_uppr = nzone_exit;		//same as the exit zone number
	furn.nzone = 1 + furn.nzone_burn + furn.nzone_nose + furn.nzone_uppr;		//total number of zones
	icount = 11;
	for (i=0; i<furn.nzone; i++)
		furn.pncell_zone[i] = piconf[icount++];	//number of cells in each zone
	furn.psh[0].npanel = piconf[icount++];
	if (furn.psh[0].npanel <=0)		//use number of panels in input file as flag to indicate if there is a superheater
		furn.nsh = 0;
	else
	{
		furn.nsh = 1;				//default=1
		furn.psh[0].npoints = piconf[icount++];
	}
	nfuels = piconf[icount++];
	for (i=0; i<nfuels; i++)
		furn.ppsd[i].nbin = piconf[icount++];	//number of size bins for particle size distribution
	furn.iprop_zone = piconf[icount++];		//flag to indicate if enclosure wall properties depend on zone
	furn.iconvect = piconf[icount++];
	if (icount!=niinput)
	{
		printf("Number of integer inputs is incorrect!\n");
		return -1;
	}
	//end of integer configuration read

	//read double variables related to geometry
	icount = 0;
	furn.zwidth = pdinput[icount++];		//furnace width in z direction
	furn.xdepth = pdinput[icount++];		//furnace depth in x direction
	furn.xhop_bot0 = pdinput[icount++];		//x coordinate of the point of hopper bottom at front wall
	furn.xhop_bot1 = pdinput[icount++];		//x coordinate of the point of hopper bottom at rear wall
	furn.xnose_frt = pdinput[icount++];		//x coordinate of the point of front wall nose tip
	furn.xnose_tip = pdinput[icount++];		//x coordinate of the point of rear wall nose tip
	furn.yheight = pdinput[icount++];		//furnace height in y direction
	furn.yhop_knuckle = pdinput[icount++];			//y coordinate of hopper knuckle
	furn.ynose_bot = pdinput[icount++];		//y coordinate of the bottom of the nose slope
	furn.ynose_tip = pdinput[icount++];		//y coordinate of the nose tip
	//y coordinate at the top of each zone
	//the first zone is always at the hopper knuckle
	for (i=1; i<furn.nzone_burn; i++)
		furn.pyzone[i] = pdinput[icount++];
	//the top of the last burner zone is always at the nose bottom
	for (i=0; i<furn.nzone_nose-1; i++)
		furn.pyzone[i+furn.nzone_burn+1] = pdinput[icount++];
	//the top of the last nose zone is always at the nose tip
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
		for (j=1; j<=nfuels; j++)
		furn.pyinlet[nlevel_burner*j+i] = furn.pyinlet[i];
	}
	//location of OFA ports
	for (i=0; i<nlevel_ofa; i++)
		furn.pyinlet[(nfuels+1)*nlevel_burner+i] = pdinput[icount++];

	//size distribution and properties for each solid fuel
	for (j=0; j<nfuels; j++)
	{
		//particle size distribution
		for (i=0; i<furn.ppsd[j].nbin; i++)
		{
			furn.ppsd[j].adp[i] = pdinput[icount++];
			furn.ppsd[j].mf[i] = pdinput[icount++];
		}
		//coal particle density
		furn.ppsd[j].denp = pdinput[icount++];
		//coal reactivity
		furn.pck[j].fswell = pdinput[icount++];
		furn.pck[j].alpha = pdinput[icount++];
		furn.pck[j].ecoco2 = pdinput[icount++];
		furn.pck[j].acoco2 = pdinput[icount++];
		furn.pck[j].echar_o2 = pdinput[icount++];
		furn.pck[j].achar_o2 = pdinput[icount++];
		furn.pck[j].nchar_o2 = pdinput[icount++];
		furn.pck[j].bchar_o2 = furn.pck[j].nchar_o2;
		furn.pck[j].echar_h2o = pdinput[icount++];
		furn.pck[j].achar_h2o = pdinput[icount++];
		furn.pck[j].nchar_h2o = pdinput[icount++];
		furn.pck[j].bchar_h2o = furn.pck[j].nchar_h2o;
		furn.pck[j].echar_co2 = pdinput[icount++];
		furn.pck[j].achar_co2 = pdinput[icount++];
		furn.pck[j].nchar_co2 = pdinput[icount++];
		furn.pck[j].bchar_co2 = furn.pck[j].nchar_co2;
	}
	
	//set stream data
	CMaterialStream	pa;
	CMaterialStream pa_extra_fuel[3];
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
	pa.coal.fcmass[0][8] = pdinput[icount++];
	pa.coal.hhv[0] = pdinput[icount++];
	pa.coal.CalcAllProperties();
	pa.UpdateProperties();
	//set additional primary streams if more than 1 solid fuel
	for (j=0; j<nfuels-1; j++)
	{
		pa_extra_fuel[j].name = "Primary_extra_fuel";
		pa_extra_fuel[j].bgas = false;
		pa_extra_fuel[j].bsolid = true;
		pa_extra_fuel[j].bliquid = false;
		pa_extra_fuel[j].temp = pa.temp;
		pa_extra_fuel[j].pres = pa.pres;
		pa_extra_fuel[j].coal.mflow = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][0] = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][1] = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][2] = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][3] = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][4] = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][5] = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][6] = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][7] = pdinput[icount++];
		pa_extra_fuel[j].coal.fcmass[0][8] = pdinput[icount++];
		pa_extra_fuel[j].coal.hhv[0] = pdinput[icount++];
		pa_extra_fuel[j].coal.CalcAllProperties();
		pa_extra_fuel[j].UpdateProperties();
	}
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
	for (j=0; j<nfuels-1; j++)
	{
		for (i=0; i<nlevel_burner; i++)
			furn.pstm_comb[istream++] = &pa_extra_fuel[j];
	}
	for (i=0; i<nlevel_burner; i++)
		furn.pstm_comb[istream++] = &sa;
	for (i=0; i<nlevel_ofa; i++)
		furn.pstm_comb[istream++] = &ofa[i];
	furn.ns_comb = istream;
	//assume ns_watr always 2 (1 for feed water and 1 for platen superheater)
	furn.ns_watr = 2;
	if (furn.nsh == 0)	//if there is no superheater
		furn.ns_watr = 1;
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
	if (furn.nsh>0)	//has one superheater, superheater and exit properties
	{
		furn.pemis_sh[0] = pdinput[icount++];
		furn.prwall_sh[0] = pdinput[icount++];
		furn.ptback_sh[0] = pdinput[icount++];
		furn.pemis_sh[1] = pdinput[icount++];
		furn.ptback_sh[1] = pdinput[icount++];
	}
	else //no superheater, exit properties
	{
		furn.pemis_sh[0] = pdinput[icount++];
		furn.ptback_sh[0] = pdinput[icount++];
	}
	//effectiveness factors
	furn.ef_kag = pdinput[icount++];
	furn.ef_kap = pdinput[icount++];
	for (i=0; i<furn.nzone; i++)
		furn.ef_rxn[i] = pdinput[icount++];
	if (icount!=ndinput)
	{
		printf("Number of double inputs is incorrect!\n");
		return -1;
	}

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