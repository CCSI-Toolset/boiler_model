// RadiantFurnace.cpp: implementation of the CRadiantFurnace class.
//
//////////////////////////////////////////////////////////////////////
#include <cmath>
#include "Constants.h"
#include "ParticleRadiationProperty.h"
#include "RadiantFurnace.h"
#include "UtilityFunctions.h"

#if defined WIN32 || defined WIN64
#include "f2c.h"
#endif

extern "C" void gasemissivity_(double* t, double* l, double* p, double f[6], double* fs, double* e);

CRadiantFurnace::CRadiantFurnace()
{
	int i;
	//effectiveness factors as model tuning parameters, currently as user inputs
	ef_kap = 1;
	ef_kag = 1;
	for (i=0; i<RF_NZONE; i++)
		ef_rxn[i] = 1;
	//hard-wired under-relaxation factors for quick and stable convergence
	urf_temp = 0.5;
	urf_kap = 0.5;
	urf_hl	= 0.99;
	urf_daf = 0.1;
	//species selected for the gas phase
	nspecies = 12;		//currently hard-wired, debug, used to be 11
	ispecies[0] = 0;	//C(S)
	ispecies[1] = 1;	//O2
	ispecies[2] = 2;	//N2
	ispecies[3] = 3;	//H2
	ispecies[4] = 4;	//CO
	ispecies[5] = 5;	//CO2
	ispecies[6] = 6;	//H2O
	ispecies[7] = 8;	//SO2
	ispecies[8] = 9;	//H2S
	ispecies[9] = 11;	//CH4
	ispecies[10] = 27;	//Ar
	ispecies[11] = 28;	//HCl
	//other member variables
	pres = 101325;
	bmodified = true;
	iconvect = 0;
	ifurn_type = 0;
	nzone_burn = 2;
	nzone_nose = 1;
	nzone_uppr = 1;
	nzone = 1+nzone_burn+nzone_nose+nzone_uppr;
	nx = 10;
	nz = 10;
	ns_comb = 1;
	ns_watr = 1;
	nsh = 0;
	npsd = 1;
	iprop_zone = 0;
	for (i=0; i<RF_NZONE; i++)
		pncell_zone[i]= 2;
	xdepth = 5;
	xhop_bot0 = 2;
	xhop_bot1 = 3;
	xnose_frt = 0;
	xnose_tip = 4;
	yheight = 20;
	yhop_knuckle = 5;
	ynose_bot = 14;
	ynose_tip = 16;
	zwidth = 5;
	for (i=0; i<RF_NSH+1; i++)
	{
		pemis_sh[i] = 0.7;
		prwall_sh[i] = 0.001;
		ptback_sh[i] = 500;
	}
	for (i=0; i<RF_NZONE; i++)
	{
		pyzone[i] = 10;
		pemis_zone[i] = 0.7;
		prwall_zone[i] = 0.001;
		ptback_zone[i] = 500;
	}
	pyzone[0] = yhop_knuckle;
	pyzone[1] = (yhop_knuckle+ynose_bot)/2;
	pyzone[2] = ynose_bot;
	pyzone[3] = ynose_tip;
	pyzone[4] = yheight;
	for (i=0; i<RF_NS_COMB; i++)
	{
		pyinlet[i] = 10;
		pstm_comb[i] = NULL;
	}
}

CRadiantFurnace::~CRadiantFurnace()
{
	//do nothing.
	//cleanup of instreams and outstreams are performed by driver code
}

CRadiantFurnace::CRadiantFurnace(const CRadiantFurnace &t)
{
	//does not include hard-wired variables
	int i, j;
	//user input variables
	bmodified = t.bmodified;
	iconvect = t.iconvect;
	ifurn_type = t.ifurn_type;
	iprop_zone = t.iprop_zone;
	nzone_burn = t.nzone_burn;
	nzone_nose = t.nzone_nose;
	nzone_uppr = t.nzone_uppr;
	nx = t.nx;
	nz = t.nz;
	ns_comb = t.ns_comb;
	ns_watr = t.ns_watr;
	nsh = t.nsh;
	npsd = t.npsd;
	nck = t.nck;
	for (i=0; i<nzone; i++)
		pncell_zone[i]= t.pncell_zone[i];
	pres = t.pres;
	xdepth = t.xdepth;
	xhop_bot0 = t.xhop_bot0;
	xhop_bot1 = t.xhop_bot1;
	xnose_frt = t.xnose_frt;
	xnose_tip = t.xnose_tip;
	yheight = t.yheight;
	yhop_knuckle = t.yhop_knuckle;
	ynose_bot = t.ynose_bot;
	ynose_tip = t.ynose_tip;
	zwidth = t.zwidth;
	for (i=0; i<nsh+1; i++)
	{
		pemis_sh[i] = t.pemis_sh[i];
		prwall_sh[i] = t.prwall_sh[i];
		ptback_sh[i] = t.ptback_sh[i];
	}
	for (i=0; i<nzone; i++)
	{
		pemis_zone[i] = t.pemis_zone[i];
		prwall_zone[i] = t.prwall_zone[i];
		ptback_zone[i] = t.ptback_zone[i];
		pyzone[i] = t.pyzone[i];
	}
	for (i=0; i<ns_comb; i++)
		pyinlet[i] = t.pyinlet[i];
	for (i=0; i<nsh; i++)
		psh[i] = t.psh[i];
	for (i=0; i<npsd; i++)
		ppsd[i] = t.ppsd[i];
	for (i=0; i<nck; i++)
		pck[i] = t.pck[i];
	//calculated variables
	nzone = t.nzone;
	ny = t.ny;
	izone_1st = t.izone_1st;
	for (i=0; i<nzone; i++)
	{
		pbzone_dead[i] = t.pbzone_dead[i];
		for (j=0; j<nzone; j++)
			ppbzone[i][j] = t.ppbzone[i][j];
	}
	for (i=0; i<ns_comb; i++)
		pizone_in[i] = t.pizone_in[i];
	for (i=0; i<nzone; i++)
	{
		phflow_zone_ad[i] = t.phflow_zone_ad[i];
		for (j=0; j<GS_NE; j++)
		{
			ppeflow_zone[i][j] = t.ppeflow_zone[i][j];
			ppdafflow_zone[i][j] = t.ppdafflow_zone[i][j];
		}
		pashflow_zone[i] = t.pashflow_zone[i];
		phhvflow_zone[i] = t.phhvflow_zone[i];
		pgasemis[i] = t.pgasemis[i];
		pkag_zone[i] = t.pkag_zone[i];
		pkat4g_zone[i] = t.pkat4g_zone[i];
		pkap_zone[i] = t.pkap_zone[i];
		pkat4p_zone[i] = t.pkat4p_zone[i];
		phl_zone[i] = t.phl_zone[i];
		phlrad_zone[i] = t.phlrad_zone[i];
		phlconv_zone[i] = t.phlconv_zone[i];
		ptemp_zone[i] = t.ptemp_zone[i];
		ptime_zone[i] = t.ptime_zone[i];
		pvol_zone[i] = t.pvol_zone[i];
		pareaw_zone[i] = t.pareaw_zone[i];
		pdiamh_zone[i] = t.pdiamh_zone[i];
		ptwall_zone[i] = t.ptwall_zone[i];
		phconv_zone[i] = t.phconv_zone[i];
		pqwall_zone[i] = t.pqwall_zone[i];
		pstm_zone[i] = t.pstm_zone[i];
	}
	for (i=0; i<RF_NRXN; i++)
		pmrxn_char[i] = t.pmrxn_char[i];
	qexit = t.qexit;
	for (i=0; i<ns_watr; i++)
		pqwater[i] = t.pqwater[i];
	for (i=0; i<nsh+1; i++)
		pqwall[i] = t.pqwall[i];
	//not include mesh, radiation equation, and material stream pointers
}

CRadiantFurnace& CRadiantFurnace::operator=(const CRadiantFurnace& t)
{
	//does not include hard-wired variables
	if (this==&t)
		return *this;
	int i, j;
	//user input variables
	bmodified = t.bmodified;
	ifurn_type = t.ifurn_type;
	iprop_zone = t.iprop_zone;
	nzone_burn = t.nzone_burn;
	nzone_nose = t.nzone_nose;
	nzone_uppr = t.nzone_uppr;
	nx = t.nx;
	nz = t.nz;
	ns_comb = t.ns_comb;
	ns_watr = t.ns_watr;
	nsh = t.nsh;
	npsd = t.npsd;
	nck = t.nck;
	for (i=0; i<nzone; i++)
		pncell_zone[i]= t.pncell_zone[i];
	pres = t.pres;
	xdepth = t.xdepth;
	xhop_bot0 = t.xhop_bot0;
	xhop_bot1 = t.xhop_bot1;
	xnose_frt = t.xnose_frt;
	xnose_tip = t.xnose_tip;
	yheight = t.yheight;
	yhop_knuckle = t.yhop_knuckle;
	ynose_bot = t.ynose_bot;
	ynose_tip = t.ynose_tip;
	zwidth = t.zwidth;
	for (i=0; i<nsh+1; i++)
	{
		pemis_sh[i] = t.pemis_sh[i];
		prwall_sh[i] = t.prwall_sh[i];
		ptback_sh[i] = t.ptback_sh[i];
	}
	for (i=0; i<nzone; i++)
	{
		pemis_zone[i] = t.pemis_zone[i];
		prwall_zone[i] = t.prwall_zone[i];
		ptback_zone[i] = t.ptback_zone[i];
		pyzone[i] = t.pyzone[i];
	}
	for (i=0; i<ns_comb; i++)
		pyinlet[i] = t.pyinlet[i];
	for (i=0; i<nsh; i++)
		psh[i] = t.psh[i];
	for (i=0; i<npsd; i++)
		ppsd[i] = t.ppsd[i];
	for (i=0; i<nck; i++)
		pck[i] = t.pck[i];
	//calculated variables
	nzone = t.nzone;
	ny = t.ny;
	izone_1st = t.izone_1st;
	for (i=0; i<nzone; i++)
	{
		pbzone_dead[i] = t.pbzone_dead[i];
		for (j=0; j<nzone; j++)
			ppbzone[i][j] = t.ppbzone[i][j];
	}
	for (i=0; i<ns_comb; i++)
		pizone_in[i] = t.pizone_in[i];
	for (i=0; i<nzone; i++)
	{
		phflow_zone_ad[i] = t.phflow_zone_ad[i];
		for (j=0; j<GS_NE; j++)
		{
			ppeflow_zone[i][j] = t.ppeflow_zone[i][j];
			ppdafflow_zone[i][j] = t.ppdafflow_zone[i][j];
		}
		pashflow_zone[i] = t.pashflow_zone[i];
		phhvflow_zone[i] = t.phhvflow_zone[i];
		pgasemis[i] = t.pgasemis[i];
		pkag_zone[i] = t.pkag_zone[i];
		pkat4g_zone[i] = t.pkat4g_zone[i];
		pkap_zone[i] = t.pkap_zone[i];
		pkat4p_zone[i] = t.pkat4p_zone[i];
		phl_zone[i] = t.phl_zone[i];
		phlrad_zone[i] = t.phlrad_zone[i];
		phlconv_zone[i] = t.phlconv_zone[i];
		ptemp_zone[i] = t.ptemp_zone[i];
		ptime_zone[i] = t.ptime_zone[i];
		pvol_zone[i] = t.pvol_zone[i];
		pareaw_zone[i] = t.pareaw_zone[i];
		pdiamh_zone[i] = t.pdiamh_zone[i];
		ptwall_zone[i] = t.ptwall_zone[i];
		phconv_zone[i] = t.phconv_zone[i];
		pqwall_zone[i] = t.pqwall_zone[i];
		pstm_zone[i] = t.pstm_zone[i];
	}
	for (i=0; i<RF_NRXN; i++)
		pmrxn_char[i] = t.pmrxn_char[i];
	qexit = t.qexit;
	for (i=0; i<ns_watr; i++)
		pqwater[i] = t.pqwater[i];
	for (i=0; i<nsh+1; i++)
		pqwall[i] = t.pqwall[i];
	//not include mesh, radiation equation, and material stream pointers
	return *this;
}

void CRadiantFurnace::Write(FILE* pf)
{
	int i;
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	//user inputs and array size variables
	fwrite(&bmodified,sizeof(bool),1,pf);
	fwrite(&ifurn_type,sizeof(int),1,pf);
	fwrite(&iprop_zone,sizeof(int),1,pf);
	fwrite(&nzone_burn,sizeof(int),1,pf);
	fwrite(&nzone_nose,sizeof(int),1,pf);
	fwrite(&nzone_uppr,sizeof(int),1,pf);
	fwrite(&nzone,sizeof(int),1,pf);
	fwrite(&izone_1st,sizeof(int),1,pf);
	fwrite(&nx,sizeof(int),1,pf);
	fwrite(&ny,sizeof(int),1,pf);
	fwrite(&nz,sizeof(int),1,pf);
	fwrite(&ns_comb,sizeof(int),1,pf);
	fwrite(&ns_watr,sizeof(int),1,pf);
	fwrite(&nsh,sizeof(int),1,pf);
	fwrite(&npsd,sizeof(int),1,pf);
	fwrite(&nck,sizeof(int),1,pf);
	fwrite(pncell_zone,sizeof(int),nzone,pf);
	fwrite(&pres,sizeof(double),1,pf);
	fwrite(&xdepth,sizeof(double),1,pf);
	fwrite(&xhop_bot0,sizeof(double),1,pf);
	fwrite(&xhop_bot1,sizeof(double),1,pf);
	fwrite(&xnose_frt,sizeof(double),1,pf);
	fwrite(&xnose_tip,sizeof(double),1,pf);
	fwrite(&yheight,sizeof(double),1,pf);
	fwrite(&yhop_knuckle,sizeof(double),1,pf);
	fwrite(&ynose_bot,sizeof(double),1,pf);
	fwrite(&ynose_tip,sizeof(double),1,pf);
	fwrite(&zwidth,sizeof(double),1,pf);
	fwrite(pemis_sh,sizeof(double),nsh+1,pf);
	fwrite(prwall_sh,sizeof(double),nsh+1,pf);
	fwrite(ptback_sh,sizeof(double),nsh+1,pf);
	fwrite(pemis_zone,sizeof(double),nzone,pf);
	fwrite(prwall_zone,sizeof(double),nzone,pf);
	fwrite(ptback_zone,sizeof(double),nzone,pf);
	fwrite(pyzone,sizeof(double),nzone,pf);
	fwrite(pyinlet,sizeof(double),ns_comb,pf);
	//calculated variables
	fwrite(pbzone_dead,sizeof(bool),nzone,pf);
	for(i=0; i<nzone; i++)
		fwrite(ppbzone[i],sizeof(bool),nzone,pf);
	fwrite(pizone_in,sizeof(int),ns_comb,pf);
	fwrite(phflow_zone_ad,sizeof(double),nzone,pf);
	for (i=0; i<nzone; i++)
	{
		fwrite(ppeflow_zone[i],sizeof(double),GS_NE,pf);
		fwrite(ppdafflow_zone[i],sizeof(double),GS_NE,pf);
	}
	fwrite(pashflow_zone,sizeof(double),nzone,pf);
	fwrite(phhvflow_zone,sizeof(double),nzone,pf);
	fwrite(pgasemis,sizeof(double),nzone,pf);
	fwrite(pkag_zone,sizeof(double),nzone,pf);
	fwrite(pkat4g_zone,sizeof(double),nzone,pf);
	fwrite(pkap_zone,sizeof(double),nzone,pf);
	fwrite(pkat4p_zone,sizeof(double),nzone,pf);
	fwrite(phl_zone,sizeof(double),nzone,pf);
	fwrite(phlrad_zone,sizeof(double),nzone,pf);
	fwrite(phlconv_zone,sizeof(double),nzone,pf);
	fwrite(ptemp_zone,sizeof(double),nzone,pf);
	fwrite(ptime_zone,sizeof(double),nzone,pf);
	fwrite(pvol_zone,sizeof(double),nzone,pf);
	fwrite(pareaw_zone,sizeof(double),nzone,pf);
	fwrite(pdiamh_zone,sizeof(double),nzone,pf);
	fwrite(ptwall_zone,sizeof(double),nzone,pf);
	fwrite(phconv_zone,sizeof(double),nzone,pf);
	fwrite(pmrxn_char,sizeof(double),RF_NRXN,pf);
	fwrite(&qexit,sizeof(double),1,pf);
	fwrite(pqwater,sizeof(double),ns_watr,pf);
	fwrite(pqwall,sizeof(double),nsh+1,pf);
	fwrite(pqwall_zone,sizeof(double),nzone,pf);
	for (i=0; i<nsh; i++)
		psh[i].Write(pf);
	for (i=0; i<npsd; i++)
		ppsd[i].Write(pf);
	for (i=0; i<nck; i++)
		pck[i].Write(pf);
	for (i=0; i<nzone; i++)
		pstm_zone[i].Write(pf);
	mesh.Write(pf);
	//inlet streams
	for (i=0; i<ns_comb; i++)
		pstm_comb[i]->Write(pf);
}

void CRadiantFurnace::Read(FILE* pf)
{
	int i;
	int iversion;
	CMaterialStream* ps;
	fread(&iversion,sizeof(int),1,pf);
	//user inputs and array size variables
	fread(&bmodified,sizeof(bool),1,pf);
	fread(&ifurn_type,sizeof(int),1,pf);
	fread(&iprop_zone,sizeof(int),1,pf);
	fread(&nzone_burn,sizeof(int),1,pf);
	fread(&nzone_nose,sizeof(int),1,pf);
	fread(&nzone_uppr,sizeof(int),1,pf);
	fread(&nzone,sizeof(int),1,pf);
	fread(&izone_1st,sizeof(int),1,pf);
	fread(&nx,sizeof(int),1,pf);
	fread(&ny,sizeof(int),1,pf);
	fread(&nz,sizeof(int),1,pf);
	fread(&ns_comb,sizeof(int),1,pf);
	fread(&ns_watr,sizeof(int),1,pf);
	fread(&nsh,sizeof(int),1,pf);
	fread(&npsd,sizeof(int),1,pf);
	fread(&nck,sizeof(int),1,pf);
	fread(pncell_zone,sizeof(int),nzone,pf);
	fread(&pres,sizeof(double),1,pf);
	fread(&xdepth,sizeof(double),1,pf);
	fread(&xhop_bot0,sizeof(double),1,pf);
	fread(&xhop_bot1,sizeof(double),1,pf);
	fread(&xnose_frt,sizeof(double),1,pf);
	fread(&xnose_tip,sizeof(double),1,pf);
	fread(&yheight,sizeof(double),1,pf);
	fread(&yhop_knuckle,sizeof(double),1,pf);
	fread(&ynose_bot,sizeof(double),1,pf);
	fread(&ynose_tip,sizeof(double),1,pf);
	fread(&zwidth,sizeof(double),1,pf);
	fread(pemis_sh,sizeof(double),nsh+1,pf);
	fread(prwall_sh,sizeof(double),nsh+1,pf);
	fread(ptback_sh,sizeof(double),nsh+1,pf);
	fread(pemis_zone,sizeof(double),nzone,pf);
	fread(prwall_zone,sizeof(double),nzone,pf);
	fread(ptback_zone,sizeof(double),nzone,pf);
	fread(pyzone,sizeof(double),nzone,pf);
	fread(pyinlet,sizeof(double),ns_comb,pf);
	//calculated variables
	fread(pbzone_dead,sizeof(bool),nzone,pf);
	for(i=0; i<nzone; i++)
		fread(ppbzone[i],sizeof(bool),nzone,pf);
	fread(pizone_in,sizeof(int),ns_comb,pf);
	fread(phflow_zone_ad,sizeof(double),nzone,pf);
	for (i=0; i<nzone; i++)
	{
		fread(ppeflow_zone[i],sizeof(double),GS_NE,pf);
		fread(ppdafflow_zone[i],sizeof(double),GS_NE,pf);
	}
	fread(pashflow_zone,sizeof(double),nzone,pf);
	fread(phhvflow_zone,sizeof(double),nzone,pf);
	fread(pgasemis,sizeof(double),nzone,pf);
	fread(pkag_zone,sizeof(double),nzone,pf);
	fread(pkat4g_zone,sizeof(double),nzone,pf);
	fread(pkap_zone,sizeof(double),nzone,pf);
	fread(pkat4p_zone,sizeof(double),nzone,pf);
	fread(phl_zone,sizeof(double),nzone,pf);
	fread(phlrad_zone,sizeof(double),nzone,pf);
	fread(phlconv_zone,sizeof(double),nzone,pf);
	fread(ptemp_zone,sizeof(double),nzone,pf);
	fread(ptime_zone,sizeof(double),nzone,pf);
	fread(pvol_zone,sizeof(double),nzone,pf);
	fread(pareaw_zone,sizeof(double),nzone,pf);
	fread(pdiamh_zone,sizeof(double),nzone,pf);
	fread(ptwall_zone,sizeof(double),nzone,pf);
	fread(phconv_zone,sizeof(double),nzone,pf);
	fread(pmrxn_char,sizeof(double),RF_NRXN,pf);
	fread(&qexit,sizeof(double),1,pf);
	fread(pqwater,sizeof(double),ns_watr,pf);
	fread(pqwall,sizeof(double),nsh+1,pf);
	fread(pqwall_zone,sizeof(double),nzone,pf);
	for (i=0; i<nsh; i++)
		psh[i].Read(pf);
	for (i=0; i<npsd; i++)
		ppsd[i].Read(pf);
	for (i=0; i<nck; i++)
		pck[i].Read(pf);
	for (i=0; i<nzone; i++)
		pstm_zone[i].Read(pf);
	mesh.Read(pf);
	//inlet streams
	for (i=0; i<ns_comb; i++)
	{
		ps = new CMaterialStream();
		ps->Read(pf);
		pstm_comb[i] = ps;
	}
	if (mesh.IsGeometrySet())
	{
		//set up radiation eqn since it is not serialized
		doeqn.pmesh = &mesh;
		doeqn.AllocateMemory();		//this call will delete existing memory first
		doeqn.Init();
	}
}

int CRadiantFurnace::CheckInputs()
{
	//return nonzero if error
	int i, j;
	double sum;
	CSolidFuel* pcoal;
	//check if input streams are set, need at least 1 fuel/oxidizer inlet and 1 water/steam inlet
	if (ns_comb<1)
	{
		printf("number of fuel/oxidizer inlets should be at least 1!");
		return 1;
	}
	if (ns_watr<1)
	{
		printf("number of water/steam inlets should be at least 1!");
		return 1;
	}
	//check individual inlet streams
	for (i=0; i<ns_comb; i++)
	{
		//check if NULL
		if (pstm_comb[i]==NULL)
		{
			printf("Inlet fuel/oxidizer stream not set yet!");
			return 1;
		}
		//check if inlet elevation is valid
		if (pyinlet[i]<0 || pyinlet[i]>yheight)
		{
			printf("Fuel/oxidizer inlet elevation invalid!");
			return 3;
		}
		//check if solid property valid
		if (pstm_comb[i]->bsolid)
		{
			pcoal = &pstm_comb[i]->coal;
			sum = 0;
			for (j=0; j<6; j++)
				sum += pcoal->fcmass[0][j];
			if (pcoal->fcmass[0][8]>sum)
			{
				printf("Volatile content of solid fuel is invalid!");
				return 4;
			}
			sum += pcoal->fcmass[0][6] + pcoal->fcmass[0][7];
			if (sum<0.999 || sum>1.001)
			{
				printf("mass fractions of solid fuel are not sum up to 1!");
				return 4;
			}
		}
	}
	//check number of superheater panels
	for (i=0; i<nsh; i++)
	{
		if (nz<(psh[i].npanel+1)*2)
		{
			printf("Number of cells between two superheater panels is less than 2!");
			return 3;
		}
	}
	return 0;
}

int CRadiantFurnace::Init()
{
	if (CheckInputs())
		return 1;
	//return nonzero if error
	//reset mesh and discrete ordinate equation members
	bool bdead;					//flag for dead zone (no fuel and oxidizer flow)
	int i, j, k, l, m, n;		//temporary integer variables
	int nx1, ny1, nz1;			//nx+1, ny+1, nz+1
	int nxface, nyface, nzface;	//number of faces in x, y, and z direction
	int npanel;					//number of panels for each superheater group, limited to nz-1
	int pizone[RF_NY];			//zone index for each cell plane along height
	int pizpanel[RF_NZ];		//face index in z direction for superheater panels
	double x, y, z;				//local variable for vertex's x, y, z
	double dx, dy, dz;			//local variable for delta x, delta y, delta z
	double pyface[RF_NY+1];		//face y
	double pxstart[RF_NY+1];	//starting point of x at each y face
	double pxend[RF_NY+1];		//ending point of x at each y face
	//first remove current mesh if any
	mesh.DeleteArrays();
	//total number of zones with 1 hopper zone only
	nzone = nzone_burn + nzone_nose + nzone_uppr + 1;
	//calculate number of cell planes in y direction
	ny = 0;
	for (i=0; i<nzone; i++)
		ny += pncell_zone[i];
	nx1 = nx+1;
	ny1 = ny+1;
	nz1 = nz+1;
	nxface = nx1*ny*nz;
	nyface = ny1*nz*nx;
	nzface = nz1*nx*ny;
	mesh.ncell = nx*ny*nz;
	mesh.nvertex = nx1*ny1*nz1;
	mesh.nface =  nxface + nyface + nzface;			//excluding shadow superheater boundary wall faces
	mesh.nbface = 2*(ny*nz + nz*nx + nx*ny);		//excluding shadow superheater boundary wall faces
	mesh.AllocateArrays();
	int* pibface = mesh.pibface;
	CVertex* pvertex = mesh.pvertex;
	CQuadFace* pface = mesh.pface;
	CHexCell* pcell = mesh.pcell;
	CQuadFace* pface_sh;			//superheater faces
	int* piface_sh_conv;			//indices of interior faces converted to superheater faces
	//calculate pyzone[] at hopper knuckle, nose bottom, nose tip and roof
	pyzone[0] = yhop_knuckle;
	pyzone[nzone_burn] = ynose_bot;
	pyzone[nzone_burn+nzone_nose] = ynose_tip;
	pyzone[nzone-1] = yheight;
	//calculate pyface[], pxstart[], and pxend[]
	pyface[0] = 0;
	pxstart[0] = xhop_bot0;
	pxend[0] = xhop_bot1;
	//hopper zone
	j = 1;
	dy = pyzone[0]/pncell_zone[0];
	for (i=0; i<pncell_zone[0]; i++)
	{
		pizone[j-1] = 0;
		y = dy*(i+1);
		pyface[j] = y;
		pxstart[j] = xhop_bot0*(yhop_knuckle-y)/yhop_knuckle;
		pxend[j] = xhop_bot1+(xdepth-xhop_bot1)*y/yhop_knuckle;
		j++;
	}
	//burner zone
	for (k=1; k<nzone_burn+1; k++)
	{
		dy = (pyzone[k]-pyzone[k-1])/pncell_zone[k];
		for (i=0; i<pncell_zone[k]; i++)
		{
			pizone[j-1] = k;
			y = pyzone[k-1] + dy*(i+1);
			pyface[j] = y;
			pxstart[j] = 0;
			pxend[j] = xdepth;
			j++;
		}
	}
	//nose zone
	for (k=nzone_burn+1; k<nzone_burn+nzone_nose+1; k++)
	{
		dy = (pyzone[k]-pyzone[k-1])/pncell_zone[k];
		for (i=0; i<pncell_zone[k]; i++)
		{
			pizone[j-1] = k;
			y = pyzone[k-1] + dy*(i+1);
			pyface[j] = y;
			pxstart[j] = xnose_frt*(y-ynose_bot)/(ynose_tip-ynose_bot);
			pxend[j] = xnose_tip+(y-ynose_tip)*(xnose_tip-xdepth)/(ynose_tip-ynose_bot);
			j++;
		}
	}
	//upper zone
	for (k=nzone_burn+nzone_nose+1; k<nzone; k++)
	{
		dy = (pyzone[k]-pyzone[k-1])/pncell_zone[k];
		for (i=0; i<pncell_zone[k]; i++)
		{
			pizone[j-1] = k;
			y = pyzone[k-1] + dy*(i+1);
			pyface[j] = y;
			pxstart[j] = xnose_frt;
			pxend[j] = xnose_tip;
			j++;
		}
	}
	//set up vertex array, calculate the vertex's x, y, z
	dz = zwidth/nz;
	for (j=0; j<ny1; j++)
	{
		y = pyface[j];
		dx = (pxend[j]-pxstart[j])/nx;
		for (i=0; i<nx1; i++)
		{
			x = dx*i + pxstart[j];
			for (k=0; k<nz1; k++)
			{
				z = dz*k;
				n = i*ny1*nz1 + j*nz1 + k;
				pvertex[n].x[0] = x;
				pvertex[n].x[1] = y;
				pvertex[n].x[2] = z;
			}
		}
	}
	//set up face array
	n = 0;		//face index
	m = 0;		//for boundary face
	//x faces first, followed by y faces and z faces
	for (i=0; i<nx1; i++)
	{
		for (j=0; j<ny; j++)
		{
			for (k=0; k<nz; k++)
			{
				pface[n].ivertex[0] = i*ny1*nz1 + j*nz1 + k;
				pface[n].ivertex[1] = i*ny1*nz1 + j*nz1 + k + 1;
				pface[n].ivertex[2] = i*ny1*nz1 + (j+1)*nz1 + k + 1;
				pface[n].ivertex[3] = i*ny1*nz1 + (j+1)*nz1 + k;
				if (i==0)
				{
					pibface[m] = n;
					m++;
					pface[n].izone = pizone[j];
					pface[n].itype = 1;
					pface[n].icell[0] = i*ny*nz + j*nz + k;
					pface[n].icell[1] = -1;
					pface[n].pbfp = new CBoundaryFaceProperty;
					pface[n].pbfp->iwall = 0;
					pface[n].pbfp->iwater = 0;
				}
				else
				{
					if (i==nx)
					{
						pibface[m] = n;
						m++;
						pface[n].izone = pizone[j];
						pface[n].itype = 1;
						pface[n].icell[0] = (i-1)*ny*nz + j*nz + k;
						pface[n].icell[1] = -1;
						pface[n].pbfp = new CBoundaryFaceProperty;
						pface[n].pbfp->iwall = 0;
						pface[n].pbfp->iwater = 0;
						if (ifurn_type==0)		//rear wall exit
						{
							if ((pyface[j]+pyface[j+1])/2>ynose_tip)	//outlet cell
							{
								pface[n].itype = 3;
								pface[n].pbfp->iwall = nsh+1;			//use outlet property
								pface[n].pbfp->iwater = 0;				//heat loss to outlet still assigned to enclosure wall
							}
						}
					}
					else
					{
						pface[n].itype = 0;
						pface[n].icell[0] = (i-1)*ny*nz + j*nz + k;
						pface[n].icell[1] = i*ny*nz + j*nz + k;
					}
				}
				n++;
			}
		}
	}
	//y faces second
	for (j=0; j<ny1; j++)
	{
		for (k=0; k<nz; k++)
		{
			for (i=0; i<nx; i++)
			{
				pface[n].ivertex[0] = i*ny1*nz1 + j*nz1 + k;
				pface[n].ivertex[1] = (i+1)*ny1*nz1 + j*nz1 + k;
				pface[n].ivertex[2] = (i+1)*ny1*nz1 + j*nz1 + k + 1;
				pface[n].ivertex[3] = i*ny1*nz1 + j*nz1 + k + 1;
				if (j==0)
				{
					pibface[m] = n;
					m++;
					pface[n].izone = 0;
					pface[n].itype = 1;
					pface[n].icell[0] = i*ny*nz + j*nz + k;
					pface[n].icell[1] = -1;
					pface[n].pbfp = new CBoundaryFaceProperty;
					pface[n].pbfp->iwall = 0;
					pface[n].pbfp->iwater = 0;
				}
				else
				{
					if (j==ny)
					{
						pibface[m] = n;
						m++;
						pface[n].izone = nzone-1;
						pface[n].itype = 1;
						pface[n].icell[0] = i*ny*nz + (j-1)*nz + k;
						pface[n].icell[1] = -1;
						pface[n].pbfp = new CBoundaryFaceProperty;
						pface[n].pbfp->iwall = 0;
						pface[n].pbfp->iwater = 0;
						if (ifurn_type==1)		//top exit, tower unit
						{
							pface[n].itype = 3;
							pface[n].pbfp->iwall = nsh+1;			//use outlet property
							pface[n].pbfp->iwater = 0;				//heat loss to outlet still assigned to enclosure wall
						}
					}
					else
					{
						pface[n].itype = 0;
						pface[n].icell[0] = i*ny*nz + (j-1)*nz + k;
						pface[n].icell[1] = i*ny*nz + j*nz + k;
					}
				}
				n++;
			}
		}
	}
	//z faces third
	for (k=0; k<nz1; k++)
	{
		for (i=0; i<nx; i++)
		{
			for (j=0; j<ny; j++)
			{
				pface[n].ivertex[0] = i*ny1*nz1 + j*nz1 + k;
				pface[n].ivertex[1] = i*ny1*nz1 + (j+1)*nz1 + k;
				pface[n].ivertex[2] = (i+1)*ny1*nz1 + (j+1)*nz1 + k;
				pface[n].ivertex[3] = (i+1)*ny1*nz1 + j*nz1 + k;
				if (k==0)
				{
					pibface[m] = n;
					m++;
					pface[n].izone = pizone[j];
					pface[n].itype = 1;
					pface[n].icell[0] = i*ny*nz + j*nz + k;
					pface[n].icell[1] = -1;
					pface[n].pbfp = new CBoundaryFaceProperty;
					pface[n].pbfp->iwall = 0;
					pface[n].pbfp->iwater = 0;
				}
				else
				{
					if (k==nz)
					{
						pibface[m] = n;
						m++;
						pface[n].izone = pizone[j];
						pface[n].itype = 1;
						pface[n].icell[0] = i*ny*nz + j*nz + k - 1;
						pface[n].icell[1] = -1;
						pface[n].pbfp = new CBoundaryFaceProperty;
						pface[n].pbfp->iwall = 0;
						pface[n].pbfp->iwater = 0;
					}
					else
					{
						pface[n].itype = 0;
						pface[n].icell[0] = i*ny*nz + j*nz + k - 1;
						pface[n].icell[1] = i*ny*nz + j*nz + k;
					}
				}
				n++;
			}
		}
	}
	//debug
	if (n!=mesh.nface)
		printf("number of face does not match!");
	if (m!=mesh.nbface)
		printf("number of boundary face does not match!");
	//set up cell array
	n = 0;		//cell index
	for (i=0; i<nx; i++)
	{
		for (j=0; j<ny; j++)
		{
			for (k=0; k<nz; k++)
			{
				pcell[n].pivertex[0] = i*ny1*nz1 + j*nz1 + k;
				pcell[n].pivertex[1] = (i+1)*ny1*nz1 + j*nz1 + k;
				pcell[n].pivertex[2] = i*ny1*nz1 + (j+1)*nz1 + k;
				pcell[n].pivertex[3] = (i+1)*ny1*nz1 + (j+1)*nz1 + k;
				pcell[n].pivertex[4] = i*ny1*nz1 + j*nz1 + k + 1;
				pcell[n].pivertex[5] = (i+1)*ny1*nz1 + j*nz1 + k + 1;
				pcell[n].pivertex[6] = i*ny1*nz1 + (j+1)*nz1 + k + 1;
				pcell[n].pivertex[7] = (i+1)*ny1*nz1 + (j+1)*nz1 + k + 1;
				pcell[n].izone = pizone[j];
				n++;
			}
		}
	}
	//debug
	if (n!=mesh.ncell)
		printf("Incorrect number of cells");
	//update cell geometry (centroid and volumne)
	mesh.UpdateCellGeometry();
	//calculate volume of each zone
	mesh.CalcZoneVolumes(nzone,pvol_zone);
	//map superheater wall cells
	pface_sh = new CQuadFace [nx*ny*nz];
	piface_sh_conv = new int [nx*ny*nz];
	z = zwidth/nz;		//cell width
	n = 0;
	for (k=0; k<nsh; k++)
	{
		//requres psh[k].npanel <= nz-1 from inputs, and no overlaping of multiple super heaters
		npanel = psh[k].npanel;
		if (npanel>nz-1)
		{
			printf("Number of superheater panels exceeded maximum allowed. Fewer panels are modeled");
			npanel = nz-1;
		}
		dz = zwidth/(npanel+1);
		for (m=0; m<psh[k].npanel; m++)
			pizpanel[m] = (int)(dz*(m+1)/z+0.4999);
		for (i=0; i<nx; i++)
		{
			for (j=0; j<ny; j++)
			{
				l = i*ny*nz + j*nz;
				if (psh[k].IsInsidePolygon(pcell[l].x[0],pcell[l].x[1]))
				{
					for (m=0; m<psh[k].npanel; m++)
					{
						l = nxface + nyface + pizpanel[m]*nx*ny + i*ny + j;
						piface_sh_conv[n] = l;
						pface[l].itype = 1;
						pface[l].izone = pizone[j];		//added on 11/6/2015 otherwise izone=-1
						pface_sh[n].CopyPartialData(&pface[l]);		//copy itype, izone and ivertex[]
						pface_sh[n].icell[0] = pface[l].icell[1];
						pface_sh[n].icell[1] = -1;
						pface_sh[n].pbfp = new CBoundaryFaceProperty;
						pface_sh[n].pbfp->iwall = psh[k].iwall;
						pface_sh[n].pbfp->iwater = psh[k].iwater;
						pface[l].icell[1] = -1;
						pface[l].pbfp = new CBoundaryFaceProperty;
						pface[l].pbfp->iwall = psh[k].iwall;
						pface[l].pbfp->iwater = psh[k].iwater;
						n++;
					}
				}
			}
		}
	}
	//add pface_sh to face and boundary face index array, memory copy is handled by AddBoundaryFaces
	//pface_sh is deleted in AddBoundaryFaces
	mesh.AddBoundaryFaces(n,pface_sh,piface_sh_conv);
	pface = mesh.pface;
	pibface = mesh.pibface;
	//update the face and neigboring cell array for each cell
	//add interior face first followed by boundary faces
	for (i=0; i<mesh.nface; i++)
	{
		if (!pface[i].itype)	//interior face
		{
			j = pface[i].icell[0];
			k = pface[i].icell[1];
			pcell[j].AddFaceAndNeighborCell(i,k);
			pcell[k].AddFaceAndNeighborCell(i,j);
		}
	}
	for (i=0; i<mesh.nface; i++)
	{
		if (pface[i].itype)		//boundary face
		{
			j = pface[i].icell[0];
			pcell[j].AddFaceAndNeighborCell(i,-1);
		}
	}
	//assign face boundary properties
	for (i=0; i<mesh.nbface; i++)
	{
		j = pibface[i];
		k = pface[j].pbfp->iwall;
		if (k==0)		//enclosure wall
		{
			if (iprop_zone)			//properties based on zone
				m = pface[j].izone;
			else					//uniform properties, use first element in array
				m = 0;
			pface[j].pbfp->emis = pemis_zone[m];
			pface[j].pbfp->rwall = prwall_zone[m];
			pface[j].pbfp->tback = ptback_zone[m];
			//initially set wall temperature to backside temperature
			pface[j].pbfp->temp = ptback_zone[m];
		}
		else			//superheater or outlet wall
		{
			pface[j].pbfp->emis = pemis_sh[k-1];
			pface[j].pbfp->rwall = prwall_sh[k-1];
			pface[j].pbfp->tback = ptback_sh[k-1];
			//for outlet face, wall temperature always equal to backside temperature (not updated)
			pface[j].pbfp->temp = ptback_sh[k-1];
		}
	}
	//update geometry data
	mesh.UpdateFaceGeomerty();
	mesh.CalcMBL();
	//set up radiation eqn
	doeqn.pmesh = &mesh;
	doeqn.AllocateMemory();		//this call will delete existing memory first
	doeqn.Init();
	//assign pizone_in[]
	for (i=0; i<ns_comb; i++)
	{
		if (pyinlet[i]<=pyzone[0])		//below hopper knuckle, allowing inlet flow in hopper
			pizone_in[i] = 0;
		else
		{
			if (pyinlet[i]>=yheight)
				pizone_in[i] = nzone-1;
			else
			{
				for (j=1; j<nzone; j++)
				{
					if (pyinlet[i]<pyzone[j])
					{
						pizone_in[i] = j;
						break;
					}
				}
			}
		}
	}
	//assign ppbzone[i][j] where i is the current zone id and j is another zone id
	//ppbzone[i][j] is true if material from zone j flow to zone i
	for (i=0; i<nzone; i++)
	{
		for (j=0; j<nzone; j++)
			ppbzone[i][j] = false;
	}
	//set to true for its own zone if there is an inlet flow to the zone
	//allows inlet to hopper
	for (i=0; i<nzone; i++)
	{
		for (j=0; j<ns_comb; j++)
		{
			if (pizone_in[j]==i)
			{
				ppbzone[i][i] = true;
				break;
			}
		}
	}
	//check individual zone based on itself and the zone below it
	pbzone_dead[0] = !ppbzone[0][0];
	for (i=1; i<nzone; i++)
	{
		bdead = true;
		for (k=0; k<nzone; k++)
		{
			ppbzone[i][k] = ppbzone[i][k] || ppbzone[i-1][k];
			if (ppbzone[i][k])
				bdead = false;
		}
		if (!bdead)
			ppbzone[i][i] = true;
		pbzone_dead[i] = bdead;
	}
	//find the lowest middle zone with flow
	izone_1st = -1;
	for (i=0; i<nzone; i++)
	{
		for (k=0; k<nzone; k++)
		{
			if (ppbzone[i][k])
			{
				izone_1st = i;
				break;
			}
		}
		if (izone_1st!=-1)
			break;
	}
	if (izone_1st==-1)
	{
		printf("Cannot find first flow zone!");
		return 1;
	}
	//calculate geometry related parameters for convective heat transfer
	mesh.CalcZoneWallArea(nzone, pareaw_zone);
	//calculate hydralic diameter and initialize phconv_zone and phlconv_zone
	for (i=0; i<nzone; i++)
	{
		pdiamh_zone[i] = 4*pvol_zone[i]/pareaw_zone[i];
		phconv_zone[i] = 0;
		phlconv_zone[i] = 0;
	}
	return 0;
}

void CRadiantFurnace::UpdateZoneResidenceTime()
{
	//calculate residence time based on zone exit gas mass flow and total gas mass in the zone
	int i;
	for (i=0; i<nzone; i++)
	{
		if (pbzone_dead[i])
			ptime_zone[i] = 0;
		else
		{
			if (!pstm_zone[i].bgas)
				ptime_zone[i] = 0;
			else
			{
				ptime_zone[i] = pvol_zone[i]*pstm_zone[i].gas.den/pstm_zone[i].gas.mflow;
				if (ifurn_type==0 && i>=nzone-nzone_uppr)		//upper zone, use 1/3 of the residence time
					ptime_zone[i] /= 3;
			}
		}
	}
}

void CRadiantFurnace::TrackSolidParticlesAndUpdateZoneStreams(double* pkap)
{
	//this function will update the unreacted solid phase in each zone and update the particle radiation absorption coefficient
	//it will also set the stream in each zone before equilibrium calculation
	//it does not update the hflow of stream in each zone
	//under-relax the mass source term based on solid daf flow
	//first calculate the residence time in each zone
	UpdateZoneResidenceTime();
	int i, j, k;
	int izone;				//zone index
	int izone_in;			//index of the zone an inlet stream enters
	int nbin;				//number of size bins of particle size distribution
	double sum;				//used to calculate sum
	double nflow;			//particle number flow rate
	double ndenp;			//particle number density
	double temp;			//temperature
	double tim;				//residence time
	double dt;				//time step for particle tracking
	double yo2;				//mole fraction of O2 in gas phase
	double yh2o;			//mole fraction of H2O in gas phase
	double yco2;			//mole fraction of CO2 in gas phase
	double rchar;			//total char mass reacted per particle with the time step [kg]
	double dp0;				//initial particle diameter
	double dp_0;			//particle diameter after swelling during devolatilization and before char oxidation
	double dp;				//current particle diameter
	double mp0;				//initial particle mass
	double mp_0;			//particle mass after devolatilization and before char reaction
	double mp;				//current particle mass
	double mp_daf;			//mass of dry ash free part of particle
	double mp_ash;			//mass of ash which does not change in the process
	double denp0;			//initial particle density
	double areap;			//external area of particle
	double loi;				//loss on ignition or mass fraction of daf coal in current particle mass
	double qa;				//particle radiation asorption efficiency
	double qa_coal;			//coal particle radiation asorption efficiency
	double qa_ash;			//ash particle radiation asorption efficiency
	double dafflow_zone[RF_NZONE][GS_NE] = {0};	//C H O N S Cl mass flow of particle phase leaving each zone
	double hhvflow_zone[RF_NZONE] = {0};		//daf HHV flow of particle phase leaving each zone
	double rates_char[RF_NRXN];					//reaction rates by O2, H2O, and CO2
	CMaterialStream* ps;
	CSolidFuel* pcoal;
	CIdealGasMixture* pgas;
	CPSD* psd;
	CCoalKinetics* pk;
	CParticleRadiationProperty* pprp = CParticleRadiationProperty::GetInstance();
	//clear current particle radiation absorption coefficient in each zone
	for (i=0; i<nzone; i++)
		pkap[i] = 0;
	//clear pmrxn_char
	for (i=0; i<RF_NRXN; i++)
		pmrxn_char[i] = 0;
	//track each inlet stream for each particle size
	for (i=0; i<ns_comb; i++)
	{
		ps = pstm_comb[i];
		if (ps->bsolid)
		{
			pcoal = &ps->coal;
			//the solid phase of a inlet stream should always have a size distribution and a set of kinetic data
			//for uniform distribution, nbin could be 1, default ipsd=0, ick=0 in CSolidFuel constructor
			psd = &ppsd[pcoal->ipsd];
			pk = &pck[pcoal->ick];
			izone_in = pizone_in[i];
			nbin = psd->nbin;
			for (j=0; j<nbin; j++)
			{
				dp0 = psd->adp[j];
				denp0 = psd->denp;
				mp0 = CT_PI/6*dp0*dp0*dp0*denp0;
				mp_ash = mp0*pcoal->fcmass[0][7];
				//number flow rate does not change if assuming no fragmentation
				nflow = pcoal->mflow*psd->mf[j]/mp0;
				//assuming moisture vaporization and devolatilization complete immediately upon entering the zone
				dp = dp0*(1+pk->fswell);		//diameter swelling after devolatilization
				mp = mp0*(1-pcoal->fcmass[0][6]-pcoal->fcmass[0][8]);
				//daf before char reactions start
				mp_daf = mp - mp_ash;
				areap = CT_PI*dp*dp;
				mp_0 = mp;
				dp_0 = dp;
				//now start char particle reactions along furnace zones
				for (izone=izone_in; izone<nzone; izone++)
				{
					temp = ptemp_zone[izone];
					tim = ptime_zone[izone];
					if (izone==izone_in)		//residence time could be a fraction of entire volume of the zone
					{
						if (izone==0)
							tim *= (pyzone[0] - pyinlet[i])/pyzone[0];
						else
							tim *= (pyzone[izone]-pyinlet[i])/(pyzone[izone]-pyzone[izone-1]);
					}
					ndenp = nflow*tim/pvol_zone[izone];
					if (mp_daf>0)		//need calculate reaction rates with multiple time steps
					{
						dt = tim/(double)RF_NSTEP;
						//apply effectiveness factor to consider non-uniform reactant gas bulk mole fractions
						yo2 = ef_rxn[izone]*pstm_zone[izone].gas.GetMoleFraction("O2");
						yh2o = ef_rxn[izone]*pstm_zone[izone].gas.GetMoleFraction("H2O");
						yco2 = ef_rxn[izone]*pstm_zone[izone].gas.GetMoleFraction("CO2");
						for (k=0; k<RF_NSTEP; k++)
						{
							pk->CalcCharReactionRate(yo2, yh2o, yco2, temp, pres, dp, rates_char);
							sum = rates_char[0] + rates_char[1] + rates_char[2];
							rchar = sum*dt;
							if (rchar>mp_daf)
							{
								rchar = mp_daf;
								mp_daf = 0;
							}
							else
								mp_daf -= rchar;
							//assign to pmrxn_char[]
							if (sum>0)
							{
								pmrxn_char[0] += rchar*rates_char[0]/sum*nflow;
								pmrxn_char[1] += rchar*rates_char[1]/sum*nflow;
								pmrxn_char[2] += rchar*rates_char[2]/sum*nflow;
							}
							mp = mp_ash + mp_daf;
							loi = mp_daf/mp;
							//update diameter and density based on Hurt's burning model parameter alpha
							//denp/denp_0 = (mp/mp_0)^alpha or dp/dp_0 = (mp/mp_0)^[(1-alpha)/3]
							dp = dp_0*pow(mp/mp_0,(1-pck->alpha)/3);
							areap = CT_PI*dp*dp;
							if (mp_daf>0)		//still have daf char left
							{
								qa_coal = pprp->GetCoalAbsorptionEfficiency(dp,temp);
								qa_ash = pprp->GetAshAbsorptionEfficiency(dp,temp);
								qa = loi*qa_coal + (1-loi)*qa_ash;
								pkap[izone] += ndenp*areap/4*qa/RF_NSTEP;
							}
							else	//ash only
							{
								qa = pprp->GetAshAbsorptionEfficiency(dp,temp);
								pkap[izone] += ndenp*areap/4*qa/RF_NSTEP*(RF_NSTEP-k);
								break;
							}
						}	//end of for k loop
						if (mp_daf>0)	//still has daf char left after traveling through the zone
						{
							//calculate solid mass leaving the zone
							for (k=0; k<6; k++)		//C H O N S Cl
								dafflow_zone[izone][k] += mp_daf*nflow*pcoal->fcmass[2][k];
							//HHV
							hhvflow_zone[izone] += mp_daf*nflow*pcoal->hhv[2];
						}
					}
					else	//ash only, only update pkap
					{
						qa = pprp->GetAshAbsorptionEfficiency(dp,temp);
						pkap[izone] += ndenp*areap/4*qa;
					}	//end of if mp_daf>0
				}	//end of izone loop
			}	//end of for j loop
		}	//end of if bsolid loop
	}	//end of for i loop
	//prepare stream for each zone for equilibrium calculation
	for (izone=izone_1st; izone<nzone; izone++)
	{
		//apply effectiveness factor for particle radiation absorption coefficient
		pkap[izone] *= ef_kap;
		//under-relax ppdafflow_zone[][] and hhvflow_zone[];
		for (i=0; i<6; i++)
			ppdafflow_zone[izone][i] = (1-urf_daf)*ppdafflow_zone[izone][i] + urf_daf*dafflow_zone[izone][i];
		phhvflow_zone[izone] = (1-urf_daf)*phhvflow_zone[izone] + urf_daf*hhvflow_zone[izone];
		ps = &pstm_zone[izone];
		pgas = &ps->gas;
		pcoal = &ps->coal;
		//first update the coal
		sum = 0;
		for (i=0; i<6; i++)
			sum += ppdafflow_zone[izone][i];
		sum += pashflow_zone[izone];
		if (sum>0)		//has solid flow
		{
			for (i=0; i<6; i++)		//C H O N S Cl
				pcoal->fcmass[0][i] = ppdafflow_zone[izone][i]/sum;
			pcoal->fcmass[0][6] = 0;	//no moisture
			pcoal->fcmass[0][7] = pashflow_zone[izone]/sum;		//ash
			pcoal->fcmass[0][8] = 0;	//no volatile
			pcoal->hhv[0] = phhvflow_zone[izone]/sum;
			pcoal->ibasis = 0;
			pcoal->temp = 298.15;
			pcoal->mflow = sum;
			pcoal->CalcAllProperties();
			ps->bsolid = true;
		}
		else
		{
			pcoal->mflow = 0;
			ps->bsolid = false;
		}
		//now update gas phase
		for (j=0; j<GS_NE; j++)
			pgas->eflow[j] = ppeflow_zone[izone][j];
		if (ps->bsolid)
		{
			for (j=0; j<GS_NE; j++)
				pgas->eflow[j] -= pcoal->eflow[0][j];
		}
		InitGasStreamForEquilibriumCalculation(pgas);
		ps->bgas = pgas->mflow>0;
		if (!ps->bgas)
			printf("Gas phase is empty in zone %d\n", izone);
		//no liquid phase
		ps->oil.mflow = 0;
		ps->bliquid = false;
		ps->mflow = pgas->mflow + pcoal->mflow;
	}
}

void CRadiantFurnace::InitGasStreamForEquilibriumCalculation(CIdealGasMixture* pigm)
{
	//initialize the gas phase array for equilibrium calculation
	int i;
	std::vector<CIdealGas*> vsl;		//vector of pointers to any gas that could be formed from the elements in pstm
	pigm->ictp = 1;					//always solve at const T, P
	pigm->InitArraysFromEflow();	//make sure the element is set
	if (pigm->mflow==0)
		return;
	//assign species
	pigm->ns = 0;
	pigm->GetValidSpeciesList(vsl);
	for (i=0; i<nspecies; i++)
		pigm->AddValidSpecies(ispecies[i],vsl);		//add if ispecies[i] is in the valid species list vsl
	if (!pigm->AreAllElementsInGasArray())
		printf("Not all elements in gas species array!");
}

int CRadiantFurnace::CalcMaterialStreamEquilibrium(CMaterialStream* pstm)
{
	//solve gas phase equilibrium combined with solid phase given eflow[] in the gas phase, composition in solid phase and total enthalpy flow
	CIdealGasMixture* pgm = &pstm->gas;
	if (!pstm->bliquid && !pstm->bsolid)	//gas phase only
	{
		pgm->h = pstm->hflow/pstm->mflow;
		pgm->ictp = 0;		//constant T and H
		pgm->Solve();
		pstm->temp = pgm->t;
		pstm->UpdateProperties();
		return 0;
	}
	double h = pstm->hflow;
	int nit;			//number of iterations
	double t0;			//lower end temperature
	double t1;			//higher end temperature
	double tn;			//new guess
	double fn;			//new (hguess - h)
	double err;			//relative error of t
	tn = 800;			//set initial guess to 800 K
	t0 = tn;
	t1 = tn;
	pgm->t = tn;
	pgm->Solve();
	pstm->temp = tn;
	pstm->UpdateProperties();
	fn = pstm->hflow-h;
	if (fn>0)
	{
		nit = 0;
		while (fn>0)		//lower end temperature too high
		{
			if (nit>=5)
			{
				return -1;	//not going to converge 800*0.8^5=262 K
			}
			nit++;
			t1 = t0;
			t0 = t0*0.8;
			pgm->t = t0;
			pgm->Solve();
			pstm->temp = t0;
			pstm->UpdateProperties();
			fn = pstm->hflow-h;
		}
	}
	else
	{
		nit = 0;
		while (fn<0)		//higher end temperature too low
		{
			if (nit>=5)
			{
				printf("Enthalpy too high!");
				return -1;	//not going to converge (800*1.45^7=5127 K)
			}
			nit++;
			t0 = t1;
			t1 = t1*1.45;
			pgm->t = t1;
			pgm->Solve();
			pstm->temp = t1;
			pstm->UpdateProperties();
			fn = pstm->hflow-h;
		}
	}
	//bisection method
	nit = 0;
	do
	{
		nit++;
		tn = (t0+t1)/2;
		pgm->t = tn;
		pgm->Solve();
		pstm->temp = tn;
		pstm->UpdateProperties();
		fn = pstm->hflow-h;
		err = fabs((tn-t0)/tn);
		if (fn<0)
		{
			t0 = tn;
		}
		else
		{
			t1 = tn;
		}
	} while (err>0.00002 && nit<20);
	//this code usually converges within 15 iterations
	if (nit>=20)
	{
		printf("Equilibrium reactor does not converge after 20 iterations");
		return -1;		//not converged within 20 iterations
	}
	return 0;
}

double CRadiantFurnace::CalcGasEmissivity(CMaterialStream* pstm)
{
	int i, j;
	double gasemis;
	CIdealGasMixture* pgm = &pstm->gas;
	CIdealGas* pgas;
	char radgases[6][5]={"CO2","H2O","CH4","CO","O2","N2"};
	double f[6] = {0};
	double den_soot = 1950;		//soot density
	double c1 = 0.1;			//mass fraction of carbon that forms soot if there is no oxidation, Glacier is 0.1
	double c2;					//fraction of soot formed that is not oxidated
	double bc;					//mass fraction of C element in gas mixture
	double er;					//gas phase ER
	double er_cr = 1.0;			//ER at which soot exist in gas phase, Glacier is 1.0
	double er_inf = 2*er_cr;	//ER above which soot does not oxidize, Glacier is also 2*er_cr
	double fvsoot;				//soot volume fraction
	double temp = pstm->temp;
	double pressure = pstm->pres;
	double mbl = mesh.mbl;
	//calculate soot volume fraction
	er = pgm->er;
	if (er>er_inf)
		c2 = 1;
	else
	{
		if (er<er_cr)
			c2 = 0;
		else
			c2 = (er-er_cr)/(er_inf-er_cr);
	}
	bc = pgm->eflow[0]*CT_AWC/pgm->mflow;
	fvsoot = c1*c2*bc*pgm->den/den_soot;
	for (i=0; i<5; i++)
	{
		pgas = CIdealGas::GetGasByName(radgases[i]);
		for (j=0; j<pgm->ns; j++)
		{
			if (pgm->pgas[j]==pgas)
			{
				f[i] = pgm->fsp_mole[j];
				break;
			}
		}
	}
	//assume all other species is equivalent to N2
	f[5] = 1-f[0]-f[1]-f[2]-f[3]-f[4];
	//call external function
	gasemissivity_(&temp, &mbl, &pressure, f, &fvsoot, &gasemis);
	return gasemis;
}

double CRadiantFurnace::CalcGasAbsorptionCoefficient(double ge)
{
	//apply effectiveness factor for gas radiation absorption coefficient
	return -log(1-ge)/mesh.mbl*ef_kag;
}

int CRadiantFurnace::CalculateOutStreams()
{
	//flue gas is the 1st in the outstreams, water wall water/steam is the 2nd in the outstreams
	int i, j;
	int nbface = mesh.nbface;
	int* pibface = mesh.pibface;
	double qface;		//heat loss to a boundary face
	CQuadFace* pf;
	CQuadFace* pface = mesh.pface;
	qexit = 0;
	for (i=0; i<ns_watr; i++)
		pqwater[i] = 0;
	for (i=0; i<nsh+1; i++)
		pqwall[i] = 0;
	for (i=0; i<RF_NZONE; i++)
		pqwall_zone[i] = 0;
	CMaterialStream* pin = &pstm_zone[nzone-1];
	pstm_flue->CopyProperties(pin);
	//calculate water streams
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		qface = (pf->pbfp->qnet + pf->pbfp->qconv)*pf->area;		//qnet and qconv are generally negative
		if (pf->itype==3)		//outlet faces
		{
			qexit -= qface;
			//add exit heat to water wall stream, assuming screen tubes are part of water wall
			pqwater[pf->pbfp->iwater] -= qface;
		}
		else
		{
			pqwall_zone[pf->izone] -= qface;
			pqwater[pf->pbfp->iwater] -= qface;
			pqwall[pf->pbfp->iwall] -= qface;
		}
	}
	return 0;
}

void CRadiantFurnace::CalcZoneElementAshAndEnthalpyFlows()
{
	int i, j;
	int izone;
	CMaterialStream* ps;
	for (izone=0; izone<nzone; izone++)
	{
		//initialize to zero
		pashflow_zone[izone] = 0;
		phflow_zone_ad[izone] = 0;
		for (j=0; j<GS_NE; j++)
			ppeflow_zone[izone][j] = 0;
		if (pbzone_dead[izone])		//skip dead zones
			continue;
		for (i=0; i<ns_comb; i++)
		{
			ps = pstm_comb[i];
			if (!ppbzone[izone][pizone_in[i]])		//inlet stream does not flow to the current zone
				continue;
			phflow_zone_ad[izone] += ps->hflow;
			if (ps->bgas)
			{
				for (j=0; j<GS_NE; j++)
					ppeflow_zone[izone][j] += ps->gas.eflow[j];
			}
			if (ps->bliquid)
			{
				for (j=0; j<GS_NE; j++)
					ppeflow_zone[izone][j] += ps->oil.eflow[j];
			}
			if (ps->bsolid)
			{
				//add molar flow of all element including moisture
				for (j=0; j<GS_NE; j++)
					ppeflow_zone[izone][j] += ps->coal.eflow[0][j];
				//add ash mass flow (or molar flow at ash MW of 1
				pashflow_zone[izone] += ps->coal.mflow*ps->coal.fcmass[0][7];
			}
		}
	}
}

void CRadiantFurnace::CalcFullyReactedAdiabaticZoneTemperatures()
{
	//also initialize ppdafflow_zone[][] and phhvflow_zone[]
	int i;
	int izone;
	CMaterialStream* ps;
	CIdealGasMixture* pgas;
	CSolidFuel* pcoal;
	for (izone=izone_1st; izone<nzone; izone++)
	{
		phhvflow_zone[izone] = 0;
		for (i=0; i<GS_NE; i++)
			ppdafflow_zone[izone][i] = 0;
		ps = &pstm_zone[izone];
		pgas = &ps->gas;
		pcoal = &ps->coal;
		//first update the ash only solid phase
		if (pashflow_zone[izone]>0)		//has ash flow
		{
			for (i=0; i<7; i++)			//C H O N S Cl moisture
				pcoal->fcmass[0][i] = 0;
			pcoal->fcmass[0][7] = 1;	//ash
			pcoal->fcmass[0][8] = 0;	//no volatile
			pcoal->hhv[0] = 0;
			pcoal->ibasis = 0;
			pcoal->temp = 298.15;
			pcoal->mflow = pashflow_zone[izone];
			pcoal->CalcAllProperties();
			ps->bsolid = true;
		}
		else
		{
			pcoal->mflow = 0;
			ps->bsolid = false;
		}
		//now update gas phase
		for (i=0; i<GS_NE; i++)
			pgas->eflow[i] = ppeflow_zone[izone][i];
		InitGasStreamForEquilibriumCalculation(pgas);
		ps->bgas = pgas->mflow>0;
		if (!ps->bgas)
			printf("Gas phase is empty in zone %d\n", izone);
		//no liquid phase
		ps->oil.mflow = 0;
		ps->bliquid = false;
		ps->mflow = pgas->mflow + pcoal->mflow;
		ps->hflow = phflow_zone_ad[izone];
		//now perform equilibrium calculation
		CalcMaterialStreamEquilibrium(ps);
		//also set ptemp_zone[] to adiabatic temperature
		ptemp_zone[izone] = ps->temp;
	}
}

void CRadiantFurnace::CalcConvectiveHeatTransferCoefficient()
{
	//calculate convective heat transfer coefficient based on fully developed turbulent pipe flow (Bird Page 399)
	//and fully developed laminar flow with uniform surface temperature (Incropera and DeWitt, Fundamentals of Heat and Mass Transfer (5th edition, pp 486-487)
	int i;
	double dy_zone;				//zone height
	double area;				//area
	double vel;					//velocity
	double visc_bulk;			//viscosity using bulk gas temperature
	double visc_wall;			//viscosity using wall temperature
	double cond_bulk;			//thermal conductivity
	double den_bulk;			//density
	double Re;					//Re number
	double Pr;					//Pr number
	double Nu;					//Nu number
	mesh.CalcZoneWallAreaTemperatureProduct(nzone, ptwall_zone);
	for (i=0; i<nzone; i++)		//calculate the average wall temperature and film temperature
		ptwall_zone[i] /= pareaw_zone[i];
	for (i=0; i<nzone; i++)
	{
		if (pbzone_dead[i])
			phconv_zone[i] = 0;
		else
		{
			if (!pstm_zone[i].bgas)
				phconv_zone[i] = 0;
			else
			{
				if (i==0)
					dy_zone = pyzone[i];
				else
					dy_zone = pyzone[i] - pyzone[i-1];
				area = pvol_zone[i]/dy_zone;
				vel = pstm_zone[i].gas.mflow/pstm_zone[i].gas.den/area;
				visc_bulk = pstm_zone[i].gas.CalcViscosity(ptemp_zone[i]);
				visc_wall = pstm_zone[i].gas.CalcViscosity(ptwall_zone[i]);
				cond_bulk = pstm_zone[i].gas.CalcConductivity(ptemp_zone[i]);
				den_bulk = pstm_zone[i].gas.den;
				Re = pdiamh_zone[i]*vel*den_bulk/visc_bulk;
				Pr = pstm_zone[i].gas.CalcCp(ptemp_zone[i])*visc_bulk/cond_bulk;
				if (Re>1.0e4)
					Nu = 0.026*pow(Re, 0.8)*pow(Pr, 0.33333)*pow(visc_bulk/visc_wall, 0.14);
				else
					Nu = 3.66;
				phconv_zone[i] = Nu*cond_bulk/pdiamh_zone[i];
				if (ifurn_type==0 && i>=nzone-nzone_uppr)		//upper zone, use 1/3 heat transfer coefficient
					phconv_zone[i] /= 3;
			}
		}
	}
}

int CRadiantFurnace::Solve()
{
	//return non-zero integer if failed
	bool bconv;						//flag for convergence
	int i, j;
	int iter = 0;					//iteration count
	double heatloss;				//accumulated heat loss
	double dtemp0;					//Zone 0 temperature adjustment
	double dtemp;					//tempreature difference
	double dtemp_max;				//maximum temperature difference
	double ptemp_prev[RF_NZONE];	//previous gas temperature array, used for checking convergence
	double hl_zone[RF_NZONE];		//current radiation and convection heat loss, used for under-relaxation
	double kap_zone[RF_NZONE];		//current particle radiation absorption coefficient, used for under-relaxation
	//skip Init() if no user inputs are modified related to the size of arrays
	if (bmodified)
	{
		if (Init())					//check inputs and reset all mesh and radiation equation members
			return 1;
	}
	//calculate module pressure as the lowest of the inlet streams
	pres = pstm_comb[0]->pres;
	for (i=1; i<ns_comb; i++)
	{
		if (pstm_comb[i]->pres<pres)
			pres = pstm_comb[i]->pres;
	}
	for (i=0; i<nzone; i++)
	{
		pstm_zone[i].pres = pres;
		pstm_zone[i].gas.p = pres;
	}
	//calculate ppeflow_zone[][], pashflow_zone[], phflow_zone_ad[] which are fixed
	CalcZoneElementAshAndEnthalpyFlows();
	if (bmodified)			//need initial guess before iterations
	{
		//use fully reacted adiabatic case as initial guess
		//this function also initializes ppdafflow_zone[][], phhvflow_zone[], and ptemp_zone[]
		CalcFullyReactedAdiabaticZoneTemperatures();
		//adjust temperature by a factor of 0.85
		for (i=izone_1st; i<nzone; i++)
		{
			if (ptemp_zone[i]>1000)
				ptemp_zone[i] *= 0.85;
		}
		//traking particles to calculate kap_zone[] and char reaction rates and update under-relaxed ppdafflow_zone[][], phhvflow_zone[] and material streams in each zone
		TrackSolidParticlesAndUpdateZoneStreams(kap_zone);
		//assign to pkap_zone[] and pkat4p_zone[] as initial guess
		for (i=izone_1st; i<nzone; i++)
		{
			pkap_zone[i] = kap_zone[i];
			pkat4p_zone[i] = kap_zone[i]*pow(ptemp_zone[i],4);
		}
		//calculate 1st non-dead zone
		pgasemis[izone_1st] = CalcGasEmissivity(&pstm_zone[izone_1st]);
		pkag_zone[izone_1st] = CalcGasAbsorptionCoefficient(pgasemis[izone_1st]);
		pkat4g_zone[izone_1st] = pkag_zone[izone_1st]*pow(ptemp_zone[izone_1st],4);
		for (i=0; i<nzone; i++)
		{
			if (i==izone_1st)		//already calculated
				continue;
			if (pbzone_dead[i])		//estimate properties for dead zones
			{
				//assign 1st non-dead zone properties to all dead zones
				ptemp_zone[i] = ptemp_zone[izone_1st];
				pgasemis[i] = pgasemis[izone_1st];
				pkag_zone[i] = pkag_zone[izone_1st];
				pkat4g_zone[i] = pkat4g_zone[izone_1st];
				//ignore particle radiation in dead zones
				pkap_zone[i] = 0;
				pkat4p_zone[i] = 0;
			}
			else		//other reaction zones
			{
				pgasemis[i] = CalcGasEmissivity(&pstm_zone[i]);
				pkag_zone[i] = CalcGasAbsorptionCoefficient(pgasemis[i]);
				pkat4g_zone[i] = pkag_zone[i]*pow(ptemp_zone[i],4);
			}
		}
		//update cell radiation properties in each zone
		doeqn.UpdateRadiationProperties(pkag_zone,pkap_zone,pkat4g_zone,pkat4p_zone);
		//initially phconv_zone[] is 0
		doeqn.UpdateConvectionBoundaryFaceProperties(phconv_zone, ptemp_zone);
		doeqn.SolvePDE();
		//calculate initial heat loss in each zone phl_zone[]
		doeqn.GetZoneRadiationHeatLoss(nzone, phlrad_zone);
		if (iconvect)		//consider convective heat transfer term
			doeqn.GetZoneConvectionHeatLoss(nzone, phlconv_zone);
		for (i=0; i<nzone; i++)
			phl_zone[i] = phlrad_zone[i] + phlconv_zone[i];
	}
	//start interations now
	do
	{
		//calculate char reaction rate and particle absorption coefficient in each zone and update material stream in each zone
		TrackSolidParticlesAndUpdateZoneStreams(kap_zone);
		//under-relax pkap_zone[] and pkat4p_zone[]
		for (i=izone_1st; i<nzone; i++)
		{
			pkap_zone[i] = (1-urf_kap)*pkap_zone[i] + urf_kap*kap_zone[i];
			pkat4p_zone[i] = (1-urf_kap)*pkat4p_zone[i] + urf_kap*kap_zone[i]*pow(ptemp_zone[i],4);
		}
		//save zone temperatures for convergence check
		for (i=0; i<nzone; i++)
			ptemp_prev[i] = ptemp_zone[i];
		doeqn.GetZoneRadiationHeatLoss(nzone, phlrad_zone);
		if (iconvect)		//consider convective heat transfer
			doeqn.GetZoneConvectionHeatLoss(nzone, phlconv_zone);
		for (i=0; i<nzone; i++)
			hl_zone[i] = phlrad_zone[i] + phlconv_zone[i];
		//under-relax phl_zone[] for non-dead zones
		for (i=0; i<nzone; i++)
		{
			if (pbzone_dead[i])
				phl_zone[i] = hl_zone[i];
			else
				phl_zone[i] = (1-urf_hl)*phl_zone[i] + urf_hl*hl_zone[i];
		}
		//update properties for each zone
		for (i=0; i<nzone; i++)
		{
			if (!pbzone_dead[i])	//hopper zone could be the first zone
			{
				//update hflow for each stream
				heatloss = 0;
				for (j=0; j<nzone; j++)
				{
					if (ppbzone[i][j])	//zone j is upstream of zone i
						heatloss += phl_zone[j];
				}
				pstm_zone[i].hflow = phflow_zone_ad[i] - heatloss;
				CalcMaterialStreamEquilibrium(&pstm_zone[i]);
				//under-relax temperature
				ptemp_zone[i] = urf_temp*pstm_zone[i].temp + (1-urf_temp)*ptemp_zone[i];
				pstm_zone[i].temp = ptemp_zone[i];
				pgasemis[i] = CalcGasEmissivity(&pstm_zone[i]);
				pkag_zone[i] = CalcGasAbsorptionCoefficient(pgasemis[i]);
				pkat4g_zone[i] = pkag_zone[i]*pow(ptemp_zone[i],4);
			}
		}
		if (izone_1st>0)	//update dead zone properties,
		{
			//calculate properties for dead zone (Zone 0) and use them for all dead zones
			//temperature determined by setting heat loss to zero
			pstm_zone[0].CopyProperties(&pstm_zone[izone_1st]);
			pstm_zone[0].temp = ptemp_zone[0];
			pstm_zone[0].UpdateProperties();
			heatloss = 0;
			for (i=0; i<nzone; i++)
			{
				if (pbzone_dead[i])
					heatloss += phl_zone[i];
			}
			//use the mass inside the dead zones to adjust the temperature
			dtemp0 = heatloss/pvol_zone[0]/pstm_zone[0].gas.den/pstm_zone[0].CalcCp();
			pstm_zone[0].temp -= dtemp0;
			if (pstm_zone[0].temp < 273)		//avoid extremely low temperature
				pstm_zone[0].temp = 273;
			ptemp_zone[0] = 0.5*pstm_zone[0].temp + 0.5*ptemp_zone[0];
			pstm_zone[0].temp = ptemp_zone[0];
			pgasemis[0] = CalcGasEmissivity(&pstm_zone[0]);
			pkag_zone[0] = CalcGasAbsorptionCoefficient(pgasemis[0]);
			pkat4g_zone[0] = pkag_zone[0]*pow(ptemp_zone[0],4);
			pkap_zone[0] = 0;		//currently set to zero for dead zones
			pkat4p_zone[0] = 0;		//currently set to zero for dead zones
			//Zone 0 properties to all dead zones
			for (i=1; i<nzone; i++)
			{
				if (pbzone_dead[i])
				{
					pstm_zone[i].CopyProperties(&pstm_zone[0]);
					ptemp_zone[i] = ptemp_zone[0];
					pgasemis[i] = pgasemis[0];
					pkag_zone[i] = pkag_zone[0];
					pkat4g_zone[i] = pkat4g_zone[0];
					pkap_zone[i] = pkap_zone[0];
					pkat4p_zone[i] = pkat4p_zone[0];
				}
			}
		}
		doeqn.UpdateRadiationProperties(pkag_zone,pkap_zone,pkat4g_zone,pkat4p_zone);
		if (iconvect)
		{
			CalcConvectiveHeatTransferCoefficient();
			doeqn.UpdateConvectionBoundaryFaceProperties(phconv_zone, ptemp_zone);
		}
		doeqn.SolvePDE();
		dtemp_max = 0;
		for (i=izone_1st; i<nzone; i++)	//does not check dead zone
		{
			if (!pbzone_dead[i])
			{
				dtemp = fabs(ptemp_zone[i]-ptemp_prev[i]);
				if (dtemp>dtemp_max)
				{
					j = i;
					dtemp_max = dtemp;
				}
			}
		}
		bconv = dtemp_max<0.01;
		//alse make dead zone convergence
		if (izone_1st)
			bconv = bconv && fabs(ptemp_zone[0]-ptemp_prev[0])<0.1;
		iter++;
		printf("Iteration %d: max temp change = %lg\n",iter,dtemp_max);
	} while (!bconv && iter<200);
	if (iter>=200)
		printf("Boiler model not converged after 200 iterations.\n");
	else
		printf("Boiler model solved in %d iterations.\n",iter);
	doeqn.CalcBoundaryNetHeatFlux();
	//update tube side water stream based on heat absorption
	CalculateOutStreams();
	bmodified = false;
	return 0;
}

void CRadiantFurnace::GetVertexCount(int* pcount)
{
	//used for export fieldview file
	int i, j;
	int ncell = mesh.ncell;
	int nvertex = mesh.nvertex;
	CHexCell* pcell = mesh.pcell;
	for (i=0; i<nvertex; i++)
		pcount[i] = 0;
	for (i=0; i<ncell; i++)
	{
		for (j=0; j<8; j++)
			pcount[pcell[i].pivertex[j]]++;
	}
}

void CRadiantFurnace::GetVertexValues(int* pcount, double* pvar_zone, float* pvar_vertex)
{
	//used for export fieldview file
	int i, j;
	int izone;
	int ncell = mesh.ncell;
	int nvertex = mesh.nvertex;
	CHexCell* pcell = mesh.pcell;
	for (i=0; i<nvertex; i++)
		pvar_vertex[i] = 0;
	for (i=0; i<ncell; i++)
	{
		izone = pcell[i].izone;
		for (j=0; j<8; j++)
			pvar_vertex[pcell[i].pivertex[j]] += (float)pvar_zone[izone];
	}
	for (i=0; i<nvertex; i++)
		pvar_vertex[i] /= (float)pcount[i];
}

void CRadiantFurnace::GetCellValues(double* pvar_zone, float* pvar_cell)
{
	//used for export Paraview vtk file
	int i;
	int izone;
	int ncell = mesh.ncell;
	CHexCell* pcell = mesh.pcell;
	for (i=0; i<ncell; i++)
	{
		izone = pcell[i].izone;
		pvar_cell[i] = (float)pvar_zone[izone];
	}
}

void CRadiantFurnace::WriteVtkFiles(char* filename)
{
	int i, j, ii;
	//add extension if needed
	std::string filename_bnd;
	std::string filename_vtk = filename;
	i = filename_vtk.rfind('.');
	if (i>-1)
		filename_vtk = filename_vtk.substr(0,i);
	filename_bnd = filename_vtk + "_bnd.vtk";
	filename_vtk = filename_vtk + ".vtk";
	int nvertex = mesh.nvertex;
	int ncell = mesh.ncell;
	int nbface = mesh.nbface;
	int* pibface = mesh.pibface;
	CVertex* pvertex = mesh.pvertex;
	CQuadFace* pface = mesh.pface;
	CHexCell* pcell = mesh.pcell;
	CHexCell* pc;
	CIdealGas* pgas;
	CIdealGasMixture* pgm;
	int itmp;
	int ncell_vertex;
	int itmpv[8];
	int* pivertex;
	double pvar_zone[RF_NZONE];
	double* pxv;
	float ftmp;
	float* pvar = new float [ncell];
	FILE* fout;
	if ((fout = fopen(filename_vtk.c_str(),"w"))==NULL)
	{
		printf("Unable to open Paraview vtk file for write!");
		return;
	}
	fprintf(fout,"# vtk DataFile Version 2.0\n");
	fprintf(fout,"Radiant Furnace Model 3-D Results\n");
	fprintf(fout,"ASCII\n");
	fprintf(fout,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fout,"POINTS %d float\n",nvertex);
	for (i=0; i<nvertex; i++)
	{
		pxv = pvertex[i].x;
		fprintf(fout,"%e %e %e\n",pxv[0], pxv[1], pxv[2]);
	}
	itmp = 0;
	for (i=0; i<ncell; i++)
		itmp += pcell[i].nvertex;
	itmp += ncell;
	fprintf(fout,"CELLS %d %d\n",ncell,itmp);
	for (i=0; i<ncell; i++)
	{
		pc = &pcell[i];
		ncell_vertex = pc->nvertex;			//ncell_vertex should be 8
		fprintf(fout,"%d",ncell_vertex);
		itmpv[0] = (int)pc->pivertex[0];
		itmpv[1] = (int)pc->pivertex[1];
		itmpv[2] = (int)pc->pivertex[3];
		itmpv[3] = (int)pc->pivertex[2];
		itmpv[4] = (int)pc->pivertex[4];
		itmpv[5] = (int)pc->pivertex[5];
		itmpv[6] = (int)pc->pivertex[7];
		itmpv[7] = (int)pc->pivertex[6];
		for (ii=0; ii<ncell_vertex; ii++)
			fprintf(fout," %d", itmpv[ii]);
		fprintf(fout,"\n");
	}
	fprintf(fout,"CELL_TYPES %d\n",ncell);
	itmp = 12;			//VTK_HEXAHEDRON
	for (i=0; i<ncell; i++)
		fprintf(fout,"%d\n",itmp);
	fprintf(fout,"CELL_DATA %d\n",ncell);
	//zone id
	fprintf(fout,"SCALARS Zone_ID float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<ncell; i++)
	{
		ftmp = (float)pcell[i].izone + 1;
		fprintf(fout,"%g ",ftmp);
	}
	//temperature
	fprintf(fout,"\nSCALARS Temperature float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	GetCellValues(ptemp_zone, pvar);
	for (i=0; i<ncell; i++)
		fprintf(fout,"%g ", pvar[i]);
	//gas emissivity
	fprintf(fout,"\nSCALARS Gas_Emissivity float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	GetCellValues(pgasemis, pvar);
	for (i=0; i<ncell; i++)
		fprintf(fout,"%g ", pvar[i]);
	//CO
	fprintf(fout,"\nSCALARS CO float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	pgas = CIdealGas::GetGasByName("CO");
	for (i=0; i<nzone; i++)
	{
		pgm = &pstm_zone[i].gas;
		for (j=0; j<pgm->ns; j++)
		{
			if (pgm->pgas[j]==pgas)
			{
				pvar_zone[i] = pgm->fsp_mole[j];
				break;
			}
		}
	}
	GetCellValues(pvar_zone, pvar);
	for (i=0; i<ncell; i++)
		fprintf(fout,"%g ", pvar[i]);
	//O2
	fprintf(fout,"\nSCALARS O2 float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	pgas = CIdealGas::GetGasByName("O2");
	for (i=0; i<nzone; i++)
	{
		pgm = &pstm_zone[i].gas;
		for (j=0; j<pgm->ns; j++)
		{
			if (pgm->pgas[j]==pgas)
			{
				pvar_zone[i] = pgm->fsp_mole[j];
				break;
			}
		}
	}
	GetCellValues(pvar_zone, pvar);
	for (i=0; i<ncell; i++)
		fprintf(fout,"%g ", pvar[i]);
	fclose(fout);
	//write bournary file
	if ((fout = fopen(filename_bnd.c_str(),"w"))==NULL)
	{
		printf("Unable to open Paraview vtk file for write!");
		return;
	}
	fprintf(fout,"# vtk DataFile Version 2.0\n");
	fprintf(fout,"Radiant Furnace Model Boundary Results\n");
	fprintf(fout,"ASCII\n");
	fprintf(fout,"DATASET POLYDATA\n");
	fprintf(fout,"POINTS %d float\n",nvertex);
	for (i=0; i<nvertex; i++)
	{
		pxv = pvertex[i].x;
		fprintf(fout,"%e %e %e\n",pxv[0], pxv[1], pxv[2]);
	}
	itmp = nbface*5;
	fprintf(fout,"POLYGONS %d %d\n",nbface,itmp);
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pivertex = pface[j].ivertex;
		fprintf(fout,"4 %d %d %d %d\n",pivertex[0],pivertex[1],pivertex[2],pivertex[3]);
	}
	fprintf(fout,"CELL_DATA %d\n",nbface);
	//wall group ID
	fprintf(fout,"SCALARS Wall_ID float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		ftmp = (float)pface[j].pbfp->iwall;
		fprintf(fout,"%g ",ftmp);
	}
	//water stream ID
	fprintf(fout,"\nSCALARS Water_Stream_ID float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		ftmp = (float)pface[j].pbfp->iwater;
		fprintf(fout,"%g ",ftmp);
	}
	//wall temperature
	fprintf(fout,"\nSCALARS Wall_Temperature float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		ftmp = (float)pface[j].pbfp->temp;
		fprintf(fout,"%g ",ftmp);
	}
	//fluid temperature
	fprintf(fout,"\nSCALARS Fluid_Temperature float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		ftmp = (float)pface[j].pbfp->tback;
		fprintf(fout,"%g ",ftmp);
	}
	//incident radiation flux
	fprintf(fout,"\nSCALARS Incident_Heat_Flux float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		ftmp = (float)pface[j].pbfp->qinc;
		fprintf(fout,"%g ",ftmp);
	}
	//net radiation flux, reverse the sign
	fprintf(fout,"\nSCALARS Net_Heat_Flux float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		ftmp = (float)(-pface[j].pbfp->qnet);
		fprintf(fout,"%g ",ftmp);
	}
	//wall emissivity
	fprintf(fout,"\nSCALARS Wall_Emissivity float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		ftmp = (float)pface[j].pbfp->emis;
		fprintf(fout,"%g ",ftmp);
	}
	//wall resistance
	fprintf(fout,"\nSCALARS Wall_Resistance float 1\n");
	fprintf(fout,"LOOKUP_TABLE default\n");
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		ftmp = (float)pface[j].pbfp->rwall;
		fprintf(fout,"%g ",ftmp);
	}
	fclose(fout);
	//delete local pointer
	delete [] pvar;
}

void CRadiantFurnace::WriteResultDescriptionFile(char* filename)
{
	int i, ncount;
	FILE* fout;
	if ((fout = fopen(filename,"w"))==NULL)
	{
		printf("Unable to open result description file for write!");
		return;
	}
	std::vector<std::string> result_descs;
	PrepareResultDescriptionStrings(result_descs);
	ncount = result_descs.size();
	for (i=0; i<ncount; i++)
		fprintf(fout,"%s\n",result_descs[i].c_str());
	fclose(fout);
}

double CRadiantFurnace::CalcWaterWallArea()
{
	//currently consider the exit plane (screen tubes) as part of waterwall
	double tmp1, tmp2;
	double len_hop_fw = sqrt(yhop_knuckle*yhop_knuckle + xhop_bot0*xhop_bot0);
	tmp1 = xdepth - xhop_bot1;
	double len_hop_rw = sqrt(yhop_knuckle*yhop_knuckle + tmp1*tmp1);
	tmp1 = xdepth - xnose_tip;
	tmp2 = ynose_tip - ynose_bot;
	double len_nose_rw = sqrt(tmp1*tmp1 + tmp2*tmp2);
	double len_nose_fw = sqrt(xnose_frt*xnose_frt + tmp2*tmp2);
	//perimeter of side wall
	tmp1 = (xhop_bot1-xhop_bot0)+len_hop_fw+len_hop_rw+(ynose_bot-yhop_knuckle)*2+len_nose_fw+len_nose_rw+(yheight-ynose_tip)*2+(xnose_tip-xnose_frt);
	tmp1 *= zwidth;
	tmp2 = xdepth-xnose_tip;
	tmp2 = yheight*xdepth - (xhop_bot0+xdepth-xhop_bot1)*yhop_knuckle/2 - (tmp2+xnose_frt)*(yheight-ynose_tip) - (tmp2+xnose_frt)*(ynose_tip-ynose_bot)/2;
	tmp1 += 2*tmp2;
	return tmp1;
}

double CRadiantFurnace::CalcSuperHeaterWallArea()
{
	double area = 0;
	for (int i=0; i<nsh; i++)
		area += psh[i].CalcArea();
	return area;
}

void CRadiantFurnace::PrepareResults(std::vector<double>& results)
{
	int i;
	double tmp, tmp1, tmp2;
	results.clear();
	//water wall area
	tmp = CalcWaterWallArea();
	results.push_back(tmp);
	//superheater wall area
	tmp = CalcSuperHeaterWallArea();
	results.push_back(tmp);
	//modeled enclosure wall area from mesh
	tmp = mesh.area_encl;
	results.push_back(tmp);
	//modeled superheater wall area from mesh
	tmp = mesh.area_sh;
	results.push_back(tmp);
	//modeled exit area from mesh
	tmp = mesh.area_exit;
	results.push_back(tmp);
	//modeled total area from mesh
	tmp = mesh.area_total;
	results.push_back(tmp);
	//modeled total volume from mesh
	tmp = mesh.volume;
	results.push_back(tmp);
	//modeled mean beam length from mesh
	tmp = mesh.mbl;
	results.push_back(tmp);
	//calculate heat losses from zones
	tmp = 0;
	tmp1 = 0;
	tmp2 = 0;
	for (i=0; i<nzone; i++)
	{
		tmp += phl_zone[i];
		tmp1 += phlrad_zone[i];
		tmp2 += phlconv_zone[i];
	}
	//total heat loss
	results.push_back(tmp);
	//heat loss to enclosure wall and supper heaters
	for (i=0; i<nsh+1; i++)
		results.push_back(pqwall[i]);
	//heat loss to exit plane
	results.push_back(qexit);
	//heat loss to individual water/steam streams,  heat loss to exit plane is assigned to feed water stream
	for (i=0; i<ns_watr; i++)
		results.push_back(pqwater[i]);
	//total radiation heat loss
	results.push_back(tmp1);
	//radiation heat loss through gas phase
	results.push_back(-doeqn.qcellg);
	//radiation heat loss through particle phase
	results.push_back(-doeqn.qcellp);
	//total convective heat loss
	results.push_back(tmp2);
	//flue gas stream data
	//entire stream data
	//flue gas pressure
	results.push_back(pstm_flue->pres);
	//flow gas temperature
	results.push_back(pstm_flue->temp);
	//gas phase data
	//gas phase mass flow rate
	results.push_back(pstm_flue->gas.mflow);
	//gas phase species mass fractions
	for (i=0; i<pstm_flue->gas.ns; i++)
		results.push_back(pstm_flue->gas.fsp_mass[i]);
	//solid phas data
	//solid phase mass flow rate
	results.push_back(pstm_flue->coal.mflow);
	//solid phase mass fractions in the order of C, H, O, N, S, Cl ash
	for (i=0; i<6; i++)
		results.push_back(pstm_flue->coal.fcmass[0][i]);
	results.push_back(pstm_flue->coal.fcmass[0][7]);
	//additional detailed model data for each zone
	//gas/solid temperature in each zone
	for (i=0; i<nzone; i++)
		results.push_back(ptemp_zone[i]);
	//average wall temperature in each zone
	//calculate ptwall_zone in case iconvect==0
	mesh.CalcZoneWallAreaTemperatureProduct(nzone, ptwall_zone);
	for (i=0; i<nzone; i++)		//calculate the average wall temperature and film temperature
		ptwall_zone[i] /= pareaw_zone[i];
	for (i=0; i<nzone; i++)
		results.push_back(ptwall_zone[i]);
	//total heat loss in each zone
	for (i=0; i<nzone; i++)
		results.push_back(phl_zone[i]);
	//radiation heat loss in each zone
	for (i=0; i<nzone; i++)
		results.push_back(phlrad_zone[i]);
	//convection heat loss in each zone
	for (i=0; i<nzone; i++)
		results.push_back(phlconv_zone[i]);
	//gas emissivity in each zone
	for (i=0; i<nzone; i++)
		results.push_back(pgasemis[i]);
	//pkag in each zone
	for (i=0; i<nzone; i++)
		results.push_back(pkag_zone[i]);
	//pkap in each zone
	for (i=0; i<nzone; i++)
		results.push_back(pkap_zone[i]);
	//unburned carbon mass fraction (LOI) in solid phase of each non-dead zone
	for (i=izone_1st; i<nzone; i++)
	{
		if (pstm_zone[i].bsolid)
			results.push_back(1-pstm_zone[i].coal.fcmass[0][7]);
	}
	//avarage wall incident heat flux
	double pqinc_avg[RF_NZONE], pqinc_min[RF_NZONE], pqinc_max[RF_NZONE];
	mesh.CalcZoneWallIncidentFlux(nzone, pqinc_avg, pqinc_min, pqinc_max);
	for (i=0; i<nzone; i++)		//average
		results.push_back(pqinc_avg[i]);
	for (i=0; i<nzone; i++)		//min
		results.push_back(pqinc_min[i]);
	for (i=0; i<nzone; i++)		//max
		results.push_back(pqinc_max[i]);
	for (i=izone_1st; i<nzone; i++)		//max
		results.push_back(ptime_zone[i]);
	//mass flow of char reacted by O2, H2O, and CO2
	for (i=0; i<RF_NRXN; i++)
		results.push_back(pmrxn_char[i]);
}

void CRadiantFurnace::PrepareResultDescriptionStrings(std::vector<std::string>& results)
{
	int i;
	char cbuf[500];
	std::string str;
	results.clear();
	results.push_back("water wall area from specified geometry");
	results.push_back("superheater wall area from specified geometry");
	results.push_back("modeled enclosure wall area from mesh");
	results.push_back("modeled superheater wall area accounting for both sides from mesh");
	results.push_back("modeled exit area");
	results.push_back("modeled total boundary area");
	results.push_back("modeled total volume");
	results.push_back("mean beam length of radiant furnace");
	results.push_back("total heat loss");
	results.push_back("heat loss to enclosure wall");
	for (i=0; i<nsh; i++)
	{
		sprintf(cbuf,"heat loss to wall of superheater %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	results.push_back("heat loss to exit plane due to radiation");
	for (i=0; i<ns_watr; i++)
	{
		sprintf(cbuf,"heat added to water stream %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	results.push_back("total radiation heat loss");
	results.push_back("radiation heat loss through gas phase");
	results.push_back("radiation heat loss through particle phase");
	results.push_back("total convection heat loss");
	//flue gas stream data
	results.push_back("flue gas pressure");
	results.push_back("flue gas temperature");
	results.push_back("gas phase mass flow rate of flue gas stream");
	for (i=0; i<pstm_flue->gas.ns; i++)
	{
		str = "mass fraction of " + pstm_flue->gas.pgas[i]->name;
		results.push_back(str);
	}
	results.push_back("solid phase mass flow rate of flue gas stream");
	results.push_back("C mass fraction in solid phase of flue gas stream");
	results.push_back("H mass fraction in solid phase of flue gas stream");
	results.push_back("O mass fraction in solid phase of flue gas stream");
	results.push_back("N mass fraction in solid phase of flue gas stream");
	results.push_back("S mass fraction in solid phase of flue gas stream");
	results.push_back("Cl mass fraction in solid phase of flue gas stream");
	results.push_back("ash mass fraction in solid phase of flue gas stream");
	//additional detailed model data
	for (i=0; i<nzone; i++)
	{
		sprintf(cbuf,"gas/solid temperature in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)
	{
		sprintf(cbuf,"average wall temperature in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)
	{
		sprintf(cbuf,"total heat loss in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)
	{
		sprintf(cbuf,"radiation heat loss in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)
	{
		sprintf(cbuf,"convection heat loss in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)
	{
		sprintf(cbuf,"gas emissivity in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)
	{
		sprintf(cbuf,"gas absorption coefficient in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)
	{
		sprintf(cbuf,"particle absorption coefficient in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	//unburned carbon mass fraction (LOI) in solid phase of each non-dead zone
	for (i=izone_1st; i<nzone; i++)
	{
		if (pstm_zone[i].bsolid)
		{
			sprintf(cbuf,"unburned carbon mass fraction in solid of zone %d",i+1);
			str = cbuf;
			results.push_back(str);
		}
	}
	for (i=0; i<nzone; i++)		//average
	{
		sprintf(cbuf,"avarage incident heat flux in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)		//min
	{
		sprintf(cbuf,"minimum incident heat flux in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=0; i<nzone; i++)		//max
	{
		sprintf(cbuf,"maximum incident heat flux in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	for (i=izone_1st; i<nzone; i++)		//residence time
	{
		sprintf(cbuf,"residence time in zone %d",i+1);
		str = cbuf;
		results.push_back(str);
	}
	results.push_back("converted char mass by O2 oxidation");
	results.push_back("converted char mass by H2O gasification");
	results.push_back("converted char mass by CO2 gasification");
}