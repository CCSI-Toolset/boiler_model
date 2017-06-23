// MaterialStream.cpp: implementation of the CMaterialStream class.
//
//////////////////////////////////////////////////////////////////////

#include <cmath>
#include <sstream>
#include "Constants.h"
#include "UtilityFunctions.h"
#include "MaterialStream.h"

///////////////////////////////////////////

CMaterialStream::CMaterialStream()
{
	name = "Stream";
	bgas = true;
	bliquid = false;
	bsolid = false;
	temp = 300;
	pres = 101325;
	er = 1;
	mflow = 1;
	hflow = 1;
	UpdateProperties();
}

CMaterialStream::~CMaterialStream()
{
}

CMaterialStream::CMaterialStream(const CMaterialStream &t)
{
	name = t.name;
	bgas = t.bgas;
	bliquid = t.bliquid;
	bsolid = t.bsolid;
	temp = t.temp;
	pres = t.pres;
	er = t.er;
	mflow = t.mflow;
	hflow = t.hflow;
	gas = t.gas;
	oil = t.oil;
	coal = t.coal;
}

CMaterialStream CMaterialStream::operator=(CMaterialStream t)
{
	name = t.name;
	bgas = t.bgas;
	bliquid = t.bliquid;
	bsolid = t.bsolid;
	temp = t.temp;
	pres = t.pres;
	er = t.er;
	mflow = t.mflow;
	hflow = t.hflow;
	gas = t.gas;
	oil = t.oil;
	coal = t.coal;
	return *this;
}

void CMaterialStream::Write(FILE* pf)
{
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::WriteString(name,pf);
	fwrite(&bgas,sizeof(bool),1,pf);
	fwrite(&bliquid,sizeof(bool),1,pf);
	fwrite(&bsolid,sizeof(bool),1,pf);
	fwrite(&temp,sizeof(double),1,pf);
	fwrite(&pres,sizeof(double),1,pf);
	fwrite(&er,sizeof(double),1,pf);
	fwrite(&mflow,sizeof(double),1,pf);
	fwrite(&hflow,sizeof(double),1,pf);
	gas.Write(pf);
	oil.Write(pf);
	coal.Write(pf);
}

void CMaterialStream::Read(FILE* pf)
{
	int iversion;
	fread(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::WriteString(name,pf);
	fread(&bgas,sizeof(bool),1,pf);
	fread(&bliquid,sizeof(bool),1,pf);
	fread(&bsolid,sizeof(bool),1,pf);
	fread(&temp,sizeof(double),1,pf);
	fread(&pres,sizeof(double),1,pf);
	fread(&er,sizeof(double),1,pf);
	fread(&mflow,sizeof(double),1,pf);
	fread(&hflow,sizeof(double),1,pf);
	gas.Read(pf);
	oil.Read(pf);
	coal.Read(pf);
}

void CMaterialStream::CopyProperties(CMaterialStream* ps)
{
	//does not copy name, locations, connectivity, color and flags (bhide, bhide_label and bknown)
	//this function is not called by CMaterialStreamDlg functions (called only when solving flowsheet)
	bgas = ps->bgas;
	bliquid = ps->bliquid;
	bsolid = ps->bsolid;
	temp = ps->temp;
	pres = ps->pres;
	er = ps->er;
	mflow = ps->mflow;
	hflow = ps->hflow;
	//save names, icalc and velocity or area of gas
	std::string name_gas = gas.name;
	std::string name_oil = oil.name;
	std::string name_coal = coal.name;
	int icalc_gas = gas.icalc;
	double vel_gas = gas.vel;
	double area_gas = gas.area;
	gas = ps->gas;
	oil = ps->oil;
	coal = ps->coal;
	//retrieve saved data
	gas.name = name_gas;
	oil.name = name_oil;
	coal.name = name_coal;
	gas.icalc = icalc_gas;
	if (icalc_gas==0)
		gas.area = area_gas;
	else
		gas.vel = vel_gas;
}

void CMaterialStream::UpdateProperties()
{
	//calculating mass, energy, loading, er, etc.
	int i;
	double massflow;
	double covfuel;
	double covoxy;
	double eflow[GS_NE] = {0};
	mflow = 0;
	hflow = 0;
	gas.t = temp;
	oil.temp = temp;
	coal.temp = temp;
	gas.p = pres;
	if (bgas)
	{
		gas.InitArraysFromFsp();
		gas.CalcAllProperties();
		massflow = gas.mflow;
		mflow += massflow;
		hflow += massflow*gas.h;
		for (i=0; i<GS_NE; i++)
			eflow[i] += gas.eflow[i];
	}
	else
	{
		gas.mflow = 0;
		gas.vflow = 0;
		for (i=0; i<GS_NE; i++)
			gas.eflow[i] = 0;
	}
	if (bliquid)
	{
		oil.CalcAllProperties();
		massflow = oil.mflow;
		mflow += massflow;
		hflow += massflow*oil.h;
		for (i=0; i<GS_NE; i++)
			eflow[i] += oil.eflow[i];
	}
	else
	{
		oil.mflow = 0;
		for (i=0; i<GS_NE; i++)
			oil.eflow[i] = 0;
	}
	if (bsolid)
	{
		coal.CalcAllProperties();
		massflow = coal.mflow;
		mflow += massflow;
		hflow += massflow*coal.h[0];
		for (i=0; i<GS_NE; i++)
			eflow[i] += coal.eflow[0][i];
	}
	else
	{
		coal.mflow = 0;
		for (i=0; i<GS_NE; i++)
			coal.eflow[0][i] = 0;
	}
	covfuel = 4*eflow[0] + eflow[1] + 4*eflow[4];
	covoxy = 2*eflow[2];
	if (covoxy>0)
		er = covfuel/covoxy;
	else
		er = -1;
}

void CMaterialStream::ModifyMassFlow(double mf, int iphase)
{
	switch (iphase)
	{
	case 0:
		if (bgas)
		{
			gas.mflow = mf;
			gas.InitArraysFromFsp();
			gas.vflow = gas.mflow / (101325*gas.mw/8314.3/273.15);
		}
		break;
	case 1:
		if (bliquid)
			oil.mflow = mf;
		break;
	case 2:
		if (bsolid)
			coal.mflow = mf;
		break;
	}
	UpdateProperties();
}

double CMaterialStream::CalcSFlow()
{
	//currently ignore liquid and solid phases
	double sflow = 0;
	gas.t = temp;
	oil.temp = temp;
	coal.temp = temp;
	gas.p = pres;
	if (bgas)
		sflow += gas.mflow*gas.CalcS();
	return sflow;
}

double CMaterialStream::CalcCp()
{
	double cp = 0;
	if (bgas)
		cp += gas.mflow/mflow*gas.CalcCp();
	if (bliquid)
	{
		oil.UpdateHCp();
		cp += oil.mflow/mflow*oil.cp;
	}
	if (bsolid)
	{
		coal.UpdateHCp();
		cp += coal.mflow/mflow*coal.cp[0];
	}
	return cp;
}

double CMaterialStream::GetHHVRate()
{
	double hhvrate = 0;
	if (bgas)
		hhvrate += gas.mflow*gas.hhv;
	if (bliquid)
		hhvrate += oil.mflow*oil.hhv;
	if (bsolid)
		hhvrate += coal.mflow*coal.hhv[0];
	return hhvrate;
}
