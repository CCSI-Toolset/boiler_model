// LiquidFuel.cpp: implementation of the CLiquidFuel class.
//
//////////////////////////////////////////////////////////////////////

#include "LiquidFuel.h"
#include "UtilityFunctions.h"
#include "Constants.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CLiquidFuel::CLiquidFuel()
{
	name = "Oil";
	mflow = 1;
	temp = 298.15;
	hhv = 4.4263e7;
	lhv = 0;
	hf = 0;
	hs = 0;
	h = 0;
	cp = 2090;
	fcmass[0] = 0.8706;
	fcmass[1] = 0.124;
	fcmass[2] = 0.0001;
	fcmass[3] = 0.0033;
	fcmass[4] = 0.002;
	CalcAllProperties();			//calculate moles, enthalpy etc.
}

CLiquidFuel::~CLiquidFuel()
{

}

void CLiquidFuel::Write(FILE* pf)
{
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::WriteString(name,pf);
	fwrite(&mflow,sizeof(double),1,pf);
	fwrite(&temp,sizeof(double),1,pf);
	fwrite(&hhv,sizeof(double),1,pf);
	fwrite(&lhv,sizeof(double),1,pf);
	fwrite(&hf,sizeof(double),1,pf);
	fwrite(&hs,sizeof(double),1,pf);
	fwrite(&h,sizeof(double),1,pf);
	fwrite(&cp,sizeof(double),1,pf);
	fwrite(fcmass,sizeof(double),5,pf);
	fwrite(eflow,sizeof(double),GS_NE,pf);
}

void CLiquidFuel::Read(FILE* pf)
{
	int iversion;
	fread(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::ReadString(name,pf);
	fread(&mflow,sizeof(double),1,pf);
	fread(&temp,sizeof(double),1,pf);
	fread(&hhv,sizeof(double),1,pf);
	fread(&lhv,sizeof(double),1,pf);
	fread(&hf,sizeof(double),1,pf);
	fread(&hs,sizeof(double),1,pf);
	fread(&h,sizeof(double),1,pf);
	fread(&cp,sizeof(double),1,pf);
	fread(fcmass,sizeof(double),5,pf);
	fread(eflow,sizeof(double),GS_NE,pf);
}

CLiquidFuel::CLiquidFuel(const CLiquidFuel &t)
{
	int i;
	name = t.name;
	mflow = t.mflow;
	temp = t.temp;
	hhv = t.hhv;
	lhv = t.lhv;
	hf = t.hf;
	hs = t.hs;
	h = t.h;
	cp = t.cp;
	for (i=0; i<5; i++)
		fcmass[i] = t.fcmass[i];
	for (i=0; i<GS_NE; i++)
		eflow[i] = t.eflow[i];
}

CLiquidFuel& CLiquidFuel::operator=(const CLiquidFuel& t)
{
	if (this==&t)
		return *this;
	int i;
	name = t.name;
	mflow = t.mflow;
	temp = t.temp;
	hhv = t.hhv;
	lhv = t.lhv;
	hf = t.hf;
	hs = t.hs;
	h = t.h;
	cp = t.cp;
	for (i=0; i<5; i++)
		fcmass[i] = t.fcmass[i];
	for (i=0; i<GS_NE; i++)
		eflow[i] = t.eflow[i];
	return *this;
}

void CLiquidFuel::UpdateMole()
{
    eflow[0] = mflow*fcmass[0]/CT_AWC;
	eflow[1] = mflow*fcmass[1]/CT_AWH;
	eflow[2] = mflow*fcmass[2]/CT_AWO;
	eflow[3] = mflow*fcmass[3]/CT_AWN;
	eflow[4] = mflow*fcmass[4]/CT_AWS;
	eflow[5] = 0;
	eflow[6] = 0;
}

void CLiquidFuel::UpdateLHV()
{
	//assume heat of vaporization=2393 kJ/kg based on Perry's formula
	lhv = hhv-fcmass[1]*8.93645*2393000;
}

void CLiquidFuel::UpdateHf()
{
	double dhLiquidFuel;	//standard heat of combustion of LiquidFuel (based on constant pressure)
	double dhc;		//heat of combustion for element C
	double dhh;		//heat of combustion for element H
	double dhs;		//heat of combustion for element S
	double hfh2o;	//heat of formation of H2O
	hfh2o= -68317.4*4.184*1000/(CT_AWH*2+CT_AWO);
	dhLiquidFuel = -hhv + 8314.3*298.15/2*(-fcmass[1]/2/CT_AWH+fcmass[2]/CT_AWO+fcmass[3]/CT_AWN);
	dhc = -94052*4.184*1000/CT_AWC*fcmass[0];
	dhh = -68317.4*4.184*1000/CT_AWH/2*fcmass[1];
	dhs = -70940*4.184*1000/CT_AWS*fcmass[4];
	hf = dhc+dhh+dhs-dhLiquidFuel;
}

void CLiquidFuel::UpdateHCp()
{
	//assume the constant Cp
	cp = 2090;
	hs = cp*(temp-298.15);
	//calculate total enthalpy
	h = hs+hf;
}

void CLiquidFuel::CalcAllProperties()
{
	UpdateMole();
	UpdateLHV();
	UpdateHf();
	UpdateHCp();
}

double CLiquidFuel::GetSum()
{
	int i;
	double sum = 0;
	for (i=0; i<5; i++)
		sum += fcmass[i];
	return sum;
}

void CLiquidFuel::Normalize()
{
	int i;
	double sum = GetSum();
	for (i=0; i<5; i++)
		fcmass[i] /= sum;
	UpdateMole();
}

