// SolidFuel.cpp: implementation of the CSolidFuel class.
//
//////////////////////////////////////////////////////////////////////

#include <cmath>
#include "Constants.h"
#include "UtilityFunctions.h"

#include "SolidFuel.h"

///////////////////////////////////////////

CSolidFuel::CSolidFuel()
{
	int i;
	name = "Coal";
	mflow = 1;
	ibasis = 0;
	ipsd = 0;
	ick = 0;
	temp = 298.15;
	hhv[0] = 2.7e7;
	hhv[1] = 0;
	hhv[2] = 0;
	for (i=0; i<3; i++)
	{
		lhv[i] = 0;
		hf[i] = 0;
		hs[i] = 0;
		h[i] = 0;
		cp[i] = 0;
		fcmass[i][0] = 0.65;
		fcmass[i][1] = 0.045;
		fcmass[i][2] = 0.09;
		fcmass[i][3] = 0.015;
		fcmass[i][4] = 0.01;
		fcmass[i][5] = 0.0;
		fcmass[i][6] = 0.09;
		fcmass[i][7] = 0.1;
		fcmass[i][8] = 0.3;
	}
	for (i=0; i<GS_NE; i++)
	{
		eflow[0][i] = 0;
		eflow[1][i] = 0;
	}
	CalcAllProperties();
}

CSolidFuel::~CSolidFuel()
{
}

void CSolidFuel::Write(FILE* pf)
{
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::WriteString(name,pf);
	fwrite(&ibasis,sizeof(int),1,pf);
	fwrite(&ipsd,sizeof(int),1,pf);
	fwrite(&ick,sizeof(int),1,pf);
	fwrite(&mflow,sizeof(double),1,pf);
	fwrite(&temp,sizeof(double),1,pf);
	fwrite(hhv,sizeof(double),3,pf);
	fwrite(lhv,sizeof(double),3,pf);
	fwrite(hf,sizeof(double),3,pf);
	fwrite(hs,sizeof(double),3,pf);
	fwrite(h,sizeof(double),3,pf);
	fwrite(cp,sizeof(double),3,pf);
	fwrite(fcmass,sizeof(double),27,pf);
	fwrite(eflow,sizeof(double),GS_NE,pf);
}

void CSolidFuel::Read(FILE* pf)
{
	int iversion;
	fread(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::ReadString(name,pf);
	fread(&ibasis,sizeof(int),1,pf);
	fread(&ipsd,sizeof(int),1,pf);
	fread(&ick,sizeof(int),1,pf);
	fread(&mflow,sizeof(double),1,pf);
	fread(&temp,sizeof(double),1,pf);
	fread(hhv,sizeof(double),3,pf);
	fread(lhv,sizeof(double),3,pf);
	fread(hf,sizeof(double),3,pf);
	fread(hs,sizeof(double),3,pf);
	fread(h,sizeof(double),3,pf);
	fread(cp,sizeof(double),3,pf);
	fread(fcmass,sizeof(double),27,pf);
	fread(eflow,sizeof(double),GS_NE,pf);
}

CSolidFuel::CSolidFuel(const CSolidFuel &t)
{
	int i, j;
	name = t.name;
	ibasis = t.ibasis;
	ipsd = t.ipsd;
	ick = t.ick;
	mflow = t.mflow;
	temp = t.temp;
	for (i=0; i<3; i++)
	{
		hhv[i] = t.hhv[i];
		lhv[i] = t.lhv[i];
		hf[i] = t.hf[i];
		hs[i] = t.hs[i];
		h[i] = t.h[i];
		cp[i] = t.cp[i];
		for (j=0; j<9; j++)
			fcmass[i][j] = t.fcmass[i][j];
	}
	for (i=0; i<GS_NE; i++)
	{
		eflow[0][i] = t.eflow[0][i];
		eflow[1][i] = t.eflow[1][i];
	}
}

CSolidFuel& CSolidFuel::operator=(const CSolidFuel& t)
{
	if (this==&t)
		return *this;
	int i, j;
	name = t.name;
	ibasis = t.ibasis;
	ipsd = t.ipsd;
	ick = t.ick;
	mflow = t.mflow;
	temp = t.temp;
	for (i=0; i<3; i++)
	{
		hhv[i] = t.hhv[i];
		lhv[i] = t.lhv[i];
		hf[i] = t.hf[i];
		hs[i] = t.hs[i];
		h[i] = t.h[i];
		cp[i] = t.cp[i];
		for (j=0; j<9; j++)
			fcmass[i][j] = t.fcmass[i][j];
	}
	for (i=0; i<GS_NE; i++)
	{
		eflow[0][i] = t.eflow[0][i];
		eflow[1][i] = t.eflow[1][i];
	}
	return *this;
}

void CSolidFuel::UpdateMass()
{
	int i;
	double d1, d2;
	//calculate fcmass[3][9]
	switch (ibasis)
	{
	case 0:			//as re'd known
		d1 = 1-fcmass[0][6];
		d2 = 1-fcmass[0][6]-fcmass[0][7];
		for (i=0; i<9; i++)
		{
			if (d1>0)
				fcmass[1][i] = fcmass[0][i]/(1-fcmass[0][6]);
			else
				fcmass[1][i] = 0;
			if (d2>0)
				fcmass[2][i] = fcmass[0][i]/(1-fcmass[0][6]-fcmass[0][7]);
			else
				fcmass[2][i] = 0;
		}
		fcmass[1][6] = 0;
		fcmass[2][6] = 0;
		fcmass[2][7] = 0;
		break;
	case 1:			//dry known, fcmass[0][6] known
		d1 = 1-fcmass[0][6];
		d2 = 1-fcmass[1][7];
		for (i=0; i<9; i++)
		{
			if (i==6) continue;
			if (d1>0)
				fcmass[0][i] = fcmass[1][i]*(1-fcmass[0][6]);
			else
				fcmass[0][i] = 0;
			if (d2>0)
				fcmass[2][i] = fcmass[1][i]/(1-fcmass[1][7]);
			else
				fcmass[2][i] = 0;
		}
		fcmass[1][6] = 0;
		fcmass[2][6] = 0;
		fcmass[2][7] = 0;
		break;
	case 2:			//daf known, fcmass[0][6] and fcmass[0][7] known
		d1 = 1-fcmass[0][6];
		d2 = 1-fcmass[0][6]-fcmass[0][7];
		for (i=0; i<9; i++)
		{
			if (i==6 || i==7) continue;
			if (d2>0)
				fcmass[0][i] = fcmass[2][i]*(1-fcmass[0][6]-fcmass[0][7]);
			else
				fcmass[0][i] = 0;
		}
		for (i=0; i<9; i++)
		{
			if (d1>0)
				fcmass[1][i] = fcmass[0][i]/(1-fcmass[0][6]);
			else
				fcmass[1][i] = 0;
		}
		fcmass[1][6] = 0;
		fcmass[2][6] = 0;
		fcmass[2][7] = 0;
		break;
	}
	//calculate eflow[6] (these elements will go to the gas phase
	eflow[0][0] = mflow*fcmass[0][0]/CT_AWC;
	eflow[0][1] = mflow*fcmass[0][1]/CT_AWH + mflow*fcmass[0][6]/(2*CT_AWH+CT_AWO)*2;
	eflow[0][2] = mflow*fcmass[0][2]/CT_AWO + mflow*fcmass[0][6]/(2*CT_AWH+CT_AWO);
	eflow[0][3] = mflow*fcmass[0][3]/CT_AWN;
	eflow[0][4] = mflow*fcmass[0][4]/CT_AWS;
	eflow[0][5] = mflow*fcmass[0][5]/CT_AWCL;
	eflow[0][6] = 0;
	eflow[1][0] = 0;
	eflow[1][1] = mflow*fcmass[0][6]/(2*CT_AWH+CT_AWO)*2;
	eflow[1][2] = mflow*fcmass[0][6]/(2*CT_AWH+CT_AWO);
	eflow[1][3] = 0;
	eflow[1][4] = 0;
	eflow[1][5] = 0;
	eflow[1][6] = 0;
}

void CSolidFuel::CalcHHVByDulong()
{
	//0=as received, 1=dry, 2=daf
	int i;
	for (i=0; i<3; i++) hhv[i] = 232367.224*
		(145.44*fcmass[i][0]+620*(fcmass[i][1]-fcmass[i][2]/8)+41*fcmass[i][4]);
}

void CSolidFuel::CalcHHVAsElements()
{
	//0=as received, 1=dry, 2=daf
	//ignore heat of combustion of Cl  (2Cl->Cl2)
	int i;
	double dhc;
	double dhh;
	double dhs;
	for (i=0; i<3; i++)
	{
		dhc = 94052*4.184*1000/CT_AWC*fcmass[i][0];
		dhh = 68317.4*4.184*1000/CT_AWH/2*fcmass[i][1];
		dhs = 70940*4.184*1000/CT_AWS*fcmass[i][4];
		hhv[i] = dhc + dhh + dhs;
	}
}

void CSolidFuel::UpdateHHV()
{
	//0=as received, 1=dry, 2=daf
	double d1 = 1-fcmass[0][6];
	double d2 = 1-fcmass[0][6]-fcmass[0][7];
	switch (ibasis)
	{
	case 0:		//hhv[0] already known
		if (d1>0)
			hhv[1] = hhv[0]/(1-fcmass[0][6]);
		else
			hhv[1] = 0;
		if (d2>0)
			hhv[2] = hhv[0]/(1-fcmass[0][6]-fcmass[0][7]);
		else
			hhv[2] = 0;
		break;
	case 1:		//hhv[1] already known
		if (d1>0)
			hhv[0] = hhv[1]*(1-fcmass[0][6]);
		else
			hhv[0] = 0;
		if (d2>0)
			hhv[2] = hhv[0]/(1-fcmass[0][6]-fcmass[0][7]);
		else
			hhv[2] = 0;
		break;
	case 2:		//hhv[2] already known
		if (d2>0)
			hhv[0] = hhv[2]*(1-fcmass[0][6]-fcmass[0][7]);
		else
			hhv[0] = 0;
		if (d1>0)
			hhv[1] = hhv[0]/(1-fcmass[0][6]);
		else
			hhv[1] = 0;
		break;
	}
}

void CSolidFuel::UpdateLHV()
{
	//assume heat of vaporization=2393 kJ/kg based on Perry's formula
	lhv[0] = hhv[0]-(fcmass[0][6]+fcmass[0][1]*8.93645)*2393000;
	lhv[1] = hhv[1]-fcmass[1][1]*8.93645*2393000;
	lhv[2] = hhv[2]-fcmass[2][1]*8.93645*2393000;
}

void CSolidFuel::UpdateHf()
{
	//assume heat of formation of ash is zero
	double dhcoal;	//standard heat of combustion of coal (based on constant pressure)
	double dhc;		//heat of combustion for element C
	double dhh;		//heat of combustion for element H
	double dhs;		//heat of combustion for element S
	double hfh2o;	//heat of formation of H2O
	hfh2o= -68317.4*4.184*1000/(CT_AWH*2+CT_AWO);
	dhcoal = -hhv[2] + 8314.3*298.15/2*(-fcmass[2][1]/2/CT_AWH+fcmass[2][2]/CT_AWO+fcmass[2][3]/CT_AWN+fcmass[2][5]/CT_AWCL);
	dhc = -94052*4.184*1000/CT_AWC*fcmass[2][0];
	dhh = -68317.4*4.184*1000/CT_AWH/2*fcmass[2][1];
	dhs = -70940*4.184*1000/CT_AWS*fcmass[2][4];
	hf[2] = dhc+dhh+dhs-dhcoal;
	hf[1] = hf[2]*(1-fcmass[1][7]);
	hf[0] = hf[2]*(1-fcmass[0][6]-fcmass[0][7]) + hfh2o*fcmass[0][6];
}

void CSolidFuel::UpdateHCp()
{
	double a;					//average atomic weight
	double hsorg;				//sensible heat of organic
	double hsash;				//sensible heat of ash
	double hsmoist;				//sensible heat of moisture
	double cporg;				//heat capacity of organic
	double cpash;				//heat capacity of ash
	double cpmoist;				//heat capacity of moisture
	double tt1,tt2,tto1,tto2;	//local varible used for organic part calculation
	double gt1,gt2,go1,go2;		//local varible used for organic part calculation
	double gcp1,gcp2;			//local varible used for organic part calculation
	double exp1,exp2;			//local varible used for organic part calculation
	double denom1,denom2;		//local varible used for organic part calculation
	tt1 = 380/temp;
	tt2 = 1800/temp;
	tto1 = 380/298.15;
	tto2 = 1800/298.15;
	exp1 = exp(tt1);
	exp2 = exp(tt2);
	denom1 = (exp1-1)/tt1;
	denom2 = (exp2-1)/tt2;
	gcp1 = exp1/denom1/denom1;
	gcp2 = exp2/denom2/denom2;
	gt1 = 1/(exp1-1);
	gt2 = 1/(exp2-1);
	go1 = 1/(exp(tto1)-1);
	go2 = 1/(exp(tto2)-1);	
	a = 1/(fcmass[2][0]/CT_AWC+fcmass[2][1]/CT_AWH+fcmass[2][2]/CT_AWO+fcmass[2][3]/CT_AWN+fcmass[2][4]/CT_AWS+fcmass[2][5]/CT_AWCL);
	hsorg = 8314.3/a*(380*(gt1-go1)+3600*(gt2-go2));
	cporg = 8314.3/a*(gcp1+2*gcp2);
	hsash = 593*(temp-298.15)+0.293*(temp*temp-298.15*298.15);
	cpash = 593 + 0.586*temp;
	hsmoist = 4184*(temp-298.15);
	cpmoist = 4184;
	//calculate sensible heat
	hs[0] = (1-fcmass[0][6]-fcmass[0][7])*hsorg+fcmass[0][6]*hsmoist+fcmass[0][7]*hsash;
	hs[1] = (1-fcmass[1][7])*hsorg+fcmass[1][7]*hsash;
	hs[2] = hsorg;
	//calculate heat capacity
	cp[0] = (1-fcmass[0][6]-fcmass[0][7])*cporg+fcmass[0][6]*cpmoist+fcmass[0][7]*cpash;
	cp[1] = (1-fcmass[1][7])*cporg+fcmass[1][7]*cpash;
	cp[2] = cporg;
	//calculate total enthalpy
	h[0] = hs[0]+hf[0];
	h[1] = hs[1]+hf[1];
	h[2] = hs[2]+hf[2];
}

void CSolidFuel::CalcAllProperties()
{
	UpdateMass();
	UpdateHHV();
	UpdateLHV();
	UpdateHf();
	UpdateHCp();
}

double CSolidFuel::GetSum()
{
	int i;
	double sum = 0;
	for (i=0; i<8; i++)
		sum += fcmass[ibasis][i];
	return sum;
}

void CSolidFuel::Normalize()
{
	int i;
	double sum = GetSum();
	for (i=0; i<8; i++)
		fcmass[ibasis][i] /= sum;
	UpdateMass();
}

