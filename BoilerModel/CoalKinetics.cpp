// CoalKinetics.cpp: implementation of the CCoalKinetics class.
//
//////////////////////////////////////////////////////////////////////

#include <cmath>
#include "Constants.h"

#include "CoalKinetics.h"

///////////////////////////////////////////

CCoalKinetics::CCoalKinetics()
{
	//using Pittsburgh # 8 coal as default
	fswell = 0.1;
	alpha = 0.1;
	ecoco2 = 1.25e8;
	echar_o2 = 1e8;
	echar_h2o = 2.4e8;
	echar_co2 = 2.51e8;			//from REI report for oxycombustion
	acoco2 = 4e4;
	achar_o2 = 82.97/12;		//converted from kg-C to kmol-C
	achar_h2o = 208;
	achar_co2 = 440;			//from REI report for oxycombustion
	bchar_o2 = 0.5;
	bchar_h2o = 1;
	bchar_co2 = 1;
	nchar_o2 = 0.5;
	nchar_h2o = 1;
	nchar_co2 = 1;
	//assign gas objects
	CIdealGas* pgas = CIdealGas::GetGasByName("N2");
	gas_n2 = *pgas;
	pgas = CIdealGas::GetGasByName("O2");
	gas_o2 = *pgas;
	pgas = CIdealGas::GetGasByName("H2O");
	gas_h2o = *pgas;
	pgas = CIdealGas::GetGasByName("CO2");
	gas_co2 = *pgas;
}

CCoalKinetics::~CCoalKinetics()
{
}

void CCoalKinetics::Write(FILE* pf)
{
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	fwrite(&fswell,sizeof(double),1,pf);
	fwrite(&alpha,sizeof(double),1,pf);
	fwrite(&ecoco2,sizeof(double),1,pf);
	fwrite(&echar_o2,sizeof(double),1,pf);
	fwrite(&echar_h2o,sizeof(double),1,pf);
	fwrite(&echar_co2,sizeof(double),1,pf);
	fwrite(&acoco2,sizeof(double),1,pf);
	fwrite(&achar_o2,sizeof(double),1,pf);
	fwrite(&achar_h2o,sizeof(double),1,pf);
	fwrite(&achar_co2,sizeof(double),1,pf);
	fwrite(&bchar_o2,sizeof(double),1,pf);
	fwrite(&bchar_h2o,sizeof(double),1,pf);
	fwrite(&bchar_co2,sizeof(double),1,pf);
	fwrite(&nchar_o2,sizeof(double),1,pf);
	fwrite(&nchar_h2o,sizeof(double),1,pf);
	fwrite(&nchar_co2,sizeof(double),1,pf);
}

void CCoalKinetics::Read(FILE* pf)
{
	int iversion;
	fread(&iversion,sizeof(int),1,pf);
	fread(&fswell,sizeof(double),1,pf);
	fread(&alpha,sizeof(double),1,pf);
	fread(&echar_o2,sizeof(double),1,pf);
	fread(&echar_o2,sizeof(double),1,pf);
	fread(&echar_h2o,sizeof(double),1,pf);
	fread(&echar_co2,sizeof(double),1,pf);
	fread(&acoco2,sizeof(double),1,pf);
	fread(&achar_o2,sizeof(double),1,pf);
	fread(&achar_h2o,sizeof(double),1,pf);
	fread(&achar_co2,sizeof(double),1,pf);
	fread(&bchar_o2,sizeof(double),1,pf);
	fread(&bchar_h2o,sizeof(double),1,pf);
	fread(&bchar_co2,sizeof(double),1,pf);
	fread(&nchar_o2,sizeof(double),1,pf);
	fread(&nchar_h2o,sizeof(double),1,pf);
	fread(&nchar_co2,sizeof(double),1,pf);
}

CCoalKinetics::CCoalKinetics(const CCoalKinetics &t)
{
	fswell = t.fswell;
	alpha = t.alpha;
	echar_o2 = t.echar_co2;
	echar_h2o = t.echar_h2o;
	echar_co2 = t.echar_co2;
	achar_o2 = t.achar_o2;
	achar_h2o = t.achar_h2o;
	achar_co2 = t.achar_co2;
	bchar_o2 = t.bchar_o2;
	bchar_h2o = t.bchar_h2o;
	bchar_co2 = t.bchar_co2;
	nchar_o2 = t.nchar_o2;
	nchar_h2o = t.nchar_h2o;
	nchar_co2 = t.nchar_co2;
	gas_n2 = t.gas_n2;
	gas_o2 = t.gas_o2;
	gas_h2o = t.gas_h2o;
	gas_co2 = t.gas_co2;
}

CCoalKinetics& CCoalKinetics::operator=(const CCoalKinetics& t)
{
	if (this==&t)
		return *this;
	fswell = t.fswell;
	alpha = t.alpha;
	echar_o2 = t.echar_co2;
	echar_h2o = t.echar_h2o;
	echar_co2 = t.echar_co2;
	achar_o2 = t.achar_o2;
	achar_h2o = t.achar_h2o;
	achar_co2 = t.achar_co2;
	bchar_o2 = t.bchar_o2;
	bchar_h2o = t.bchar_h2o;
	bchar_co2 = t.bchar_co2;
	nchar_o2 = t.nchar_o2;
	nchar_h2o = t.nchar_h2o;
	nchar_co2 = t.nchar_co2;
	gas_n2 = t.gas_n2;
	gas_o2 = t.gas_o2;
	gas_h2o = t.gas_h2o;
	gas_co2 = t.gas_co2;
	return *this;
}

void CCoalKinetics::CalcCharReactionRate(double yo2, double yh2o, double yco2, double temp, double pres, double dp, double* rates)
{
	//calculate total heterogeneous rate in kg-char/s plus the fraction due to O2, H2O, and CO2 reactions
	//rates: array of rates due to O2, H2O, and CO2 heterogeneous reactions
	//yo2: mole fraction of O2
	//yh2o: mole fraction of H2O
	//yco2: mole fraction of CO2
	//char oxidation: C + gamma O2 -> (2*gamma-1) CO2 + (2-2*gamma)CO,  net: (1-gamma) mole per mole C, where gamma = (2+phi)/(2+2*phi) and phi = CO/CO2 molar ratio
	//char gasification by H2O: C + H2O -> CO + H2, net 1 mole per mole C
	//char gasification by CO2: C + CO2 -> 2 CO, net 1 mole per mole C
	//net molar flow needs to be solved iteratively, so is the blowing factor
	//currently use half order oxidation kinetics and 1st order gasification kinetics
	int i;
	double ctotal = pres/CT_UR/temp;				//total concentration in kmol/m3
	double coco2 = acoco2*exp(-ecoco2/CT_UR/temp);	//CO to CO2 molar ratio phi
	double gamma = (2+coco2)/2/(1+coco2);			//gamma
	double ji;										//flux of reactant gas i
	double jnet;									//net molar flux kmol/m2-s, positive for outward flow
	double phi;										//parameter for blowing factor of mass transfer, phi = jnet/kxi
	double fblowing;								//blowing factor for high mass transfer
	double diff_o2;									//diffusivity of O2 in gas mixture (use N2 to represent gas mixture)
	double diff_h2o;								//diffusivity of H2O in gas mixture (use N2 to represent gas mixture)
	double diff_co2;								//diffusivity of CO2 in gas mixture (use N2 to represent gas mixture)
	double kxi;										//mass transfer coefficient based on mole fraction of species i
	double kxi_loc;									//local kxi, before blowing factor correction, calculated assumming Sh=2
	double kri;										//rate coefficient due to reaction of species i
	double xis;										//surface mole fraction of species i
	double sum_net;									//sum of net flux
	diff_o2 = gas_o2.CalcDiffusivity(&gas_n2,temp,pres);
	diff_h2o = gas_h2o.CalcDiffusivity(&gas_n2,temp,pres);
	diff_co2 = gas_co2.CalcDiffusivity(&gas_n2,temp,pres);
	jnet = 0;		//use 0 as initial net flux guess
	for (i=0; i<3; i++)		//iterate three times to converge jnet and blowing factors for mass transfer
	{
		sum_net = 0;
		//calculate char oxidation by O2 first, half order kinetics
		//C + gamma O2 -> (2*gamma-1) CO2 + (2-2*gamma)CO
		kxi_loc = ctotal*2*diff_o2/dp;
		phi = jnet/kxi_loc;
		if (phi<0.001)
			fblowing = 1;
		else
			fblowing = phi/(exp(phi)-1);
		kxi = fblowing*kxi_loc;
		kri = achar_o2*pow(temp,bchar_o2)*exp(-echar_o2/CT_UR/temp)*pow(ctotal,nchar_o2)*gamma;
		//solve surface mole fraction
		if (nchar_o2>0.499 && nchar_o2<0.501)
		{
			xis = (sqrt(kri*kri+4*kxi*yo2*(kxi+jnet))-kri)/(kxi+jnet)/2;
			xis *= xis;
		}
		else
		{
			if (nchar_o2>0.999 && nchar_o2<1.001)
			{
				xis = kxi*yo2/(kri+kxi+jnet);
			}
			else
			{
				printf("invalid reaction order for char oxidation!\n");
				return;
			}
		}
		//calculate flux of O2
		if (xis<0.1*yo2)		//use mass transfer formula
			ji = kxi*(yo2-xis) - jnet*xis;
		else
			ji = kri*pow(xis,nchar_o2);
		sum_net += (1-gamma)/gamma*ji;
		rates[0] = ji/gamma;
		
		//calculate char gasification by H2O, C + H2O -> CO + H2
		kxi_loc = ctotal*2*diff_h2o/dp;
		phi = jnet/kxi_loc;
		if (phi<0.001)
			fblowing = 1;
		else
			fblowing = phi/(exp(phi)-1);
		kxi = fblowing*kxi_loc;
		kri = achar_h2o*pow(temp,bchar_h2o)*exp(-echar_h2o/CT_UR/temp)*pow(ctotal,nchar_h2o);
		if (nchar_h2o>0.999 && nchar_h2o<1.001)
			xis = kxi*yh2o/(kri+kxi+jnet);
		else
		{
			printf("invalid reaction order for char gasification by H2O!\n");
			return;
		}
		//calculate flux of H2O
		if (xis<0.1*yh2o)		//use mass transfer formula
			ji = kxi*(yh2o-xis) - jnet*xis;
		else
			ji = kri*pow(xis,nchar_h2o);
		sum_net += ji;
		rates[1] = ji;

		//calculate char gasification by CO2, C + CO2 -> 2CO
		kxi_loc = ctotal*2*diff_co2/dp;
		phi = jnet/kxi_loc;
		if (phi<0.001)
			fblowing = 1;
		else
			fblowing = phi/(exp(phi)-1);
		kxi = fblowing*kxi_loc;
		kri = achar_co2*pow(temp,bchar_co2)*exp(-echar_co2/CT_UR/temp)*pow(ctotal,nchar_co2);
		if (nchar_co2>0.999 && nchar_co2<1.001)
			xis = kxi*yco2/(kri+kxi+jnet);
		else
		{
			printf("invalid reaction order for char gasification by CO2!\n");
			return;
		}
		//calculate flux of CO2
		if (xis<0.1*yco2)		//use mass transfer formula
			ji = kxi*(yco2-xis) - jnet*xis;
		else
			ji = kri*pow(xis,nchar_co2);
		sum_net += ji;
		rates[2] = ji;
		//update jnet
		jnet = sum_net;
	}
	//calculate mass reaction rate of a particle
	for (i=0; i<3; i++)
		rates[i] *= CT_PI*dp*dp*CT_AWC;
}
