// LiquidFuel.h: interface for the CLiquidFuel class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __LIQUIDFUEL_H__
#define __LIQUIDFUEL_H__

#include <cstdio>
#include "IdealGas.h"

class CLiquidFuel
{
public:
	//data
	std::string name;		//name of LiquidFuel
	double	mflow;			//mass flow rate [kg/s]
	double	temp;			//temperature [K]
	double	hhv;			//high heating value [J/kg]
	double	lhv;			//low heating value [J/kg]
	double	hf;				//heat of formation [J/kg]
	double	hs;				//sensible heat [J/kg]
	double	h;				//total enthalpy [J/kg]
	double	cp;				//heat capacity [J/kgK] (currently assume constant)
	double	fcmass[5];		//mass fraction of each component 0=C,1=H,2=O,3=N,4=S, no ash
	double	eflow[GS_NE];	//kmole flow rate of each element (0=C, 1=H, 2=O, 3=N, 4=S, 5=Cl, 6=Ar) [kmole/s]

	//functions
	CLiquidFuel();
	virtual ~CLiquidFuel();
	CLiquidFuel(const CLiquidFuel &t);
	CLiquidFuel& operator=(const CLiquidFuel& t);
	void UpdateMole();
	void UpdateLHV();
	void UpdateHf();
	void UpdateHCp();
	void CalcAllProperties();
	double GetSum();
	void Normalize();
	void Write(FILE* pf);
	void Read(FILE* pf);
};

#endif