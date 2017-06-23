// SolidFuel.h: interface for the CSolidFuel class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __SOLIDFUEL_H__
#define __SOLIDFUEL_H__

#include <cstdio>
#include "IdealGas.h"			//GS_NE defined in CIdealGas

class CSolidFuel
{
public:
	//member data
	//index of array: 0=as re'd, 1=dry, 2=daf
	std::string name;			//name of coal
	int		ibasis;				//basis for input. 0=as re'd, 1=dry, 2=daf
	int		ipsd;				//particle size distribution index, -1=N/A
	int		ick;				//coal kinetics index, -1=N/A
	double	mflow;				//mass flow rate [kg/s]
	double	temp;				//temperature of the coal [K]
	double	hhv[3];				//high heating value [J/kg], 0=as re'd, 1=dry, 2=daf
	double	lhv[3];				//low heating value [J/kg]
	double	hf[3];				//heat of formation [J/kg]
	double	hs[3];				//sensible heat [J/kg]
	double	h[3];				//total enthalpy [J/kg]
	double	cp[3];				//heat capacity of the coal [J/kgK]
	double	fcmass[3][9];		//mass fraction of each component in the coal 0=C,1=H,2=O,3=N,4=S,5=Cl, 6=moisture,7=ash,8=volatile
	double	eflow[2][GS_NE];	//kmole flow rate of each element (0=C, 1=H, 2=O, 3=N, 4=S, 5=Cl, 6=Ar) [kmole/s], eflow[0][] for all elements, eflow[1][] for moisture only

	//member functions
	CSolidFuel();
	virtual ~CSolidFuel();
	CSolidFuel(const CSolidFuel &t);
	CSolidFuel& operator=(const CSolidFuel& t);
	void UpdateMass();
	void CalcHHVByDulong();
	void CalcHHVAsElements();
	void UpdateHHV();
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