// CoalKinetics.h: interface for the CCoalKinetics class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __COALKINETICS_H__
#define __COALKINETICS_H__

#include <cstdio>
#include "IdealGas.h"			//GS_NE defined in CIdealGas

class CCoalKinetics
{
public:
	double	fswell;				//swelling factor after devolatilization
	double	alpha;				//burning mode parameter, denp/denp0 = (mp/mp0)^alpha, alpha<1
	double  ecoco2;				//activation energy of CO/CO2 molar ratio for char oxidation
	double	echar_o2;			//activation energy of char oxidation
	double	echar_h2o;			//activation energy of char gasification by H2O
	double	echar_co2;			//activation energy of char gasification by CO2
	double	acoco2;				//pre-exponential factor of CO/CO2 molar ratio for char oxidation
	double	achar_o2;			//pre-exponential factor of char oxidation, in kmol-C/m2-s-K^b-Ci^n
	double	achar_h2o;			//pre-exponential factor of char gasification by H2O, in kmol-C/m2-s-K^b-Ci^n
	double	achar_co2;			//pre-exponential factor of char gasification by CO2, in kmol-C/m2-s-K^b-Ci^n
	double	bchar_o2;			//temperature exponent of char oxdiation
	double	bchar_h2o;			//temperature exponent of char gasification by H2O
	double	bchar_co2;			//temperature exponent of char gasification by CO2
	double	nchar_o2;			//reaction order of char combustion
	double	nchar_h2o;			//reaction order of char gasification by H2O
	double	nchar_co2;			//reaction order of char gasification by CO2

	CIdealGas gas_n2;			//N2 gas, assigned in constructor, not saved to file
	CIdealGas gas_o2;			//O2 gas, assigned in constructor, not saved to file
	CIdealGas gas_h2o;			//H2O gas, assigned in constructor, not saved to file
	CIdealGas gas_co2;			//CO2 gas, assigned in constructor, not saved to file

	//member functions
	CCoalKinetics();
	virtual ~CCoalKinetics();
	CCoalKinetics(const CCoalKinetics &t);
	CCoalKinetics& operator=(const CCoalKinetics& t);
	void CalcCharReactionRate(double yo2, double yh2o, double yco2, double temp, double pres, double dp, double* rates);
	void Write(FILE* pf);
	void Read(FILE* pf);
};

#endif