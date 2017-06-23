// IdealGas.h: interface for the CIdealGas class.
//////////////////////////////////////////////////////////////////////

#ifndef __IDEALGAS_H__
#define __IDEALGAS_H__

#include <string>

const int GS_NE = 7;					//number of elements
const int GS_NSP = 29;					//number of availabe species

class CIdealGas
{
public:
	static double ljkte[79];				//L-J kte points
	static double ljomega_mk[79];			//L-J omega for viscosity and thermal conductivity
	static double ljomega_d[79];			//L-J omega for diffussivity

	std::string name;	//name of the gas
	double tempth;		//threshold of the temperature [K]
	double mw;			//molecular weight of the species [kg/kmol]
	double cp;			//molar heat capacity [J/kmol-K]
	double h;			//molar enthalpy [J/kmol]
	double s;			//molar entropy [J/kmol-K]
	double g;			//molar Gibbs-free energy [J/kmol]
	double vis;			//viscosity [kg/m-sec]
	double cond;		//thermal conductivity [J/m-sec-K]
	double x0[2][2];	//initial ln(ni) for equilibrium calculation
	int iformat;		//format of the property coefficient, 0=JANAF, 1=NIST
	int nelem;			//number of elements
	int matom[GS_NE];	//atom number, order based on molecular formula
	int natom[GS_NE];	//atom number for each elements 0=C, 1=H, 2=O, 3=N, 4=S, 5=Cl, 6=Ar
	char elem[GS_NE];	//atom name, order based on molecular formula
	double z[2][7];		//upper and lower temperature data
	double	ljsigma;	//L-J parameter sigma in A (1e-10 m)
	double	ljek;		//L-J parameter e/k

private:
	static CIdealGas** pGasList;

public:
	CIdealGas();
	CIdealGas(std::string _name);
	virtual ~CIdealGas();
	CIdealGas(const CIdealGas &tt);
	CIdealGas& operator=(const CIdealGas& tt);
	void SetName(std::string str);
	void AssignData();
	double GetX0(double t, double er);
	double CalcMW();
	double CalcCp(double t);
	double CalcH(double t);
	double CalcS(double t, double p=1.013e5);
	double CalcG(double t, double p=1.013e5);
	double CalcViscosity(double t);
	double CalcConductivity(double t);
	double CalcDiffusivity(CIdealGas* pg, double t, double p);
	int GetIndex();

	static CIdealGas** GetGasList();
	static void DeleteGasList();
	static CIdealGas* GetGasByName(std::string str);
	static CIdealGas* GetGasByIndex(int index);
	static int GetIndexByName(std::string str);
};

#endif