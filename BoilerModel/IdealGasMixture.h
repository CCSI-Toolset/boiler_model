//IdealGasMixture.h: interface for the CIdealGasMixture class.
/////////////////////////////////////////////////
#ifndef __IDEALGASMIXTURE_H__
#define __IDEALGASMIXTURE_H__

#include <cstdio>
#include <vector>
#include "IdealGas.h"

class CIdealGasMixture
{
public:
	std::string name;		//name of the mixture
	bool bdry;				//dry basis for mol/mass fractions
	int ivol;				//specified by volumetric flow
	int imol;				//show mole fraction if 0 and mass fraction if 1 for species
	int ictp;				//const T and P if 1 and const H and P if 0
	int icalc;				//calculation flag 0=velocity, 1=area, 2=mass flow
	int ne;					//number of elements
	int ns;					//number of species
	double mflow;			//mass flow rate [kg/s]
	double vflow;			//volumetric flow rate [Nm3/s]
	double area;			//area [m2]
	double vel;				//velocity [m/s]
	double p;				//total pressure [Pa]
	double t;				//temperature [K]
	double h;				//total enthalpy in 1 kg mixture [J/kg]
	double g;				//total Gibbs free energy in 1 kg mixture [J/kg]
	double hhv;				//high heating value [J/kg]
	double lhv;				//low heating value [J/kg]
	double hsh;				//sensible heat [J/kg]
	double cp;				//heat capacity [J/kg-K]
	double mw;				//average molecular weight of mixture [kg/kmol]
	double den;				//density [kg/m3]
	double er;				//equivalence ratio
	double eflow[GS_NE];	//kmole flow rate of each element (0=C, 1=H, 2=O, 3=N, 4=S, 5=Cl, 6=Ar) [kmole/s]
	double fsp_mole[50];	//species mole fractions
	double fsp_mass[50];	//species mass fractions
	double bj[GS_NE];		//kmole of atoms in 1 kg mixture [kmol/kg], (not necesarily in C,H,O,N,S,Cl order)
	double ni[50];			//specific mole, ni for each species

	CIdealGas* pgas[50];			//gas array

protected:
	char elem[GS_NE];			//element name array, maximum size is 5 (not necesarily in C,H,O,N,S,Cl,Ar order)
	int ielem[GS_NE];			//index of each element (0=C,1=H,2=O,3=N,4=S,5=Cl,6=Ar)
	double aw[GS_NE];			//atomic weight of each element (not necesarily in C,H,O,N,S,Cl,Ar order)
	static double x[50+GS_NE];			//ln(ni) and Lagrangian coefficiets
	static double dx[50+GS_NE];			//correction vector for x in Newton-Raphson method
	static double aij[50][GS_NE];		//number of atoms of each element in a molecule
	static double jacob[50+GS_NE][51+GS_NE];	//augmented Jacobian matrix

	//protected member functions
	void CalcAij();
	void CalcNiFromX();
	void CalcJacobian();
	int SolveLAE(int n, double a[50+GS_NE][51+GS_NE], double x[50+GS_NE]);
	int SolveXAtTP();
	int SolveXAtHP();

public:
	//constructor and destructor
	CIdealGasMixture();
	~CIdealGasMixture();
	//public member functions
	CIdealGasMixture(const CIdealGasMixture &tt);
	CIdealGasMixture operator=(CIdealGasMixture tt);
	void Write(FILE* pf);
	void Read(FILE* pf);
	void InitArraysFromFsp();					//used by gas substream input dialog
	void InitArraysFromEflow();					//used to calculate equilibrium based on eflow[] and hflow
	void GetDryFsp(double f[50]);				//used by gas substream input dialog
	double GetFspSum();							//used by gas substream input dialog
	void NormalizeFsp();						//used by gas substream input dialog
	void CalcFspFromNi();						//used for equilibrium dialog
	bool AreAllElementsInGasArray();			//check if any element not associated with a gas species
	void GetValidSpeciesList(std::vector<CIdealGas*>& plist);	//used for equilibrium dialg
	void AddValidSpecies(std::string str, std::vector<CIdealGas*>& pvl);			//used for equilibrium dialog
	void AddValidSpecies(int isp, std::vector<CIdealGas*>& pvl);				//used for equilibrium reactor
	void DeleteSpecies(int i);
	void DeleteSpecies(std::string str);
	void SetT(double tt) {t = tt;}
	void SetH(double hh) {h = hh;}
	void SetP(double pp=1.013e5) {p = pp;}
	int Solve();
	double CalcH();
	double CalcS();
	double CalcCp();
	double CalcCp(double temp);
	double CalcEr();
	double CalcViscosity(double temp);
	double CalcConductivity(double temp);
	void CalcAllProperties();
	int GetNs() {return ns;}
	int GetNe() {return ne;}
	double GetT() {return t;}
	double GetP() {return p;}
	double GetH() {return h;}
	double GetG() {return g;}
	double GetHHV() {return hhv;}
	double GetLHV() {return lhv;}
	double GetHSH() {return hsh;}
	double GetCp() {return cp;}
	double GetDen() {return den;}
	double GetMw() {return mw;}
	double GetEr() {return er;}
	int GetSpeciesIndex(std::string sp);
	double GetMoleFraction(std::string sp);
	double GetMoleFraction(int i);
	double GetSpeciesMoleFlow(std::string sp);
	double GetValenceO2Fraction(int imm, int iwd);
};

#endif