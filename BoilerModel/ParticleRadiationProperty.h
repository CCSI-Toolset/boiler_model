// ParticleRadiationProperty.h: interface for the CParticleRadiationProperty class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __PAETICLERADIATIONPROPERTY_H__
#define __PAETICLERADIATIONPROPERTY_H__

const int RP_NWL = 27;						//number of wavelengths
const int RP_NPD = 19;						//number of diameters
const int RP_NPT = 14;						//number of temperatures

class CParticleRadiationProperty  
{
private:
	static double refre_coal;				//real part of refrective index of coal
	static double refim_coal;				//imaginary part of refrective index of coal
	static double refre_ash[RP_NWL];		//real part of refrective index of ash, dependent on wavelength
	static double refim_ash[3][RP_NWL];		//imaginary part of refrective index of ash, dependent on wavelength and temperature
	static double pdiam[RP_NPD];			//discretized particle diameter
	static double ptemp[RP_NPT];			//discretized particle temperature
	static CParticleRadiationProperty* pinstance;
	//instance variables
	double qa_coal[RP_NPD][RP_NPT];	//absorption efficiency table for coal
	double qa_ash[RP_NPD][RP_NPT];	//absorption efficiency table for ash

	
	CParticleRadiationProperty();
	virtual ~CParticleRadiationProperty();
	
	void Mie(double x, double refrelr, double refreli, int nang, double& qext, double& qsca, double& qback);
	double Planck(double temp, double* q);
	double Rosseland(double temp, double* q);
	void SetUpCoalEfficiencyTable();
	void SetUpAshEfficiencyTable();
	int FindLowerIndexAndInterpolationFactor(int nsize, double x, double* px, double& f);

public:	
	static CParticleRadiationProperty* GetInstance();
	static void DeleteInstance();
	double GetCoalAbsorptionEfficiency(double d, double t);
	double GetAshAbsorptionEfficiency(double d, double t);

};

#endif
