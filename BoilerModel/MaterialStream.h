// MaterialStream.h: interface for the CMaterialStream class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __MATERIALSTREAM_H__
#define __MATERIALSTREAM_H__

#include "IdealGasMixture.h"
#include "LiquidFuel.h"
#include "SolidFuel.h"
#include <vector>
#include <map>

class CMaterialStream
{
public:
	std::string name;		//name of the stream
	bool bgas;				//contains gas phase
	bool bliquid;			//contains liquid phase
	bool bsolid;			//contains solid phase
	double temp;			//temperature of the stream [K]
	double pres;			//pressure of the stream [Pa]
	double er;				//equivalence ratio of the stream
	double mflow;			//total mass flow of the stream [kg/s]
	double hflow;			//total enthalpy flow of the stream [W]

	CIdealGasMixture gas;	//gas phase object
	CLiquidFuel oil;		//liquid phase object
	CSolidFuel coal;		//solid phase object
	
	//functions
	CMaterialStream();
	virtual ~CMaterialStream();
	CMaterialStream(const CMaterialStream &t);
	CMaterialStream operator=(CMaterialStream t);
	void Write(FILE* pf);
	void Read(FILE* pf);

	void CopyProperties(CMaterialStream* ps);
	void UpdateProperties();
	void ModifyMassFlow(double mf, int iphase);
	double CalcSFlow();
	double CalcCp();
	double GetHHVRate();

};

#endif