// PSD.h: header file for the CPSD class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __PSD_H__
#define __PSD_H__
#include <cstdio>

const int PSD_NBIN = 50;	//maximum number of size bins

class CPSD
{
public:
	int nbin;				//number of size bins
	double adp[PSD_NBIN];	//average diameter in for each size bin [m]
	double mf[PSD_NBIN];	//mass fraction for each size bin
	double denp;			//particle density, could be moved to coal property

	//functions
	CPSD();
	~CPSD() {}
	CPSD (const CPSD &t);
	CPSD& operator=(const CPSD& t);
	double GetMfSum();
	void Write(FILE* pf);
	void Read(FILE* pf);
};

#endif