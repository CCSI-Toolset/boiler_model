// PSD.cpp: implementation of the CPSD class.
//
//////////////////////////////////////////////////////////////////////

#include "PSD.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPSD::CPSD()
{
	int i;
	nbin = 11;
	for (i=0; i<PSD_NBIN; i++)
	{
		mf[i]=0;
		adp[i]=5e-5;
	}
	mf[0] = 0.025;
	mf[1] = 0.05;
	mf[2] = 0.075;
	mf[3] = 0.1;
	mf[4] = 0.15;
	mf[5] = 0.2;
	mf[6] = 0.15;
	mf[7] = 0.1;
	mf[8] = 0.075;
	mf[9] = 0.05;
	mf[10] = 0.025;
	denp = 1350;
}

void CPSD::Write(FILE* pf)
{
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	fwrite(&nbin,sizeof(int),1,pf);
	fwrite(adp,sizeof(double),nbin,pf);
	fwrite(mf,sizeof(double),nbin,pf);
	fwrite(&denp,sizeof(double),1,pf);
}

void CPSD::Read(FILE* pf)
{
	int iversion;
	fread(&iversion,sizeof(int),1,pf);
	fread(&nbin,sizeof(int),1,pf);
	fread(adp,sizeof(double),nbin,pf);
	fread(mf,sizeof(double),nbin,pf);
	fread(&denp,sizeof(double),1,pf);
}

CPSD::CPSD(const CPSD &t)
{
	int i;
	nbin = t.nbin;	
	for (i=0; i<nbin; i++)
	{
		adp[i]=t.adp[i];
		mf[i]=t.mf[i];
	}
}


CPSD& CPSD::operator=(const CPSD& t)
{
	if (this==&t)
		return *this;
	int i;
	nbin = t.nbin;	
	for (i=0; i<nbin; i++)
	{
		adp[i]=t.adp[i];
		mf[i]=t.mf[i];
	}
	return *this;
}

double CPSD::GetMfSum()
{
	int i;
	double sum = 0;
	for (i=0; i<nbin; i++)
		sum += mf[i];
	return sum;
}
