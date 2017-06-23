// SuperHeater.h: interface for the CSuperHeater class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __SUPERHEATER_H__
#define __SUPERHEATER_H__

const int SH_NP = 50;		//maximum number of polygon points

#include <cstdio>
#include <string>

class CSuperHeater
{
public:
	std::string name;		//name of the Shape
	int iwall;				//wall index
	int	iwater;				//water stream index
	int npanel;				//number of panels along the furnace width
	int npoints;			//number of points for polygon shape
	double xp[2*SH_NP];		//50 2-D x points for general shape
	double yp[2*SH_NP];		//50 2-D y points for general shape

public:
	CSuperHeater();
	virtual ~CSuperHeater();
	CSuperHeater(const CSuperHeater &tt);
	CSuperHeater& operator=(const CSuperHeater& tt);
	void UpdateExtraPoint();
	bool IsInsidePolygon(double x, double y);
	double CalcArea();
	void Write(FILE* pf);
	void Read(FILE* pf);
};

#endif