// SuperHeater.cpp: implementation of the CSuperHeater class.
//
//////////////////////////////////////////////////////////////////////

#include <cmath>
#include "SuperHeater.h"
#include "UtilityFunctions.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CSuperHeater::CSuperHeater()
{
	int i;
	name = "SH";
	iwall = 1;
	iwater = 1;
	npanel = 1;
	npoints = 4;
	for (i=0; i<SH_NP; i++)
	{
		xp[i] = 0;
		yp[i] = 0;
	}
	xp[0] = 0;
	yp[0] = 0;
	xp[1] = 5;
	yp[1] = 0;
	xp[2] = 5;
	yp[2] = 5;
	xp[3] = 0;
	yp[3] = 5;
	UpdateExtraPoint();
}

CSuperHeater::~CSuperHeater()
{
}

double CSuperHeater::CalcArea()
{
	//use cross product to calculate area from the first point
	int i;
	double x1, x2, y1, y2;
	double area = 0;
	for (i=2; i<npoints; i++)
	{
		x1 = xp[i-1] - xp[0];
		x2 = xp[i] - xp[0];
		y1 = yp[i-1] - yp[0];
		y2 = yp[i] - yp[0];
		area += x1*y2 - x2*y1;
	}
	area = fabs(area)*npanel/2;
	return area;
}

void CSuperHeater::Write(FILE* pf)
{
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::WriteString(name,pf);
	fwrite(&iwall,sizeof(int),1,pf);
	fwrite(&iwater,sizeof(int),1,pf);
	fwrite(&npanel,sizeof(int),1,pf);
	fwrite(&npoints,sizeof(int),1,pf);
	fwrite(xp,sizeof(double),npoints,pf);
	fwrite(yp,sizeof(double),npoints,pf);
}

void CSuperHeater::Read(FILE* pf)
{
	int iversion;
	fread(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::ReadString(name,pf);
	fread(&iwall,sizeof(int),1,pf);
	fread(&iwater,sizeof(int),1,pf);
	fread(&npanel,sizeof(int),1,pf);
	fread(&npoints,sizeof(int),1,pf);
	fread(xp,sizeof(double),npoints,pf);
	fread(yp,sizeof(double),npoints,pf);
}

CSuperHeater::CSuperHeater(const CSuperHeater &tt)
{
	int i;
	name = tt.name;
	iwall = tt.iwall;
	iwater = tt.iwater;
	npanel = tt.npanel;
	npoints = tt.npoints;
	for (i=0; i<2*npoints; i++)
	{
	  xp[i] = tt.xp[i];
	  yp[i] = tt.yp[i];
	}
}


CSuperHeater& CSuperHeater::operator=(const CSuperHeater& tt)
{
	if (this==&tt)
		return *this;
	int i;
	name = tt.name;
	iwall = tt.iwall;
	iwater = tt.iwater;
	npanel = tt.npanel;
	npoints = tt.npoints;
	for (i=0; i<2*npoints; i++)
	{
	  xp[i] = tt.xp[i];
	  yp[i] = tt.yp[i];
	}
	return *this;
}

void CSuperHeater::UpdateExtraPoint()
{
	int i;
	for (i=0; i<npoints; i++)
	{
		xp[i+npoints] = xp[i];
		yp[i+npoints] = yp[i];
	}
}

bool CSuperHeater::IsInsidePolygon(double x, double y)
{
	int i,j,k;
	int istart;		//starting array index
	int nabove;		//number of points above
	bool bleft;		//is from left side
	bool bhit;		//hit the point
	double yi;		//intersection y
	nabove = 0;
	//first find a point that is not collide with x
	istart = -1;
	for (i=0; i<npoints; i++)
	{
		if (x!=xp[i])
		{
			istart = i;
			break;
		}
	}
	if (istart==-1) return false;
	bhit = false;
	if (x>xp[istart])
		bleft = true;
	else
		bleft = false;
	for (i=0; i<npoints; i++)
	{
		j = i+istart;
		k = j+1;
		if (bhit)
		{
			if (x==xp[k])		//check if point on line
			{
				if ((y-yp[j])*(y-yp[k])<=0) return true;
			}
			else
			{
				bhit = false;
				if (bleft)
				{
					if (x<xp[k])
					{
						bleft = false;
						if (y<yp[j]) nabove ++;
					}
				}
				else
				{
					if (x>xp[k])
					{
						bleft = true;
						if (y<yp[j]) nabove++;
					}
				}
			}
		}
		else
		{
			if (x==xp[k])
				bhit = true;
			else
			{
				if (bleft)
				{
					if (x<xp[k])
					{
						bleft = false;
						yi = yp[j] + (x-xp[j])/(xp[k]-xp[j])*(yp[k]-yp[j]);
						if (y==yi) return true;
						if (y<yi) nabove++;
					}
				}
				else
				{
					if (x>xp[k])
					{
						bleft = true;
						yi = yp[j] + (x-xp[j])/(xp[k]-xp[j])*(yp[k]-yp[j]);
						if (y==yi) return true;
						if (y<yi) nabove++;
					}
				}
			}
		}
	}
	if (nabove%2) return true;
	return false;
}
