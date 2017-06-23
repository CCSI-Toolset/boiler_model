// QuadFace.cpp: implementation of the CQuadFace class.
//
//////////////////////////////////////////////////////////////////////
#include <cmath>
#include "QuadFace.h"

CQuadFace::CQuadFace()
{
	itype = 0;
	nvertex = 4;
	icell[0] = icell[1] = -1;
	izone = -1;			//for interior face zone index is -1
	pbfp = NULL;
}

CQuadFace::~CQuadFace()
{
	if (pbfp!=NULL)
	{
		delete pbfp;
		pbfp = NULL;
	}
}

void CQuadFace::Write(FILE* pf)
{
	int i;
	fwrite(&itype,sizeof(char),1,pf);
	fwrite(&nvertex,sizeof(short),1,pf);
	fwrite(&izone,sizeof(int),1,pf);
	fwrite(ivertex,sizeof(int),4,pf);
	fwrite(icell,sizeof(int),2,pf);
	fwrite(&area,sizeof(double),1,pf);
	fwrite(x,sizeof(double),3,pf);
	fwrite(normal,sizeof(double),3,pf);
	if (pbfp!=NULL)
	{
		i = 1;
		fwrite(&i,sizeof(int),1,pf);
		pbfp->Write(pf);
	}
	else
	{
		i = 0;
		fwrite(&i,sizeof(int),1,pf);
	}
}

void CQuadFace::Read(FILE* pf)
{
	int i;
	fread(&itype,sizeof(char),1,pf);
	fread(&nvertex,sizeof(short),1,pf);
	fread(&izone,sizeof(int),1,pf);
	fread(ivertex,sizeof(int),4,pf);
	fread(icell,sizeof(int),2,pf);
	fread(&area,sizeof(double),1,pf);
	fread(x,sizeof(double),3,pf);
	fread(normal,sizeof(double),3,pf);
	fread(&i,sizeof(int),1,pf);
	if (pbfp)		//delete in case the current case contains object
			delete pbfp;
	if (i)
	{
		
		pbfp = new CBoundaryFaceProperty();
		pbfp->Read(pf);
	}
	else
		pbfp = NULL;
}

void CQuadFace::CopyPartialData(CQuadFace* pf)
{
	int i;
	itype = pf->itype;
	izone = pf->izone;
	for (i=0; i<nvertex; i++)
		ivertex[i] = pf->ivertex[i];
}

void CQuadFace::CalcCentroidX(CVertex* pvertex)
{
	//calculate face centroid from vertices
	int i;
	short ii, jj;
	x[0] = x[1] = x[2] = 0;
	for (ii=0; ii<nvertex; ii++)
	{
		i = ivertex[ii];
		for (jj=0; jj<3; jj++)
			x[jj] += pvertex[i].x[jj];
	}
	for (jj=0; jj<3; jj++)
		x[jj] /= (double)nvertex;
}

void CQuadFace::CalcNormalVectorAndArea(CVertex* pvertex)
{
	//calculate normal vecter using the first 3 points and face area
	short ii;
	int i0, i1, i2;
	double a[3], b[3];
	double c0, c1, c2, c;
	i0 = ivertex[0];
	i1 = ivertex[1];
	i2 = ivertex[2];
	for (ii=0; ii<3; ii++)
	{
		a[ii] = pvertex[i1].x[ii]-pvertex[i0].x[ii];
		b[ii] = pvertex[i2].x[ii]-pvertex[i0].x[ii];
	}
	c0 = a[1]*b[2] - b[1]*a[2];
	c1 = a[2]*b[0] - b[2]*a[0];
	c2 = a[0]*b[1] - b[0]*a[1];
	c = (double)sqrt(c0*c0+c1*c1+c2*c2);
	normal[0] = c0/c;
	normal[1] = c1/c;
	normal[2] = c2/c;
	area = c/2;
	if (nvertex>3)			//quard
	{
		i1 = ivertex[2];
		i2 = ivertex[3];
		for (ii=0; ii<3; ii++)
		{
			a[ii] = pvertex[i1].x[ii]-pvertex[i0].x[ii];
			b[ii] = pvertex[i2].x[ii]-pvertex[i0].x[ii];
		}
		c0 = a[1]*b[2] - b[1]*a[2];
		c1 = a[2]*b[0] - b[2]*a[0];
		c2 = a[0]*b[1] - b[0]*a[1];
		area += (double)sqrt(c0*c0+c1*c1+c2*c2)/2;
	}
}

void CQuadFace::UpdateNormal(CVertex* pvertex, CHexCell* pcell)
{
	//this function switch the normal direction if it is opposite to (C1-C0) or (F-C0)
	//this ensures that the normal face vector is pointing to outside of the boundary if it is a boundary face
	short ii;
	int i0, i1;
	double l = 0;
	if (icell[1]<0)		//boundary face
	{
		i0 = icell[0];
		for (ii=0; ii<3; ii++)
			l += (x[ii] - pcell[i0].x[ii])*normal[ii];
	}
	else
	{
		i0 = icell[0];
		i1 = icell[1];
		for (ii=0; ii<3; ii++)
			l += (pcell[i1].x[ii] - pcell[i0].x[ii])*normal[ii];
	}
	if (l<0)	//if normal direction opposite of (C1-C0), switch normal direction
	{
		for (ii=0; ii<3; ii++)
			normal[ii] = -normal[ii];
	}
}

double CQuadFace::GetWallTemperature()
{
	if (itype==1) return pbfp->temp;
	return 0;
}

double CQuadFace::GetConvectiveHeatFlux()
{
	if (itype==1) return pbfp->qconv;
	return 0;
}

double CQuadFace::GetRadiationIncidentHeatFlux()
{
	if (itype) return pbfp->qinc;
	return 0;
}

double CQuadFace::GetRadiationNetHeatFlux()
{
	if (itype) return pbfp->qnet;
	return 0;
}

double CQuadFace::GetRadiationBoundaryFactor()
{
	if (itype) return pbfp->fradadj;
	return 0;
}
