// Mesh.cpp: implementation of the CMesh class.
//
//////////////////////////////////////////////////////////////////////
#include <cstring>		//required for memcpy() function
#include "Mesh.h"

CMesh::CMesh()
{
	nvertex = 8;
	nface = 6;
	ncell = 1;
	nbface = 0;
	pvertex = NULL;
	pface = NULL;
	pcell = NULL;
	pibface = NULL;
}

CMesh::~CMesh()
{
	DeleteArrays();
}

void CMesh::Write(FILE* pf)
{
	int i;
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	fwrite(&nvertex,sizeof(int),1,pf);
	fwrite(&nface,sizeof(int),1,pf);
	fwrite(&ncell,sizeof(int),1,pf);
	fwrite(&nbface,sizeof(int),1,pf);
	fwrite(&mbl,sizeof(double),1,pf);
	if (pvertex && pface && pcell && pibface)	//not null pointers, geometry data available
	{
		i = 1;
		fwrite(&i,sizeof(int),1,pf);
		fwrite(pibface,sizeof(int),nbface,pf);
		for (i=0; i<nvertex; i++)
			pvertex[i].Write(pf);
		for (i=0; i<nface; i++)
			pface[i].Write(pf);
		for (i=0; i<ncell; i++)
			pcell[i].Write(pf);
	}
	else
	{
		i = 0;
		fwrite(&i,sizeof(int),1,pf);
	}
}

void CMesh::Read(FILE* pf)
{
	int i;
	int iversion;
	fread(&iversion,sizeof(int),1,pf);
	fread(&nvertex,sizeof(int),1,pf);
	fread(&nface,sizeof(int),1,pf);
	fread(&ncell,sizeof(int),1,pf);
	fread(&nbface,sizeof(int),1,pf);
	fread(&mbl,sizeof(double),1,pf);
	fread(&i,sizeof(int),1,pf);
	if (i)	//not null pointers, geometry data available
	{
		AllocateArrays();
		fread(pibface,sizeof(int),nbface,pf);
		for (i=0; i<nvertex; i++)
			pvertex[i].Read(pf);
		for (i=0; i<nface; i++)
			pface[i].Read(pf);
		for (i=0; i<ncell; i++)
			pcell[i].Read(pf);
	}
	else	//null pointer
		DeleteArrays();
}

void CMesh::AllocateArrays()
{
	DeleteArrays();
	pvertex = new CVertex [nvertex];
	pface = new CQuadFace [nface];
	pcell = new CHexCell [ncell];
	pibface = new int [nbface];
}

void CMesh::DeleteArrays()
{
	if (pvertex!=NULL)
	{
		delete [] pvertex;
		pvertex = NULL;
	}
	if (pface!=NULL)
	{
		delete [] pface;
		pface = NULL;
	}
	if (pcell!=NULL)
	{
		delete [] pcell;
		pcell = NULL;
	}
	if (pibface!=NULL)
	{
		delete [] pibface;
		pibface = NULL;
	}
}

void CMesh::AddBoundaryFaces(int n, CQuadFace* pf, int* pif_conv)
{
	//n is the number of faces to be added
	//pf is the array of faces to be added
	//pif_conv is the array of indices of faces converted from interior to wall faces
	int i;
	int nface_old = nface;
	int nbface_old = nbface;
	int* pibface_old = pibface;
	CQuadFace* pface_old = pface;
	nface += n;
	nbface += 2*n;
	pibface = new int [nbface];
	pface = new CQuadFace [nface];
	memcpy(pibface, pibface_old, nbface_old*sizeof(int));
	memcpy(pface, pface_old, nface_old*sizeof(CQuadFace));
	delete [] pibface_old;
	for (i=0; i<nface_old; i++)
		pface_old[i].pbfp = NULL;
	delete [] pface_old;
	for (i=0; i<n; i++)
	{
		pibface[nbface_old+i] = nface_old + i;
		pibface[nbface_old+n+i] = pif_conv[i];
	}
	memcpy(&pface[nface_old], pf, n*sizeof(CQuadFace));
	for (i=0; i<n; i++)
		pf[i].pbfp = NULL;
	delete [] pf;
	delete [] pif_conv;
}

void CMesh::UpdateCellGeometry()
{
	int i;
	for (i=0; i<ncell; i++)
	{
		pcell[i].CalcCentroidX(pvertex);
		pcell[i].CalcVolume(pvertex);
	}
}

void CMesh::UpdateFaceGeomerty()
{
	//calculate face normal, face centroid and cell centriod
	int i;
	for (i=0; i<nface; i++)
	{
		pface[i].CalcCentroidX(pvertex);
		pface[i].CalcNormalVectorAndArea(pvertex);
		pface[i].UpdateNormal(pvertex,pcell);
	}
}

void CMesh::CalcZoneVolumes(int nz, double* pvz)
{
	int i;
	for (i=0; i<nz; i++)
		pvz[i] = 0;
	for (i=0; i<ncell; i++)
		pvz[pcell[i].izone] += pcell[i].vol;
}

void CMesh::CalcZoneWallArea(int nz, double* pawz)
{
	int i, j;
	CQuadFace* pf;
	for (i=0; i<nz; i++)
		pawz[i] = 0;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		if (pf->itype==1)		//wall face
			pawz[pf->izone] += pf->area;
	}
}

void CMesh::CalcZoneWallAreaTemperatureProduct(int nz, double* patwz)
{
	int i, j;
	CQuadFace* pf;
	for (i=0; i<nz; i++)
		patwz[i] = 0;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		if (pf->itype==1)		//wall face
			patwz[pf->izone] += pf->area*pf->pbfp->temp;
	}
}

void CMesh::CalcZoneWallIncidentFlux(int nz, double* pavg, double* pmin, double* pmax)
{
	//added to compare the UK test furnace experimental data
	//calculate the average, minimum, and maximum incident heat flux for each zone
	int i, j;
	int iz;
	int izone;
	double x, xmin, xmax, xavg;
	double area, sum_area;
	CQuadFace* pf;
	for (iz=0; iz<nz; iz++)
	{
		xmin = 1e50;
		xmax = 0;
		xavg = 0;
		sum_area = 0;
		for (i=0; i<nbface; i++)
		{
			j = pibface[i];
			pf = &pface[j];
			izone = pf->izone;
			if (izone==iz)
			{
				if (pf->itype==1)	//wall boundary
				{
					if (pf->normal[1]<0.999 && pf->normal[1]>-0.999)	//not bottom or roof
					{
						area = pf->area;
						x = pf->pbfp->qinc;
						sum_area += area;
						xavg += area*x;
						if (x<xmin)
							xmin = x;
						if (x>xmax)
							xmax = x;
					}
				}
			}
		}
		xavg /= sum_area;
		pavg[iz] = xavg;
		pmin[iz] = xmin;
		pmax[iz] = xmax;
	}
}

double CMesh::CalcMBL()
{
	int i, j;
	volume = 0;
	area_encl = 0;
	area_sh = 0;
	area_exit = 0;
	for (i=0; i<ncell; i++)
		volume += pcell[i].vol;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		if (pface[j].itype == 3)	//exit
			area_exit += pface[j].area;
		else	//wall
		{
			if (pface[j].pbfp->iwall>0)		//superheater wall
				area_sh += pface[j].area;
			else		//enclosure wall
				area_encl += pface[j].area;
		}
	}
	area_total = area_encl + area_sh + area_exit;
	mbl = 3.6*volume/area_total;
	return mbl;
}

bool CMesh::IsGeometrySet()
{
	bool bset = pvertex && pface && pcell && pibface;
	return bset;
}