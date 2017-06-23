// HexCell.cpp: implementation of the CHexCell class.
//
//////////////////////////////////////////////////////////////////////

//my include
#include <cmath>
//end of my include

#include "HexCell.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CHexCell::CHexCell()
{
	nvertex = 8;
	nface = 0;			//used for count
	ncell_nbr = 0;		//used for count
	izone = 0;
}

CHexCell::~CHexCell()
{
}

void CHexCell::Write(FILE* pf)
{
	fwrite(&nvertex,sizeof(short),1,pf);
	fwrite(&nface,sizeof(short),1,pf);
	fwrite(&ncell_nbr,sizeof(short),1,pf);
	fwrite(&izone,sizeof(int),1,pf);
	fwrite(pivertex,sizeof(int),8,pf);
	fwrite(piface,sizeof(int),6,pf);
	fwrite(picell,sizeof(int),6,pf);
	fwrite(x,sizeof(double),3,pf);
	fwrite(&vol,sizeof(double),1,pf);
}

void CHexCell::Read(FILE* pf)
{
	fread(&nvertex,sizeof(short),1,pf);
	fread(&nface,sizeof(short),1,pf);
	fread(&ncell_nbr,sizeof(short),1,pf);
	fread(&izone,sizeof(int),1,pf);
	fread(pivertex,sizeof(int),8,pf);
	fread(piface,sizeof(int),6,pf);
	fread(picell,sizeof(int),6,pf);
	fread(x,sizeof(double),3,pf);
	fread(&vol,sizeof(double),1,pf);
}

void CHexCell::AddFaceAndNeighborCell(int iface, int icell)
{
	//note that the interior face should be added first followed by boundary faces
	piface[nface] = iface;
	nface++;
	if (icell>-1)
	{
		picell[ncell_nbr] = icell;
		ncell_nbr++;
	}
}

void CHexCell::CalcCentroidX(CVertex* pvertex)
{
	short ii, jj;
	int i;
	x[0] = x[1] = x[2] = 0;
	for (ii=0; ii<nvertex; ii++)
	{
		i = pivertex[ii];
		for (jj=0; jj<3; jj++)
			x[jj] += pvertex[i].x[jj];
	}
	x[0] /= (double)nvertex;
	x[1] /= (double)nvertex;
	x[2] /= (double)nvertex;
}

double CHexCell::GetTetrahedronVolume(short iv[4], CVertex* pvertex)
{
	short ii;
	int i0, i1, i2, i3;
	double vol4;
	double a[3], b[3], c[3];
	i0 = pivertex[iv[0]];
	i1 = pivertex[iv[1]];
	i2 = pivertex[iv[2]];
	i3 = pivertex[iv[3]];
	for (ii=0; ii<3; ii++)
	{
		a[ii] = pvertex[i1].x[ii]-pvertex[i0].x[ii];
		b[ii] = pvertex[i2].x[ii]-pvertex[i0].x[ii];
		c[ii] = pvertex[i3].x[ii]-pvertex[i0].x[ii];
	}
	vol4 = (a[0]*b[1]*c[2]+a[1]*b[2]*c[0]+a[2]*b[0]*c[1]-a[2]*b[1]*c[0]-a[1]*b[0]*c[2]-a[0]*b[2]*c[1])/6;
	vol4 = fabs(vol4);
	return vol4;
}

void CHexCell::CalcVolume(CVertex* pvertex)
{
	//vertex numbering order based on FieldView specifications
	short iv[4];
	iv[0] = 0;
	iv[1] = 1;
	iv[2] = 2;
	iv[3] = 4;
	vol = GetTetrahedronVolume(iv,pvertex);
	iv[0] = 3;
	iv[1] = 1;
	iv[2] = 2;
	iv[3] = 7;
	vol += GetTetrahedronVolume(iv,pvertex);
	iv[0] = 5;
	iv[1] = 1;
	iv[2] = 4;
	iv[3] = 7;
	vol += GetTetrahedronVolume(iv,pvertex);
	iv[0] = 6;
	iv[1] = 2;
	iv[2] = 4;
	iv[3] = 7;
	vol += GetTetrahedronVolume(iv,pvertex);
	iv[0] = 1;
	iv[1] = 2;
	iv[2] = 4;
	iv[3] = 7;
	vol += GetTetrahedronVolume(iv,pvertex);
	//debug
	if (vol<=0)
		printf("negative volume");
}
