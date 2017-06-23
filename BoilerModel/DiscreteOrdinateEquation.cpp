// DiscreteOrdinateEquation.cpp: implementation of the CDiscreteOrdinateEquation class.
//
//////////////////////////////////////////////////////////////////////

#include <cmath>
#include "Constants.h"
#include "UtilityFunctions.h"
#include "DiscreteOrdinateEquation.h"

CDiscreteOrdinateEquation::CDiscreteOrdinateEquation()
{
	ndos = 20;				//based on uniform dodecahedron
	nite_conv = 0;
	urf_twall = 0.95;
	urf_srcrad = 0.4;
	errmax = 0;
	errconv = 5e-6;			//old value is 1e-5;
	fqinc = 1;
	pwt = NULL;
	pdcos = NULL;
	pint = NULL;
	parea = NULL;
	pkag = NULL;
	pkat4g = NULL;
	pkap = NULL;
	pkat4p = NULL;
	pgir = NULL;
	pqdivg = NULL;
	pqdivp = NULL;
	pisc = NULL;
	pmesh = NULL;
}

CDiscreteOrdinateEquation::~CDiscreteOrdinateEquation()
{
	DeleteMemory();
}

void CDiscreteOrdinateEquation::AllocateMemory()
{
	int i;
	int ndos2;				//half of number of discrete ordinates
	DeleteMemory();			//delete first in case of updating from previous data
	ndos2 = ndos/2;
	int ncell = pmesh->ncell;
	int nface = pmesh->nface;
	pwt = new double [ndos];
	pdcos = new double [ndos][3];
	pint = new double* [ndos];
	for (i=0; i<ndos; i++)
		pint[i] = new double [ncell];
	parea = new double [nface];
	pkag = new double [ncell];
	pkat4g = new double [ncell];
	pkap = new double [ncell];
	pkat4p = new double [ncell];
	pgir = new double [ncell];
	pqdivg = new double [ncell];
	pqdivp = new double [ncell];
	pisc = new int* [ndos2];
	for (i=0; i<ndos2; i++)
		pisc[i] = new int [ncell];
}

void CDiscreteOrdinateEquation::DeleteMemory()
{
	int i;
	int ndos2 = ndos/2;
	if (pwt!=NULL)
	{
		delete [] pwt;
		pwt = NULL;
	}
	if (pdcos!=NULL)
	{
		delete [] pdcos;
		pdcos = NULL;
	}
	if (pint!=NULL)
	{
		for (i=0; i<ndos; i++)
			delete [] pint[i];
		delete [] pint;
		pint = NULL;
	}
	if (parea!=NULL)
	{
		delete [] parea;
		parea = NULL;
	}
	if (pkag!=NULL)
	{
		delete [] pkag;
		pkag = NULL;
	}
	if (pkat4g!=NULL)
	{
		delete [] pkat4g;
		pkat4g = NULL;
	}
	if (pkap!=NULL)
	{
		delete [] pkap;
		pkap = NULL;
	}
	if (pkat4p!=NULL)
	{
		delete [] pkat4p;
		pkat4p = NULL;
	}
	if (pgir!=NULL)
	{
		delete [] pgir;
		pgir = NULL;
	}
	if (pqdivg!=NULL)
	{
		delete [] pqdivg;
		pqdivg = NULL;
	}
	if (pqdivp!=NULL)
	{
		delete [] pqdivp;
		pqdivp = NULL;
	}
	if (pisc!=NULL)
	{
		for (i=0; i<ndos2; i++)
			delete [] pisc[i];
		delete [] pisc;
		pisc = NULL;
	}
}

void CDiscreteOrdinateEquation::UpdateDirectionCosines()
{
	//ndos has to be an even number
	//direction i should be opposite to direction i+ndos/2
	//specify +z direction first, -z direction is calculated from +z direction through a for loop
	int i, j;
	int ndos2 = ndos/2;
	double r = sqrt(3.0);
	double phi = (1+sqrt(5.0))/2;
	double r1 = 1/r;
	double r2 = phi/r;
	double r3 = 1/phi/r;
	double wt1;
	//20 vertices of a dodecahedron
	wt1 = CT_PI/5;
	for (i=0; i<ndos; i++)
		pwt[i] = wt1;
	pdcos[0][0] = r1;
	pdcos[0][1] = r1;
	pdcos[0][2] = r1;
	pdcos[1][0] = -r1;
	pdcos[1][1] = r1;
	pdcos[1][2] = r1;
	pdcos[2][0] = r1;
	pdcos[2][1] = -r1;
	pdcos[2][2] = r1;
	pdcos[3][0] = -r1;
	pdcos[3][1] = -r1;
	pdcos[3][2] = r1;
	pdcos[4][0] = 0;
	pdcos[4][1] = r3;
	pdcos[4][2] = r2;
	pdcos[5][0] = 0;
	pdcos[5][1] = -r3;
	pdcos[5][2] = r2;
	pdcos[6][0] = r2;
	pdcos[6][1] = 0;
	pdcos[6][2] = r3;
	pdcos[7][0] = -r2;
	pdcos[7][1] = 0;
	pdcos[7][2] = r3;
	pdcos[8][0] = r3;
	pdcos[8][1] = r2;
	pdcos[8][2] = 0;
	pdcos[9][0] = -r3;
	pdcos[9][1] = r2;
	pdcos[9][2] = 0;
	for (i=0; i<ndos2; i++)
	{
		j = i+ndos2;
		pdcos[j][0] = -pdcos[i][0];
		pdcos[j][1] = -pdcos[i][1];
		pdcos[j][2] = -pdcos[i][2];
	}
}

void CDiscreteOrdinateEquation::UpdateSortedCellIndices()
{
	int i;
	int ii;
	int ndos2 = ndos/2;
	double* pxc;								//cell centroid location
	double* pdc;								//direction cosine
	double* pxdo = new double [pmesh->ncell];	//projected location in discrete ordinate direction
	CHexCell* pcell = pmesh->pcell;
	int ncell = pmesh->ncell;
	for (ii=0; ii<ndos2; ii++)
	{
		pdc = pdcos[ii];
		for (i=0; i<ncell; i++)
		{
			pxc = pcell[i].x;
			pxdo[i] = pxc[0]*pdc[0] + pxc[1]*pdc[1] + pxc[2]*pdc[2];
			pisc[ii][i] = i;
		}
		CUtilityFunctions::QuickSortWithInsertion(pxdo,pisc[ii],0,ncell-1);
	}
	delete [] pxdo;
}

void CDiscreteOrdinateEquation::CalcProjectedFaceArea(int ido)
{
	//ido is the index of the discrete ordinate
	//parea is array that contains the calculated projected area of each face in ido direction (contains sign)
	int i;
	double* pnormal;
	CQuadFace* pf;
	CQuadFace* pface = pmesh->pface;
	int nface = pmesh->nface;
	for (i=0; i<nface; i++)
	{
		pf = &pface[i];
		pnormal = pf->normal;
		parea[i] = pf->area*(pnormal[0]*pdcos[ido][0] + pnormal[1]*pdcos[ido][1] +pnormal[2]*pdcos[ido][2]);
	}
}

void CDiscreteOrdinateEquation::UpdateBoundaryFaceFactor()
{
	int ii;
	int i, j;
	double proj;
	double sum;
	double* pnormal;
	CQuadFace* pf;
	int nbface = pmesh->nbface;
	int*	pibface = pmesh->pibface;
	CQuadFace* pface = pmesh->pface;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		pnormal = pf->normal;
		sum = 0;
		for (ii=0; ii<ndos; ii++)
		{
			proj = pnormal[0]*pdcos[ii][0] + pnormal[1]*pdcos[ii][1] + pnormal[2]*pdcos[ii][2];
			if (proj>0)
				sum += proj*pwt[ii];
		}
		pf->pbfp->fradadj = CT_PI/sum;
	}
}

bool CDiscreteOrdinateEquation::CalcBoundaryIncidentHeatFlux()
{
	int ii;
	int i, j, k;
	double proj;
	double sum;
	double err;
	double* pnormal;
	CQuadFace* pf;
	int nbface = pmesh->nbface;
	int*	pibface = pmesh->pibface;
	CQuadFace* pface = pmesh->pface;
	errmax = 0;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		k = pf->icell[0];
		pnormal = pf->normal;
		sum = 0;
		for (ii=0; ii<ndos; ii++)
		{
			proj = pnormal[0]*pdcos[ii][0] + pnormal[1]*pdcos[ii][1] + pnormal[2]*pdcos[ii][2];
			if (proj>0)
				sum += proj*pwt[ii]*pint[ii][k];
		}
		sum *= pf->pbfp->fradadj;
		err = fabs(pf->pbfp->qinc - sum);
		if (sum>0)
			err /= sum;
		if (err>errmax)
			errmax = err;
		pf->pbfp->qinc = sum;
	}
	return errmax<errconv;
}

void CDiscreteOrdinateEquation::Init()
{
	int ii;
	int i;
	int ncell = pmesh->ncell;
	//initialize pint to zero, pisc initialized in UpdateSortedCellIndices()
	for (ii=0; ii<ndos; ii++)
	{
		for (i=0; i<ncell; i++)
			pint[ii][i] = 1e-20;		//seting to a very small positive number
	}
	for (i=0; i<ncell; i++)
	{
		pkag[i] = 0;
		pkat4g[i] = 0;
		pkap[i] = 0;
		pkat4p[i] = 0;
		pgir[i] = 0;
		pqdivg[i] = 0;
		pqdivp[i] = 0;
	}
	UpdateDirectionCosines();
	UpdateSortedCellIndices();
	UpdateBoundaryFaceFactor();
	bfirst = true;
}

void CDiscreteOrdinateEquation::SolvePDE()
{
	bool bconv;				//convergence flag
	int nite = 0;
	int ii, ii1;
	int i0, i0max;			//direction index and direction index with maximum dot product
	int ndos2 = ndos/2;
	short jj;
	short ncell_nbr;		//number of neighbors in a cell
	short nface;			//number of faces in a cell
	int i, j, k;
	double area;			//projected area
	double sum_area;		//sum of projected area
	double sum_int;			//sum of indensity
	double dotp;			//dot product of two vectors
	double dotp_max;		//maximum of dotp
	double vect[3];			//vector of incoming direction for symmetry face
	double* pfnormal;		//face normal vector
	CHexCell* pc;
	CQuadFace* pf;
	CBoundaryFaceProperty* pbfp;
	int ncell = pmesh->ncell;
	CHexCell* pcell = pmesh->pcell;
	CQuadFace* pface = pmesh->pface;
	do
	{
		nite++;
		for (ii=0; ii<ndos2; ii++)
		{
			ii1 = ndos2+ii;
			CalcProjectedFaceArea(ii);
			//positive z direction
			for (i=0; i<ncell; i++)
			{
				j = pisc[ii][i];
				pc = &pcell[j];
				ncell_nbr = pc->ncell_nbr;
				nface = pc->nface;
				sum_area = 0;
				sum_int = 0;
				for (jj=0; jj<ncell_nbr; jj++)		//interior neighbors
				{
					k = pc->piface[jj];
					pf = &pface[k];
					if (pf->icell[0]==j)			//from cell j to outside
						area = parea[k];
					else
						area = -parea[k];
					if (area<0)						//upstream face
					{
						sum_area -= area;
						sum_int -= area*pint[ii][pc->picell[jj]];
					}
				}
				for (jj=ncell_nbr; jj<nface; jj++)	//boundary faces
				{
					k = pc->piface[jj];
					area = parea[k];
					if (area<0)		//upstream face
					{
						sum_area -= area;
						pf = &pface[k];
						pbfp = pf->pbfp;
						switch(pf->itype)
						{
						case 1:		//wall boundary
							sum_int -= area*(pbfp->qinc*(1-pbfp->emis)+pbfp->emis*CT_RSIGMA*pow(pbfp->temp,4))/CT_PI;
							break;
						case 2:		//inlet boundary
							sum_int -= area*(pbfp->qinc*(1-pbfp->emis)+pbfp->emis*CT_RSIGMA*pow(pbfp->temp,4))/CT_PI;
							break;
						case 3:		//outlet boundary
							sum_int -= area*(pbfp->qinc*(1-pbfp->emis)+pbfp->emis*CT_RSIGMA*pow(pbfp->temp,4))/CT_PI;
							break;
						case 4:		//symmetry boundary
							//better to interpolate intensity, currently use the closest ray
							//calculate the vector of incoming ray before reflection
							pfnormal = pf->normal;
							dotp = pfnormal[0]*pdcos[ii][0] + pfnormal[1]*pdcos[ii][1] + pfnormal[2]*pdcos[ii][2];
							vect[0] = pdcos[ii][0] - 2*dotp*pfnormal[0];
							vect[1] = pdcos[ii][1] - 2*dotp*pfnormal[1];
							vect[2] = pdcos[ii][2] - 2*dotp*pfnormal[2];
							//find the best match of the direction cosine (maximum dot product)
							dotp_max = 0;
							i0max = 0;
							for (i0=0; i0<ndos; i0++)
							{
								dotp = vect[0]*pdcos[i0][0] + vect[1]*pdcos[i0][1] + vect[2]*pdcos[i0][2];
								if (dotp>dotp_max)
								{
									dotp_max = dotp;
									i0max = i0;
								}
							}
							//debug
							if (dotp_max<0.95)
								printf("Failed to find incoming direction");
							sum_int -= area*pint[i0max][j];
							break;
						}
					}
				}
				//need in-phase scattering term if particles are scattering, currently assume no scattering or 100% foward scattering
				pint[ii][j] = (sum_int+CT_RSIGMA*(pkat4g[j]+pkat4p[j])/CT_PI*pc->vol)/((pkag[j]+pkap[j])*pc->vol+sum_area);
			}
			//negative z direction
			for (i=ncell-1; i>=0; i--)
			{
				j = pisc[ii][i];
				pc = &pcell[j];
				ncell_nbr = pc->ncell_nbr;
				nface = pc->nface;
				sum_area = 0;
				sum_int = 0;
				for (jj=0; jj<ncell_nbr; jj++)		//interior neighbors
				{
					k = pc->piface[jj];
					pf = &pface[k];
					if (pf->icell[0]==j)			//from cell j to outside
						area = parea[k];
					else
						area = -parea[k];
					if (area>0)						//upstream face
					{
						sum_area += area;
						sum_int += area*pint[ii1][pc->picell[jj]];
					}
				}
				for (jj=ncell_nbr; jj<nface; jj++)	//boundary faces
				{
					k = pc->piface[jj];
					area = parea[k];
					if (area>0)		//upstream face
					{
						sum_area += area;
						pf = &pface[k];
						pbfp = pf->pbfp;
						switch(pf->itype)
						{
						case 1:		//wall boundary
							sum_int += area*(pbfp->qinc*(1-pbfp->emis)+pbfp->emis*CT_RSIGMA*pow(pbfp->temp,4))/CT_PI;
							break;
						case 2:		//inlet boundary
							sum_int += area*(pbfp->qinc*(1-pbfp->emis)+pbfp->emis*CT_RSIGMA*pow(pbfp->temp,4))/CT_PI;
							break;
						case 3:		//outlet boundary
							sum_int += area*(pbfp->qinc*(1-pbfp->emis)+pbfp->emis*CT_RSIGMA*pow(pbfp->temp,4))/CT_PI;
							break;
						case 4:		//symmetry boundary
							//better to interpolate intensity, currently use the closest ray
							//calculate the vector of incoming ray before reflection
							pfnormal = pf->normal;
							dotp = pfnormal[0]*pdcos[ii1][0] + pfnormal[1]*pdcos[ii1][1] + pfnormal[2]*pdcos[ii1][2];
							vect[0] = pdcos[ii1][0] - 2*dotp*pfnormal[0];
							vect[1] = pdcos[ii1][1] - 2*dotp*pfnormal[1];
							vect[2] = pdcos[ii1][2] - 2*dotp*pfnormal[2];
							//find the best match of the direction cosine (maximum dot product)
							dotp_max = 0;
							i0max = 0;
							for (i0=0; i0<ndos; i0++)
							{
								dotp = vect[0]*pdcos[i0][0] + vect[1]*pdcos[i0][1] + vect[2]*pdcos[i0][2];
								if (dotp>dotp_max)
								{
									dotp_max = dotp;
									i0max = i0;
								}
							}
							//debug
							if (dotp_max<0.95)
								printf("Failed to find incoming direction");
							sum_int += area*pint[i0max][j];
							break;
						}
					}
				}
				//need in-phase scattering term if particles are involved
				pint[ii1][j] = (sum_int+CT_RSIGMA*(pkat4g[j]+pkat4p[j])/CT_PI*pc->vol)/((pkag[j]+pkap[j])*pc->vol+sum_area);
			}
		}
		bconv = CalcBoundaryIncidentHeatFlux();
	}while (!bconv && nite<15);
	nite_conv = nite;
	//finally prepare source term
	CalcFluxDivergence();
	UpdateWallTemperature();
}

void CDiscreteOrdinateEquation::CalcFluxDivergence()
{
	//pqdiv*vol should be the heat loss term due to radiation
	//under-relax if not first call after Init()
	int ii;
	int i;
	double girrad;			//irradiance
	double qdivg;			//gas divergence
	double qdivp;			//particle divergence
	double sigma4 = 4*CT_RSIGMA;
	double urf_srcrad1 = 1 - urf_srcrad;
	int ncell = pmesh->ncell;
	CHexCell* pcell = pmesh->pcell;
	qcellg = 0;
	qcellp = 0;
	for (i=0; i<ncell; i++)
	{
		girrad = 0;
		for (ii=0; ii<ndos; ii++)
			girrad += pint[ii][i]*pwt[ii];
		pgir[i] = girrad;
		qdivg = girrad*pkag[i] - sigma4*pkat4g[i];
		qdivp = girrad*pkap[i] - sigma4*pkat4p[i];
		qcellg += qdivg*pcell[i].vol;
		qcellp += qdivp*pcell[i].vol;
		if (bfirst)		//first call, not under-relax
		{
			pqdivg[i] = qdivg;
			pqdivp[i] = qdivp;
		}
		else
		{
			pqdivg[i] = pqdivg[i]*urf_srcrad1 + qdivg*urf_srcrad;
			pqdivp[i] = pqdivp[i]*urf_srcrad1 + qdivp*urf_srcrad;
		}
	}
	bfirst = false;
}

void CDiscreteOrdinateEquation::GetZoneRadiationHeatLoss(int nzone, double* phlrad)
{
	//nzone is the number of zones
	int i;
	int ncell = pmesh->ncell;
	CHexCell* pcell = pmesh->pcell;
	for (i=0; i<nzone; i++)
		phlrad[i] = 0;
	for (i=0; i<ncell; i++)
		phlrad[pcell[i].izone] -= (pqdivg[i] + pqdivp[i])*pcell[i].vol;
}

void CDiscreteOrdinateEquation::GetZoneConvectionHeatLoss(int nzone, double* phlconv)
{
	//nzone is the number of zones
	int i, j;
	int nbface = pmesh->nbface;
	int* pibface = pmesh->pibface;
	CQuadFace* pface = pmesh->pface;
	CQuadFace* pf;
	CBoundaryFaceProperty* pbfp;
	for (i=0; i<nzone; i++)
		phlconv[i] = 0;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		pbfp = pf->pbfp;
		if (pf->itype==1)	//wall face
		{
			pbfp->qconv = pbfp->hconv*(pbfp->temp - pbfp->tgas);	//qconv is generally negative
			phlconv[pf->izone] -= pbfp->qconv*pf->area;
		}
	}
}

void CDiscreteOrdinateEquation::CalcBoundaryNetHeatFlux()
{
	//net heat flux is defined as heat that goes to fluid
	//also adjust the incident heat flux to meet overall radiative energy balance.  The discrepancy is caused by face normal not aligned with predefined discrete ordinates
	//note: pbfp->qinc is not modified since it is used by next Solve(). when calculate qnet, use fqinc*pbfp->qinc
	int i, j;
	double qnet;
	double area_emis;
	double qsum_em;
	double qsum_ab;
	CBoundaryFaceProperty* pbfp;
	CQuadFace* pf;
	CQuadFace* pface = pmesh->pface;
	int nbface = pmesh->nbface;
	int* pibface = pmesh->pibface;
	qsum_em = 0;
	qsum_ab = 0;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		pbfp = pf->pbfp;
		switch (pf->itype)
		{
		case 1:		//wall
			area_emis = pf->area*pbfp->emis;
			qsum_em += area_emis*CT_RSIGMA*pow(pbfp->temp,4);
			qsum_ab += area_emis*pbfp->qinc;
			break;
		case 2:		//inlet
			area_emis = pf->area*pbfp->emis;
			qsum_em += area_emis*CT_RSIGMA*pow(pbfp->temp,4);
			qsum_ab += area_emis*pbfp->qinc;
			break;
		case 3:		//outlet
			area_emis = pf->area*pbfp->emis;
			qsum_em += area_emis*CT_RSIGMA*pow(pbfp->temp,4);
			qsum_ab += area_emis*pbfp->qinc;
			break;
		case 4:		//symmetry
			break;
		}
	}
	fqinc = (qsum_em - qcellg - qcellp)/qsum_ab;
	//debug, print fqinc
	printf("fqinc=%lg\n",fqinc);
	qface = 0;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		pbfp = pf->pbfp;
		switch (pf->itype)
		{
		case 1:		//wall
			qnet = pbfp->emis*(CT_RSIGMA*pow(pbfp->temp,4)-fqinc*pbfp->qinc);
			pbfp->qnet = qnet;
			qface += qnet*pf->area;
			break;
		case 2:		//inlet
			qnet = pbfp->emis*(CT_RSIGMA*pow(pbfp->temp,4)-fqinc*pbfp->qinc);
			pbfp->qnet = qnet;
			qface += qnet*pf->area;
			break;
		case 3:		//outlet
			qnet = pbfp->emis*(CT_RSIGMA*pow(pbfp->temp,4)-fqinc*pbfp->qinc);
			pbfp->qnet = qnet;
			qface += qnet*pf->area;
			break;
		case 4:		//symmetry
			pbfp->qnet = 0;
			break;
		}
	}
	//note qface should be the same as (qcellg+qcellp). It is not used anywhere else
	//debug, print fqinc
	printf("qcellg=%lg, qcellp=%lg, qface=%lg\n",qcellg, qcellp, qface);
}

void CDiscreteOrdinateEquation::UpdateWallTemperature()
{
	int i, j;
	CBoundaryFaceProperty* pbfp;
	CQuadFace* pf;
	CQuadFace* pface = pmesh->pface;
	int nbface = pmesh->nbface;
	int* pibface = pmesh->pibface;
	double urf_twall1 = 1 - urf_twall;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		if (pf->itype==1)		//for wall face only, outlet face temperature is fixed
		{
			pbfp = pf->pbfp;
			pbfp->temp = urf_twall1*pbfp->temp + urf_twall*SolveWallTemperature(pbfp->emis, pbfp->rwall, pbfp->hconv, pbfp->qinc, pbfp->tgas, pbfp->tback);
		}
	}
}

void CDiscreteOrdinateEquation::UpdateRadiationProperties(double* pkagz, double* pkapz, double* pkat4gz, double* pkat4pz)
{
	int i;
	int izone;
	int ncell = pmesh->ncell;
	CHexCell* pcell = pmesh->pcell;
	for (i=0; i<ncell; i++)
	{
		izone = pcell[i].izone;
		pkag[i] = pkagz[izone];
		pkap[i] = pkapz[izone];
		pkat4g[i] = pkat4gz[izone];
		pkat4p[i] = pkat4pz[izone];
	}
}

void CDiscreteOrdinateEquation::UpdateConvectionBoundaryFaceProperties(double* phconv, double* ptgas)
{
	int i, j;
	int izone;
	int	nbface = pmesh->nbface;
	int* pibface = pmesh->pibface;
	CQuadFace*	pface = pmesh->pface;
	CQuadFace*	pf;
	CBoundaryFaceProperty* pbfp;
	for (i=0; i<nbface; i++)
	{
		j = pibface[i];
		pf = &pface[j];
		if (pf->itype==1)	//wall boundary
		{
			pbfp = pf->pbfp;
			izone = pf->izone;
			pbfp->hconv = phconv[izone];
			pbfp->tgas = ptgas[izone];
		}
	}
}

double CDiscreteOrdinateEquation::SolveWallTemperature(double em, double rw, double hc, double qi, double tg, double tb)
{
	short i;
	double t, t1, t2, t3;
	double func;
	double ers = em*rw*CT_RSIGMA;
	double hr1 = hc*rw+1;
	double term = em*rw*qi + hc*rw*tg + tb;
	t = 500;
	func = t*(ers*t*t*t+hr1);
	i = 0;
	if (func>term)
	{
		do
		{
			i++;
			if (i>5)
			{
				printf("Cannot solve wall temperature!");
				return t;
			}
			t2 = t;
			t *= 0.75;
			func = t*(ers*t*t*t+hr1);
		}while (func>term);
		t1 = t;
	}
	else
	{
		do
		{
			i++;
			if (i>5)
			{
				printf("Cannot solve wall temperature!");
				return t;
			}
			t1 = t;
			t *= 1.5;
			func = t*(ers*t*t*t+hr1);
		}while (func<term);
		t2 = t;
	}
	i = 0;
	t = (t1+t2)/2;
	do
	{
		i++;
		if (i>20)
		{
			printf("wall temperature does not converge after 20 iterations!");
			return t;
		}
		func = t*(ers*t*t*t+hr1);
		if (func>term)
			t2 = t;
		else
			t1 = t;
		t = (t1+t2)/2;
		t3 = t-t1;
	}while (t3<-0.01 || t3>0.01);
	return t;
}