//IdealGasMixture.cpp: implementation of the CIdealGasMixture class.
//
///////////////////////////////////////////////////

//my incude
#include <cmath>
#include "Constants.h"
#include "IdealGasMixture.h"
#include "UtilityFunctions.h"

///////////////////////////////////////////////////////
//Constuctor and Destructor
///////////////////////////////////////////////////////
double CIdealGasMixture::x[50+GS_NE]={0};
double CIdealGasMixture::dx[50+GS_NE]={0};
double CIdealGasMixture::aij[50][GS_NE]={0};
double CIdealGasMixture::jacob[50+GS_NE][51+GS_NE]={0};

CIdealGasMixture::CIdealGasMixture()
{
	int i;
	name = "Gas Mixture";
	mflow = 1;
	area = 1;
	vel = 1;
	bdry = false;
	ivol = 0;
	imol = 0;
	ictp = 1;			//for calculations not including equilibrium, should set to 0
	icalc = 0;
	ne = 4;
	p = 101325;
	t = 298.15;
	//elements
	elem[0] = 'C';
	elem[1] = 'H';
	elem[2] = 'O';
	elem[3] = 'N';
	elem[4] = 'S';
	elem[5] = 'L';
	elem[6] = 'A';
	//default species
	pgas[0] = CIdealGas::GetGasByName("O2");
	pgas[1] = CIdealGas::GetGasByName("N2");
	pgas[2] = CIdealGas::GetGasByName("H2O");
	pgas[3] = CIdealGas::GetGasByName("CO2");
	ns = 4;
	//default species fraction
	for (i=0; i<50; i++)
	{
		fsp_mole[i] = 0;
		fsp_mass[i] = 0;
	}
	fsp_mole[0] = 0.206201;
	fsp_mole[1] = 0.777811;
	fsp_mole[2] = 0.0156534;
	fsp_mole[3] = 0.000334678;
	InitArraysFromFsp();
	CalcAllProperties();
}

CIdealGasMixture::~CIdealGasMixture()
{
}

CIdealGasMixture::CIdealGasMixture(const CIdealGasMixture &tt)
{
	int i, j;
	name = tt.name;
	mflow = tt.mflow;
	vflow = tt.vflow;
	area = tt.area;
	vel = tt.vel;
	bdry = tt.bdry;
	ivol = tt.ivol;
	imol = tt.imol;
	ictp = tt.ictp;
	icalc = tt.icalc;
	ne = tt.ne;	
	ns = tt.ns;
	p = tt.p;
	t = tt.t;
	h = tt.h;
	g = tt.g;
	hhv = tt.hhv;
	lhv = tt.lhv;
	hsh = tt.hsh;
	cp = tt.cp;
	mw = tt.mw;
	den = tt.den;
	er = tt.er;
	for (j=0; j<GS_NE; j++)
		eflow[j] = tt.eflow[j];
	for (j=0; j<ne; j++)
	{
		elem[j] = tt.elem[j];
		ielem[j] = tt.ielem[j];
		aw[j] = tt.aw[j];
		bj[j] = tt.bj[j];
	}
	for (i=0; i<ns; i++)
	{
		ni[i] = tt.ni[i];
		fsp_mole[i] = tt.fsp_mole[i];
		fsp_mass[i] = tt.fsp_mass[i];
		pgas[i] = tt.pgas[i];
	}
	for(i = ns; i < 50; i++)
	{
		fsp_mole[i] = fsp_mass[i] = 0.0;
		pgas[i] = NULL;
	}
}

CIdealGasMixture CIdealGasMixture::operator=(CIdealGasMixture tt)
{
	int i, j;
	name = tt.name;
	mflow = tt.mflow;
	vflow = tt.vflow;
	area = tt.area;
	vel = tt.vel;
	bdry = tt.bdry;
	ivol = tt.ivol;
	imol = tt.imol;
	ictp = tt.ictp;
	icalc = tt.icalc;
	ne = tt.ne;	
	ns = tt.ns;
	p = tt.p;
	t = tt.t;
	h = tt.h;
	g = tt.g;
	hhv = tt.hhv;
	lhv = tt.lhv;
	hsh = tt.hsh;
	cp = tt.cp;
	mw = tt.mw;
	den = tt.den;
	er = tt.er;
	for (j=0; j<GS_NE; j++)
		eflow[j] = tt.eflow[j];
	for (j=0; j<ne; j++)
	{
		elem[j] = tt.elem[j];
		ielem[j] = tt.ielem[j];
		aw[j] = tt.aw[j];
		bj[j] = tt.bj[j];
	}
	for (i=0; i<ns; i++)
	{
		ni[i] = tt.ni[i];
		fsp_mole[i] = tt.fsp_mole[i];
		fsp_mass[i] = tt.fsp_mass[i];
		pgas[i] = tt.pgas[i];
	}	
	for(i = ns; i < 50; i++)
	{
		fsp_mole[i] = fsp_mass[i] = 0.0;
		pgas[i] = NULL;
	}
	return *this;
}

void CIdealGasMixture::Write(FILE* pf)
{
	int i, j;
	int iversion = 0;
	fwrite(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::WriteString(name,pf);
	fwrite(&bdry,sizeof(bool),1,pf);
	fwrite(&ivol,sizeof(int),1,pf);
	fwrite(&imol,sizeof(int),1,pf);
	fwrite(&ictp,sizeof(int),1,pf);
	fwrite(&icalc,sizeof(int),1,pf);
	fwrite(&ne,sizeof(int),1,pf);
	fwrite(&ns,sizeof(int),1,pf);
	fwrite(&mflow,sizeof(double),1,pf);
	fwrite(&vflow,sizeof(double),1,pf);
	fwrite(&area,sizeof(double),1,pf);
	fwrite(&vel,sizeof(double),1,pf);
	fwrite(&p,sizeof(double),1,pf);
	fwrite(&t,sizeof(double),1,pf);
	fwrite(&h,sizeof(double),1,pf);
	fwrite(&g,sizeof(double),1,pf);
	fwrite(&hhv,sizeof(double),1,pf);
	fwrite(&lhv,sizeof(double),1,pf);
	fwrite(&hsh,sizeof(double),1,pf);
	fwrite(&cp,sizeof(double),1,pf);
	fwrite(&mw,sizeof(double),1,pf);
	fwrite(&den,sizeof(double),1,pf);
	fwrite(&er,sizeof(double),1,pf);
	fwrite(eflow,sizeof(double),GS_NE,pf);
	fwrite(ielem,sizeof(int),ne,pf);
	fwrite(elem,sizeof(double),ne,pf);
	fwrite(aw,sizeof(double),ne,pf);
	fwrite(bj,sizeof(double),ne,pf);
	fwrite(ni,sizeof(double),ns,pf);
	fwrite(fsp_mole,sizeof(double),ns,pf);
	fwrite(fsp_mass,sizeof(double),ns,pf);
	for (i=0; i<ns; i++)
	{
		j = pgas[i]->GetIndex();
		fwrite(&j,sizeof(int),1,pf);
	}
}

void CIdealGasMixture::Read(FILE* pf)
{
	int i, j;
	int iversion = 0;
	fread(&iversion,sizeof(int),1,pf);
	CUtilityFunctions::ReadString(name,pf);
	fread(&bdry,sizeof(bool),1,pf);
	fread(&ivol,sizeof(int),1,pf);
	fread(&imol,sizeof(int),1,pf);
	fread(&ictp,sizeof(int),1,pf);
	fread(&icalc,sizeof(int),1,pf);
	fread(&ne,sizeof(int),1,pf);
	fread(&ns,sizeof(int),1,pf);
	fread(&mflow,sizeof(double),1,pf);
	fread(&vflow,sizeof(double),1,pf);
	fread(&area,sizeof(double),1,pf);
	fread(&vel,sizeof(double),1,pf);
	fread(&p,sizeof(double),1,pf);
	fread(&t,sizeof(double),1,pf);
	fread(&h,sizeof(double),1,pf);
	fread(&g,sizeof(double),1,pf);
	fread(&hhv,sizeof(double),1,pf);
	fread(&lhv,sizeof(double),1,pf);
	fread(&hsh,sizeof(double),1,pf);
	fread(&cp,sizeof(double),1,pf);
	fread(&mw,sizeof(double),1,pf);
	fread(&den,sizeof(double),1,pf);
	fread(&er,sizeof(double),1,pf);
	fread(eflow,sizeof(double),GS_NE,pf);
	fread(ielem,sizeof(int),ne,pf);
	fread(elem,sizeof(double),ne,pf);
	fread(aw,sizeof(double),ne,pf);
	fread(bj,sizeof(double),ne,pf);
	fread(ni,sizeof(double),ns,pf);
	fread(fsp_mole,sizeof(double),ns,pf);
	fread(fsp_mass,sizeof(double),ns,pf);
	for (i=0; i<ns; i++)
	{
		fread(&j,sizeof(int),1,pf);
		pgas[i] = CIdealGas::GetGasByIndex(j);
	}
}

double CIdealGasMixture::GetFspSum()
{
	int j;
	double sum = 0;
	if (imol)
	{
		for (j=0; j<ns; j++)
		{
			sum += fsp_mass[j];
		}
	}
	else
	{
		for (j=0; j<ns; j++)
		{
			sum += fsp_mole[j];
		}
	}
	return sum;
}

void CIdealGasMixture::NormalizeFsp()
{
	int i;
	double* fsp;
	double sum = GetFspSum();	
	if (sum<=0)
	{
		printf("Sum of fractions <= 0!");
		return;
	}
	if (imol)
		fsp = fsp_mass;
	else
		fsp = fsp_mole;
	for (i=0; i<ns; i++)
	{
		fsp[i] /= sum;
	}
}

void CIdealGasMixture::InitArraysFromFsp()
{
	//calculate ielem[], elem[], aw[], bj[], eflow[], ni[] and set ne, mw
	int i, j;
	double sum;
	double mol_chons[GS_NE] = {0};
	double aw_chons[GS_NE] = {CT_AWC, CT_AWH, CT_AWO, CT_AWN, CT_AWS, CT_AWCL, CT_AWAR};
	char ele_chons[GS_NE] = {'C', 'H', 'O', 'N', 'S', 'L', 'A'};
	sum = 0;
	if (imol)		//mass fraction
	{
		for (i=0; i<ns; i++)
		{
			fsp_mole[i] = fsp_mass[i]/pgas[i]->mw;
			sum += fsp_mole[i];
		}
		for (i=0; i<ns; i++)
		{
			fsp_mole[i] /= sum;
		}
	}
	else
	{
		for (i=0; i<ns; i++)
		{
			fsp_mass[i] = fsp_mole[i]*pgas[i]->mw;
			sum += fsp_mass[i];
		}
		for (i=0; i<ns; i++)
		{
			fsp_mass[i] /= sum;
		}
	}
	mw = 0;
	for (i=0; i<ns; i++)
	{
		for (j=0; j<GS_NE; j++)
		{
			mol_chons[j] += pgas[i]->natom[j]*fsp_mole[i];
		}
		mw += pgas[i]->mw*fsp_mole[i];
	}
	for (i=0; i<ns; i++)
	{
		ni[i] = fsp_mole[i]/mw;
	}
	ne = 0;
	for (i=0; i<GS_NE; i++)
	{
		if (mol_chons[i]>0)
		{
			ielem[ne] = i;
			elem[ne] = ele_chons[i];
			aw[ne] = aw_chons[i];
			bj[ne] = mol_chons[i]/mw;
			ne++;
		}
		eflow[i] = mflow*mol_chons[i]/mw;
	}
}

void CIdealGasMixture::InitArraysFromEflow()
{
	//calculate ielem[], elem[], aw[], bj[], and set ne
	int i;
	double aw_chons[GS_NE] = {CT_AWC, CT_AWH, CT_AWO, CT_AWN, CT_AWS, CT_AWCL, CT_AWAR};
	char ele_chons[GS_NE] = {'C', 'H', 'O', 'N', 'S', 'L', 'A'};
	//first calculate mflow from eflow
	mflow = 0;
	for (i=0; i<GS_NE; i++)
		mflow += eflow[i]*aw_chons[i];
	if (mflow<=0)
	{
		mflow = 0;
		printf("Invalid elemental molar flow inputs\n");
		return;
	}
	ne = 0;
	for (i=0; i<GS_NE; i++)
	{
		if (eflow[i]>0)
		{
			ielem[ne] = i;
			elem[ne] = ele_chons[i];
			aw[ne] = aw_chons[i];
			bj[ne] = eflow[i]/mflow;
			ne++;
		}
	}
}

void CIdealGasMixture::GetDryFsp(double f[50])
{
	int i;
	double sum = 0;
	double fh2o = 0;		//fraction of H2O and H2O(L)
	double fh2o1;			//sum-fh2o
	double* fsp;
	std::string str;
	if (imol)
		fsp = fsp_mass;
	else
		fsp = fsp_mole;
	for (i=0; i<ns; i++)
	{
		sum += fsp[i];
		str = pgas[i]->name;
		if (str=="H2O" || str=="H2O(L)")
		{
			fh2o += fsp[i];
		}
	}
	fh2o1 = sum-fh2o;
	for (i=0; i<ns; i++)
	{
		str = pgas[i]->name;
		if (str=="H2O" || str=="H2O(L)")
		{
			f[i] = 0;
		}
		else
		{
			f[i] = fsp[i]/fh2o1;
		}
	}
}

void CIdealGasMixture::CalcFspFromNi()
{
	int i;
	double n = 1.0/mw;
	for (i=0; i<ns; i++)
	{
		fsp_mole[i] = ni[i]/n;
		fsp_mass[i] = ni[i]*pgas[i]->mw;
	}
}

void CIdealGasMixture::CalcAij()
//calculate number of atoms of each elements in a molecule
{
	int i,j;
	for (i=0; i<ns; i++)
	{
		for (j=0; j<ne; j++)
		{
			aij[i][j] = pgas[i]->natom[ielem[j]];
		}
	}
}

void CIdealGasMixture::CalcNiFromX()
{
	int i;
	for (i=0; i<ns; i++)
	{
		ni[i] = exp(x[i]);
	}
}

void CIdealGasMixture::CalcJacobian()
//CalcAij() must be called first
//Calculate augmented Jacobian Matrix
//also calculate the Gibbs free energy of mixture and grand function
{
	CalcAij();
	int i,j;
	double rt = CT_UR*t;		//RT
	double rtdn;				//RT/N
	double n = 0;				//total kmoles in 1 kg mixture
	double sum;					//for g calculation
	double gbar;				//partial molar Gibbs free energy for each species
	//calculate total kmoles
	for (i=0; i<ns; i++)
	{
		n += ni[i];
	}
	mw = 1.0/n;				//assign molecular weight of mixture
	rtdn = rt/n;
	//assign the upper left corner of Jacobian matrix
	for (i=0; i<ns; i++)
	{
		for (j=0; j<ns; j++)
		{
			jacob[i][j] = -rtdn*ni[j];
		}
		jacob[i][i] += rt;
	}
	//assign the upper right corners of Jacobian matrix
	for (i=0; i<ns; i++)
	{
		for (j=0; j<ne; j++)
		{
			jacob[i][j+ns] = aij[i][j];			
		}
	}
	//assign the lower left corners of Jacobian matrix
	for (j=0; j<ne; j++)
	{
		for (i=0; i<ns; i++)
		{
			jacob[j+ns][i] = aij[i][j]*ni[i];		
		}
	}
	//assign the lower right corner of Jacobian matrix
	for (i=0; i<ne; i++)
	{
		for (j=0; j<ne; j++)
		{
			jacob[i+ns][j+ns] = 0;
		}
	}
	//calculate -fi (last column in the augumented Jacobian matrix
	//first calculate the partial molar Gibbs free energy
	g = 0;
	for (i=0; i<ns; i++)
	{
		sum = 0;				//sigmaj(lamitaj*aij)
		for (j=0; j<ne; j++)
		{
			sum += x[ns+j]*aij[i][j];
		}
		gbar = pgas[i]->CalcG(t,p)+rt*log(ni[i]/n);
		jacob[i][ns+ne] = -sum-gbar;
		g += gbar*ni[i];
	}
	for (j=0; j<ne; j++)
	{
		sum = 0;				//sigmai(aij*ni)
		for (i=0; i<ns; i++)
		{
			sum += aij[i][j]*ni[i];
		}
		jacob[j+ns][ns+ne] = bj[j] - sum;
	}
}

int CIdealGasMixture::SolveLAE(int n, double a[50+GS_NE][51+GS_NE], double x[50+GS_NE])
//solve linear algebraic equations using Gaussian elimination
//with largest pivot.  a is an augmented matrix
{
	int i,j,k;
	int irow;
	double pmax;
	double aabs;
	double aswap;
	double fac;
	double sum;
	for (i=0; i<n-1; i++)
	{
		pmax = 0;
		irow = i;
		for (k=i; k<n; k++)
		{
			aabs = fabs(a[k][i]);
			if (aabs>pmax)
			{
				pmax = aabs;
				irow = k;
			}
		}
		if (pmax<=0) return -1;		//singular
		if (irow!=i)		//swap
		{
			for (k=i; k<=n; k++)
			{
				aswap = a[i][k];
				a[i][k] = a[irow][k];
				a[irow][k] = aswap;
			}
		}
		//do elimination
		for (j=i+1; j<n; j++)
		{
			fac = a[j][i]/a[i][i];
			for (k=i; k<=n; k++)
			{
				a[j][k] = a[j][k]-fac*a[i][k];
			}
		}
	}
	//back calculation
	k = n-1;
	x[k] = a[k][n]/a[k][k];
	for (i=n-2; i>=0; i--)
	{
		sum = 0;
		for (j=i+1; j<n; j++)
		{
			sum += a[i][j]*x[j];
		}
		x[i] = (a[i][n]-sum)/a[i][i];
	}
	return 0;
}

bool CIdealGasMixture::AreAllElementsInGasArray()
{
	int i, j;
	int sumj;
	std::string strerr = "  is not in any species in the species list!";
	for (j=0; j<ne; j++)
	{
		sumj = 0;
		for (i=0; i<ns; i++)
		{
			sumj += pgas[i]->natom[ielem[j]];
		}
		if (sumj==0)
		{
			strerr[0] = elem[j];
			printf("%s\n",strerr.c_str());
			return false;
		}
	}
	return true;
}

void CIdealGasMixture::GetValidSpeciesList(std::vector<CIdealGas*>& plist)
{
	int i, j, k;
	bool bvalid, bfound;
	CIdealGas* pg;
	CIdealGas** gaslist = CIdealGas::GetGasList();	//build-in gas species list
	plist.clear();
	for (k=0; k<GS_NSP; k++)
	{
		pg = gaslist[k];
		bvalid = true;
		for (i=0; i<GS_NE; i++)
		{
			if (pg->natom[i]!=0)		//the element is in the species
			{
				bfound = false;
				for (j=0; j<ne; j++)
				{
					if (ielem[j]==i)
					{
						bfound = true;
						break;
					}
				}
				if (!bfound)			//the element is not found in the selection list
				{
					bvalid = false;
					break;
				}
			}
		}
		if (bvalid)
			plist.push_back(pg);
	}
}

void CIdealGasMixture::AddValidSpecies(std::string str, std::vector<CIdealGas*>& pvl)
{
	if (ns>50)
	{
		printf("Cannot exceed 50 species!");
		return;
	}
	int i;
	int nsize = pvl.size();
	for (i=0; i<nsize; i++)
	{
		if (pvl[i]->name==str)
		{
			pgas[ns] = pvl[i];
			ns++;
			break;
		}
	}
}

void CIdealGasMixture::AddValidSpecies(int isp, std::vector<CIdealGas*>& pvl)
{
	if (ns>50)
	{
		printf("Cannot exceed 50 species!");
		return;
	}
	CIdealGas* pgadd = CIdealGas::GetGasByIndex(isp);
	int i;
	int nsize = pvl.size();
	for (i=0; i<nsize; i++)
	{
		if (pvl[i]==pgadd)
		{
			pgas[ns] = pvl[i];
			ns++;
			break;
		}
	}
}

void CIdealGasMixture::DeleteSpecies(int i)
{
	int k;
	if (i>=ns || i<0) 
	{
		printf("Species index out of range");
		return;
	}
	for (k=i; k<ns+1; k++)
	{
		pgas[k] = pgas[k+1];
		fsp_mole[k] = fsp_mole[k+1];
		fsp_mass[k] = fsp_mass[k+1];
	}
	ns--;
}

void CIdealGasMixture::DeleteSpecies(std::string str)
{
	int i;
	for (i=0; i<ns; i++)
	{
		if (pgas[i]->name == str)
		{
			DeleteSpecies(i);
			return;
		}
	}
	printf("Species not found in species list!");
	return;
}

int CIdealGasMixture::SolveXAtTP()
//solve Lagrange equation set at const temperature and pressure
{
	int i,j;
	int nn = ns+ne;			//number of equations to be solved
	int nit = 0;			//number of iteration count
	bool bconv;				//converged
	//initial guess, this may be improved based on equivalence ratio
	for (i=0; i<ns; i++)
	{
		x[i] = pgas[i]->GetX0(t,er);
	}
	for (j=0; j<ne; j++)
	{
		x[j+ns] = 100000;	//Lagrangian coefficients seem to be large numbers
	}
	do
	{
		nit++;
		CalcNiFromX();
		CalcJacobian();
		SolveLAE(nn,jacob,dx);			//call linear equations solver
		for (i=0; i<ns; i++)
		{		
			if (dx[i]>5) dx[i] = 2;		//limit maximum change
			if (dx[i]<-5) dx[i] = -2;	//limit maximum change
			if (x[i]<-100) dx[i] = 0;	//ni[i]<3.7e-44, ignore it
			x[i] += dx[i];				//don't use under-relaxation
		}
		for (j=0; j<ne; j++)
		{
			x[j+ns] += dx[j+ns];		//don't use under-relaxation
		}
		//check convergence
		bconv = true;
		for (i=0; i<ns; i++)
		{
			if (fabs(dx[i])>5e-11)
			{
				bconv = false;
				break;
			}
		}
	}while (!bconv && nit<150);
	if (nit>=150)
	{
		printf("SolveXAtTP() not converged after 150 iterations");
		return -1;
	}
	return 0;					//converged
}


int CIdealGasMixture::SolveXAtHP()
//solve Lagrange equation set at const enthalpy and pressure
//based on a guessed temperature use SolveXAtTP() to match enthalpy
//member variable t is updated
{
	int nit;			//number of iterations
	double t0 = 300;	//initial min temperature
	double t1 = 2500;	//initial max temperature
	double tn;			//new guess
	double f0;			//min of (hguess - h)
	double f1;			//max of (hguess - h)
	double fn;			//new (hguess - h)
	double err;			//relative error of t
	t = t0;
	SolveXAtTP();
	f0 = CalcH()-h;
	nit = 0;
	while (f0>0)		//min temperature too high
	{
		if (nit>3) 
		{
			printf("Enthalpy too low!");
			return -1;	//not going to converge
		}
		nit++;
		t0 = t0/1.5;
		t = t0;
		SolveXAtTP();
		f0 = CalcH()-h;
	}
	t = t1;
	SolveXAtTP();
	f1 = CalcH()-h;
	nit = 0;
	while (f1<0)		//max temperature too low
	{
		if (nit>3) 
		{
			printf("Enthalpy too high!");
			return -1;	//not going to converge
		}
		nit++;
		t1 = t1*1.5;
		t = t1;
		SolveXAtTP();
		f1 = CalcH()-h;
	}
	nit = 0;
	do
	{
		nit++;
		tn = (t0+t1)/2;
		t = tn;
		SolveXAtTP();
		fn = CalcH()-h;
		err = fabs((tn-t0)/tn);
		if (fn<0)
		{
			t0 = tn;
		}
		else
		{
			t1 = tn;
		}
	} while (err>0.001 && nit<20);	
	//this code usually converges within 15 iterations
	if (nit>=20)
	{
		printf("SolveXAtHP() not converged after 20 iterations");
		return -1;		//not converged within 20 iterations
	}
	return 0;			//converged
}

int CIdealGasMixture::Solve()
{
	//CalcJacobian() must have been called
	int ierr = 0;
	if (ictp==1)
		ierr = SolveXAtTP();
	else
		ierr = SolveXAtHP();
	CalcFspFromNi();
	return ierr;
}

double CIdealGasMixture::CalcH()
//calculate using the current t
//does not update member variable h
{
	int i;
	double xh = 0;
	for (i=0; i<ns; i++)
	{
		xh += ni[i]*pgas[i]->CalcH(t);
	}
	return xh;
}

double CIdealGasMixture::CalcS()
//calculate using the current t, p
//sum of entropy of individual species
{
	int i;
	double xs = 0;
	for (i=0; i<ns; i++)
	{
		xs += ni[i]*pgas[i]->CalcS(t,p);
	}
	return xs;
}

double CIdealGasMixture::CalcCp()
//calculate Cp of mixture [J/kg-K] using the current t
//does not update member variable cp
{
	int i;
	double xcp = 0;
	for (i=0; i<ns; i++)
	{
		xcp += ni[i]*pgas[i]->CalcCp(t);
	}
	return xcp;
}

double CIdealGasMixture::CalcCp(double temp)
//calculate Cp of mixture [J/kg-K] at temp
{
	int i;
	double xcp = 0;
	for (i=0; i<ns; i++)
	{
		xcp += ni[i]*pgas[i]->CalcCp(temp);
	}
	return xcp;
}

double CIdealGasMixture::CalcEr()
{
	//does not update the current member variable er
	//assume C->CO2, H->H2O, O->O2, N->N2, S->SO2
	//equilvalence ratio is defined as ratio of fuel existed
	//to the stoichiometric fuel needed based on oxygen available
	//C:4, H:1, O:-2, N:0, S:4
	int j;
	double covfuel = 0;		//covalence of fuel
	double covoxy = 0;		//covalence of oxygen
	for (j=0; j<ne; j++)
	{
		switch (elem[j])
		{
		case 'C':
			covfuel += 4*bj[j];
			break;
		case 'H':
			covfuel += bj[j];
			break;
		case 'O':
			covoxy += 2*bj[j];
			break;
		case 'S':
			covfuel += 4*bj[j];
			break;
		}
	}
	if (covoxy>0)
		return covfuel/covoxy;
	else
		return -1;		//no oxygen
}

double CIdealGasMixture::CalcViscosity(double temp)
{
	int i, j;
	double visc;
	double sum;
	double mwt[50];			//molecular weight of species
	double vis[50];			//viscosity of individual species
	double phi[50][50];	//phi based on Bird's formula on page 24
	//update molecular weight and calulate species viscosity
	for (i=0; i<ns; i++)
	{
		mwt[i] = pgas[i]->mw;
		vis[i] = pgas[i]->CalcViscosity(temp);
	}
	for (i=0; i<ns; i++)
	{
		for (j=0; j<ns; j++)
		{
			if (i==j)
				phi[i][j] = 1;
			else
			{
				phi[i][j] = 1 + sqrt(vis[i]/vis[j]*sqrt(mwt[j]/mwt[i]));
				phi[i][j] *= phi[i][j]/sqrt(8*(1+mwt[i]/mwt[j]));
			}
		}
	}
	visc = 0;
	for (i=0; i<ns; i++)
	{
		sum = 0;
		for (j=0; j<ns; j++)
			sum += fsp_mole[j]*phi[i][j];
		visc += fsp_mole[i]*vis[i]/sum;
	}
	return visc;
}

double CIdealGasMixture::CalcConductivity(double temp)
{
	int i, j;
	double cond;
	double sum;
	double mwt[50];			//molecular weight of species
	double con[50];			//viscosity of individual species
	double phi[50][50];	//phi based on Bird's formula on page 258
	//update molecular weight and calulate species viscosity
	for (i=0; i<ns; i++)
	{
		mwt[i] = pgas[i]->mw;
		con[i] = pgas[i]->CalcConductivity(temp);
	}
	for (i=0; i<ns; i++)
	{
		for (j=0; j<ns; j++)
		{
			if (i==j)
				phi[i][j] = 1;
			else
			{
				phi[i][j] = 1 + sqrt(con[i]/con[j]*sqrt(mwt[j]/mwt[i]));
				phi[i][j] *= phi[i][j]/sqrt(8*(1+mwt[i]/mwt[j]));
			}
		}
	}
	cond = 0;
	for (i=0; i<ns; i++)
	{
		sum = 0;
		for (j=0; j<ns; j++)
			sum += fsp_mole[j]*phi[i][j];
		cond += fsp_mole[i]*con[i]/sum;
	}
	return cond;
}

void CIdealGasMixture::CalcAllProperties()
//assume equilibrium already solved
{
	int i;
	double fh2o;		//mole fraction of H2O
	double fh2ol;		//mole fraction of H2O(L)
	double fwater;		//fh2o+fh2ol
	double href;		//reference h at 298.15 K
	double hcnoh2o;		//heat of combustion without considering H2O
	if (ictp==1)		//const T P
	{
		h = CalcH();
		if (fabs(h)<0.05)
			h = 0;
	}
	cp = CalcCp();
	er = CalcEr();
	//calculate other properties
	//density
	den = p*mw/CT_UR/t;
	//sensible heat
	href = 0;
	for (i=0; i<ns; i++)
	{
		href += ni[i]*pgas[i]->CalcH(298.15);
	}
	hsh = h-href;
	if (fabs(hsh)<0.05) hsh = 0;	//eliminate round-off error
	//calculate heat of combustion (close to high heating value)
	double bjj[6] = {0,0,0,0,0,0};	//kmole of C, H, O, N, S, Cl per kg of mixture
	CIdealGas co2, h2o, h2ol, so2;
	co2.SetName("CO2");
	co2.AssignData();
	h2o.SetName("H2O");
	h2o.AssignData();
	h2ol.SetName("H2O(L)");
	h2ol.AssignData();
	so2.SetName("SO2");
	so2.AssignData();
	for (i=0; i<ne; i++)
	{
		bjj[ielem[i]] = bj[i];
	}
	//get rid of H2O and H2O(L) in reactants
	fh2o = GetMoleFraction("H2O");
	fh2ol = GetMoleFraction("H2O(L)");
	fwater = fh2o+fh2ol;
	href -= fh2o/mw*h2o.CalcH(298.15);
	href -= fh2ol/mw*h2ol.CalcH(298.15);
	hcnoh2o = bjj[0]*co2.CalcH(298.15)+bjj[4]*so2.CalcH(298.15);
	hhv = href-hcnoh2o-(bjj[1]/2-fwater/mw)*h2ol.CalcH(298.15);
	lhv = href-hcnoh2o-(bjj[1]/2-fwater/mw)*h2o.CalcH(298.15);
	if (fabs(hhv)<0.05) hhv = 0;	//get rid of round-off error
	if (fabs(lhv)<0.05) lhv = 0;	//get rid of round-off error
	//calculate velcocity, area, or mass flow
	switch (icalc)
	{
	case 0:
		if (area>0)
			vel = mflow/den/area;
		else
			vel = 0;
		break;
	case 1:
		if (vel>0)
			area = mflow/den/vel;
		else
			area = 0;
		break;
	case 2:
		mflow = vel*den*area;
		break;
	}
	vflow = mflow / (101325*mw/8314.3/273.15);
}

int CIdealGasMixture::GetSpeciesIndex(std::string sp)
{
	int i;
	for (i=0; i<ns; i++)
	{
		if (pgas[i]->name==sp)
			return i;
	}
	return -1;			//species not found
}

double CIdealGasMixture::GetMoleFraction(std::string sp)
//return mole fraction of species sp
//may be useful for char oxidation and vaporization models
{
	int i;
	for (i=0; i<ns; i++)
	{
		if (pgas[i]->name==sp)
		{
			return ni[i]*mw;
		}
	}
	return 0;		//if species name not found
}

double CIdealGasMixture::GetMoleFraction(int i)
{
	return ni[i]*mw;
}

double CIdealGasMixture::GetSpeciesMoleFlow(std::string sp)
{
	int i;
	for (i=0; i<ns; i++)
	{
		if (pgas[i]->name==sp)
		{
			return ni[i]*mflow;
		}
	}
	return 0;		//if species name not found
}

double CIdealGasMixture::GetValenceO2Fraction(int imm, int iwd)
{
	//imm=0 for mol fraction 1 for mass fraction
	//iwd=0 for wet and 1 for dry
	//mflow must be positive
	if (mflow<=0) return 0;
	double fo2mol;
	if (er<1)
	{
		fo2mol = GetMoleFraction("O2");
		if (iwd==0)		//wet
		{
			if (imm==0)	//mole
				return fo2mol;
			else
				return fo2mol*CT_AWO*2/mw;
		}
		else
		{
			double fh2o = GetMoleFraction("H2O");
			if (imm==0)	//mole
				return fo2mol/(1-fh2o);
			else
				return fo2mol*CT_AWO*2/(mw-fh2o*(CT_AWH*2+CT_AWO));
		}
	}
	//for er>1, calculated based on eflow[]
	double eflowco2 = eflow[0];
	double eflowso2 = eflow[4];
	double eflowh2o = eflow[1]/2;
	double eflown2 = eflow[3]/2;
	double eflowcl2 = eflow[5]/2;
	double eflowar = eflow[6];
	double eflowo2 = eflow[2]/2 - eflowco2 - eflowh2o/2 - eflowso2;		//negative value
	double eflowsum = eflowco2 + eflowso2 + eflowh2o + eflowo2 + eflown2 + eflowcl2 + eflowar;
	if (iwd==0)		//wet
	{
		if (imm==0)		//mole based
			return eflowo2/eflowsum;
		else
			return eflowo2*CT_AWO*2/mflow;
	}
	else				//dry
	{
		if (imm==0)		//mole based
			return eflowo2/(eflowsum-eflowh2o);
		else
			return eflowo2*CT_AWO*2/(mflow-eflowh2o*(CT_AWH*2+CT_AWO));
	}
}
