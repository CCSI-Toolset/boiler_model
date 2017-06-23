// ParticleRadiationProperty.cpp: implementation of the CParticleRadiationProperty class.
//
//////////////////////////////////////////////////////////////////////

#include <complex>
#include "Constants.h"

#include "ParticleRadiationProperty.h"

CParticleRadiationProperty* CParticleRadiationProperty::pinstance = NULL;
double CParticleRadiationProperty::refre_coal = 1.85;
double CParticleRadiationProperty::refim_coal = 0.25;
double CParticleRadiationProperty::refre_ash[RP_NWL] =	  {	1.55,   1.55,  1.55,  1.55,  1.55,   1.5,   1.5,   1.5,   1.5,   1.5,  1.45,  1.45,  1.4, 1.35,  1.3,  1.2, 1.0,0.9,1.1,1.4,1.7,2.0,2.3,2.1, 1.8, 1.7,1.7};
double CParticleRadiationProperty::refim_ash[3][RP_NWL] = {0.001,0.00003,0.0003,0.0004,0.0004,0.0003,0.0003,0.0004,0.0004,0.0008, 0.002, 0.006,0.009,0.012,0.015,0.022,0.08,0.5,0.9,0.9,0.9,0.9,0.4,0.1,0.06,0.15,0.2,
												   0.001,0.00003,0.0003,0.0004,0.0004,0.0003,0.0003,0.0004,0.0004, 0.001,0.0025,0.0075,0.012,0.016, 0.02,0.029,0.08,0.5,0.9,0.9,0.9,0.9,0.4,0.1,0.06,0.15,0.2,
												   0.001,0.00003,0.0003,0.0004,0.0004,0.0003,0.0003,0.0004,0.0004,0.0012, 0.003, 0.009,0.014,0.019,0.024,0.035,0.08,0.5,0.9,0.9,0.9,0.9,0.4,0.1,0.06,0.15,0.2};
double CParticleRadiationProperty::pdiam[RP_NPD] = {0.5,1,2,4,6,8,10,12.5,15,17.5,20,25,30,40,50,75,100,200,300};
double CParticleRadiationProperty::ptemp[RP_NPT] = {250,400,550,700,850,1000,1150,1300,1450,1600,1800,2000,2500,3000};

void CParticleRadiationProperty::Mie(double x, double refrelr, double refreli, int nang, double& qext, double& qsca, double& qback)
{
	int j, n, nn, nstop, nmx;
	double psi0, psi1, psi, dn, dx;
	double xstop;
	double ymod;
	double dang;
	double rn;
	double amu[100],theta[100],pi[100],tau[100],pi0[100],pi1[100];
	std::complex<double> d[3000], y, an, bn, s1[200], s2[200];
	std::complex<double> refrel(refrelr,refreli);
   	dx = x;
   	y = x*refrel;
	xstop = x+4.0*pow(x,0.3333)+2.0;
	nstop = (int)xstop;
	ymod = std::abs(y);
	if(xstop>ymod)
		nmx = (int)xstop + 15;
	else
		nmx = (int)ymod + 15;
	dang = CT_PI/2/(nang-1);
	for(j=0; j<nang; j++)
	{
		theta[j] = j*dang;
		amu[j] = cos(theta[j]);
	}
	std::complex<double> cmplx(0.0,0.0);
	d[nmx-1] = cmplx;
	nn = nmx-1;
	for(n=0; n<nn; n++)
	{
		rn=(double)(nmx-n);
		std::complex<double> one1(1.0,0.0), num2(rn,0.0);
		d[nmx-n-2] = (num2/y)-(one1/(d[nmx-n-1] + num2/y));
	}
	for(j=0; j<nang; j++)
	{
		pi0[j] = 0.0;
		pi1[j] = 1.0;
	}
	nn = 2*nang-1;
	for(j=0; j<nn; j++)
	{
		s1[j] = cmplx;
		s2[j] = cmplx;
	}
	psi0 = cos(dx);
	psi1 = sin(dx);
	double chi0 = -sin(x);
	double chi1 = cos(x);
	double apsi0 = psi0;
	double apsi1 = psi1;
	double fn;
	double apsi;
	double chi;
	double p;
	double t;
	int jj;
	std::complex<double> xi0(apsi0,-chi0);
	std::complex<double> xi1(apsi1,-chi1);
	qsca=0.0;
	n=1;
	do
	{
		dn = (double)n;
		rn = (double)n;
		fn=(2.0*rn+1.0)/(rn*(rn+1.0));
		psi=(2.0*dn-1.0)*psi1/dx-psi0;
		apsi=psi;
		chi=(2.0*rn-1.0)*chi1/x-chi0;
		std::complex<double> xi(apsi,-chi);
		an=(d[n-1]/refrel+rn/x)*apsi-apsi1;
		an=an/((d[n-1]/refrel+rn/x)*xi-xi1);
		bn=(refrel*d[n-1]+rn/x)*apsi-apsi1;
		bn=bn/((refrel*d[n-1]+rn/x)*xi-xi1);
		qsca=qsca+(2.0*rn+1.0)*(std::abs(an)*std::abs(an)+std::abs(bn)*std::abs(bn));
		for(j=0; j<nang; j++)
		{
			jj=2*nang-j-2;
			pi[j]=pi1[j];
			tau[j]=rn*amu[j]*pi[j]-(rn+1.0)*pi0[j];
			p=pow(-1.0,n-1);
			s1[j]=s1[j]+fn*(an*pi[j]+bn*tau[j]);
			t=pow(-1.0,n);
			s2[j]=s2[j]+fn*(an*tau[j]+bn*pi[j]);
			if(j!=jj)
			{
				s1[jj]=s1[jj]+fn*(an*pi[j]*p+bn*tau[j]*t);
				s2[jj]=s2[jj]+fn*(an*tau[j]*t+bn*pi[j]*p);
			}
		}
		psi0=psi1;
		psi1=psi;
		apsi1=psi1;
		chi0=chi1;
		chi1=chi;
		std::complex<double> cmplx1(apsi1,-chi1);
		xi1=cmplx1;
		n=n+1;
		rn=(double)n;
		for(j=0; j<nang; j++)
		{
			pi1[j]=((2.0*rn-1.0)/(rn-1.0))*amu[j]*pi[j];
			pi1[j]=pi1[j]-rn*pi0[j]/(rn-1.0);
			pi0[j] = pi[j];
		}
	}while(n-1-nstop<0);
	qsca=(2.0/(x*x))*qsca;
	qext=(4.0/(x*x))*s1[0].real();
	qback=(4.0/(x*x))*std::abs(s1[2*nang-2])*std::abs(s1[2*nang-2]);
}

double CParticleRadiationProperty::Planck(double temp, double* q)
{
	int i;
	int iw;
	int nh = 5000;
	double c2 = 14400;
	double a = 0.01;
	double b = 50;
	double c2dt = c2/temp;
	double h = (b-a)/nh;
	double fa = q[26]*a*a*a/(exp(a)-1);
	double fb = q[0]*b*b*b/(exp(b)-1);
	double xi1 = 0;
	double xi2 = 0;
	double z;
	double wl;
	double fx;
	for(i=1; i<nh; i++)
	{
		z = a + i*h;
		wl = c2dt/z;
		if(wl<=0.25) iw = 0;
		else if(wl>=12.75) iw = 26;
		else iw = (int)((wl - 0.25)/0.5) + 1;
		fx = q[iw]*z*z*z/(exp(z)-1);
		if (i%2)
			xi1 += fx;
		else
			xi2 += fx;
	}
	return 15/CT_PI4*h/3*(fa + fb + 2*xi2 + 4*xi1);
}

double CParticleRadiationProperty::Rosseland(double temp, double* q)
{
	int i;
	int iw;
	int nh = 5000;
	double z;
	double wl;
	double ez;
	double ez1;
	double fx;
	double c2 = 14400;
	double a = 0.01;
	double b = 50;
	double c2dt = c2/temp;
	double h = (b-a)/nh;
	ez = exp(a);
	ez1 = ez-1;
	double fa = 1/q[26]*a*a*a*a*ez/ez1/ez1;
	ez = exp(b);
	ez1 = ez-1;
	double fb = 1/q[0]*b*b*b*b*ez/ez1/ez1;
	double xi1 = 0;
	double xi2 = 0;
	for(i=1; i<nh; i++)
	{
		z = a + i*h;
		wl = c2dt/z;
		if(wl<=0.25) iw = 0;
		else if(wl>=12.75) iw = 26;
		else iw = (int)((wl - 0.25)/0.5) + 1;
		ez = exp(z);
		ez1 = ez-1;
		fx = 1/q[iw]*z*z*z*z*ez/ez1/ez1;
		if (i%2)
			xi1 += fx;
		else
			xi2 += fx;
	}
	return 1/(15/4/CT_PI4*h/3*(fa + fb + 2*xi2 + 4*xi1));
}

void CParticleRadiationProperty::SetUpCoalEfficiencyTable()
{
	int i, j, k;
	int nang = 2;			//number of angles for phase function calculation by calling Mie(), does not affect extinction and scattering efficiency
	double qa[RP_NWL];		//spectral absorption efficiency
	double refre;			//real part of refrective index
	double refim;			//imaginary part of refrective index
	double qap;				//Planck average absorption efficiency
	double qar;				//Rosseland average absorption efficiency
	double temp;			//particle temperature
	double wavel;			//wavelength in microm
	double x;				//size parameter
	double qext;			//extinction efficiency
	double qsca;			//scattering efficiency
	double qback;			//back efficiency
	double diam = 0.5;		//diameter in microm
	//calculate spectral efficiencies first since refractive index is independent of temperature for coal
	refre = refre_coal;
	refim = refim_coal;
	for(j=0; j<RP_NPD; j++)
	{
		diam = pdiam[j];
		wavel = 0.25;	//shortest wavelength
		for(i=0; i<RP_NWL; i++)
		{
			x = CT_PI*diam/wavel;
			if (x>1500) x = 1500;
			Mie(x,refre,refim,nang,qext,qsca,qback);
			qa[i] = qext - qsca;
			if(qa[i]<0) qa[i] = 1.0e-6;
			if (i==0)
				wavel = wavel + 0.25;
			else
				wavel = wavel + 0.5;
		}
		for(k=0; k<RP_NPT; k++)
		{
			temp = ptemp[k];
			qap = Planck(temp,qa);
			qar = Rosseland(temp,qa);
			qa_coal[j][k] = (qap+qar)/2;
		}
	}
}


void CParticleRadiationProperty::SetUpAshEfficiencyTable()
{
	int i, j, k;
	int itmp;				//temperature range index for ash refrective index lookup
	int nang = 2;			//number of angles for phase function calculation by calling Mie(), does not affect extinction and scattering efficiency
	double qa[RP_NWL];		//spectral absorption efficiency
	double refre;			//real part of refrective index
	double refim;			//imaginary part of refrective index
	double qap;				//Planck average absorption efficiency
	double qar;				//Rosseland average absorption efficiency
	double temp;			//particle temperature
	double wavel;			//wavelength in microm
	double x;				//size parameter
	double qext;			//extinction efficiency
	double qsca;			//scattering efficiency
	double qback;			//back efficiency
	double diam = 0.5;		//diameter in microm
	for(j=0; j<RP_NPD; j++)
	{
		diam = pdiam[j];
		for(k=0; k<RP_NPT; k++)
		{
			temp = ptemp[k];
			wavel = 0.25;	//shortest wavelength
			if (temp<=1290)
				itmp = 0;
			else
			{
				if(temp<=1650)
					itmp = 1;
				else
					itmp = 2;
			}
			for(i=0; i<RP_NWL; i++)
			{
				refre = refre_ash[i];
				refim = refim_ash[itmp][i];
				x = CT_PI*diam/wavel;
				if (x>1500) x = 1500;
				Mie(x,refre,refim,nang,qext,qsca,qback);
				qa[i] = qext - qsca;
				if(qa[i]<0.0) qa[i] = 1.0e-6;
				if (i==0)
					wavel = wavel + 0.25;
				else
					wavel = wavel + 0.5;
			}
			qap = Planck(temp,qa);
			qar = Rosseland(temp,qa);
			qa_ash[j][k] = (qap+qar)/2;
		}
	}
}

int CParticleRadiationProperty::FindLowerIndexAndInterpolationFactor(int nsize, double x, double* px, double& f)
{
	//using binary search to find the lower index of the interpolation region and the intepolation factor for lower point
	//the factor for the upper index should be 1-f
	//nsize is the size of the sorted ascending array
	int i;				//midpoint index
	int i0 = 0;			//lower end index
	int i1 = nsize-1;	//upper end index
	if (x<=px[0])
	{
		f = 1;
		return 0;
	}
	if (x>=px[i1])
	{
		f = 0;
		return nsize-2;
	}
	while (true)
	{
		if (i1-i0<=1)
		{
			f = (px[i1]-x)/(px[i1]-px[i0]);
			return i0;
		}
		i = (i0+i1)/2;
		if (x<px[i])
			i1 = i;
		else
			i0 = i;
	}

}

double CParticleRadiationProperty::GetCoalAbsorptionEfficiency(double d, double t)
{
	//interpolate absorption efficiency based on particle diameter d and temperature t
	int id0;			//lower index for diameter
	int it0;			//lower index for temperature
	int id1;			//upper index for diameter
	int it1;			//upper index for temperature
	double fd0;			//lower point interpolation factor for diameter
	double ft0;			//lower point interpolation factor for temperature
	double fd1;			//upper point factor for diamter
	double ft1;			//upper point factor for temperature
	double dm = 1e6*d;	//convert d from m to microm
	id0 = FindLowerIndexAndInterpolationFactor(RP_NPD,dm,pdiam,fd0);
	it0 = FindLowerIndexAndInterpolationFactor(RP_NPT,t,ptemp,ft0);
	id1 = id0+1;
	it1 = it0+1;
	fd1 = 1-fd0;
	ft1 = 1-ft0;
	return fd0*(qa_coal[id0][it0]*ft0 + qa_coal[id0][it1]*ft1) + fd1*(qa_coal[id1][it0]*ft0 + qa_coal[id1][it1]*ft1);
}

double CParticleRadiationProperty::GetAshAbsorptionEfficiency(double d, double t)
{
	//interpolate absorption efficiency based on particle diameter d and temperature t
	int id0;			//lower index for diameter
	int it0;			//lower index for temperature
	int id1;			//upper index for diameter
	int it1;			//upper index for temperature
	double fd0;			//lower point interpolation factor for diameter
	double ft0;			//lower point interpolation factor for temperature
	double fd1;			//upper point factor for diamter
	double ft1;			//upper point factor for temperature
	double dm = 1e6*d;	//convert d from m to microm
	id0 = FindLowerIndexAndInterpolationFactor(RP_NPD,dm,pdiam,fd0);
	it0 = FindLowerIndexAndInterpolationFactor(RP_NPT,t,ptemp,ft0);
	id1 = id0+1;
	it1 = it0+1;
	fd1 = 1-fd0;
	ft1 = 1-ft0;
	return fd0*(qa_ash[id0][it0]*ft0 + qa_ash[id0][it1]*ft1) + fd1*(qa_ash[id1][it0]*ft0 + qa_ash[id1][it1]*ft1);
}

CParticleRadiationProperty::CParticleRadiationProperty()
{

}

CParticleRadiationProperty::~CParticleRadiationProperty()
{

}

CParticleRadiationProperty* CParticleRadiationProperty::GetInstance()
{
	if (pinstance==NULL)
	{
		pinstance = new CParticleRadiationProperty;
		pinstance->SetUpCoalEfficiencyTable();
		pinstance->SetUpAshEfficiencyTable();
	}
	return pinstance;
}

void CParticleRadiationProperty::DeleteInstance()
{
	if (pinstance!=NULL)
	{
		delete pinstance;
		pinstance = NULL;
	}
}
