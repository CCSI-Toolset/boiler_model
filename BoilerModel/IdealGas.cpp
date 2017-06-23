// IdealGas.cpp: implementation of the CIdealGas class.
//
//////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IdealGas.h"
#include "Constants.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
double CIdealGas::ljkte[79] = {0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100};
double CIdealGas::ljomega_mk[79] = {2.785,2.628,2.492,2.368,2.257,2.156,2.065,1.982,1.908,1.841,1.78,1.725,1.675,1.629,1.587,1.549,1.514,1.482,1.452,1.424,1.399,1.375,1.353,1.333,1.314,1.296,1.279,1.264,1.248,1.234,1.221,1.209,1.197,1.186,1.175,1.156,1.138,1.122,1.107,1.093,1.081,1.069,1.058,1.048,1.039,1.03,1.022,1.014,1.007,0.9999,0.9932,0.987,0.9811,0.9755,0.97,0.9649,0.96,0.9553,0.9507,0.9464,0.9422,0.9382,0.9343,0.9305,0.9269,0.8963,0.8727,0.8538,0.8379,0.8242,0.7432,0.7005,0.6718,0.6504,0.6335,0.6194,0.6076,0.5973,0.5882};
double CIdealGas::ljomega_d[79] = {2.662,2.476,2.318,2.184,2.066,1.966,1.877,1.798,1.729,1.667,1.612,1.562,1.517,1.476,1.439,1.406,1.375,1.346,1.32,1.296,1.273,1.253,1.233,1.215,1.198,1.182,1.167,1.153,1.14,1.128,1.116,1.105,1.094,1.084,1.075,1.057,1.041,1.026,1.012,0.9996,0.9878,0.977,0.9672,0.9576,0.949,0.9406,0.9328,0.9256,0.9186,0.912,0.9058,0.8998,0.8942,0.8888,0.8836,0.8788,0.874,0.8694,0.8652,0.861,0.8568,0.853,0.8492,0.8456,0.8422,0.8124,0.7896,0.7712,0.7556,0.7424,0.664,0.6232,0.596,0.5756,0.5596,0.5464,0.5352,0.5256,0.517};

CIdealGas** CIdealGas::pGasList = NULL;

CIdealGas::CIdealGas()
{
	//currently all the JANAF table threshold temperature is 1000 K (TMID in pcgc3)
	iformat = 0;
	tempth = 1000;
	//set default as N2
	name = "N2";
	AssignData();
}

CIdealGas::CIdealGas(std::string _name)
{
	iformat = 0;
	tempth = 1000;
	name = _name;
	AssignData();
}

CIdealGas::~CIdealGas()
{
}

CIdealGas::CIdealGas(const CIdealGas &tt)
{
	int i, j;
	name = tt.name;
	tempth = tt.tempth;
	mw = tt.mw;
	cp = tt.cp;
	h = tt.h;
	s = tt.s;
	g = tt.g;
	vis = tt.vis;
	cond = tt.cond;
	ljsigma = tt.ljsigma;
	ljek = tt.ljek;
	iformat = tt.iformat;
	nelem = tt.nelem;
	for (i=0; i<GS_NE; i++)
	{
		matom[i] = tt.matom[i];
		natom[i] = tt.natom[i];
		elem[i] = tt.elem[i];
	}
	for (i=0; i<2; i++)
	{
		for (j=0; j<2; j++)
			x0[i][j] = tt.x0[i][j];
		for (j=0; j<7; j++)
		{
			z[i][j] = tt.z[i][j];
		}
	}
}


CIdealGas& CIdealGas::operator=(const CIdealGas& tt)
{
	if (this==&tt)
		return *this;
	int i, j;
	name = tt.name;
	tempth = tt.tempth;
	mw = tt.mw;
	cp = tt.cp;
	h = tt.h;
	s = tt.s;
	g = tt.g;
	vis = tt.vis;
	cond = tt.cond;
	ljsigma = tt.ljsigma;
	ljek = tt.ljek;
	iformat = tt.iformat;
	nelem = tt.nelem;
	for (i=0; i<GS_NE; i++)
	{
		matom[i] = tt.matom[i];
		natom[i] = tt.natom[i];
		elem[i] = tt.elem[i];
	}
	for (i=0; i<2; i++)
	{
		for (j=0; j<2; j++)
			x0[i][j] = tt.x0[i][j];
		for (j=0; j<7; j++)
		{
			z[i][j] = tt.z[i][j];
		}
	}
	return *this;
}

void CIdealGas::SetName(std::string str)
{
	name = str;
	AssignData();
}

void CIdealGas::AssignData()
//contains a list of existing data
{
	int i, j;
	mw = 0;
	iformat = 0;
	ljsigma = 3.681;
	ljek = 91.5;
	x0[0][0] = x0[0][1] = x0[1][0] = x0[1][1] = 1;
	for (i=0; i<2; i++)
	{
		for (j=0; j<7; j++)
			z[i][j] = 0;
	}
	for (i=0; i<GS_NE; i++)
	{
		natom[i] = 0;
		matom[i] = 0;
		elem[i] = ' ';
	}
	if (name=="C(S)")
	{
		//element 13 is modified slightly to make sure h at 298.15 is zero
		double d[2][7]={ 0.13604937E+01, 0.19182237E-02,-0.84040386E-06, 0.16448705E-09,-0.11672670E-13,
						-0.65713843E+03,-0.80070200E+01,-0.44778049E+00, 0.53690970E-02,-0.39775563E-06,
						-0.40459263E-08, 0.21134925E-11,-0.94622075E+02, 0.16840773E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=1;
		nelem=1;
		elem[0] = 'C';
		matom[0] = 1;
		CalcMW();
		x0[0][0] = -30;
		x0[0][1] = -10;
		x0[1][0] = -30;
		x0[1][1] = -10;
		return;
	}
	if (name=="O2")
	{
		//element 13 is modified slightly to make sure h at 298.15 is zero
		double d[2][7]={ 0.36219521E+01, 0.73618256E-03,-0.19652219E-06, 0.36201556E-10,-0.28945623E-14,
						-0.12019822E+04, 0.36150942E+01, 0.36255980E+01,-0.18782183E-02, 0.70554543E-05,
						-0.67635071E-08, 0.21555977E-11,-0.10474773E+04, 0.43052769E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[2]=2;
		nelem=1;
		elem[0] = 'O';
		matom[0] = 2;
		ljsigma = 3.433;
		ljek = 113.0;
		CalcMW();
		x0[0][0] = 1;
		x0[0][1] = -10;
		x0[1][0] = 1;
		x0[1][1] = -10;
		return;
	}
	if (name=="H2")
	{
		double d[2][7]={ 0.31001883E+01, 0.51119458E-03, 0.52644204E-07,-0.34909964E-10, 0.36945341E-14,
						-0.87738013E+03,-0.19629412E+01, 0.30574446E+01, 0.26765198E-02,-0.58099149E-05,
						 0.55210343E-08,-0.18122726E-11,-0.98926469E+03,-0.22997046E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[1]=2;
		nelem=1;
		elem[0] = 'H';
		matom[0] = 2;
		ljsigma = 2.915;
		ljek = 38.0;
		CalcMW();
		x0[0][0] = -10;
		x0[0][1] = 1;
		x0[1][0] = -10;
		x0[1][1] = 1;
		return;
	}
	if (name=="CO")
	{
		double d[2][7]={ 0.29840689E+01, 0.14891387E-02,-0.57899678E-06, 0.10364576E-09,-0.69353499E-14,
						-0.14245227E+05, 0.63479147E+01, 0.37100916E+01,-0.16190964E-02, 0.36923584E-05,
						-0.20319673E-08, 0.23953344E-12,-0.14356309E+05, 0.29555340E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=1;
		natom[2]=1;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 1;
		elem[1] = 'O';
		matom[1] = 1;
		ljsigma = 3.59;
		ljek = 110.0;
		CalcMW();
		x0[0][0] = -10;
		x0[0][1] = 1;
		x0[1][0] = -10;
		x0[1][1] = 1;
		return;
	}
	if (name=="CO2")
	{
		double d[2][7]={ 0.44608040E+01, 0.30981717E-02,-0.12392566E-05, 0.22741323E-09,-0.15525948E-13,
						-0.48961438E+05,-0.98635978E+00, 0.24007788E+01, 0.87350905E-02,-0.66070861E-05,
						 0.20021860E-08, 0.63274039E-15,-0.48377520E+05, 0.96951447E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=1;
		natom[2]=2;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 1;
		elem[1] = 'O';
		matom[1] = 2;
		ljsigma = 3.996;
		ljek = 190.0;
		CalcMW();
		x0[0][0] = 1;
		x0[0][1] = -10;
		x0[1][0] = 1;
		x0[1][1] = -10;
		return;
	}
	if (name=="H2O")
	{
		double d[2][7]={ 0.27167616E+01, 0.29451370E-02,-0.80224368E-06, 0.10226681E-09,-0.48472104E-14,
						-0.29905820E+05, 0.66305666E+01, 0.40701275E+01,-0.11084499E-02, 0.41521180E-05,
						-0.29637404E-08, 0.80702101E-12,-0.30279719E+05,-0.32270038E+00};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[1]=2;
		natom[2]=1;
		nelem=2;
		elem[0] = 'H';
		matom[0] = 2;
		elem[1] = 'O';
		matom[1] = 1;
		ljsigma = 2.641;
		ljek = 809.1;
		CalcMW();
		return;
	}
	if (name=="H2O(L)")
	{
		double d[2][7]={-1.00000000E+04, 0.00000000E+00, 0.0           , 0.0           , 0.0,
						 3.69780000E+06, 5.92300000E+04, 0.12712782E+02,-0.17662790E-01,-0.22556633E-04,
						 0.20820897E-06,-0.24078606E-09,-0.37483195E+05,-0.59115326E+02};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[1]=2;
		natom[2]=1;
		nelem=2;
		elem[0] = 'H';
		matom[0] = 2;
		elem[1] = 'O';
		matom[1] = 1;
		ljsigma = 2.641;
		ljek = 809.1;
		CalcMW();
		tempth = 373.15;	//different temperature threshold
		x0[0][0] = 1;
		x0[0][1] = 1;
		x0[1][0] = -70;
		x0[1][1] = -70;
		return;
	}
	if (name=="SO2")
	{
		double d[2][7]={ 0.52451364E+01, 0.19704204E-02,-0.80375769E-06, 0.15149969E-09,-0.10558004E-13,
						-0.37558227E+05,-0.10873524E+01, 0.32665338E+01, 0.53237902E-02, 0.68437552E-06,
						-0.52810047E-08, 0.25590454E-11,-0.36908148E+05, 0.96513476E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[2]=2;
		natom[4]=1;
		nelem=2;
		elem[0] = 'S';
		matom[0] = 1;
		elem[1] = 'O';
		matom[1] = 2;
		ljsigma = 4.29;
		ljek = 252.0;
		CalcMW();
		return;
	}
	if (name=="H2S")
	{
		double d[2][7]={ 0.28479103E+01, 0.38415990E-02,-0.14099367E-05, 0.24278754E-09,-0.15783283E-13,
						-0.34469788E+04, 0.74781412E+01, 0.38811293E+01,-0.13211856E-03, 0.36517726E-05,
						-0.21820445E-08, 0.28783779E-12,-0.36350917E+04, 0.25161511E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[1]=2;
		natom[4]=1;
		nelem=2;
		elem[0] = 'H';
		matom[0] = 2;
		elem[1] = 'S';
		matom[1] = 1;
		ljsigma = 3.49;
		ljek = 343.0;
		CalcMW();
		x0[0][0] = -50;
		x0[0][1] = -10;
		x0[1][0] = -50;
		x0[1][1] = -10;
		return;
	}
	if (name=="COS")
	{
		double d[2][7]={ 0.52392000E+01, 0.24100584E-02,-0.96064522E-06, 0.17778347E-09,-0.12235704E-13,
						-0.18480455E+05,-0.30910517E+01, 0.24625321E+01, 0.11947992E-01,-0.13794370E-04,
						 0.80707736E-08,-0.18327653E-11,-0.17803987E+05, 0.10792556E+02};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=1;
		natom[2]=1;
		natom[4]=1;
		nelem=3;
		elem[0] = 'C';
		matom[0] = 1;
		elem[1] = 'O';
		matom[1] = 1;
		elem[2] = 'S';
		matom[2] = 1;
		ljsigma = 4.13;
		ljek = 335.0;
		CalcMW();
		x0[0][0] = -20;
		x0[0][1] = -20;
		x0[1][0] = -20;
		x0[1][1] = -20;
		return;
	}
	if (name=="CH4")
	{
		double d[2][7]={ 0.15027056E+01, 0.10416795E-01,-0.39181514E-05, 0.67777872E-09,-0.44283706E-13,
						-0.99787031E+04, 0.10707143E+02, 0.38261929E+01,-0.39794557E-02, 0.24558321E-04,
						-0.22732920E-07, 0.69626952E-11,-0.10144945E+05, 0.86690062E+00};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=1;
		natom[1]=4;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 1;
		elem[1] = 'H';
		matom[1] = 4;
		ljsigma = 3.822;
		ljek = 137.0;
		CalcMW();
		x0[0][0] = -50;
		x0[0][1] = -10;
		x0[1][0] = -50;
		x0[1][1] = -10;
		return;
	}
	if (name=="C2H2")
	{
		double d[2][7]={ 0.45751083E+01, 0.51238358E-02,-0.17452354E-05, 0.28673065E-09,-0.17951426E-13,
						 0.25607428E+05,-0.35737940E+01, 0.14102768E+01, 0.19057275E-01,-0.24501390E-04,
						 0.16390872E-07,-0.41345447E-11, 0.26188208E+05, 0.11393827E+02};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=2;
		natom[1]=2;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 2;
		elem[1] = 'H';
		matom[1] = 2;
		ljsigma = 4.221;
		ljek = 185.0;
		CalcMW();
		x0[0][0] = -50;
		x0[0][1] = -10;
		x0[1][0] = -50;
		x0[1][1] = -10;
		return;
	}
	if (name=="C2H4")
	{
		double d[2][7]={ 0.34552152E+01, 0.11491803E-01,-0.43651750E-05, 0.76155095E-09,-0.50123200E-13,
						 0.44773119E+04, 0.26987959E+01, 0.14256821E+01, 0.11383140E-01, 0.79890006E-05,
						-0.16253679E-07, 0.67491256E-11, 0.53370755E+04, 0.14621819E+02};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=2;
		natom[1]=4;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 2;
		elem[1] = 'H';
		matom[1] = 4;
		ljsigma = 4.232;
		ljek = 205.0;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -30;
		x0[1][0] = -70;
		x0[1][1] = -30;
		return;
	}
	if (name=="C2H6")
	{
		double d[2][7]={ 1.67107058E+00, 1.88078150E-02,-6.98943156E-06, 1.16385735E-09,-7.17707692E-14,
						-1.14683543E+04, 1.26317347E+01, 1.92453270E+00, 1.68224303E-02,-2.24906498E-06,
						-3.40875417E-09, 1.49239675E-12,-1.14789269E+04, 1.16292438E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=2;
		natom[1]=6;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 2;
		elem[1] = 'H';
		matom[1] = 6;
		ljsigma = 4.418;
		ljek = 230.0;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -40;
		x0[1][0] = -70;
		x0[1][1] = -40;
		return;
	}
	if (name=="C3H6")
	{
		double d[2][7]={ 0.06732257E+02, 0.14908336E-01,-0.04949899E-04, 0.07212022E-08,-0.03766204E-12,
						-0.09235703E+04,-0.13313348E+02, 0.14933071E+01, 0.02092517E+00, 0.04486794E-04,
						-0.16689121E-07, 0.07158146E-10, 0.10748264E+04, 0.16145340E+02};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=3;
		natom[1]=6;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 3;
		elem[1] = 'H';
		matom[1] = 6;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -30;
		x0[1][0] = -70;
		x0[1][1] = -30;
		return;
	}
	if (name=="C3H8")
	{
		double d[2][7]={ 0.07525217E+02, 0.01889034E+00,-0.06283924E-04, 0.09179373E-08,-0.04812410E-12,
						-0.16464548E+05,-0.01784390E+03, 0.08969208E+01, 0.02668986E+00, 0.05431425E-04,
						-0.02126000E-06, 0.09243330E-10,-0.13954918E+05, 0.01935533E+03};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=3;
		natom[1]=8;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 3;
		elem[1] = 'H';
		matom[1] = 8;
		ljsigma = 5.061;
		ljek = 254.0;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -30;
		x0[1][0] = -70;
		x0[1][1] = -30;
		return;
	}
	if (name=="C4H6")
	{
		double d[2][7]={ 0.08046583E+02, 0.16485251E-01,-0.05522227E-04, 0.08123593E-08,-0.04295078E-12,
						 0.13701305E+05,-0.01800457E+03, 0.03197108E+02, 0.02025591E+00, 0.06510192E-04,
						-0.16584423E-07, 0.06400282E-10, 0.15715203E+05, 0.09895660E+02};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=4;
		natom[1]=6;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 4;
		elem[1] = 'H';
		matom[1] = 6;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -30;
		x0[1][0] = -70;
		x0[1][1] = -30;
		return;
	}
	if (name=="C4H8")
	{
		double d[2][7]={ 0.02053584E+02, 0.03435050E+00,-0.15883196E-04, 0.03308966E-07,-0.02536104E-11,
						-0.02139723E+05, 0.15543201E+02, 0.11811380E+01, 0.03085338E+00, 0.05086524E-04,
						-0.02465488E-06, 0.11110192E-10,-0.01790400E+05, 0.02106247E+03};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=4;
		natom[1]=8;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 4;
		elem[1] = 'H';
		matom[1] = 8;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -30;
		x0[1][0] = -70;
		x0[1][1] = -30;
		return;
	}
	if (name=="C4H10")
	{
		double d[2][7]={ 0.01998784E+03, 0.10372807E-01,-0.09610818E-05,-0.04623017E-08, 0.08202828E-12,
						-0.02625571E+06,-0.08837907E+03,-0.02256618E+02, 0.05881732E+00,-0.04525782E-03,
						 0.02037115E-06,-0.04079458E-10,-0.01760233E+06, 0.03329595E+03};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=4;
		natom[1]=10;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 4;
		elem[1] = 'H';
		matom[1] = 10;
		ljsigma = 5.341;
		ljek = 313.0;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -30;
		x0[1][0] = -70;
		x0[1][1] = -30;
		return;
	}
	if (name=="C5H12")
	{
		double d[2][7]={ 0.16677979E+02, 0.02114483E+00,-0.03533321E-04,-0.05742202E-08, 0.15159483E-12,
						-0.02553670E+06,-0.06372940E+03, 0.01877907E+02, 0.04121645E+00, 0.12532337E-04,
						-0.03701536E-06, 0.15255685E-10,-0.02003815E+06, 0.01877256E+03};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=5;
		natom[1]=12;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 5;
		elem[1] = 'H';
		matom[1] = 12;
		ljsigma = 5.769;
		ljek = 345.0;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -30;
		x0[1][0] = -70;
		x0[1][1] = -30;
		return;
	}
	if (name=="C6H6")
	{
		double d[2][7]={ 0.12910740E+02, 0.01723296E+00,-0.05024210E-04, 0.05893497E-08,-0.01947521E-12,
						 0.03664511E+05,-0.05002699E+03,-0.03138012E+02, 0.04723103E+00,-0.02962207E-04,
						-0.03262819E-06, 0.01718691E-09, 0.08890031E+05, 0.03657573E+03};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=6;
		natom[1]=6;
		nelem=2;
		elem[0] = 'C';
		matom[0] = 6;
		elem[1] = 'H';
		matom[1] = 6;
		ljsigma = 5.27;
		ljek = 440.0;
		CalcMW();
		x0[0][0] = -70;
		x0[0][1] = -30;
		x0[1][0] = -70;
		x0[1][1] = -30;
		return;
	}
	if (name=="HCN")
	{
		double d[2][7]={ 0.37068121E+01, 0.33382803E-02,-0.11913320E-05, 0.19992917E-09,-0.12826452E-13,
						 0.14962636E+05, 0.20794904E+01, 0.24513556E+01, 0.87208371E-02,-0.10094203E-04,
						 0.67255698E-08,-0.17626959E-11, 0.15213002E+05, 0.80830085E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[0]=1;
		natom[1]=1;
		natom[3]=1;
		nelem=3;
		elem[0] = 'H';
		matom[0] = 1;
		elem[1] = 'C';
		matom[1] = 1;
		elem[2] = 'N';
		matom[2] = 1;
		ljsigma = 3.63;
		ljek = 569.0;
		CalcMW();
		return;
	}
	if (name=="NH3")
	{
		double d[2][7]={ 0.24165177E+01, 0.61871211E-02,-0.21785136E-05, 0.37599090E-09,-0.24448856E-13,
						-0.64747177E+04, 0.77043482E+01, 0.35912768E+01, 0.49388668E-03, 0.83449322E-05,
						-0.83833385E-08, 0.27299092E-11,-0.66717143E+04, 0.22520966E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[1]=3;
		natom[3]=1;
		nelem=2;
		elem[0] = 'N';
		matom[0] = 1;
		elem[1] = 'H';
		matom[1] = 3;
		ljsigma = 3.15;
		ljek = 358.0;
		CalcMW();
		return;
	}
	if (name=="NO")
	{
		double d[2][7]={ 0.31890000E+01, 0.13382281E-02,-0.52899318E-06, 0.95919332E-10,-0.64847932E-14,
						 0.98283290E+04, 0.67458126E+01, 0.40459521E+01,-0.34181783E-02, 0.79819190E-05,
						-0.61139316E-08, 0.15919076E-11, 0.97453934E+04, 0.29974988E+01};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[2]=1;
		natom[3]=1;
		nelem=2;
		elem[0] = 'N';
		matom[0] = 1;
		elem[1] = 'O';
		matom[1] = 1;
		ljsigma = 3.47;
		ljek = 119.0;
		CalcMW();
		x0[0][0] = -20;
		x0[0][1] = -20;
		x0[1][0] = -8;
		x0[1][1] = -8;
		return;
	}
	if (name=="OH")
	{
		double d[2][7]={ 0.29106417E+01, 0.95931627E-03,-0.19441700E-06, 0.13756646E-10, 0.14224542E-15,
						 0.39353811E+04, 0.54423428E+01, 0.38375931E+01,-0.10778855E-02, 0.96830354E-06,
						 0.18713971E-09,-0.22571089E-12, 0.36412820E+04, 0.49370009E+00};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[1]=1;
		natom[2]=1;
		nelem=2;
		elem[0] = 'O';
		matom[0] = 1;
		elem[1] = 'H';
		matom[1] = 1;
		ljsigma = 3.1;
		ljek = 60.0;
		CalcMW();
		x0[0][0] = -50;
		x0[0][1] = -50;
		x0[1][0] = -20;
		x0[1][1] = -20;
		return;
	}
	if (name=="CL2")
	{
		double d[2][7]={ 4.68722369E+00,-2.88785951E-04, 2.05481238E-07,-4.90451057E-11, 3.77253004E-15,
						-1.53958856E+03,-1.18221739E-01, 3.59183238E+00, 2.35122034E-03,-2.34527045E-06,
						 1.12156277E-09,-2.07823110E-13,-1.15680719E+03, 5.75327148E+00};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[5]=2;
		nelem=1;
		elem[0] = 'L';
		matom[0] = 2;
		CalcMW();
		x0[0][0] = -5;
		x0[0][1] = -5;
		x0[1][0] = -5;
		x0[1][1] = -5;
		return;
	}
	if (name=="AR")
	{
		double d[2][7]={ 2.50104422E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
						-7.45686320E+02, 4.36103012E+00, 2.50104422E+00, 0.00000000E+00, 0.00000000E+00,
						 0.00000000E+00, 0.00000000E+00,-7.45686320E+02, 4.36103012E+00};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[6]=1;
		nelem=1;
		elem[0] = 'A';
		matom[0] = 1;
		CalcMW();
		x0[0][0] = -5;
		x0[0][1] = -5;
		x0[1][0] = -5;
		x0[1][1] = -5;
		return;
	}
	if (name=="HCL")	//fixed temperature range error when copied from Cantera thermal file
	{
		double d[2][7]={ 2.76658840E+00, 1.438188300E-03,-4.69930000E-07, 7.34994080E-11,-4.37311060E-15,
			            -1.19174680E+04, 6.471506290E+00, 3.52481710E+00, 2.99848620E-05,-8.62218910E-07, 
						 2.09797210E-09,-9.865819100E-13,-1.21505090E+04, 2.40892359E+00};
		for (i=0; i<2; i++)
		{
			for (j=0; j<7; j++)
				z[i][j] = d[i][j];
		}
		natom[1] = 1;
		natom[5] = 1;
		nelem=2;
		elem[0] = 'H';
		elem[1] = 'L';
		matom[0] = 1;
		matom[1] = 1;
		ljsigma = 3.1;
		ljek = 60.0;
		CalcMW();
		x0[0][0] = -5;
		x0[0][1] = -5;
		x0[1][0] = -5;
		x0[1][1] = -5;
		return;
	}
	//default nitrogen, element 13 is modified slightly to make sure h at 298.15 is zero
	double d[2][7]={ 0.28963194E+01, 0.15154863E-02,-0.57235275E-06, 0.99807385E-10,-0.65223536E-14,
					-0.90586182E+03, 0.61615143E+01, 0.36748257E+01,-0.12081496E-02, 0.23240100E-05,
					-0.63217520E-09,-0.22577253E-12,-0.10611273E+04, 0.23580418E+01};
	for (i=0; i<2; i++)
	{
		for (j=0; j<7; j++)
			z[i][j] = d[i][j];
	}
	natom[3]=2;
	nelem=1;
	elem[0] = 'N';
	matom[0] = 2;
	ljsigma = 3.681;
	ljek = 91.5;
	CalcMW();
}

double CIdealGas::GetX0(double t, double er)
{
	int i = 0;
	int j = 0;
	if (t>tempth) i = 1;
	if (er>1) j = 1;
	return x0[i][j];
}

double CIdealGas::CalcMW()
{
	mw = CT_AWC*natom[0]+CT_AWH*natom[1]+CT_AWO*natom[2]+CT_AWN*natom[3]+CT_AWS*natom[4]+CT_AWCL*natom[5]+CT_AWAR*natom[6];
	return mw;
}

double CIdealGas::CalcCp(double t)
//Cp is a function of temperature for ideal gas
{
	int i;
	if (t<tempth)
		i = 1;
	else
		i = 0;
	if (iformat==0)		//JANAF format
	{
		cp = CT_UR*(z[i][0]+t*z[i][1]+t*t*z[i][2]+t*t*t*z[i][3]+t*t*t*t*z[i][4]);
	}
	else				//NIST format
	{	//original formula in J/mol-K, IS unit is J/kmol-K
		t /= 1000;
		cp = 1000*(z[i][0]+t*z[i][1]+t*t*z[i][2]+t*t*t*z[i][3]+z[i][4]/t/t);
	}
	return cp;
}


double CIdealGas::CalcH(double t)
//enthalpy is a function of temperature for ideal gas
{
	int i;
	if (t<tempth)
		i = 1;
	else
		i = 0;
	if (iformat==0)		//JANAF format
	{
		h = CT_UR*(z[i][5]+z[i][0]*t+z[i][1]*t*t/2+z[i][2]*t*t*t/3+z[i][3]*t*t*t*t/4+z[i][4]*t*t*t*t*t/5);
	}
	else				//NIST format
	{	//original formula in KJ/mol, IS unit is J/kmol
		t /= 1000;
		h = 1e6*(z[i][5]+z[i][0]*t+z[i][1]*t*t/2+z[i][2]*t*t*t/3+z[i][3]*t*t*t*t/4-z[i][4]/t);
	}
	return h;
}

double CIdealGas::CalcS(double t, double p)
//entropy is a function of temperature and pressure
//default pressure is 1.013e5 Pa
{
	int i;
	if (t<tempth)
		i = 1;
	else
		i = 0;
	if (iformat==0)		//JANAF format
	{
		s = CT_UR*(z[i][6]+z[i][0]*log(t)+z[i][1]*t+z[i][2]*t*t/2+z[i][3]*t*t*t/3+z[i][4]*t*t*t*t/4 - log(p/101325));
	}
	else				//NIST format
	{	//original formula in J/mol-K, IS unit is J/kmol-K
		t /= 1000;
		s = 1000*(z[i][6]+z[i][0]*log(t)+z[i][1]*t+z[i][2]*t*t/2+z[i][3]*t*t*t/3-z[i][4]/t/t/2) - CT_UR*log(p/101325);
	}
	return s;
}

double CIdealGas::CalcG(double t, double p)
//Gibbs free energy is a function of temperature and pressure
//default pressure is 1.013e5 Pa
{
	//G = H - T*S + 298.15*sum(GAMMAi*Si,298.15)
	//Standard entropy from http://www.codata.org/databases/key1.html
	int i;
	double sum = 0;
	double s0[GS_NE] = {5740, 65340, 102576, 95804.5, 32054, 111500, 154800};  //entropy of C(S), 0.5H2, 0.5O2, 0.5N2, S, 0.5Cl2, Ar at 298.15 K
	CalcH(t);
	CalcS(t,p);
	g = h - t*s;
	for (i=0; i<GS_NE; i++)
	{
		sum += natom[i]*s0[i];
	}
	g = h - t*s + 298.15*sum;
	return g;
}

double CIdealGas::CalcViscosity(double t)
{
	//requires mw, which is always available since it is calculated in AssignData()
	int i;
	double dkte;
	double omega;
	double kte = t/ljek;
	if (kte<=0.3)
		omega = ljomega_mk[0];
	else
	{
		if (kte>=100.0)
			omega = ljomega_mk[78];
		else
		{
			if (kte<=2)
			{
				dkte = 0.05;
				i = (int)((kte-0.3)/0.05);
			}
			else
			{
				if (kte<=5)
				{
					dkte = 0.1;
					i = (int)((kte-2)/0.1)+34;
				}
				else
				{
					if (kte<=10)
					{
						dkte = 1;
						i = (int)(kte-5)+64;
					}
					else
					{
						dkte = 10;
						i = (int)((kte-10)/10)+69;
					}
				}
			}
			omega = ljomega_mk[i] + (kte-ljkte[i])/dkte*(ljomega_mk[i+1]-ljomega_mk[i]);
		}
	}
	vis = 2.69693e-6*sqrt(mw*t)/ljsigma/ljsigma/omega;
	return vis;
}

double CIdealGas::CalcConductivity(double t)
{
	//needs to calculate viscosity (vis) and heat capacity (cp) first
	//requires mw, which is always available since it is calculated in AssignData()
	//based on Eucken formula (Bird, p257), valid for both monoatomic and polyatomic gases at low density
	double cp_at_t = CalcCp(t);
	double vis_at_t = CalcViscosity(t);
	cond = (cp_at_t+1.25*CT_UR)/mw*vis_at_t;
	return cond;
}

double CIdealGas::CalcDiffusivity(CIdealGas* pg, double t, double p)
{
	//requires mw, which is always available since it is calculated in AssignData()
	//pg is the other gas pointer (binary diffusion), p is pressure in Pa
	int i;
	double dif;			//diffusivity to be returned
	double dkte;
	double omega;
	double sigma = (ljsigma+pg->ljsigma)/2;		//average sigma
	double kte = t/sqrt(ljek*pg->ljek);			//average kte
	if (kte<=0.3)
		omega = ljomega_d[0];
	else
	{
		if (kte>=100.0)
			omega = ljomega_d[78];
		else
		{
			if (kte<=2)
			{
				dkte = 0.05;
				i = (int)((kte-0.3)/0.05);
			}
			else
			{
				if (kte<=5)
				{
					dkte = 0.1;
					i = (int)((kte-2)/0.1)+34;
				}
				else
				{
					if (kte<=10)
					{
						dkte = 1;
						i = (int)(kte-5)+64;
					}
					else
					{
						dkte = 10;
						i = (int)((kte-10)/10)+69;
					}
				}
			}
			omega = ljomega_d[i] + (kte-ljkte[i])/dkte*(ljomega_d[i+1]-ljomega_d[i]);
		}
	}
	dif = 0.018829*sqrt(t*t*t*(1/mw+1/pg->mw))/p/sigma/sigma/omega;
	return dif;		//diffusivity in [m2/s]
}

CIdealGas** CIdealGas::GetGasList()
{
	if (pGasList!=NULL)
		return pGasList;
	int i = 0;
	pGasList = new CIdealGas* [GS_NSP];
	CIdealGas* pg;
	pg = new CIdealGas("C(S)");	//0
	pGasList[i++] = pg;
	pg = new CIdealGas("O2");	//1
	pGasList[i++] = pg;
	pg = new CIdealGas("N2");	//2
	pGasList[i++] = pg;
	pg = new CIdealGas("H2");	//3
	pGasList[i++] = pg;
	pg = new CIdealGas("CO");	//4
	pGasList[i++] = pg;
	pg = new CIdealGas("CO2");	//5
	pGasList[i++] = pg;
	pg = new CIdealGas("H2O");	//6
	pGasList[i++] = pg;
	pg = new CIdealGas("H2O(L)");//7
	pGasList[i++] = pg;
	pg = new CIdealGas("SO2");	//8
	pGasList[i++] = pg;
	pg = new CIdealGas("H2S");	//9
	pGasList[i++] = pg;
	pg = new CIdealGas("COS");	//10
	pGasList[i++] = pg;
	pg = new CIdealGas("CH4");	//11
	pGasList[i++] = pg;
	pg = new CIdealGas("C2H2");	//12
	pGasList[i++] = pg;
	pg = new CIdealGas("C2H4");	//13
	pGasList[i++] = pg;
	pg = new CIdealGas("C2H6");	//14
	pGasList[i++] = pg;
	pg = new CIdealGas("C3H6");	//15
	pGasList[i++] = pg;
	pg = new CIdealGas("C3H8");	//16
	pGasList[i++] = pg;
	pg = new CIdealGas("C4H6");	//17
	pGasList[i++] = pg;
	pg = new CIdealGas("C4H8");	//18
	pGasList[i++] = pg;
	pg = new CIdealGas("C4H10");//19
	pGasList[i++] = pg;
	pg = new CIdealGas("C5H12");//20
	pGasList[i++] = pg;
	pg = new CIdealGas("C6H6");	//21
	pGasList[i++] = pg;
	pg = new CIdealGas("HCN");	//22
	pGasList[i++] = pg;
	pg = new CIdealGas("NH3");	//23
	pGasList[i++] = pg;
	pg = new CIdealGas("NO");	//24
	pGasList[i++] = pg;
	pg = new CIdealGas("OH");	//25
	pGasList[i++] = pg;
	pg = new CIdealGas("CL2");	//26
	pGasList[i++] = pg;
	pg = new CIdealGas("AR");	//27
	pGasList[i++] = pg;
	pg = new CIdealGas("HCL");	//28
	pGasList[i++] = pg;
	return pGasList;
}

void CIdealGas::DeleteGasList()
{
	int i;
	if (pGasList!=NULL)
	{
		for (i=0; i<GS_NSP; i++)
			delete pGasList[i];
		delete pGasList;
		pGasList=NULL;
	}
}

CIdealGas* CIdealGas::GetGasByName(std::string str)
{
	int i;
	CIdealGas* pg;
	CIdealGas** plist = GetGasList();
	for (i=0; i<GS_NSP; i++)
	{
		pg = plist[i];
		if (pg->name==str)
			return pg;
	}
	pg = NULL;
	return pg;
}

CIdealGas* CIdealGas::GetGasByIndex(int index)
{
	CIdealGas** plist = GetGasList();
	if (index>=0 && index<GS_NSP)
		return plist[index];
	return NULL;
}

int CIdealGas::GetIndex()
{
	int i;
	CIdealGas** plist = GetGasList();
	for (i=0; i<GS_NSP; i++)
	{
		if (plist[i]==this)
			return i;
	}
	return -1;
}

int CIdealGas::GetIndexByName(std::string str)
{
	int i;
	CIdealGas** plist = GetGasList();
	for (i=0; i<GS_NSP; i++)
	{
		if (plist[i]->name==str)
			return i;
	}
	return -1;
}