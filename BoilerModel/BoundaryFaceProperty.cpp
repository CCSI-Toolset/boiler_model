// BoundaryFaceProperty.cpp: implementation of the CBoundaryFaceProperty class.
//
//////////////////////////////////////////////////////////////////////
#include "BoundaryFaceProperty.h"


CBoundaryFaceProperty::CBoundaryFaceProperty()
{
	iwall = 0;			//for water wall
	iwater = 0;			//for water wall
	hconv = 0;	
	emis = 0.6;	
	rwall = 0.001;	
	tback = 500;	
	temp = 500;
	tgas = 500;
	fradadj = 1;
	qinc = 1e-10;		//need small amount for solving radiation equation
	qnet = 0;	
	qconv = 0;	
}

CBoundaryFaceProperty::~CBoundaryFaceProperty()
{

}

void CBoundaryFaceProperty::Write(FILE* pf)
{
	int iversion = 0;
	fwrite(&iversion, sizeof(int), 1, pf);
	fwrite(&iwall, sizeof(int), 1, pf);
	fwrite(&iwater, sizeof(int), 1, pf);
	fwrite(&hconv, sizeof(double), 1, pf);
	fwrite(&emis, sizeof(double), 1, pf);
	fwrite(&rwall, sizeof(double), 1, pf);
	fwrite(&tback, sizeof(double), 1, pf);
	fwrite(&temp, sizeof(double), 1, pf);
	fwrite(&tgas, sizeof(double), 1, pf);
	fwrite(&fradadj, sizeof(double), 1, pf);
	fwrite(&qinc, sizeof(double), 1, pf);
	fwrite(&qnet, sizeof(double), 1, pf);
	fwrite(&qconv, sizeof(double), 1, pf);
}

void CBoundaryFaceProperty::Read(FILE* pf)
{
	int iversion;
	fread(&iversion, sizeof(int), 1, pf);
	fread(&iwall, sizeof(int), 1, pf);
	fread(&iwater, sizeof(int), 1, pf);
	fread(&hconv, sizeof(double), 1, pf);
	fread(&emis, sizeof(double), 1, pf);
	fread(&rwall, sizeof(double), 1, pf);
	fread(&tback, sizeof(double), 1, pf);
	fread(&temp, sizeof(double), 1, pf);
	fread(&tgas, sizeof(double), 1, pf);
	fread(&fradadj, sizeof(double), 1, pf);
	fread(&qinc, sizeof(double), 1, pf);
	fread(&qnet, sizeof(double), 1, pf);
	fread(&qconv, sizeof(double), 1, pf);
}
