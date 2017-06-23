// Vertex.cpp: implementation of the CVertex class.
//
//////////////////////////////////////////////////////////////////////
#include "Vertex.h"


CVertex::CVertex()
{

}

CVertex::~CVertex()
{

}

void CVertex::Write(FILE* pf)
{
	fwrite(x,sizeof(double),3,pf);
}

void CVertex::Read(FILE* pf)
{
	fread(x,sizeof(double),3,pf);
}
