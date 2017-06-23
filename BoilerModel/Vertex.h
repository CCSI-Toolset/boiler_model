// Vertex.h: interface for the CVertex class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __VERTEX_H__
#define __VERTEX_H__

#include <cstdio>

class CVertex
{
public:
	double	x[3];			//x, y, z coordinates of a vertex in an unstructured mesh

	CVertex();
	virtual ~CVertex();
	void Write(FILE* pf);
	void Read(FILE* pf);
};

#endif