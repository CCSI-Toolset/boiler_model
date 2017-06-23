// QuadFace.h: interface for the CQuadFace class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __QUADFACE_H__
#define __QUADFACE_H__

#include "Vertex.h"
#include "HexCell.h"
#include "BoundaryFaceProperty.h"

//we may need to replace CBoundaryProperty object by a property mapping method object
class CQuadFace
{
public:
	char		itype;				//face type, 0=flow(interior), 1=wall, 2=inlet, 3=outlet, 4=symmetry
	short		nvertex;			//number of vertices, always 4
	int			izone;				//zone index
	int			ivertex[4];			//indices of vertices that form the face
	int			icell[2];			//indices of 2 cells that share the face, icell[1]==-1 for boundary face
	
	double		area;				//face area
	double		x[3];				//face centroid x, y, z value
	double		normal[3];			//normal unit vector of the face, always points to the direction away from icell[0]

	CBoundaryFaceProperty* pbfp;	//pointer to a boundary face property
	
	CQuadFace();
	virtual ~CQuadFace();
	void CopyPartialData(CQuadFace* pf);
	void CalcCentroidX(CVertex* pvertex);
	void CalcNormalVectorAndArea(CVertex* pvertex);
	void UpdateNormal(CVertex* pvertex, CHexCell* pcell);

	double GetWallTemperature();
	double GetConvectiveHeatFlux();
	double GetRadiationIncidentHeatFlux();
	double GetRadiationNetHeatFlux();
	double GetRadiationBoundaryFactor();
	void Write(FILE* pf);
	void Read(FILE* pf);
};

#endif
