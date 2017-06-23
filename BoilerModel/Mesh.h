// Mesh.h: interface for the CMesh class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __MESH_H__
#define __MESH_H__

#include "Vertex.h"
#include "QuadFace.h"
#include "HexCell.h"

class CMesh
{
public:
	//data related to geometry and connectivity
	int		nvertex;		//total number of vertices
	int		nface;			//total number of faces
	int		ncell;			//total number of cells
	int		nbface;			//total number of boundary faces
	double	area_encl;		//total area of enclosure wall modeled
	double	area_sh;		//total area of superheater wall modeled
	double	area_exit;		//total area of total exit boundary modeled
	double  area_total;		//total boundary area modeled
	double	volume;			//total volume of the modeled
	double	mbl;			//meam beam length

	int*		pibface;		//array of boundary face indices
	CVertex*	pvertex;		//array of vertices
	CQuadFace*	pface;			//array of faces
	CHexCell*	pcell;			//array of cells

	CMesh();
	virtual ~CMesh();
	void Write(FILE* pf);
	void Read(FILE* pf);
	void AllocateArrays();
	void DeleteArrays();
	void AddBoundaryFaces(int n, CQuadFace* pf, int* pif_conv);
	void UpdateCellGeometry();
	void UpdateFaceGeomerty();
	void CalcZoneVolumes(int nz, double* pvz);
	void CalcZoneWallArea(int nz, double* pawz);
	void CalcZoneWallAreaTemperatureProduct(int nz, double* patwz);
	void CalcZoneWallIncidentFlux(int nz, double* pavg, double* pmin, double* pmax);
	double CalcMBL();
	bool IsGeometrySet();
};

#endif