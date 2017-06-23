// HexCell.h: interface for the CHexCell class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __HEXCELL_H__
#define __HEXCELL_H__

#include "Vertex.h"

class CHexCell
{
public:
	short	nvertex;		//number of vertices, determine polyhedron type
	short	nface;			//total number of faces including boundary faces
	short	ncell_nbr;		//number of neighboring cells (equal to nface if there is no boundary face)
	int		izone;			//zone index
	int		pivertex[8];	//vertex indices for base points
	int		piface[6];		//face indices, interior faces followed by boundary faces, order same as picell for the first ncell_nbr faces
	int		picell[6];		//cell indices, neighbor cells, -1 if the corresponding face is boundary face

	double		x[3];		//centroid x, y, z location
	double		vol;		//volume of the cell

	CHexCell();
	virtual ~CHexCell();
	void AddFaceAndNeighborCell(int iface, int icell);
	void CalcCentroidX(CVertex* pvertex);
	double GetTetrahedronVolume(short iv[4], CVertex* pvertex);
	void CalcVolume(CVertex* pvertex);
	bool HasBoundary() {return ncell_nbr<nface;}
	void Write(FILE* pf);
	void Read(FILE* pf);
};

#endif
