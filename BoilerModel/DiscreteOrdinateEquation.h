// DiscreteOrdinateEquation.h: interface for the CDiscreteOrdinateEquation class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __DISCRETEORDINATEEQUATION_H__
#define __DISCRETEORDINATEEQUATION_H__

#include "Mesh.h"

const int RE_NDOS = 20;					//maximum number of discrete ordinates

class CDiscreteOrdinateEquation  
{
public:
	bool		bfirst;					//true if first call to CalcFluxDivergence() after Init()
	int			ndos;					//number of discrete ordinate, currently set at 20
	int			nite_conv;				//number of iterations to reach convergence
	double		urf_twall;				//under-relaxation factor for wall temperature
	double		urf_srcrad;				//under-relaxation factor for radiation
	double		errmax;					//maximum error
	double		errconv;				//relative error allowed at convergence
	double		qcellg;					//total heat absorption by gas phase
	double		qcellp;					//total heat absorption by particle phase
	double		qface;					//total heat absorption from boundary faces
	double		fqinc;					//adjustment factor for qinc to force qface=qcellg+qcellp
	double*		pwt;					//weight of each discrete ordinate
	double		(*pdcos)[3];			//direction cosine of each discrete ordinate
	double**	pint;					//radiation intensity
	double*		parea;					//projected area of faces
	double*		pkag;					//gas absorption coefficient, calculated from narrow-band model
	double*		pkat4g;					//gas kag*temp^4
	double*		pkap;					//particle absorption coefficient, calculated during particle tracking
	double*		pkat4p;					//particle ka*temp^4
	double*		pgir;					//irradiance, integration of intensity around 4*pi solid angle
	double*		pqdivg;					//divergence of radiation flux for gas phase
	double*		pqdivp;					//divergence of radiation flux for particle phase
	int**		pisc;					//sorted cell index for radiation calculation in half of discrete ordinates

	CMesh*		pmesh;					//a pointer to CMesh object

	CDiscreteOrdinateEquation();
	virtual ~CDiscreteOrdinateEquation();
	void AllocateMemory();
	void DeleteMemory();
	void UpdateDirectionCosines();
	void UpdateSortedCellIndices();
	void CalcProjectedFaceArea(int ido);
	void UpdateBoundaryFaceFactor();
	bool CalcBoundaryIncidentHeatFlux();
	void Init();
	void SolvePDE();
	void CalcFluxDivergence();
	void GetZoneRadiationHeatLoss(int nzone, double* phlrad);
	void GetZoneConvectionHeatLoss(int nzone, double* phlconv);
	void CalcBoundaryNetHeatFlux();
	void UpdateWallTemperature();
	void UpdateRadiationProperties(double* pkagz, double* pkapz, double* pkat4gz, double* pkat4pz);
	void UpdateConvectionBoundaryFaceProperties(double* phconv, double* ptgas);
	double SolveWallTemperature(double em, double rw, double hc, double qi, double tg, double tb);
};

#endif
