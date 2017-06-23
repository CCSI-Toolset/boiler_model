// RadiantFurnace.h: interface for the CRadiantFurnace class.
// This is the main class for the 1D/3D hybrid model with flow variables in 1D and radiation variables in 3D.
// The 3-D mesh is built on Cartesian coordinates where x, y and z are in depth, height, and width directions, respectively.
// x = 0 at front wall, y = 0 at bottom of hopper, z=0 at left hand side wall
// A typical boiler has a hopper and nose and could exit at the top or the back. Arch-fired configuration is not included in the model.
//
//////////////////////////////////////////////////////////////////////
#ifndef __RADIANTFURNACE_H__
#define __RADIANTFURNACE_H__

#include "DiscreteOrdinateEquation.h"
#include "SuperHeater.h"
#include "PSD.h"
#include "CoalKinetics.h"
#include "MaterialStream.h"
#include "Mesh.h"

const int RF_NZONE = 20;					//maximum number of total zones in y direction
const int RF_NX = 50;						//maximum number of cells in x direction
const int RF_NY = 100;						//maximum number of cells in y direction
const int RF_NZ = 50;						//maximum number of cells in z direction
const int RF_NS_COMB = 20;					//maximum number of inlet streams for combustion
const int RF_NS_WATR = 5;					//maximum number of inlet streams on water/steam side
const int RF_NSH = 3;						//maximum number of superheaters
const int RF_NPSD = 3;						//maximum number of particle size distribution
const int RF_NSPECIES = 12;					//maximum number of species in gas phase
const int RF_NCK = 3;						//maximum number of coal kinetics
const int RF_NSTEP = 100;					//number of Lagrange iteration time steps in each zone
const int RF_NRXN = 3;						//number of char heterogeneous reactions, currently by O2, H2O and CO2 only

class CRadiantFurnace
{
public:
	//variables hard-wired in constructor
	int		nspecies;						//number of species selected for equilibrium calculation
	int		ispecies[RF_NSPECIES];			//species index array, currently 11 species are considered
	double	urf_temp;						//under-relaxation factor for solving temperature in each zone
	double	urf_kap;						//under-relaxation factor for particle radiation absorption coefficient
	double	urf_hl;							//under-relaxation factor for heat loss
	double	urf_daf;						//under-relaxation factor for char daf mass
	double	ef_kap;							//effectiveness factor for particle radiation absorption coefficient due to uniform zone property assumption
	double	ef_kag;							//effectiveness factor for gas radiation absorption coefficient due to uniform zone property assumption
	double	ef_rxn[RF_NZONE];				//effectiveness factor for char/gas heterogenous reactions due to uniform zone bulk concentration assumption
	
	//user input variables
	bool	bmodified;						//true if user input has changed
	int		iconvect;						//option for convective heat transfer, 0=off, 1=on
	int		ifurn_type;						//furnace type, 0=rear wall exit, 1=top exit
	int		iprop_zone;						//method to specify enclosure wall properties, 0 means uniform, 1 means zone dependent
	int		nzone_burn;						//number of zones in burner region between hopper knuckle and nose bottom
	int		nzone_nose;						//number of zones in nose region between nose bottom and tip of nose
	int		nzone_uppr;						//number of zones in upper region between tip of nose and roof
	int		nx;								//number of cells in x (depth direction)
	int		nz;								//number of cells in z (width direction)
	int		ns_comb;						//number of inlet streams for combustion
	int		ns_watr;						//number of inlet streams on water/steam side
	int		nsh;							//number of superheaters
	int		npsd;							//number of particle size distributions
	int		nck;							//number of coal kinetics
	int		pncell_zone[RF_NZONE];			//number of cells in y direction of each zone
	double	pres;							//pressure of outlet streams (minimum of the pressures of inlet streams)
	double	xdepth;							//depth in x direction
	double	xhop_bot0;						//left hopper bottom point in x direction
	double	xhop_bot1;						//right hopper bottom point in x direction
	double	xnose_frt;						//front wall x at nose tip, added for UK test furnace
	double	xnose_tip;						//nose tip point in x direction
	double	yheight;						//height in y direction
	double	yhop_knuckle;					//hopper knuckle point in y direction
	double	ynose_bot;						//nose bottom point in y direction
	double	ynose_tip;						//nose tip point in y direction
	double	zwidth;							//width in z direction
	double	pemis_sh[RF_NSH+1];				//superheater and outlet wall emissivity, last=outlet
	double	prwall_sh[RF_NSH+1];			//superheater and outlet wall heat transfer resistance, last=outlet
	double	ptback_sh[RF_NSH+1];			//superheater and outlet wall backside temperature, last=outlet
	double	pemis_zone[RF_NZONE];			//enclosure wall emissivity, use pemis_zone[0] if uniform
	double	prwall_zone[RF_NZONE];			//enclosure wall heat transfer resistance, use prwall_zone[0] if uniform
	double	ptback_zone[RF_NZONE];			//enclosure wall backside temperature, use ptback_zone[0] if uniform
	double	pyzone[RF_NZONE];				//upper end y for each zone, used to determine inlet zone index
	double	pyinlet[RF_NS_COMB];			//y of inlet combustion streams, if less than yhop_knuckle, put in the first burner zone
	CSuperHeater psh[RF_NSH];				//array of shapes (polygon) for SH
	CPSD ppsd[RF_NPSD];						//array of particle size distribution
	CCoalKinetics pck[RF_NCK];				//array of coal kinetics
	CMaterialStream* pstm_comb[RF_NS_COMB];	//array of inlet combustion stream pointers
	
	//calculated variables
	bool	pbzone_dead[RF_NZONE];			//dead zone flag
	bool	ppbzone[RF_NZONE][RF_NZONE];	//flag for upstream and current zones
	int		nzone;							//total number of zones
	int		ny;								//number of cells in y (height direction)
	int		izone_1st;						//zone index for the first zone (with inlet flow)
	int		pizone_in[RF_NS_COMB];			//index of zone to which each inlet stream enters
	double	phflow_zone_ad[RF_NZONE];		//adiabatic enthalpy flow of material leaving each zone, calculated once only
	double	ppeflow_zone[RF_NZONE][GS_NE];	//element molar flow leaving each zone, calculated once only
	double	pashflow_zone[RF_NZONE];		//mass flow of ash leaving each zone, calculated once only
	double	ppdafflow_zone[RF_NZONE][GS_NE];//C H O N S Cl daf element mass flow of particle phase leaving each zone, calculated by kinetics
	double	phhvflow_zone[RF_NZONE];		//daf HHV flow of particle phase leaving each zone
	double	pgasemis[RF_NZONE];				//gas emissivity for each zone
	double	pkag_zone[RF_NZONE];			//gas absorption coefficient for each zone
	double	pkat4g_zone[RF_NZONE];			//Kag*Tg^4
	double	pkap_zone[RF_NZONE];			//particle absorption coefficient for each zone
	double	pkat4p_zone[RF_NZONE];			//Kap*Tp^4
	double	phl_zone[RF_NZONE];				//radiation heat loss in each zone
	double	phlrad_zone[RF_NZONE];			//radiation heat loss in each zone
	double	phlconv_zone[RF_NZONE];			//convection heat loss in each zone
	double	ptemp_zone[RF_NZONE];			//temperature in each zone, used for under relaxation
	double	ptime_zone[RF_NZONE];			//residence time of each zone, calculated by gas flow rate and volume
	double	pvol_zone[RF_NZONE];			//volume of each zone, calculated by summing up volumes of individual cells
	double	pareaw_zone[RF_NZONE];			//wall area of each zone, calculated by summing up area of boundary wall faces
	double	pdiamh_zone[RF_NZONE];			//hydraulic diameter of each zone, formula=4*V/A
	double  ptwall_zone[RF_NZONE];			//average wall temperature of each zone, used to calculate film temperature for convective heat transfer
	double	phconv_zone[RF_NZONE];			//average convective heat transfer coefficient, 0 for dead zone, 0 if iconv==0
	double	pmrxn_char[RF_NRXN];			//mass of char reacted by individual reactions (by O2, H2O, and CO2)
	double	qexit;							//heat loss to exit
	double	pqwater[RF_NS_WATR];			//heat loss to water streams
	double	pqwall[RF_NSH+1];				//heat loss to enclosure and SH walls
	double	pqwall_zone[RF_NZONE];			//heat loss to wall/supper heater in each zone
	
	CMesh mesh;								//3-D mesh
	CDiscreteOrdinateEquation doeqn;		//radiation equation
	CMaterialStream pstm_zone[RF_NZONE];	//array of product streams for each zone
	CMaterialStream* pstm_flue;				//pointer to a flue gas outlet stream

public:
	//functions
	CRadiantFurnace();
	virtual ~CRadiantFurnace();
	CRadiantFurnace(const CRadiantFurnace &t);
	CRadiantFurnace& operator=(const CRadiantFurnace& t);
	void Write(FILE* pf);
	void Read(FILE* pf);
	int CheckInputs();
	int Init();
	void CalcZoneElementAshAndEnthalpyFlows();
	void CalcFullyReactedAdiabaticZoneTemperatures();
	void CalcConvectiveHeatTransferCoefficient();
	void UpdateZoneResidenceTime();
	void TrackSolidParticlesAndUpdateZoneStreams(double* pkap);
	void InitGasStreamForEquilibriumCalculation(CIdealGasMixture* pigm);
	int CalcMaterialStreamEquilibrium(CMaterialStream* pstm);
	double CalcGasEmissivity(CMaterialStream* pstm);
	double CalcGasAbsorptionCoefficient(double ge);
	int Solve();
	int CalculateOutStreams();
	void GetVertexCount(int* pcount);
	void GetVertexValues(int* pcount, double* pvar_zone, float* pvar_vertex);
	void GetCellValues(double* pvar_zone, float* pvar_cell);
	void WriteVtkFiles(char* filename);
	void WriteResultDescriptionFile(char* filename);
	double CalcWaterWallArea();
	double CalcSuperHeaterWallArea();
	void PrepareResults(std::vector<double>& results);
	void PrepareResultDescriptionStrings(std::vector<std::string>& results);
};

#endif