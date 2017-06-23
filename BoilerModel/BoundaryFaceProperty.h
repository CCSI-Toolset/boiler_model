// BoundaryFaceProperty.h: interface for the CBoundaryFaceProperty class.
//
//////////////////////////////////////////////////////////////////////
#ifndef __BOUNDARYFACEPROPERTY_H__
#define __BOUNDARYFACEPROPERTY_H__

#include <cstdio>

class CBoundaryFaceProperty
{
public:
	int		iwall;			//wall index, 0=water wall, >0=SH wall
	int		iwater;			//water stream index
	double	hconv;			//convective heat tranfer coefficient, enabled on 11/6/2015
	double	emis;			//boundary emisivity
	double	rwall;			//wall heat transfer resistance
	double	tback;			//backside fluid temperature
	double	temp;			//inner wall temperature, used for radiation
	double	tgas;			//temperature of gas in contact with the boundary, the same as zone temperature, added on 11/6/2015
	double	fradadj;		//adjusting factor for radiation hemisphere flux calculation
	double	qinc;			//incident radiation heat flux
	double	qnet;			//net radiation heat flux
	double	qconv;			//convective heat flux
	
	CBoundaryFaceProperty();
	virtual ~CBoundaryFaceProperty();
	void Write(FILE* pf);
	void Read(FILE* pf);
};

#endif