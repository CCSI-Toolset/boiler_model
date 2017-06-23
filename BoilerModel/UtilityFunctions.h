// UtilityFunctions.h: interface for the CUtilityFunctions class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __UTILITYFUNCTIONS_H__
#define __UTILITYFUNCTIONS_H__

#include <cstdio>
#include <string>

const int UF_CUTOFF = 10;				//cut-off for minimum partition size used in quick sort

class CUtilityFunctions
{
public:
	CUtilityFunctions();
	virtual ~CUtilityFunctions();

	static void QuickSort(double* pdata, int* pindx, int ileft, int iright);
	static void QuickSortWithInsertion(double* pdata, int* pindx, int ileft, int iright);
	static void WriteString(std::string& str, FILE* pf);
	static void ReadString(std::string& str, FILE* pf);
};

#endif
