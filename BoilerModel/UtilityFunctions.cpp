// UtilityFunctions.cpp: implementation of the CUtilityFunctions class.
//
//////////////////////////////////////////////////////////////////////
#include "UtilityFunctions.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CUtilityFunctions::CUtilityFunctions()
{

}

CUtilityFunctions::~CUtilityFunctions()
{

}

void CUtilityFunctions::QuickSort(double* pdata, int* pindx, int ileft, int iright)
{
	double dsave;
	double dmedian;
	int isave;
	int j, k;
	int ileft1 = ileft+1;
	int imid = (ileft+iright)/2;
	//swap between ileft+1 and imid
	isave = pindx[ileft1];
	dsave = pdata[ileft1];
	pindx[ileft1] = pindx[imid];
	pdata[ileft1] = pdata[imid];
	pindx[imid] = isave;
	pdata[imid] = dsave;
	if (pdata[ileft1]>pdata[iright])		//swap ileft1 and iright
	{
		isave = pindx[ileft1];
		dsave = pdata[ileft1];
		pindx[ileft1] = pindx[iright];
		pdata[ileft1] = pdata[iright];
		pindx[iright] = isave;
		pdata[iright] = dsave;
	}
	if (pdata[ileft]>pdata[iright])			//swap ileft and iright
	{
		isave = pindx[ileft];
		dsave = pdata[ileft];
		pindx[ileft] = pindx[iright];
		pdata[ileft] = pdata[iright];
		pindx[iright] = isave;
		pdata[iright] = dsave;
	}
	if (pdata[ileft1]>pdata[ileft])			//swap ileft1 and ileft
	{
		isave = pindx[ileft1];
		dsave = pdata[ileft1];
		pindx[ileft1] = pindx[ileft];
		pdata[ileft1] = pdata[ileft];
		pindx[ileft] = isave;
		pdata[ileft] = dsave;
	}
	j = ileft1;
	k = iright;
	dmedian = pdata[ileft];
	while (true)
	{
		do
		{
			j++;
		}while (pdata[j]<dmedian);			//keep looking if equal to avoid too many swap
		do
		{
			k--;
		}while (pdata[k]>dmedian);
		if (j<k)							//swap j and k
		{
			isave = pindx[j];
			dsave = pdata[j];
			pindx[j] = pindx[k];
			pdata[j] = pdata[k];
			pindx[k] = isave;
			pdata[k] = dsave;
		}
		else
			break;
	}
	//swap ileft and k
	isave = pindx[ileft];
	dsave = pdata[ileft];
	pindx[ileft] = pindx[k];
	pdata[ileft] = pdata[k];
	pindx[k] = isave;
	pdata[k] = dsave;
	if (k-ileft>UF_CUTOFF)
		QuickSort(pdata, pindx, ileft, k-1);
	if (iright-k>UF_CUTOFF)
		QuickSort(pdata, pindx, k+1, iright);	
}

void CUtilityFunctions::QuickSortWithInsertion(double* pdata, int* pindx, int ileft, int iright)
{
	int i, j;
	int isave;
	double dsave;
	//do quick sort with cutoff
	QuickSort(pdata, pindx, ileft, iright);
	//do insertion sort
	for (i=ileft+1; i<=iright; i++)
	{
		dsave = pdata[i];
		isave = pindx[i];
		for (j=i; j>0 && pdata[j-1]>dsave; j--)
		{
			pdata[j] = pdata[j-1];
			pindx[j] = pindx[j-1];
		}
		pdata[j] = dsave;
		pindx[j] = isave;
	}
}

void CUtilityFunctions::WriteString(std::string& str, FILE* pf)
{
	int strlen = str.size();
	fwrite(&strlen,sizeof(int),1,pf);
	fwrite(str.c_str(),sizeof(char),strlen,pf);
}

void CUtilityFunctions::ReadString(std::string& str, FILE* pf)
{
	int strlen;
	char* pbuffer;
	fread(&strlen,sizeof(int),1,pf);
	pbuffer = new char [strlen+1];
	pbuffer[strlen] = '\0';
	fread(pbuffer,sizeof(char),strlen,pf);
	str = pbuffer;
	delete [] pbuffer;
}