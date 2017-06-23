#include <cstdio>
#include <stdlib.h>
#include <cstring>
#include "SolveBoilerModel.h"

const int NIINPUT = 50;
const int NDINPUT = 400;
const int NDOUTPUT = 500;

int main(int argc, char* argv[])
{
	//main function to test SolveBoilerModel function
	//command line options:
	//-i for input file, if not specified, use default file name "boiler_model_input.txt"
	//-o for output file name, if not specified, use default file name "boiler_model_output.txt"
	char input_filename[500] = "boiler_model_input.txt";
	char output_filename[500] = "boiler_model_output.txt";
	char cline[501];			//a line of characters
	char strdouble[20];			//a string contain double
	int i;						//used in loop
	int nfuels;					//number of solid fuels in primary stream
	int niinput;				//number of integer inputs
	int ndinput;				//number of double inputs
	int noutput;				//number of output variables
	int nitem;					//number of items read succesfully, used to handle empty line and invalid double number
	int iinput[NIINPUT];		//integer input array
	double dinput[NDINPUT];		//double input array
	double doutput[NDOUTPUT];	//double output array
	//check argument array to find the input and output file names
	for (i=1; i<argc; i++)
	{
		if (!strcmp(argv[i],"-i"))
		{
			if (i+1<argc)
				strcpy(input_filename,argv[i+1]);
			i++;
		}
		else
		{
			if (!strcmp(argv[i],"-o"))
			{
				if (i+1<argc)
					strcpy(output_filename,argv[i+1]);
				i++;
			}
			else
			{
				if (argv[i][0]=='-')
				{
					printf("Invalid commandline option!\n");
					return 1;
				}
			}
		}
	}
	//read the input file that contains integer and double input arrays
	FILE* fin = fopen(input_filename,"r");
	if (!fin)
	{
		printf("Cannot open boiler model input file \"%s\"!\n",input_filename);
		return 1;
	}
	FILE* fout = fopen(output_filename,"w");
	if (!fout)
	{
		printf("Cannot open boiler model output file \"%s\"!\n",output_filename);
		return 1;
	}
	//read 11 integer input variables
	for (i=0; i<11; i++)
	{
		fscanf(fin,"%d",&iinput[i]);
		fgets(cline,500,fin);
	}
	//calculate total integers needed
	niinput = 13 + iinput[6] + iinput[7] + iinput[8];
	for (i=11; i<niinput; i++)
	{
		fscanf(fin,"%d",&iinput[i]);
		fgets(cline,500,fin);
	}
	if (iinput[niinput-1]>0)	//number of superheater panels > 0, read number of points for the polygon
	{
		fscanf(fin,"%d",&iinput[niinput++]);
		fgets(cline,500,fin);
	}
	//read number of solid fuels in primary stream
	fscanf(fin,"%d",&nfuels);
	fgets(cline,500,fin);
	iinput[niinput++] = nfuels;
	for (i=0; i<nfuels; i++)
	{
		fscanf(fin,"%d",&iinput[niinput++]);
		fgets(cline,500,fin);
	}
	//flag for zone based boundary condition
	fscanf(fin,"%d",&iinput[niinput++]);
	fgets(cline,500,fin);
	//flag for including convective heat transfer
	fscanf(fin,"%d",&iinput[niinput++]);
	fgets(cline,500,fin);
	//check if the last integer input is valid
	if (iinput[niinput-1]!=0 && iinput[niinput-1]!=1)
	{
		printf("Invalid last integer input!\n");
		fclose(fin);
		return 1;
	}
	//read double input varaibles
	ndinput = 0;
	for (i=0; i<NDINPUT; i++)
	{
		nitem = fscanf(fin,"%s",strdouble);
		if (nitem!=1)	//empty line
			break;
		nitem = sscanf(strdouble,"%lg",&dinput[i]);
		if (nitem!=1)	//invalid double
		{
			printf("Invalid double in input file!\n");
			fclose(fin);
			return 1;
		}
		fgets(cline,500,fin);
		ndinput++;
		if (ferror(fin) || feof(fin))
			break;
	}
	fclose(fin);
	noutput = SolveBoilerModel(niinput, ndinput, iinput, dinput, doutput);
	if (noutput==-1)
	{
		printf("Failed to solve boiler model!\n");
		return 1;
	}
	//write results to output file
	for (i=0; i<noutput; i++)
		fprintf(fout,"%lg\n",doutput[i]);
	fclose(fout);
	printf("The end of main() function.\n");
	return 0;
}

