/*******************************************************************************
 **
 **  MIToolboxMex.cpp
 **  is the MATLAB entry point for the MIToolbox functions when called from
 **  a MATLAB/OCTAVE script.
 **
 **  Copyright 2010 Adam Pocock, The University Of Manchester
 **  www.cs.manchester.ac.uk
 **
 **  This file is part of MIToolbox.
 **
 **  MIToolbox is free software: you can redistribute it and/or modify
 **  it under the terms of the GNU Lesser General Public License as published by
 **  the Free Software Foundation, either version 3 of the License, or
 **  (at your option) any later version.
 **
 **  MIToolbox is distributed in the hope that it will be useful,
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **  GNU Lesser General Public License for more details.
 **
 **  You should have received a copy of the GNU Lesser General Public License
 **  along with MIToolbox.  If not, see <http://www.gnu.org/licenses/>.
 **
 *******************************************************************************/
#include "MIToolbox.h"
#include "ArrayOperations.h"
#include "CalculateProbability.h"
#include "Entropy.h"
#include "MutualInformation.h"

#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS REAL 2D FULL DOUBLE(P) && mxGetNumberOfElements(P) == 1)
#define MM(m,n,D) (m + (D*n))
#define xMM(A,m,n,D) (A[m+(D*n)])


double tzsum(double* vec1,double* vec2,int vLen)
{
	int i = 0;
	double zout = 0;
	printf("In tzsum\n");
	for (i=0;i<vLen;++i)
	{
		printf("%d) %f\n",i,vec2[i]);
		zout = zout + vec1[i];
	}

	return zout;
}


/*******************************************************************************
 **entry point for the mex call
 **nlhs - number of outputs
 **plhs - pointer to array of outputs
 **nrhs - number of inputs
 **prhs - pointer to array of inputs
 *******************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*****************************************************************************
	 ** this function takes a flag and a variable number of arguments
	 ** depending on the value of the flag and returns either a construct 
	 ** containing probability estimates, a merged vector or a double value 
	 ** representing an entropy or mutual information
	 *****************************************************************************/

	/*  
    int flag, i,j, numberOfSamples, checkSamples, thirdCheckSamples, numberOfFeatures, checkFeatures, thirdCheckFeatures;
	  int  numArities, errorTest;
	  double *dataVector, *condVector, *targetVector, *firstVector, *secondVector, *output, *numStates;
	  double *matrix, *mergedVector, *arities;
	  int *outputIntVector, *intArities;
	  double *jointOutput, *numJointStates, *firstOutput, *numFirstStates, *secondOutput, *numSecondStates;
	
	  ProbabilityState state;
	  JointProbabilityState jointState;
	
	  if (nlhs != 1)
	  {
	    printf("Incorrect number of output arguments\n");
	  }//if not 1 output
	  */


	switch (nrhs)
	{
		case 2:
		{
			/*printf("Must be H(X), calculateProbability(X), merge(X), normaliseArray(X)\n");*/
			break;
		}
		case 3:
		{
			/*printf("Must be H(XY), H(X|Y), calculateJointProbability(XY), I(X;Y)\n");*/
			break;
		}
		case 4:
		{
			/*printf("Must be I(X;Y|Z)\n");*/
			break;
		}
		default:
		{
			printf("Incorrect number of arguments, format is MIToolbox(\"FLAG\",varargin)\n");
			break;
		}
	}

	  /* number to function map
	  ** 1 = I(X;Y)
	  ** 2 = I(X;Y|Z)
	  */
	  

	int i,j;
	int flag = *mxGetPr(prhs[0]);
	switch (flag)
	{
		case 1:
		{
			/**
			 **I(X;Y)
			 */
			int rows = mxGetM(prhs[1]); /* patients */
			int cols = mxGetN(prhs[1]); /* Genes */

			double* matData = mxGetPr(prhs[1]);


            plhs[0] = mxCreateDoubleMatrix(cols,cols,mxREAL);
 			double* output = (double *)mxGetPr(plhs[0]);

 			/*
 			i = 2; j=3;
 			double zz = tzsum(&matData[i*rows],&matData[j*rows],rows);
            double zz = calculateMutualInformation(&matData[i*rows],&matData[j*rows],rows);
            printf("%f\n\n",zz);

            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            double* output = (double *)mxGetPr(plhs[0]);

            i = 0;
            j = 3;
            printf("%d %d %u %u\n",rows,cols,&matData[i*rows],&matData[j*rows]);
            printf("i = %f\tj = %f\n",matData[i*rows],matData[1+j*rows]);



            printf("%f\n\n",zz);


			for (j=0;j<cols;++j)
			{
            	for (i=0;i<rows;++i)
				{
                    printf("\t%f",matData[(i + j*rows)]);

				}
                printf("\n");
			}

			i = 0;j=1;
			xMM(output,i,j,rows) =  zz;
*/

            for (i=0;i<cols;++i)
			{
				for (j=(i+1);j<cols;++j)
				{
					
					xMM(output,i,j,cols) = calculateMutualInformation(&matData[i*rows],&matData[j*rows],rows);
				}
			}

            /*double calculateMutualInformation(double *firstVector, double *secondVector, int vectorLength);*/
           
			break;
		}/*case 7 - I(X;Y)*/
		default:
		{
			printf("Unrecognised flag\n");
			break;
		}/*default*/
	}/*switch(flag)*/

	return;
}/*mexFunction()*/

