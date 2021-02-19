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
#include "mex.h"
#include "matrix.h"
#include "math.h"

#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS REAL 2D FULL DOUBLE(P) && mxGetNumberOfElements(P) == 1)
#define MM(m,n,D) (m + (D*n))

#define xMM(A,m,n,D) (A[m+(D*n)])

/*******************************************************************************
** calculate Chi2 residuals for a pair of binary vectors.
**
** length(dataVector) == length(targetVector) == vectorLength otherwise there
** will be a segmentation fault
** dataA / dataB not binary this will error
** outResid needds to be a 4 element double array
** 		 B
** __|_0_|_1_|___
** 0 | 0 | 2 | 0
** --|---|---|---   A
** 1 | 1 | 3 | 1
** --|---|---|---
**   | 2 | 3 |
*******************************************************************************/
int calculateChi2BinaryResid(uint16_T* dataA, uint16_T* dataB, int vectorLength, double* outResid)
{
	uint32_T countTable[] = {0,0,0,0};
	int pos = 0;
	for (int i = 0; i < vectorLength; ++i)
	{
		pos = dataA[i]+((dataB[i])*2);
		++countTable[pos];
	}

	double countTableD[] = {0.0, 0.0, 0.0, 0.0};
	countTableD[0] = (double)countTable[0];
	countTableD[1] = (double)countTable[1];
	countTableD[2] = (double)countTable[2];
	countTableD[3] = (double)countTable[3];

	//mexPrintf("%f\t%f\n%f\t%f\n",countTableD[0],countTableD[2],countTableD[1],countTableD[3]);

	double N = (double)vectorLength;
	double marginalTable[] = {0.0, 0.0, 0.0, 0.0};
	marginalTable[0] = (countTableD[0]+countTableD[2])/N;
	marginalTable[1] = (countTableD[1]+countTableD[3])/N;
	marginalTable[2] = (countTableD[0]+countTableD[1])/N;
	marginalTable[3] = (countTableD[2]+countTableD[3])/N;

	//mexPrintf("%f\t%f\n%f\t%f\n",marginalTable[0],marginalTable[2],marginalTable[1],marginalTable[3]);

	double expectedTable[] = {0.0, 0.0, 0.0, 0.0};
	expectedTable[0] = marginalTable[0]*marginalTable[2]*N;
	expectedTable[1] = marginalTable[1]*marginalTable[2]*N;
	expectedTable[2] = marginalTable[0]*marginalTable[3]*N;
	expectedTable[3] = marginalTable[1]*marginalTable[3]*N;

	outResid[0] = (countTableD[0] - expectedTable[0])/sqrt(expectedTable[0]);
	outResid[1] = (countTableD[1] - expectedTable[1])/sqrt(expectedTable[1]);
	outResid[2] = (countTableD[2] - expectedTable[2])/sqrt(expectedTable[2]);
	outResid[3] = (countTableD[3] - expectedTable[3])/sqrt(expectedTable[3]);

	return 0;
}

int findMax(uint16_T* dataIn, size_t numEl)
{
	int maxVal = 0;
	for (size_t i = 0; i<numEl; ++i) { if (maxVal < dataIn[i]) {maxVal = dataIn[i];} }

	return maxVal;
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
		case 5:
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
			 ** Chi2 of two binary vectors
			 */
			int rows = mxGetM(prhs[1]); /* patients */
			int cols = mxGetN(prhs[1]); /* Genes */

			uint16_T* matData = mxGetData(prhs[1]);

//			if (findMax(matData,mxGetNumberOfElements(prhs[1])) > 1)
//			{
//				mexPrintf("Error: this verseion of chi2 only works on binary vectors");
//				break;
//			}

//			if (nlhs < 3)
//			{
//				mexPrintf("Error: this verseion of chi2 returns 4 arrays");
//			    break;
//			}

            plhs[0] = mxCreateDoubleMatrix(cols,cols,mxREAL);
 			double* output11 = (double *)mxGetPr(plhs[0]);

// 			plhs[1] = mxCreateDoubleMatrix(cols,cols,mxREAL);
// 			double* output21 = (double *)mxGetPr(plhs[1]);
//
//            plhs[2] = mxCreateDoubleMatrix(cols,cols,mxREAL);
// 			double* output12 = (double *)mxGetPr(plhs[2]);
//
//            plhs[3] = mxCreateDoubleMatrix(cols,cols,mxREAL);
// 			double* output22 = (double *)mxGetPr(plhs[3]);

 			double nan = mxGetNaN();
 			double chi2R[] = {0.0, 0.0, 0.0, 0.0};

            for (i=0;i<cols;++i)
			{
				for (j=(i+1);j<cols;++j)
				{
					// mexPrintf("Debug: %d %d\n",i,j);
					if ( calculateChi2BinaryResid(&matData[i*rows],&matData[j*rows],rows,chi2R) == 0)
					{
						xMM(output11,i,j,cols) = pow(chi2R[0],2.0) +pow(chi2R[1],2.0) +pow(chi2R[2],2.0) + pow(chi2R[3],2.0);
//						xMM(output21,i,j,cols) = chi2R[1];
//					    xMM(output12,i,j,cols) = chi2R[2];
//					    xMM(output22,i,j,cols) = chi2R[3];
					}
					else
					{
						mexPrintf("Error %d %d\n",i,j);
						xMM(output11,i,j,cols) = nan;
//						xMM(output21,i,j,cols) = nan;
//					    xMM(output12,i,j,cols) = nan;
//					    xMM(output22,i,j,cols) = nan;
					}
				}
			}

			break;
		}/*case 1 - I(X;Y)*/
//		case 2:
//		{
//			/**
//			 **I(X;Y) - For a subset
//			 */
//			int rows = mxGetM(prhs[1]); /* patients */
//			int cols = mxGetN(prhs[1]); /* Genes */
//			double* matData = mxGetPr(prhs[1]);
//
//			int rowSelectLen = mxGetN(prhs[2]);
//			int32_T* rowSelect = (uint32_T*)mxGetData(prhs[2]);
//
//			int colSelectLen = mxGetN(prhs[3]);
//			int32_T* colSelect = (uint32_T*)mxGetData(prhs[3]);
//			/* printf("Input: %d %d %d %d\n",rows,cols,rowSelectLen,colSelectLen); */
//
//            plhs[0] = mxCreateDoubleMatrix(rowSelectLen,colSelectLen,mxREAL);
// 			double* output = (double *)mxGetPr(plhs[0]);
// 			/* printf("Created: %d %d %d %d\n",rows,cols,rowSelectLen,colSelectLen); */
//
// 			double nan = mxGetNaN();
//
// 			int32_T i,j;
// 			int32_T ix,jx;
// 			for (ix=0;ix<rowSelectLen;++ix)
//			{
// 				i = rowSelect[ix]-1;
//				for (jx=0;jx<colSelectLen;++jx)
//				{
//					j = colSelect[jx]-1;
//					/* printf("Calc: (%d %d)<- (%d %d)\n",ix,jx,i,j); */
//					if (j>i)
//					{
//						xMM(output,ix,jx,rowSelectLen) = calculateMutualInformation(&matData[i*rows],&matData[j*rows],rows);
//					}
//					else
//					{
//						xMM(output,ix,jx,rowSelectLen) = nan;
//					}
//				}
//			}
// 			/* printf("Done: %d %d %d %d\n",rows,cols,rowSelectLen,colSelectLen); */
//
//			break;
//		}/*case 2 - I(X;Y) for a subset of indices*/
//		case 3:
//		{
//			/*
//			**I(X;Y|Z)
//			*/
//			int rows = mxGetM(prhs[1]); /* patients */
//			int cols = mxGetN(prhs[1]); /* Genes */
//			double* matData = mxGetPr(prhs[1]);
//
//			int condRows = mxGetM(prhs[2]); /* patients */
//            int condCols = mxGetN(prhs[2]); /* Genes */
//    		double* condData = mxGetPr(prhs[2]);
//
//    		if (condRows != rows || condCols != cols)
//    		{
//    			mexPrintf("Error data matrix not identical in dim to conditional matrix");
//    			break;
//    		}
//
//            plhs[0] = mxCreateDoubleMatrix(cols,cols,mxREAL);
// 			double* output = (double *)mxGetPr(plhs[0]);
// 			double* mergedVector = (double *) mxCalloc(condRows,sizeof(double));
// 			double* stateVectA;
// 			double* stateVectB;
//
//            for (i=0;i<cols;++i)
//			{
//				for (j=(i+1);j<cols;++j)
//				{
//					/*int mergeMultipleArrays(double *inputMatrix, double *outputVector, int matrixWidth, int vectorLength)*/
//					stateVectA = &condData[i*rows];
//					stateVectB = &condData[j*rows];
//					mergeArrays(stateVectA,stateVectB, mergedVector,condRows);
//					xMM(output,i,j,cols) = calculateConditionalMutualInformation(&matData[i*rows],&matData[j*rows],mergedVector,rows);
//				}
//			}
//            mxFree(mergedVector);
//			break;
//
//		}/*case 8 - I(X;Y|Z)*/
		default:
		{
			printf("Unrecognised flag\n");
			break;
		}/*default*/
	}/*switch(flag)*/

	return;
}/*mexFunction()*/

