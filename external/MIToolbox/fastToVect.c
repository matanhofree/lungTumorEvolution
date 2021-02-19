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

	int flag = *mxGetPr(prhs[0]);
	switch (flag)
	{
		case 1:
		{
			/**
			 ** ToVect dense matrix
			 */
			int i,j;
			int rows = mxGetM(prhs[1]); /* patients */
			int cols = mxGetN(prhs[1]); /* Genes */

			double* matData = mxGetPr(prhs[1]);

			size_t numTriU = (cols*(cols-1))/2;
            plhs[0] = mxCreateDoubleMatrix(numTriU,1,mxREAL);
            double* output = (double *)mxGetPr(plhs[0]);

 			size_t cp = 0;
 			for (j=1;j<cols;++j)
			{
 				for (i=0;i<j;++i)
				{
					output[cp++] = xMM(matData,i,j,rows);
				}
			}

			break;
		}
		case 2:
		{
			/**
			 ** To Vect sparse matrix
			 */
			mwSize i,j;
			int rows = mxGetM(prhs[1]); /* patients */
			int cols = mxGetN(prhs[1]); /* Genes */

			double  *sr;
			mwIndex *irs;
			mwIndex *jcs;
			mwIndex startI, endI;

		    sr  = mxGetPr(prhs[1]);
		    irs = mxGetIr(prhs[1]);
		    jcs = mxGetJc(prhs[1]);

		    // mwSize nnzN = mxGetNzmax(plhs[0]);
			size_t numTriU = (cols*(cols-1))/2;

            plhs[0] = mxCreateDoubleMatrix(numTriU,1,mxREAL);
            double* output  = (double *)mxGetPr(plhs[0]);

            //mexPrintf("Start\n");
            mwIndex pCol = 0;
            for (j=1;j<cols;++j)
			{
            	//mexPrintf("J%d\n",j);
            	startI = jcs[j];
            	endI = jcs[j+1];
            	pCol += (j-1);

            	if (startI == endI)
            	{
            		continue;
            	}

				for (i=startI;i<endI;++i)
				{
					if (irs[i] < j)
					{
						// xMM(output,irs[i],jcs[j],rows-1) = sr[i];
						output[irs[i]+pCol] = sr[i];
					}
					//mexPrintf("J%d,I%d\n",i,j);
				}

            	// pCol = pCol + (rows - j - 1); For a lower triangular

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

