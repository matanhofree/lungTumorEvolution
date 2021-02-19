/*******************************************************************************
** CalculateChi2.h
** Part of the mutual information toolbox
**
** Contains functions to calculate the chi2 residuals for a pair of vectors
** 
**
*******************************************************************************/

#ifndef __CalculateChi2_H
#define __CalculateChi2_H

#ifdef __cplusplus
extern "C" {
#endif 

/*******************************************************************************
** calculate Chi2 residuals for a pair of binary vectors.
**
** length(dataVector) == length(targetVector) == vectorLength otherwise there
** will be a segmentation fault
*******************************************************************************/
void calculateChi2BinaryResid(uint16_T *dataA, uint16_T *dataB, int vectorLength, )
{

}

#ifdef __cplusplus
}
#endif

#endif

