% Compiles the MIToolbox functions

mex -v MIToolboxMex.c MutualInformation.c Entropy.c CalculateProbability.c ArrayOperations.c
mex -v RenyiMIToolboxMex.c RenyiMutualInformation.c RenyiEntropy.c CalculateProbability.c ArrayOperations.c
mex -v WeightedMIToolboxMex.c WeightedMutualInformation.c WeightedEntropy.c CalculateProbability.c ArrayOperations.c
