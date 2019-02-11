################################################################################

                    `````````Top 'n' good peaks/hits`````````````

                                 Halfon Lab
                              Date: Sept 2018		
################################################################################


INTRODUCTION

This script takes into account both the raw SCRMshaw scores and the summed scores (amplitude) of the SCRMshaw-HD peaks. The cutoff points are based on determining the “elbow” points of the score and amplitude curves (calculated as the point furthest from the line connecting the first and last points on the curve. For SCRMshaw-HD, defining the top predictions is a two-part process. First, all peaks with amplitude above the amplitude cutoff point are accepted. Then, the SCRMshaw score curve is constructed as follows: first, each peak is evaluated to determine the maximum SCRMshaw score for any sequence window within the peak. These scores are then ranked, and the “elbow” point calculated. Peaks which also pass this cutoff are accepted as the set of top predictions.

1. INPUT
2. USAGE
3. PARAMETERS
4. OUTPUT

1.INPUT

The complete SCRMshaw output, (combined or individual training set's) containing 5000+ predicted CRMs for each training set, is used as input for this script. List of training set could also be provided.
 
2.USAGE

Following basic modules are required to run this script.Please make sure these modules have already been properly installed and are recognizable. pybedtools, statistics, scipy, numpy, subprocess, shutil and itertools. 
Following is an example of command line execution.
python cutoff.py -po SCRMshawHDConcatenatedOutput.bed -listTset listOfTrainingSets

3. PARAMETERS
 
 -po  The concatenated output of post processing SCRMshaw HD script
 -listTset List of Training sets expecting in output.
 
4.OUTPUT

List of Top CRMs will be generated (again based on amplitude and SCRMhsaw score) for each of the training set.
