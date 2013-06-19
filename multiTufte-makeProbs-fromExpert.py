#!/usr/bin/env python

#
# system imports
#
import cPickle, Image
from scipy.misc import fromimage
import numpy as np
import sys

#
# user defined imports
#
from DefinitionsAndUtils        import *
from ImageProcessing            import thresholdNDArray
from GraphAndHistogramUtilities import timeToIdx, toProbs

removeNonresponders   = True
removeSaturated       = True
useLogscaleCounts     = True
applyExpertThresholds = True

threshPath = "/home/mfenner/scipy_prep/data/input/"
threshFile = threshPath + "manual-thresholds-all-conditions.csv"

outputPath = "/home/mfenner/scipy_prep/data/output/pickle-probs/"
outputFile = outputPath + "expThresh-trimmedEnds-log-probs-by-OrgTime.pck"

threshDicts = readThresholdFileAsDictionaries(threshFile)
cndKeys = ["Organelle", "Stain", "Time"] # , "Series"]

#
# takes threshold file as "canonical" list of experiments and
# files
#
storedProbs = {}
#for thisStack in breakIntoStacks(threshDicts):
#    for threshes in thisStack:
for thisCndSet in breakByConditions(threshDicts):
    stackArray = np.empty((0,3), dtype=np.uint8)
    print tuple(thisCndSet[0][k] for k in cndKeys), "...", 
    for threshes in thisCndSet:
        # print threshes["Slice"],

        # get file->image and thresholds together and apply them

        exampleFilename = fileNameFormat[threshes["Organelle"]] % threshes
        currentImage = Image.open(imageDataPath+exampleFilename)
        asArray = fromimage(currentImage).reshape((numImagePoints, 3))

        # port to thresholdNDArray
        if applyExpertThresholds:
            expertThresholds = dict(((c, threshes[c]) for c in colorNames))
            thresholdNDArray(asArray, expertThresholds)

        # asArray = fromimage(currentImage).reshape((numImagePoints, 3))
        stackArray = np.concatenate((stackArray, asArray))

    org = simplifyOrgStain(threshes["Organelle"], threshes["Stain"])
    # r = threshes["Series"] - 1
    t = timeToIdx(threshes["Time"])
    
    print org, t

    # convert image stack to counts and add to histograms
    for c1, c2 in colorPairs:
        probs = toProbs(stackArray[:,c1], stackArray[:,c2],
                        removeNonresponders = removeNonresponders,
                        removeSaturated     = removeSaturated,
                        useLogscaleCounts   = useLogscaleCounts)
        cnd = (org, c1, c2, t)
        storedProbs[cnd] = probs
    
outputFile = open(outputFile, "wb")
cPickle.dump(storedProbs, outputFile, -1)
outputFile.close()

