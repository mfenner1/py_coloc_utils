#!/usr/bin/env python

# system imports
import Image
from scipy.misc import fromimage
import numpy as np

# user defined imports
from DefinitionsAndUtils import *
from ImageProcessing     import thresholdNDArray

# from GraphAndHistogramUtilities import samplePoints

threshPath = "/home/mfenner/scipy_prep/data/input/"
threshFile = threshPath + "manual-thresholds-all-conditions.csv"
# threshFile = threshPath + "manual-thresholds-sampleof-conditions.csv"
threshDicts = readThresholdFileAsDictionaries(threshFile)
cndKeys = ["Organelle", "Stain", "Time"] # , "Series"]

allPixelCt = 0
nonZeroPixelCt = 0

#
# takes threshold file as "canonical" list of experiments and files
#
# for thisStack in breakIntoStacks(threshDicts):
#    for threshes in thisStack:
for thisCndSet in breakByConditions(threshDicts):
    cnd = tuple((k, thisCndSet[0][k]) for k in cndKeys)
    print cnd

    stackArray = np.empty((0,3), np.uint8)
    for threshes in thisCndSet:
        print threshes["Slice"],
        expertThresholds = dict(((c, threshes[c]) for c in colorNames))

        #
        # get file->image and thresholds together
        # apply threshold
        #
        exampleFilename = fileNameFormat[threshes["Organelle"]] % threshes
        currentImage = Image.open(imageDataPath+exampleFilename)

        asArray2 = fromimage(currentImage).reshape((numImagePoints,3))
        thresholdNDArray(asArray2, expertThresholds)

        stackArray = np.concatenate((stackArray, asArray2))

    # zero removal
    currAllPixelCt     = stackArray.shape[0]
    stackArray = stackArray[np.any(stackArray, 1)]

    currNonZeroPixelCt = stackArray.shape[0]
    print "(%d --> %d)" % (currAllPixelCt, currNonZeroPixelCt)
    
    allPixelCt     += currAllPixelCt
    nonZeroPixelCt += currNonZeroPixelCt


print "%15s: %12s" % ("with [0,0,0]",    allPixelCt)
print "%15s: %12s" % ("removed [0,0,0]", nonZeroPixelCt)
