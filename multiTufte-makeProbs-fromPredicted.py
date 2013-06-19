#!/usr/bin/env python

#
# system imports
#
import cPickle, Image, glob
from scipy.misc import fromimage
import numpy as np

#
# user defined imports
#
from DefinitionsAndUtils        import *
from ImageProcessing            import thresholdNDArray
from GraphAndHistogramUtilities import timeToIdx, toProbs, countQuantiles
from CurrentLM import ileNames, iles, applyCurrentLM

removeNonresponders   = True
removeSaturated       = True
useLogscaleCounts     = True
applyPredThresholds   = True


outputPath = "/home/mfenner/scipy_prep/data/output/pickle-probs/"
outputFile = outputPath + "predThresh-trimmedEnds-log-probs-by-OrgTime.pck"

#
# take allExpCnds and the file hierarchy as canonical
#
storedProbs = {}

#########################################
# forward to next big comment is necessary to get joined repetitions
#########################################
for cnd in allNonSerExpCnds:
    expArray = np.empty((0,3), np.uint8)
    for ser in series:
        cnd["Series"] = ser

        #
        # determine the number of slices in this 
        # stack by "ls'ing" the directory
        #
        cnd["Slice"] = "*"
        globpath = imageDataPath + fileNameFormat[cnd["Organelle"]] % cnd
        # ls -l blah/*/blah.tif (brittle FIXME)
        matchingFiles = glob.glob(globpath) 
        sliceValues = makeSliceValues(len(matchingFiles))

        for sli in sliceValues:
            cnd["Slice"] = sli
            exampleFilename = fileNameFormat[cnd["Organelle"]] % cnd
            # print exampleFilename,
            if exampleFilename.rsplit("/",1)[1] in \
                {"15m60xendser301.TIF",
                 "15m60xendser301.TIF",
                 "15m60xendser301.TIF",
                 "120m60xac17ser24.TIF",
                 "120m60xac17ser27.TIF"}:
                # print ".....skipping:"
                continue

            currentImage = Image.open(imageDataPath+exampleFilename)
            pixels = fromimage(currentImage).reshape((numImagePoints,3))        
            if applyPredThresholds:
                # zero removal for quantile computation
                pixels = pixels[np.any(pixels, 1)]        
                counts = {c:np.bincount(pixels[:,c], 
                                        minlength=256) for c in colors}
                qs     = {c:countQuantiles(counts[c], iles) for c in colors}

                # clunky:
                qDict = {"R8D": qs[R][0], "R9D": qs[R][1],
                         "G8D": qs[G][0], "G9D": qs[G][1],
                         "B8D": qs[B][0], "B9D": qs[B][1]}

                predThreshes = dict((c, applyCurrentLM(qDict, 
                                                       c)) for c in colorNames)
                # print ",".join("%5.4f" % t for t in predThresholds.values())
                thresholdNDArray(pixels, predThreshes, dropSaturated=True)

            expArray = np.concatenate((expArray, pixels))
#########################################
# all to here is necessary to get joined series
#########################################

    org = simplifyOrgStain(cnd["Organelle"], cnd["Stain"])
    t = timeToIdx(cnd["Time"])

    print org, t

    # convert image stack to counts and add to histograms
    for c1, c2 in colorPairs:
        probs = toProbs(expArray[:,c1], expArray[:,c2],
                        removeNonresponders = removeNonresponders,
                        removeSaturated     = removeSaturated,
                        useLogscaleCounts   = useLogscaleCounts)
        cnd = (org, c1, c2, t)
        storedProbs[cnd] = probs
    
outputFile = open(outputFile, "wb")
cPickle.dump(storedProbs, outputFile, -1)
outputFile.close()

