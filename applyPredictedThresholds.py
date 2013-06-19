#!/usr/bin/env python

# system imports
import Image, glob
from scipy.misc import fromimage
import numpy as np

# user defined imports
from DefinitionsAndUtils import *
from ImageProcessing     import thresholdNDArray
from GraphAndHistogramUtilities import countQuantiles
from CurrentLM import ileNames, iles, applyCurrentLM

REMOVE_ZEROS = True

for cnd in allExpCnds:
    # determine the number (really, just need the number of digits)
    # of slices in this stack by "ls'ing" the directory
    cnd["Slice"] = "*"
    globpath = imageDataPath + fileNameFormat[cnd["Organelle"]] % cnd
    matchingFiles = glob.glob(globpath) # ls -l blah/*/blah.tif (brittle FIXME)
    sliceValues = makeSliceValues(len(matchingFiles))
    
    stackArray = np.empty((0,3), np.uint8)
    for sli in sliceValues:
        cnd["Slice"] = sli
        exampleFilename = fileNameFormat[cnd["Organelle"]] % cnd
        print exampleFilename,
        if exampleFilename.rsplit("/",1)[1] in \
            {"15m60xendser301.TIF",
             "15m60xendser301.TIF",
             "15m60xendser301.TIF",
             "120m60xac17ser24.TIF",
             "120m60xac17ser27.TIF"}:
            print ".....skipping:"
            continue
        
        currentImage = Image.open(imageDataPath+exampleFilename)
        pixels = fromimage(currentImage).reshape((numImagePoints,3))        
        # zero removal
        if REMOVE_ZEROS:
            pixels = pixels[np.any(pixels, 1)]        
        counts = {c:np.bincount(pixels[:,c], minlength=256) for c in colors}
        qs     = {c : 
                  countQuantiles(counts[c], iles) for c in colors}
        # import pdb; pdb.set_trace()
        
        # clunky:
        qDict = {"R8D": qs[R][0], "R9D": qs[R][1],
                 "G8D": qs[G][0], "G9D": qs[G][1],
                 "B8D": qs[B][0], "B9D": qs[B][1]}

        predictedThresholds = \
                 dict((c, applyCurrentLM(qDict, c)) for c in colorNames)
        print ",".join("%5.4f" % t for t in predictedThresholds.values())
        thresholdNDArray(pixels, predictedThresholds)

        stackArray = np.concatenate((stackArray, pixels))
        
