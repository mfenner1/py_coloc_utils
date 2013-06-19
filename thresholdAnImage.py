import os.path as osp

import Image
from scipy.misc import fromimage
import numpy as np

from ImageProcessing import thresholdNDArray
from DefinitionsAndUtils import *
from GraphAndHistogramUtilities import countQuantiles
from CurrentLM import applyCurrentLM, iles, ileNames

def applyPredThresh(pixels):
    # zero removal for quantile computation
    nonzp = pixels[np.any(pixels, 1)]        
    counts = {c:np.bincount(nonzp[:,c], minlength=256) for c in colors}
    qs     = {c : 
              countQuantiles(counts[c], iles) for c in colors}

    # clunky:
    qDict = {"R8D": qs[R][0], "R9D": qs[R][1],
             "G8D": qs[G][0], "G9D": qs[G][1],
             "B8D": qs[B][0], "B9D": qs[B][1]}

    predictedThreshes = \
             dict((c, applyCurrentLM(qDict, c)) for c in colorNames)
    # print ",".join("%5.4f" % t for t in predictedThresholds.values())
    thresholdNDArray(pixels, predictedThreshes, dropSaturated=True)

origFilename = osp.join(imageDataPath, 
                        "Endosomes/10minend/ser3/10m60xendser31.TIF")

expThreshes = {"R":17, "G":41, "B":34}

currentImage = Image.open(origFilename)

# get base arrays
asArray = {}
asArray['exp'] = fromimage(currentImage).reshape((numImagePoints,3))
asArray['pred'] = asArray['exp'].copy()

# apply thresholds
thresholdNDArray(asArray['exp'], expThreshes, dropSaturated=True)
applyPredThresh(asArray['pred'])

# reconstruct images and write out
outputPath = "/home/mfenner/scipy_prep/final-images"

for name, arr in asArray.items():
    #     outImage = Image.merge(currentImage.mode, (arr[R], arr[G], arr[B]))
    outImage = Image.fromarray(arr.astype(np.uint8).reshape((1024,1024,3)))
    outImage.save(osp.join(outputPath, 
                           "10minEndosomeSeries3Slice1-"+name+".tif"))
