import Image, math
#from scipy.misc import fromimage
# from scipy import dot
from numpy import (mean, float_, dot, interp,
                   uint8, uint16, uint64, log10, any as np_any,
                   all as np_all)
from operator import itemgetter as ig


from DefinitionsAndUtils import *
from NumericalStatistics import quickQuantiles

#
# Return desired quantiles (R,G,B) on an image
#
def channelQuantilesFromImage(img, iles = (.5), useSampling=False):
    imageAsList = list(img.getdata())

#    if useSample:
#        print "*** !!!UNTESTED CODE. UNTESTED CODE!!! ***"
#        indices = random.sample(xrange(numPoints), int(numPoints*samplePct))
#        sampledlst = ig(*indices)(currentAsList) # a tuple
#        evensize = int(math.sqrt(samplePct) * currentImage.size[0])
#        
#        imageAsList = sampledlst[:evensize*evensize]

    
    pixelCount = len(imageAsList)

    # quantiles
    redq   =  quickQuantiles((ig(R)(c) for c in imageAsList),
                             pixelCount, quantiles=iles)
    greenq =  quickQuantiles((ig(G)(c) for c in imageAsList),
                             pixelCount, quantiles=iles)
    blueq  =  quickQuantiles((ig(B)(c) for c in imageAsList),
                             pixelCount, quantiles=iles)

    return redq, greenq, blueq
    

###################################################
###################################################
###################################################
#
# Apply a threshold {R:rthresh, etc.} to an image
#
###################################################
###################################################
###################################################

# def makePointMapper(threshold):
#     def newValue(v, T=threshold):
#         slope = (255 / (255.0 - T))
#         intercept = -slope * T
#         newV = int(slope * v + intercept)
#         return max(0, newV)
#     return newValue


# thresholds["R"] = 27
# thresholds["B"] = 35
# thresholds["G"] = 35
def thresholdNDArray(array, thresholds, dropSaturated=False):
    maxValue = 254 if dropSaturated else 255
    
    # threshold[c] --> 0 and 255 --> maxValue
    for c in colors:
        array[:,c] = interp(array[:,c],
                            [thresholds[colorToName[c]], maxValue],
                            [0,255]
                            ).astype(uint8)


# drop 255 b/c it's saturated
# def thresholdAndDrop255NDArray(array, thresholds):
#     # threshold[c] --> 0 and 254 --> 255
#     clean = array<255
#     array = array[np_all(clean, 1)]

#     for c in colors:
#         array[:,c] = interp(array[:,c],
#                             [thresholds[colorToName[c]], 254],
#                             [0,255]
#                             ).astype(uint8)

# def thresholdImage(image, thresholds):
#     sourceChannels = image.split()

#     # 0,1,2 are RGB for this image format
#     colorChannelLists = {R : list(sourceChannels[0].getdata()),
#                          G : list(sourceChannels[1].getdata()),
#                          B : list(sourceChannels[2].getdata())}

#     newColorChannelLists = {}
#     for c in colorChannelLists:
#         ptMapFunc = makePointMapper(thresholds[colorToName[c]])
#         newColorChannelLists[c] = sourceChannels[c].point(ptMapFunc)

#     threshedImage = Image.merge(image.mode,
#                                 (newColorChannelLists[R],
#                                  newColorChannelLists[G],
#                                  newColorChannelLists[B]))
#     return threshedImage


###################################################
###################################################
###################################################
#
# Image processing with numpy functions
#
###################################################
###################################################
###################################################


# safe for uints/ints
def mydot(a,b):
    prods = float_(a)*float_(b)
    return prods.sum()

#def mydot(a,b):
#    return dot(a.astype(float_), b.astype(float_))

#def mydot(a,b):
#    return float_(dot(a, b))

#
# this expects a 1d of pixels, with 2nd dimension as colors
#
# ccc => compute colocalization coefficients 
def cccOnFlatArray(ia):
    # in flat form, 2nd dimension holds channel;
    # cca -> color channel arrays
    cca = {}
    cca[R] = ia[:,R] #.flatten()  I don't think flatten does any additional
    cca[G] = ia[:,G] #.flatten()  work here
    cca[B] = ia[:,B] #.flatten()

    #
    # compute pieces
    #

    ### FIXME  consider replacing **2 pow(,2) with multiplication
    ###        there is a sum of squares function in scipy.stats
    #means          = dict((c, mean(cca[c])              for c in colors)
    # can replace mean(cca[c]) with means[c] but check it!
    meanDiffs      = dict((c, cca[c]-mean(cca[c]))      for c in colors)
    sqMeanDiffs    = dict((c, meanDiffs[c]**2)          for c in colors) # ok for float
    sumSqMeanDiffs = dict((c, sqMeanDiffs[c].sum())     for c in colors)
    sumSqrs        = dict((c, mydot(cca[c], cca[c]))    for c in colors)
    indicator      = dict((c, cca[c]>0)                 for c in colors)
    sums           = dict((c, cca[c].sum(dtype=float_)) for c in colors) # .sum() is ok

    crossDot = {}
    for c1, c2 in ((R,G), (R,B), (G,B)):
        crossDot[(c1,c2)] = mydot(cca[c1], cca[c2])

    result = {}
    for c1, c2 in ((R,G), (R,B), (G,B)):
        theseCoeffs = {}
        # scipy.dot() is ok for float 
        theseCoeffs["Pearson"]     = dot(meanDiffs[c1],meanDiffs[c2]) / \
                                     math.sqrt(dot(sumSqMeanDiffs[c1],
                                                   sumSqMeanDiffs[c2]))

        theseCoeffs["Manders"]     = crossDot[(c1,c2)] / \
                                     math.sqrt(sumSqrs[c1]*sumSqrs[c2])

        theseCoeffs["Coloc(m)1"]   = cca[c1][indicator[c2]].sum() / sums[c1]
        theseCoeffs["Coloc(m)2"]   = cca[c2][indicator[c1]].sum() / sums[c2]

        theseCoeffs["Overlap(k)1"] = crossDot[(c1,c2)] / sumSqrs[c1]
        theseCoeffs["Overlap(k)2"] = crossDot[(c1,c2)] / sumSqrs[c2]

        result[(c1,c2)] = theseCoeffs

    return result


def toyingCoeffs(ia, **auxmsrs):
    # in flat form, 2nd dimension holds channel;
    # cca -> color channel arrays
    cca = {}
    cca[R] = ia[:,R]
    cca[G] = ia[:,G]
    cca[B] = ia[:,B]

    # mean gives float64 back
    means        = dict((c, mean(cca[c]))                   for c in colors)
    indicator    = dict((c, cca[c]>0)                       for c in colors)
    indicatorSum = dict((c, indicator[c].sum(dtype=float_)) for c in colors)
    selfSum      = dict((c, cca[c].sum(dtype=float_))       for c in colors)

    bigN = ia.shape[0]

    # crossDot = {}
    # for c1, c2 in ((R,G), (R,B), (G,B)):
    #      crossDot[(c1,c2)] = mydot(cca[c1], cca[c2])

    results = {}
    for c in colors:
        myCoeffs = {}
        myCoeffs["Mean"]      = means[c]
        myCoeffs["Sum"]       = log10(selfSum[c])
        myCoeffs["NumOn"]     = log10(indicatorSum[c])
        myCoeffs["SumToOn"]   = selfSum[c] / indicatorSum[c]
        myCoeffs["SumToBigN"] = selfSum[c] / bigN
        myCoeffs["OnToBigN"]  = indicatorSum[c] / bigN
        myCoeffs["OnToMasterN"]  = indicatorSum[c] / auxmsrs["masterN"]

        results[c] = myCoeffs

    return results


# good for products of uint8s:  convert to 16s
def safedot(a,b):
    return (a.astype(uint16) * b.astype(uint16)).sum().astype(float_)

# it appears that dot uses the incoming datatype as the accumulator
# type as well (unlike "arr.sum()" which will bump up to machine
# int size)
#def safedot(a,b):
#    return dot(a.astype(uint16), b.astype(uint16)).astype(float_)

def cccOnFlatArray64bit(ia):
    # in flat form, 2nd dimension holds channel;
    # cca -> color channel arrays
    cca = {}
    cca[R] = ia[:,R] # should come in as uint8
    cca[G] = ia[:,G]
    cca[B] = ia[:,B]

    #promoted = {}
    #for c in colors:
    #    promoted[c] = uint16(cca[c]) # .astype(uint64)

    # accumulating to a uint64 should be ok (approx 2^27 values with max
    # value of 256=2^8 gives max sum of 2^35)

    # mean() products a float64
    #means          = dict((c, mean(cca[c])              for c in colors)
    meanDiffs      = dict((c, cca[c]-mean(cca[c]))      for c in colors) 
    sqMeanDiffs    = dict((c, meanDiffs[c]**2)          for c in colors)
    sumSqMeanDiffs = dict((c, sqMeanDiffs[c].sum())     for c in colors)

    sumSqrs        = dict((c, safedot(cca[c], cca[c]))
                                                        for c in colors)
    indicator      = dict((c, cca[c]>0)                 for c in colors)
    sums           = dict((c, float_(cca[c].sum()))     for c in colors)

    crossDot = {}
    for c1, c2 in ((R,G), (R,B), (G,B)):
        crossDot[(c1,c2)] = safedot(cca[c1], cca[c2])

    result = {}
    for c1, c2 in ((R,G), (R,B), (G,B)):
        theseCoeffs = {}
        theseCoeffs["Pearson"]     = dot(meanDiffs[c1],meanDiffs[c2]) / \
                                     math.sqrt(dot(sumSqMeanDiffs[c1],
                                                   sumSqMeanDiffs[c2]))

        theseCoeffs["Manders"]     = crossDot[(c1,c2)] / \
                                     math.sqrt(sumSqrs[c1]*sumSqrs[c2])

        theseCoeffs["Coloc(m)1"]   = cca[c1][indicator[c2]].sum() / sums[c1]
        theseCoeffs["Coloc(m)2"]   = cca[c2][indicator[c1]].sum() / sums[c2]

        theseCoeffs["Overlap(k)1"] = crossDot[(c1,c2)] / sumSqrs[c1]
        theseCoeffs["Overlap(k)2"] = crossDot[(c1,c2)] / sumSqrs[c2]

        result[(c1,c2)] = theseCoeffs

    return result



# def cccOnFlatZeroSensitive(ia):
#     # in flat form, 2nd dimension holds channel;
#     # cca -> color channel arrays
#     cca = {}
#     cca[R] = ia[:,R] # should come in as uint8
#     cca[G] = ia[:,G]
#     cca[B] = ia[:,B]

#     meanDiffs      = dict((c, cca[c]-mean(cca[c]))      for c in colors) 
#     sqMeanDiffs    = dict((c, meanDiffs[c]**2)          for c in colors)
#     sumSqMeanDiffs = dict((c, sqMeanDiffs[c].sum())     for c in colors)

#     result = {}
#     for c1, c2 in ((R,G), (R,B), (G,B)):
#         theseCoeffs = {}
#         theseCoeffs["Pearson"]     = dot(meanDiffs[c1],meanDiffs[c2]) / \
#                                      math.sqrt(dot(sumSqMeanDiffs[c1],
#                                                    sumSqMeanDiffs[c2]))    
#         result[(c1,c2)] = theseCoeffs
#     return result

#
# safe for 64-bit; unsafe for 32-bit
#
# def cccOnFlatZeroInsensitive(ina):
#     # in flat form, 2nd dimension holds channel;
#     # cca -> color channel arrays
#     cca = {}
#     # should come in as uint8
#     ia = ina[np_any(ina, 1)]
#     cca[R] = ia[:,R]; cca[G] = ia[:,G]; cca[B] = ia[:,B]

#     # accumulating to a uint64 should be ok (approx 2^27 values with max
#     # value of 256=2^8 gives max sum of 2^35)
#     sumSqrs        = dict((c, safedot(cca[c], cca[c]))
#                                                         for c in colors)
#     indicator      = dict((c, cca[c]>0)                 for c in colors)
#     sums           = dict((c, float_(cca[c].sum()))     for c in colors)

#     crossDot = {}
#     for c1, c2 in ((R,G), (R,B), (G,B)):
#         crossDot[(c1,c2)] = safedot(cca[c1], cca[c2])

#     result = {}
#     for c1, c2 in ((R,G), (R,B), (G,B)):
#         theseCoeffs = {}
#         theseCoeffs["Manders"]     = crossDot[(c1,c2)] / \
#                                      math.sqrt(sumSqrs[c1]*sumSqrs[c2])

#         theseCoeffs["Coloc(m)1"]   = cca[c1][indicator[c2]].sum() / sums[c1]
#         theseCoeffs["Coloc(m)2"]   = cca[c2][indicator[c1]].sum() / sums[c2]

#         theseCoeffs["Overlap(k)1"] = crossDot[(c1,c2)] / sumSqrs[c1]
#         theseCoeffs["Overlap(k)2"] = crossDot[(c1,c2)] / sumSqrs[c2]

#         result[(c1,c2)] = theseCoeffs

#     return result
