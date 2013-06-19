#!/usr/bin/env python

#
# Imports (duh)
#
# python
import random

# numpy, scipy, matplotlib, PIL
import numpy as np
import matplotlib.pyplot as pp
from matplotlib.colors import LinearSegmentedColormap as LSC

# mine
from DefinitionsAndUtils import *

colorOrder = "bgrcmy"

myColorMaps = []
for c in colorOrder:
    thisMap = LSC.from_list("mine", ["w", c])
    myColorMaps.append(thisMap)

# # for testing
# # from rpy import r as externR
# def rquant(vals, q):
#     return externR.quantile(vals, q, type=1)

# # check quantiles using R (on original values)
# def checkQuantiles(arr):
#     neededQs = (.1, .25, .5, .75, .9)
#     print "x axis linear scale (red)?: "
#     rqs = rquant(arr[:,R], neededQs)
#     for k in sorted(rqs.keys()):
#         print k, "%3.3f"%rqs[k],
#     print "\n" + "*"*50

#     print "y axis linear scale (green)?: "
#     rqs = rquant(arr[:,G], neededQs)
#     for k in sorted(rqs.keys()):
#         print k, "%3.3f"%rqs[k],
#     print "\n" + "*"*50


# numpy based zero removal
def removeAllZeroCases(arr):
    # np.any(arr,1) returns true only if one or more along the 2nd
    # dimension (aka, dimension 1) is not False (non-zero)
    return arr[np.any(arr,1)]


# remove zeros on specified channels
def removeZeroes(arr):
    nonz = []
    for p in arr:
        # if either is nonzero, keep it
        if p[R] != 0 or p[G] != 0:
            nonz.append(p)
    nonz = np.array(nonz)
    return nonz

# sample to given percent or total # values
def samplePoints(arr, size):
    numPoints = arr.size / 3
    if size <= 1.0:
        size = int(numPoints * size)
    # indices = random.sample(xrange(numPoints), int(numPoints * .3))
    indices = random.sample(xrange(numPoints), size)
    return arr[indices,:]

# verified against r.quantile
def countQuantiles(freqs, qs, asprobs=False):
    if not asprobs: # must convert counts to probs first
        cumCounts = freqs.cumsum()
        N = cumCounts[-1]
        cumProbs = cumCounts / float(N)
    else:
        cumProbs = freqs.cumsum()
        
    cumIdx = 0; quantIdx = 0

    foundScores = []
    # while we haven't hit the top quantile
    while cumProbs[cumIdx] < qs[-1]:
        while quantIdx < len(qs) and cumProbs[cumIdx] >= qs[quantIdx]:
            foundScores.append(cumIdx)
            quantIdx += 1
        cumIdx += 1

    # clean up any quantiles left over if cumSum "jumps"
    # above needed[-1] and skips other quantiles along the way
    scoresLeft = len(qs) - len(foundScores)
    if scoresLeft > 0:
        foundScores += [cumIdx]*scoresLeft

    return foundScores

#
# to make tufte grids
#step 0
#from mpl_toolkits.axes_grid import AxesGrid
#
#fig = pp.figure(1, (20, 24))
# step 1
#theGrid = AxesGrid(fig, 111, nrows_ncols = (12, 10), axes_pad=0.1)
# step 2
#addToTufteAxesGrid(theGrid, (r, t), graphprobs)
# step 3
#finalizeTufteAxesGrid(theGrid)
#pp.savefig("tufte-real.pdf")
#pp.show()

# scatterCoord is relative to other scatter grams
# not to everything on the AxesGrid
def addToTufteAxesGrid(thisGrid, scatterCoords, data):
    nrows = thisGrid._nrows
    ncols = thisGrid._ncols
    def quickb(axes, data, orientation="hor"):
        if orientation == "hor":
            drawer = axes.hlines
        else:
            drawer = axes.vlines

        neededQs = (.1, .25, .5, .75, .9)
        iles = countQuantiles(data, qs=neededQs, asprobs=True)

        # thin line deciles
        drawer(0, iles[0], iles[4], linewidth=1)
        # thick line quartiles
        drawer(0, iles[1], iles[3], linewidth=4)
        # vertical break line at median
        delta = 1
        drawer(0, iles[2]-delta, iles[2]+delta, linewidth=4, color='w')


    # r, t  repetition, time
    # rep is row     (6 total reps)
    # time is col    (5 total times .... 5 pairs across a row)
    def rtToGridScatt(r,t):
        rowValue = r * (ncols * 2) # 2 Axes per item
                                   # and advance two total rows
        colValue = (t * 2) + 1
        return rowValue + colValue

    def rtToGridBP(r,t):
        # boxplot - 1       for vertical (y)   (back one spot)
        # boxplot + rowsize for horizontal (x) (forward ncol spots)
        return (rtToGridScatt(r,t) - 1 ,
                rtToGridScatt(r,t) + ncols) 

    rep, time = scatterCoords
    gridLoc = rtToGridScatt(rep, time)

    thisCM = myColorMaps[rep]


    # use 1.0-data to get black points on white background
    # use interpolation="nearest" for "discrete" bins

    # SEE:  ***NOTE*** axis orientation below
    # ndarray versus image flips row/col so, we have to take the
    # transpose

    # natural (untested)
    thisGrid[gridLoc].imshow(data.T, interpolation="bilinear",
                             cmap=thisCM, extent=(0,255,0,255),
                             origin="lower")

    # flipped
    #thisGrid[gridLoc].imshow(1.0-data, interpolation="bilinear",
    #                         cmap=pp.cm.gray, extent=(0,255,0,255),
    #                         origin="lower")

    # note sure why this is still here:
    #thisGrid[gridLoc].imshow(data, cmap=pp.cm.gray,
    #                         interpolation='nearest', extent=(0,255,0,255))


    #
    # this works for the data as probabilities
    # (sum over rows to get marginal cols /
    #  sum over cols to get marginal rows)
    #
    bpLocs = rtToGridBP(rep,time)
    #
    # ***NOTE*** axis orientation
    #
    # note these must correspond to the "data" or "data.T" in imshow
    # above ... imshow operates opposite of numpy ...
    # numpy[x,y] is imshow[y,x] ... however, axis=0 ... the numpy x ...
    # adds over the OTHER stuff ... so numpy.sumover(x,:) gives totals
    # for a particular y value ... putting that on the bottom corresponds
    # to data
    #
    # in short:  with data, use data.sum(axis=0) on horizontal
    #                       use data.sum(axis=1) on vertical
    #            (but if data[R,G] then R on horizontal, G on vertical
    #
    # with data.T, use data.sum(axis=1) on vertical;
    #              use data.sum(axis=0) on horizontal
    #

    # natural (untested)
    quickb(thisGrid[bpLocs[1]], data.sum(axis=1), orientation="hor")  # xs
    quickb(thisGrid[bpLocs[0]], data.sum(axis=0), orientation="vert") # ys
    
    # flipped
    #quickb(thisGrid[bpLocs[1]], data.sum(axis=0), orientation="hor")  # xs
    #quickb(thisGrid[bpLocs[0]], data.sum(axis=1), orientation="vert") # ys

def finalizeTufteAxesGrid(thisGrid):
    nrows = thisGrid._nrows
    ncols = thisGrid._ncols
    ngrids = thisGrid.ngrids

    # turn off all axes
    for g in range(ngrids):
        for side in ("top", "bottom", "left", "right"):
            thisGrid[g].axis[side].set_visible(False)

    # set y-lims on 0,2,4,6... everything
    # the mf'ing .imshow() insists on padding -50, 300 on data 0,255
    # even with the "extent" parameter
    # could it be from the vlines?  doubt it, but ...
    for row in range(0,nrows, 2):
        for col in range(ncols):
            # thisGrid[row * 10 + col].set_ylim(0,255)
            thisGrid[row * ncols + col].set_ylim(0,255)

    # set thin verticals (all in column 0, 2, 4, ..., ncols-2)
    # 10 columns ... 0-9 ... thin on 0,2,4,6,8
    # in each even column
    for row in range(nrows):
        for col in range(0,ncols,2):
            # thisGrid[row * 10 + col].set_xlim(-3,3)            
            thisGrid[row * ncols + col].set_xlim(-3,3)

    # set thin horizontals (all in row 1, 3, 5, ... nrows)
    # 10 rows ... 0-9 ... thin on 1,3,5,7,9
    # in each odd row
    for col in range(ncols):
        for row in range(1, nrows+1, 2): # +1 to get the last one
            # thisGrid[row*10 + col].set_ylim(-3,3)
            thisGrid[row*ncols + col].set_ylim(-3,3)

def convert(num, base):
    if num < base:
        return (num,)
    else:
        return convert(num/base, base) + (num%base,)

# l is length of number in base-b digits
def wconvert(num, base, l):
    result = convert(num, base)
    return (0,)* (l-len(result)) + result

def timeToIdx(time):
    times = [10, 15, 30, 60, 120]
    return times.index(time)


def fullJoint(arr):
    b3val = (numIntensityValues * numIntensityValues * arr[:,R].astype(np.uint32)) + \
            (numIntensityValues * arr[:,G]) + \
            (1                  * arr[:,B])
    qcounts = np.bincount(b3val)
    qcounts.resize(numIntensityTriples)   # pad out with zeros

    import scipy.sparse as sps
    qcounts = sps.lil_matrix(np.vstack((qcounts,
                                        np.zeros(numIntensityTriples))))

    N = qcounts.sum()
    probabilities = qcounts / float(N)
    return probabilities

def toProbs(channel1, channel2,
            removeNonresponders=False,
            removeSaturated=False,
            useLogscaleCounts=False):
    # could incorporate sampling here
    #stackArray = samplePoints(stackArray, .3)
    #stackArray = samplePoints(stackArray, 100)

    b3val = (numIntensityValues * channel1.astype(np.uint32)) + \
            (1                  * channel2)
    qcounts = np.bincount(b3val)
    qcounts.resize(numIntensityPairs)   # pad out with zeros
    qcounts = qcounts.reshape((numIntensityValues,  
                               numIntensityValues))

    if removeNonresponders:
        qcounts[0,0] = 0

    if removeSaturated:
        qcounts[255,   0] = 0
        qcounts[0,   255] = 0
        qcounts[255, 255] = 0
        
    if useLogscaleCounts:
        # logcts = np.log(np.log(qcounts+1))
        logcts = np.log(qcounts+1)
        N = logcts.sum()
        probabilities = logcts / float(N)
    else:
        N = qcounts.sum()
        probabilities = qcounts / float(N)

    return probabilities



# to use for single histograms:
#
# #histotype=["mine", "npheat", "hexbin-raw",
# #           "hexbin-counts", "simplescatter"]
# #histotype=["mine"]
# histotype = ["hexbin-counts"]
# if histotype:
#     makeSingleHistogram(histotype, stackArray, graphprobs)
#     raw_input("Enter to continue:")

def makeSingleHistogram(histo, stackArray, graphprobs):
    figCt = 1
    forHM=samplePoints(stackArray, .3)
    if "mine" in histo:
        # from probs
        figCt += 1; pp.figure(figCt)
        pp.title("Mine: probs -> image (sketchy)")
        # invert the values for black on white
        invprobs = 1.0-graphprobs
        pp.imshow(invprobs, interpolation="nearest",
                  cmap=pp.cm.gray, extent=(0,256,0,256),
                  origin="lower")
    if "npheat" in histo:
        # from sampled raw data
        hm, xedg, yedg = np.histogram2d(forHM[:,R], forHM[:,G],
                                        bins=64)
        # remove pesky (overbearing) 0,0 count
        hm[0,0] = 0
        # use log scale counts
        hm = np.log(hm+1)
        extent = [xedg[0], xedg[-1], yedg[0], yedg[-1]]
        figCt += 1; pp.figure(figCt)
        pp.title("NP 2D-Histogram (Fixed) -- no (0,0)")
        # pp.imshow(hm, extent=extent, interpolation="nearest", origin="lower")
        # bilinear is default
        pp.imshow(hm, extent=extent, interpolation="bilinear", origin="lower")
    if "hexbin-raw" in histo:
        # from sampled raw data
        figCt += 1; pp.figure(figCt)
        pp.title("Hexbin from Raw Data (sampled)")
        pp.hexbin(forHM[:,R], forHM[:,G], extent=(0, 255, -10, 255),
                  gridsize=64, cmap=pp.cm.jet, bins="log")
    if "hexbin-counts" in histo:
        # from probs

        # simple technique to unroll implicit indices to explicit
        # x,y values
        tx = ty = np.linspace(0,255, 256)
        xs, ys = np.meshgrid(tx, ty)
        xs = xs.ravel(); ys = ys.ravel(); cs = graphprobs.ravel()

        figCt += 1; pp.figure(figCt)
        pp.title("Hexbin from computed probs")
        # note, must ravel (unravel?) the 2-d graphprobs
        pp.hexbin(xs, ys, C=cs, extent=(-10,255,-10,255),
                  gridsize=64, cmap=pp.cm.jet, marginals=True)
    if "simplescatter" in histo:
        # from raw data
        figCt += 1; pp.figure(figCt)
        pp.title("Simple scatter plot of Raw Data (sampled) -- on or off")
        pp.scatter(forHM[:,R], forHM[:,G], s=1)
    pp.show()
