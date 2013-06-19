#
# This file is pure python.
#

#
# Two utilities that should disappear by making a class
# for experiments that can be used as a dictionary key
# (do I want to be a be able to specify repetition and slice
#  as well? ... ie specify an image with that key?)
# 
# see expSetupKey below
#
def dictToKeyTuple(d):
    return tuple(d.items())
    # modified APR 20
    # return tuple( (k,d[k]) for k in d )

def keyTupleToDict(t):
    return dict(t)
##################################################
##################################################
# Color names for RGB Image in PIL
##################################################
##################################################

R, G, B = 0, 1, 2
colors  = (R,G,B)
colorNames  = ("R", "G", "B")
colorToName = {R:"R", G:"G", B:"B"}
colorPairs = [(R,G), (R,B), (G,B)]
colorPairNames = ["RG", "RB", "GB"]
colorToExp = {R:"BDNF", G:"TrkB.tc", B:"Organelle"}

numImagePoints = 1024 * 1024
numIntensityValues = 256
numIntensityPairs = numIntensityValues ** 2
numIntensityTriples = numIntensityValues ** 3

##################################################
##################################################
# Getting to files from identifying information
##################################################
##################################################
imageDataPath = "/home/mfenner/scipy_prep/data/input/All_stacks_60x/"


fileNameFormat = {}
fileNameFormat["Endosomes"] = \
    "%(Organelle)s/%(Time)dminend/ser%(Series)d/" + \
    "%(Time)dm60xendser%(Series)d%(Slice)s.TIF"

fileNameFormat["Vesicles"] = \
    "%(Organelle)s/%(Stain)s/%(Time)dmin/ser%(Series)d/" + \
    "%(Time)dm60x%(ReducedStain)sser%(Series)d%(Slice)s.TIF"

fileNameFormat["Lysosomes"] = \
    "%(Organelle)s/%(Time)dmin/ser%(Series)d/" + \
    "%(Time)dm60xac17ser%(Series)d%(Slice)s.TIF"


# this is here b/c it works with the filename formats
# convert a number to the string form needed for the
# file path (1-3 needs 1,2,3; 1-12 needs 01, 02, ..., 11, 12)
def makeSliceValues(n):
    if 10 <= n <99:
        formatString = "%02d"
    elif 1 <= n <= 9:
        formatString = "%1d"
    return [formatString%x for x in range(1,n+1)]

# FIXME
# should probably make empty stain values
# "" (empty string) instead of None
#
def reduceStain(sta):
    try:
        return sta[:3].lower()
    except TypeError:
        if sta == None:
            return None
        else:
            raise TypeError
##################################################
##################################################
# Setting up experimental conditions
##################################################
##################################################

organelles = ["Endosomes", "Vesicles", "Lysosomes"]
# FIXME None --> "" ??? might be better
stains     = {"Endosomes" : [None],
              "Vesicles"   : ["Arf", "Rab4"],
              "Lysosomes"  : [None]}
stainStrs  = {"Endosomes" : ["NA"],
              "Vesicles"   : ["Arf", "Rab4"],
              "Lysosomes"  : ["NA"]}
series     = [x for x in range(1,7)]
times      = [10, 15, 30, 60, 120]


def simplifyOrgStain(org, stain):
    simpleName = org
    if simpleName == "Vesicles":
        simpleName += "(" + stain + ")"
    return simpleName

organelleStainStrings = ["Endosomes", "Lysosomes",
                         "Vesicles(Arf)", "Vesicles(Rab4)"]

#
# ser is really not part of the condition, but it simplifies things to
# specify it here
#
allExpCnds = [ {"Organelle":o, "Stain":sta, "ReducedStain":reduceStain(sta),
                "Series":ser, "Time":t}
               for o in organelles
               for sta in stains[o]
               for t in times
               for ser in series]

allNonSerExpCnds = [ {"Organelle":o, "Stain":sta, 
                      "ReducedStain":reduceStain(sta), "Time":t}
                   for o in organelles
                   for sta in stains[o]
                   for t in times]

coefficients = ["Pearson", "Manders",
                "Coloc(m)1", "Coloc(m)2",
                "Overlap(k)1", "Overlap(k)2"]

#
# FIXME
# perhaps the allExpCnds should be instances of expSetupKey
#
class expSetupKey:
    def __init__(self, org, sta, time):
        self.organelle = org
        self.stain = sta
        self.time = time
    def __hash__(self):
        return hash(self.organelle) ^ hash(self.stain) ^ hash(self.time)
    def __cmp__(self, other):
        return cmp(self.organelle, other.organelle) and \
               cmp(self.stain, other.stain) and \
               cmp(self.time, other.time)
 

##################################################
##################################################
# Computing the predicted thresholds for background
# from the linear model learned on expert values
##################################################
##################################################

#
# these should be determined from the LM formula
# brittle if LM changes FIXME
#
# iles = (.8, .9)

# def applyTM2(example, color):
#     # even worse:  ex["Slice"] is a string
#     prediction = 17.9453 # intercept
#     prediction += float(example["Slice"]) * 0.2175
#     if color == "G":
#         prediction += -14.1770
#     elif color == "R":
#         prediction += -15.3520
#     prediction += example[color+"8D"] * -0.3105
#     prediction += example[color+"9D"] * 0.7656

#     return prediction


# def applyLM(example, color)
#      prediction = 17.9453 # intercept
#      if color == "G":
#          prediction += -14.1770
#      elif color == "R":
#          prediction += -15.3520
#      prediction += example[color+"8D"] * -0.3105
#      prediction += example[color+"9D"] * 0.7656
#      return prediction

##################################################
##################################################
# Read expert thresholds from file into dictionary
##################################################
##################################################

#
# this is used to work with the expert thresholded files
# (1) develop LM; (2) compare expert against predicted; (3) other???
#
def readThresholdFileAsDictionaries(thresholdFilename):
    tfile = open(thresholdFilename)

    # not pretty, needed for fileNameFormats b/c they expect certain things
    # (can they be made "agnostic"?  probably if everything is a string
    #  but some have to be used as number ... for example in computing
    #  predicted thresholds)
    forceString = ("Slice",)
    forceInteger = ("Series", "Time")

    #
    # get dictionary keys from header line in file
    #
    orderedKeys = tfile.readline().strip().split(",")
    orderedKeys = [k.strip("\"") for k in orderedKeys]

    #
    # build a list of examples:  each example is a dict
    # with keys from orderedKeys
    #
    thresholdExamples = []

    for line in tfile.readlines():
        if line.startswith("#"):
            continue
        exampleDict = {}

        # walk in lock step across the ordered keys
        # and the values on this line ... put them
        # together as exampleDict[key] = value
        # with a couple special cases
        for value, key in zip(line.split(","), orderedKeys):
            value = value.strip()

            try:
                result = eval(value)
            except NameError: # variable symbol is unquoted string
                result = str(value)

            #
            # forced keys enforce a conversion
            #
            if key in forceString:
                result = str(result)
            elif key in forceInteger:
                result = int(result)

            if key == "Slice" and result == "10":
                # back pad previous Slice values (9 of them)
                # add 10

                # FIXME:  April 14, 2010
                # ACK!  not necessarily ... have to walk back until
                # we get a previous Repetition or back to 0
                # ugh!!!

                # hopefully fixed, apr 14, 2010

                #for prevEx in thresholdExamples[-9:]:
                #    prevEx["Slice"] = "0" + prevEx["Slice"]

                thisSeries = exampleDict["Series"]
                prevIdx = -1
                prevEx = thresholdExamples[prevIdx]
                while prevEx["Series"] == thisSeries:
                    prevEx["Slice"] = "0" + prevEx["Slice"]
                    prevIdx = prevIdx - 1
                    try:
                        prevEx = thresholdExamples[prevIdx]
                    except IndexError:
                        break

                    
            elif key == "Stain":
                # add key for reducedStain value
                exampleDict["ReducedStain"] = result[:3].lower()

            exampleDict[key] = result

        thresholdExamples.append(exampleDict)
    tfile.close()
    return thresholdExamples


# from a list of dictionaries of (expert) thresholded examples
# return groups of stacks (found by new slice >= new slice + 1 and
# series == series)
#
# added >= slice to allow for skipped slices due to bad data/image
#
# we do this b/c we process a stack at a time for
# coloc. coefficients
def breakIntoStacks(all):
    startIndex = 0

    while startIndex < len(all):
        initial=all[startIndex]
        nextSlice = int(initial["Slice"]) + 1
        thisSeries = initial["Series"]

        currentImgStack = [initial]

        i = startIndex + 1
        while i < len(all) and \
                  int(all[i]["Slice"]) >= nextSlice and \
                  all[i]["Series"] == thisSeries:
            currentImgStack.append(all[i])
            nextSlice += 1
            i += 1
        yield currentImgStack
        startIndex = i

def breakByConditions(all):
    def makeCnd(d):
        return (d["Organelle"], d["Stain"], d["Time"])
    
    currentStartIndex = 0
    currentCnd = makeCnd(all[currentStartIndex])

    while currentStartIndex < len(all):
        currentCndImgList = [all[currentStartIndex]]

        nextIndex = currentStartIndex + 1
        nextCnd = makeCnd(all[nextIndex])
        while nextCnd == currentCnd:
            currentCndImgList.append(all[nextIndex])
            nextIndex += 1
            try:
                nextCnd = makeCnd(all[nextIndex])
            except IndexError:  # nextIndex == len(all)
                break           # so exit loop, yield this stack,
                                # nextIdx -> currIdx and end outer loop
                                # possibly could: yield current
                                # then:           raise StopIteration
        yield currentCndImgList

        currentStartIndex = nextIndex
        currentCnd = nextCnd
        

def breakConditionIntoStacks(cnd):
    startIndex = 0

    while startIndex < len(cnd):
        initial=cnd[startIndex]
        nextSlice = int(initial["Slice"]) + 1
        thisSeries = initial["Series"]

        currentImgStack = [initial]

        i = startIndex + 1
        while i < len(cnd) and \
                  int(cnd[i]["Slice"]) >= nextSlice and \
                  cnd[i]["Series"] == thisSeries:
            currentImgStack.append(cnd[i])
            nextSlice += 1
            i += 1
        yield currentImgStack
        startIndex = i
