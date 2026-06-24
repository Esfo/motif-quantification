from generalfunctions import intersection_merge, radius_neighbors_hard_tolerance
from elementalcomponents import proton

from collections import defaultdict
from itertools import product
from time import time
import numpy as np
import pickle

def charge_handling(processingdirectory, nprocs):
    
    regionfile = ''.join((processingdirectory, 'regions.pickle'))
    with open(regionfile, 'rb') as pick:
        regions = pickle.load(pick)
    
    trackedgroupsfile = ''.join((processingdirectory, 'trackedgroups.pickle'))
    with open(trackedgroupsfile, 'rb') as pick:
        trackedgroups = pickle.load(pick)
    
    nodistfile = ''.join((processingdirectory, 'nodists.pickle'))
    with open(nodistfile, 'rb') as pick:
        nodists = pickle.load(pick)
    
    solodistfile = ''.join((processingdirectory, 'solodists.pickle'))
    with open(solodistfile, 'rb') as pick:
        solodists = pickle.load(pick)
    
    finalindfile = ''.join((processingdirectory, 'finalindex.pickle'))
    with open(finalindfile, 'rb') as pick:
        finaldefiniteind = pickle.load(pick)
    
    t1 = time()
    
    distributionmasses = {} #distid: ordered masses
    distributioncharges = {} #distid: charge
    distributionsoflines = {} #line: distid
    linesofdistributions = {} #distid: mass-ordered linedkeys
    distributiontimelimits = {} #distid [starting rt, ending rt]
    distributionintensities = {} #distid: mass-ordered intensities
    #distributionsbycharge = defaultdict(dict) #charge: dists: mass-ordered linekeys
    distributionsbycharge = defaultdict(list) #charge: dists: mass-ordered linekeys
    for charge, dists in solodists.items():
        for dist, lines in dists.items():
            dmasses = regions[lines,7]
            lineorder = dmasses.argsort().tolist()
            sortedlines = [lines[i] for i in lineorder]
            sortedmasses = regions[sortedlines,7]
            dintensities = regions[sortedlines,5]
            rtlimits = regions[sortedlines,2:4]
            minrt = rtlimits.min()
            maxrt = rtlimits.max()
            distributionmasses[dist] = sortedmasses
            distributioncharges[dist] = charge
            for line in lines:
                distributionsoflines[line] = dist
            linesofdistributions[dist] = sortedlines
            distributiontimelimits[dist] = [minrt, maxrt]
            distributionintensities[dist] = dintensities
            distributionsbycharge[charge].append(dist)
    
    nodistmasses = regions[nodists,7]
    
    nodistkeys = []
    dr = finaldefiniteind
    for line in nodists.tolist():
        dreg = regions[line]
        dmass = dreg[7]
        dintensity = dreg[5]
        rtlimit = dreg[2:4]
        minrt = rtlimit.min()
        maxrt = rtlimit.max()
        distributionmasses[dr] = np.array([dmass])
        distributioncharges[dr] = 0
        distributionsoflines[line] = dr
        linesofdistributions[dr] = [line]
        distributiontimelimits[dr] = [minrt, maxrt]
        distributionintensities[dr] = np.array([dintensity])
        distributionsbycharge[0].append(dr)
        nodistkeys.append(dr)
        dr += 1 #continuing from solodists count
    
    massranges = regions[:,:2]
    minmass = massranges.min() - 1
    maxmass = massranges.max() + 1
    
    sortednodistkeys = np.array(nodistkeys)[nodistmasses.argsort()].tolist()
    sortednodistmasses = np.sort(nodistmasses)
    
    if nodists.size > 0:
        #keeping the 0-charge nodists out
        chargeiterations = sorted(distributionsbycharge)[1:]
    else:
        chargeiterations = sorted(distributionsbycharge)
    
    nodistcharges = defaultdict(dict) #matched dist: nodistkey: would-be charge of nodist distribution
    chargegroups = []
    oldcharge = -1
    for charge in chargeiterations: #ignoring 0-charge nodists here
        if charge - oldcharge > 1: #adjacent charge searches only
            moving = False
        dists = distributionsbycharge[charge]
        distkeys = []
        massdiffs = []
        distmasses = []
        for d in dists:
            dmasses = distributionmasses[d]
            basemasses = dmasses * charge - proton * charge
            distkeys.extend(d for _ in dmasses)
            distmasses.extend(basemasses)
            massdiff = np.diff(basemasses)
            massdiffs.extend(massdiff.tolist())
        massdiffs = np.array(massdiffs)
        diffcuts = np.abs(proton - massdiffs)
        ctol = diffcuts.mean() #max for for exploration I suppose
        distmasses = np.array(distmasses)
        distkeys = np.array(distkeys)
        distkeys = distkeys[distmasses.argsort()].tolist()
        distmasses = np.sort(distmasses).tolist()
        if moving:
            #train it on the new and find old masses + radii
            roundmatches = set()
            oldmatches = set()
            matches = radius_neighbors_hard_tolerance(distmasses, oldmasses, oldradii)
            for ind, fm in matches.items():
                o = oldkeys[ind]
                matchkeys = []
                for m in fm:
                    ominrt, omaxrt = distributiontimelimits[o]
                    omasses = distributionmasses[o]
                    obmasses = omasses * oldcharge - proton * oldcharge
                    intensities = distributionintensities[o]
                    intensityranks = intensities.argsort()[::-1] #i want to see if substituting the ranked intensities for the ranked intensities diffs (ordered by mass) would be a better mechanism, this currently messes with the mass alignment [being the basis for the rank comparison] below and throws a valueerror
                    mk = distkeys[m]
                    nminrt, nmaxrt = distributiontimelimits[mk]
                    nmasses = distributionmasses[mk]
                    nbmasses = nmasses * charge - proton * charge
                    #below is requiring a majority of the matchable masses have sufficiently overlapping retention times
                    if ominrt < nmaxrt and omaxrt > nminrt: #overlap exists
                        basemasses = [obmasses, nbmasses]
                        dlines = [[regions[j,2:4].tolist() for j in linesofdistributions[i]] for i in [o, mk]]
                        sizes = [i.size for i in basemasses]
                        maxind = sizes.index(max(sizes))
                        lineup = basemasses[maxind]
                        retentionboundaries = defaultdict(list)
                        for n, (sm, dlims) in enumerate(zip(basemasses, dlines)):
                            sdiff = np.abs(lineup - sm[:,None])
                            alignmentloc = np.argwhere(sdiff == sdiff.min())[0]
                            alignmentloc -= alignmentloc.min()
                            luind = alignmentloc[1] - alignmentloc[0]
                            outinds = [luind, luind + sm.size]
                            for ind, dlim in zip(range(*outinds), dlims):
                                retentionboundaries[ind].append(dlim)
                        overpass = 0
                        matchables = 0
                        for ind, lims in retentionboundaries.items():
                            if len(lims) > 1:
                                matchables += 1
                                (lminrt, lmaxrt), (rminrt, rmaxrt) = lims
                                if lminrt > rminrt and lmaxrt < rmaxrt: #old encompassed
                                    overpass += 1
                                elif rminrt > lminrt and rmaxrt < lmaxrt: #new encompassed
                                    overpass += 1
                                else:
                                    overlap = min(rmaxrt, lmaxrt) - max(rminrt, lminrt)
                                    fullrange = max(rmaxrt, lmaxrt) - min(rminrt, lminrt)
                                    if overlap / fullrange > 0.5:
                                        overpass += 1
                        #i think that rather than this, i'd prefer to check if the stronger ions are matched here, or maybe if they have the best matches?
                        if overpass > matchables / 2:
                            diffmatrix = np.abs(obmasses - nbmasses[:,None])
                            matchmatrix = diffmatrix < ctol
                            mmshape = matchmatrix.shape
                            majoraxis = np.argmax(mmshape)
                            minoraxis = np.argmin(mmshape)
                            matrixmatches = matchmatrix.any(axis=majoraxis)
                            if matrixmatches.sum() > matrixmatches.size / 2: #matching at least half of the smaller with >=, make this a 'majority' with >
                                alignmentloc = np.argwhere(diffmatrix == diffmatrix.min())[0]
                                minindex = min(alignmentloc)
                                alignmentloc = [i-minindex for i in alignmentloc] #for when the closest match isn't the earliest matching index
                                rmax = nbmasses.size - alignmentloc[0]
                                lmax = obmasses.size - alignmentloc[1]
                                maxsize = min(rmax, lmax)
                                matchintensities = distributionintensities[mk]
                                matchintensityranks = matchintensities.argsort()[::-1]
                                scoutrankslice = intensityranks[alignmentloc[1]:alignmentloc[1]+maxsize]
                                matchrankslice = matchintensityranks[alignmentloc[0]:alignmentloc[0]+maxsize]
                                slicesubtraction = np.abs(scoutrankslice - matchrankslice)
                                if slicesubtraction.sum() <= slicesubtraction.size - 1: #close order matching passes
                                    matchkeys.append(mk)
                                    roundmatches.add(mk)
                                    oldmatches.add(o)
                                #else:
                                    #record non-matches i suppose
                if matchkeys:
                    matchkeys.append(o)
                    chargegroups.append(matchkeys)
                #else:
                    #nodist match against old, nodists as higher charge
                    #nonoldmatches.add(o)
            ndchargemasses = (sortednodistmasses * oldcharge - proton * oldcharge).tolist()
            #nodist match against new, nodists as lower charge
            #there are other nonmatches that don't pass the above, add them here I think
            nonmatched = set(distkeys).difference(roundmatches)
            nonmatchedmasses = np.array([distmasses[n] for n in range(len(distkeys)) if distkeys[n] in nonmatched])
            nonmatchedkeys = [distkeys[n] for n in range(len(distkeys)) if distkeys[n] in nonmatched]
            nonmatchedkeys = np.array(nonmatchedkeys)[nonmatchedmasses.argsort()].tolist()
            nonmatchedmasses = np.sort(nonmatchedmasses)
            if nonmatchedmasses.size > 0: #this is a new addition, not fully tested
                downmatches = radius_neighbors_hard_tolerance(ndchargemasses, nonmatchedmasses.tolist(), oldradii)
                for ind, fm in downmatches.items():
                    o = nonmatchedkeys[ind]
                    matchkeys = []
                    for m in fm:
                        ominrt, omaxrt = distributiontimelimits[o]
                        omasses = distributionmasses[o]
                        obmasses = omasses * charge - proton * charge
                        intensities = distributionintensities[o]
                        intensityranks = intensities.argsort()[::-1]
                        mk = sortednodistkeys[m]
                        matchintensities = distributionintensities[mk]
                        if matchintensities.max() < intensities.max():
                            nminrt, nmaxrt = distributiontimelimits[mk]
                            nmasses = distributionmasses[mk]
                            nbmasses = nmasses * oldcharge - proton * oldcharge
                            #below is requiring a majority of the matchable masses have sufficiently overlapping retention times
                            if ominrt < nmaxrt and omaxrt > nminrt: #overlap exists
                                basemasses = [obmasses, nbmasses]
                                dlines = [[regions[j,2:4].tolist() for j in linesofdistributions[i]] for i in [o, mk]]
                                sizes = [i.size for i in basemasses]
                                maxind = sizes.index(max(sizes))
                                lineup = basemasses[maxind]
                                retentionboundaries = defaultdict(list)
                                for n, (sm, dlims) in enumerate(zip(basemasses, dlines)):
                                    sdiff = np.abs(lineup - sm[:,None])
                                    alignmentloc = np.argwhere(sdiff == sdiff.min())[0]
                                    alignmentloc -= alignmentloc.min()
                                    luind = alignmentloc[1] - alignmentloc[0]
                                    outinds = [luind, luind + sm.size]
                                    for ind, dlim in zip(range(*outinds), dlims):
                                        retentionboundaries[ind].append(dlim)
                                overpass = 0
                                matchables = 0
                                for ind, lims in retentionboundaries.items():
                                    if len(lims) > 1:
                                        matchables += 1
                                        (lminrt, lmaxrt), (rminrt, rmaxrt) = lims
                                        if lminrt > rminrt and lmaxrt < rmaxrt: #old encompassed
                                            overpass += 1
                                        elif rminrt > lminrt and rmaxrt < lmaxrt: #new encompassed
                                            overpass += 1
                                        else:
                                            overlap = min(rmaxrt, lmaxrt) - max(rminrt, lminrt)
                                            fullrange = max(rmaxrt, lmaxrt) - min(rminrt, lminrt)
                                            if overlap / fullrange > 0.5:
                                                overpass += 1
                                if overpass > matchables / 2:
                                    diffmatrix = np.abs(obmasses - nbmasses[:,None])
                                    matchmatrix = diffmatrix < ctol
                                    mmshape = matchmatrix.shape
                                    majoraxis = np.argmax(mmshape)
                                    minoraxis = np.argmin(mmshape)
                                    matrixmatches = matchmatrix.any(axis=majoraxis)
                                    if matrixmatches.sum() > matrixmatches.size / 2: #matching at least half of the smaller with >=, make this a 'majority' with >
                                        alignmentloc = np.argwhere(diffmatrix == diffmatrix.min())[0]
                                        minindex = min(alignmentloc)
                                        alignmentloc = [i-minindex for i in alignmentloc] #for when the closest match isn't the earliest matching index
                                        rmax = nbmasses.size - alignmentloc[0]
                                        lmax = obmasses.size - alignmentloc[1]
                                        maxsize = min(rmax, lmax)
                                        matchintensityranks = matchintensities.argsort()[::-1]
                                        scoutrankslice = intensityranks[alignmentloc[1]:alignmentloc[1]+maxsize]
                                        matchrankslice = matchintensityranks[alignmentloc[0]:alignmentloc[0]+maxsize]
                                        slicesubtraction = np.abs(scoutrankslice - matchrankslice)
                                        if slicesubtraction.sum() <= slicesubtraction.size - 1: #close order matching passes
                                            nodistcharges[o][mk] = oldcharge
                    if matchkeys:
                        matchkeys.append(o)
                        chargegroups.append(matchkeys)
                    #else:
                        #nodist match against old, nodists as higher charge
                        #nonoldmatches.add(o)
            #nodist match against old, nodists as higher charge
            nonoldmatches = set(oldkeys).difference(oldmatches) #if you end up needing more competition in the mix, you can just match everything to nodists and you'll get competing matches, I've left them out via routes like this so as to lessen the complications atm, competing nodist vs matches of the same charge can be annoying when determining if a distribution should extend beyond a nodist match
            nonmatchedmasses = np.array([oldmasses[n] for n in range(len(oldkeys)) if oldkeys[n] in nonoldmatches])
            nonmatchedkeys = [oldkeys[n] for n in range(len(oldkeys)) if oldkeys[n] in nonoldmatches]
            nonmatchedkeys = np.array(nonmatchedkeys)[nonmatchedmasses.argsort()].tolist()
            nonmatchedmasses = np.sort(nonmatchedmasses)
            if nonmatchedmasses.size > 0: #this is a new addition, not fully tested
                upmatches = radius_neighbors_hard_tolerance(ndchargemasses, nonmatchedmasses.tolist(), oldradii)
                for ind, fm in upmatches.items():
                    o = nonmatchedkeys[ind]
                    matchkeys = []
                    for m in fm:
                        ominrt, omaxrt = distributiontimelimits[o]
                        omasses = distributionmasses[o]
                        obmasses = omasses * oldcharge - proton * oldcharge
                        intensities = distributionintensities[o]
                        intensityranks = intensities.argsort()[::-1]
                        mk = sortednodistkeys[m]
                        matchintensities = distributionintensities[mk]
                        if matchintensities.max() < intensities.max():
                            nminrt, nmaxrt = distributiontimelimits[mk]
                            nmasses = distributionmasses[mk]
                            nbmasses = nmasses * charge - proton * charge
                            #below is requiring a majority of the matchable masses have sufficiently overlapping retention times
                            if ominrt < nmaxrt and omaxrt > nminrt: #overlap exists
                                basemasses = [obmasses, nbmasses]
                                dlines = [[regions[j,2:4].tolist() for j in linesofdistributions[i]] for i in [o, mk]]
                                sizes = [i.size for i in basemasses]
                                maxind = sizes.index(max(sizes))
                                lineup = basemasses[maxind]
                                retentionboundaries = defaultdict(list)
                                for n, (sm, dlims) in enumerate(zip(basemasses, dlines)):
                                    sdiff = np.abs(lineup - sm[:,None])
                                    alignmentloc = np.argwhere(sdiff == sdiff.min())[0]
                                    alignmentloc -= alignmentloc.min()
                                    luind = alignmentloc[1] - alignmentloc[0]
                                    outinds = [luind, luind + sm.size]
                                    for ind, dlim in zip(range(*outinds), dlims):
                                        retentionboundaries[ind].append(dlim)
                                overpass = 0
                                matchables = 0
                                for ind, lims in retentionboundaries.items():
                                    if len(lims) > 1:
                                        matchables += 1
                                        (lminrt, lmaxrt), (rminrt, rmaxrt) = lims
                                        if lminrt > rminrt and lmaxrt < rmaxrt: #old encompassed
                                            overpass += 1
                                        elif rminrt > lminrt and rmaxrt < lmaxrt: #new encompassed
                                            overpass += 1
                                        else:
                                            overlap = min(rmaxrt, lmaxrt) - max(rminrt, lminrt)
                                            fullrange = max(rmaxrt, lmaxrt) - min(rminrt, lminrt)
                                            if overlap / fullrange > 0.5:
                                                overpass += 1
                                if overpass > matchables / 2:
                                    diffmatrix = np.abs(obmasses - nbmasses[:,None])
                                    matchmatrix = diffmatrix < ctol
                                    mmshape = matchmatrix.shape
                                    majoraxis = np.argmax(mmshape)
                                    minoraxis = np.argmin(mmshape)
                                    matrixmatches = matchmatrix.any(axis=majoraxis)
                                    if matrixmatches.sum() > matrixmatches.size / 2: #matching at least half of the smaller with >=, make this a 'majority' with >
                                        alignmentloc = np.argwhere(diffmatrix == diffmatrix.min())[0]
                                        minindex = min(alignmentloc)
                                        alignmentloc = [i-minindex for i in alignmentloc] #for when the closest match isn't the earliest matching index
                                        rmax = nbmasses.size - alignmentloc[0]
                                        lmax = obmasses.size - alignmentloc[1]
                                        maxsize = min(rmax, lmax)
                                        matchintensityranks = matchintensities.argsort()[::-1]
                                        scoutrankslice = intensityranks[alignmentloc[1]:alignmentloc[1]+maxsize]
                                        matchrankslice = matchintensityranks[alignmentloc[0]:alignmentloc[0]+maxsize]
                                        slicesubtraction = np.abs(scoutrankslice - matchrankslice)
                                        if slicesubtraction.sum() <= slicesubtraction.size - 1: #close order matching passes
                                            nodistcharges[o][mk] = charge
                    if matchkeys:
                        matchkeys.append(o)
                        chargegroups.append(matchkeys)
        oldmasses = distmasses.copy()
        oldradii = ctol
        oldkeys = distkeys.copy()
        oldcharge = charge
        moving = True
    
    print(time() - t1, 'initial charge matching')
    t2 = time()
    
    #combining redundant matches
    chargesets = intersection_merge(chargegroups)
    
    chargelayers = defaultdict(lambda: defaultdict(set)) #charge groupid: charge: [connections]
    extralayercharges = defaultdict(dict) #charge groupid: nodist: charge
    distchargesofnodists = defaultdict(dict) #nodist line: nodist dist (this is currently the dist that isnt the nodist dist???): charge
    for n, cs in enumerate(chargesets):
        for c in cs:
            chargelayers[n][distributioncharges[c]].add(c)
            if c in nodistcharges:
                for nd, ndc in nodistcharges[c].items():
                    chargelayers[n][ndc].add(nd)
                    extralayercharges[n][nd] = ndc
                    #distchargesofnodists[linesofdistributions[nd][0]][c] = ndc #i believe this was an error previously, it raised errors later on
                    distchargesofnodists[linesofdistributions[nd][0]][nd] = ndc
    
    cid = 0
    conpairs = {} #cid: adjacent charge-state matches as pairs
    chargeconsets = set() #finalized charge-state groups for lookup
    chargecongroups = defaultdict(dict) #connection: [all chargegroups its involved in]
    #distmatchcount = defaultdict(int) #it's for plotting, I guess
    for ck, cl in chargelayers.items():
        for connections in product(*cl.values()):
            connectioncharges = {}
            for con in connections:
                if distributioncharges[con] > 0:
                    bc = distributioncharges[con]
                else:
                    bc = extralayercharges[ck][con]
                connectioncharges[con] = bc
            connections = sorted(connections, key=lambda x: connectioncharges[x])
            clen = len(connections)
            startind = 0
            endind = 2
            coniters = []
            restart = True
            while endind <= clen:
                if restart:
                    tempind = 0
                    restart = False
                coniters.append(tuple(connections[startind+tempind:endind+tempind]))
                tempind += 1
                if endind + tempind > clen:
                    restart = True
                    endind += 1
            for cons in coniters:
                basemasses = []
                for con in cons:
                    if distributioncharges[con] > 0:
                        bc = distributioncharges[con]
                    else:
                        bc = extralayercharges[ck][con]
                    bm = distributionmasses[con] * bc - proton * bc
                    basemasses.append(bm)
                intensities = [distributionintensities[i] for i in cons] 
                retentiontimes = [distributiontimelimits[i] for i in cons]
                #^this wouldn't be so great, I want individual isotopomer overlaps to be assessed here
                dlines = [[regions[j,2:4].tolist() for j in linesofdistributions[i]] for i in cons]
                sizes = [i.size for i in intensities]
                maxind = sizes.index(max(sizes))
                lineup = basemasses[maxind]
                masslines = defaultdict(list)
                intensitylines = defaultdict(list)
                retentionboundaries = defaultdict(list)
                lineinds = defaultdict(list)
                for n, (sm, sints, dlims) in enumerate(zip(basemasses, intensities, dlines)):
                    sdiff = np.abs(lineup - sm[:,None])
                    alignmentloc = np.argwhere(sdiff == sdiff.min())[0]
                    alignmentloc -= alignmentloc.min()
                    luind = alignmentloc[1] - alignmentloc[0]
                    outinds = [luind, luind + sm.size]
                    #below allows for mean diffs to be made on unequal sized arrays
                    for ind, mass, sint, dlim in zip(range(*outinds), sm.tolist(), sints, dlims):
                        masslines[ind].append(mass)
                        intensitylines[ind].append(sint)
                        retentionboundaries[ind].append(dlim)
                        lineinds[ind].append(n)
                #all isotopomers of a basemass must overlap -> the whole group must overlap with everything in that group, and > 1/2 of the groups must fulfill this
                #negative values [that form via non-overlaps] will get removed from mainoverlaps below
                matchcount = 0
                passedmatches = 0
                for k, v in retentionboundaries.items():
                    if len(v) > 1:
                        matchcount += 1
                        v = np.array(v)
                        if np.logical_and(v[:,0,None] <= v[:,1], v[:,1,None] >= v[:,0]).all():
                            passedmatches += 1
                if passedmatches > matchcount / 2:
                    massmeans = []
                    intensitylocations = defaultdict(list)
                    rtoverlaps = []
                    for n, inds in lineinds.items():
                        if len(inds) > 1:
                            massline = np.array(masslines[n])
                            massmeans.extend(np.abs(massline.mean() - massline).tolist())
                            intensityline = np.array(intensitylines[n])
                            intensityperc = intensityline / intensityline.sum()
                            for p, i in zip(intensityperc.tolist(), inds):
                                intensitylocations[i].append(p)
                            rbounds = np.array(retentionboundaries[n])
                            #not doing the fully vectorized approach so that selfoverlaps can be excluded while not infringing upon any actual 100% overlaps, although i could probably just take out a diagonal somewhere
                            for ni in range(len(rbounds)-1):
                                mainrb = rbounds[ni].copy()
                                nonmainrbs = rbounds[ni+1:].copy()
                                leftinners = nonmainrbs[:,0].copy()
                                leftinners[leftinners < mainrb[0]] = mainrb[0]
                                rightinners = nonmainrbs[:,1].copy()
                                rightinners[rightinners > mainrb[1]] = mainrb[1]
                                leftouters = nonmainrbs[:,0].copy()
                                leftouters[leftouters > mainrb[0]] = mainrb[0]
                                rightouters = nonmainrbs[:,1].copy()
                                rightouters[rightouters < mainrb[1]] = mainrb[1]
                                mainoverlaps = (rightinners - leftinners) / (rightouters - leftouters)
                                mainoverlaps = mainoverlaps[mainoverlaps > 0]
                                if mainoverlaps.size > 0:
                                    rtoverlaps.extend(mainoverlaps.tolist())
                    intensitymeans = []
                    for k, v in intensitylocations.items():
                        v = np.array(v)
                        intensitymeans.extend(np.abs(v.mean() - v).tolist())
                    massmeandiff = np.mean(massmeans)
                    intensitymeandiff = np.mean(intensitymeans)
                    generalizedoverlap = 1 / np.mean(rtoverlaps)
                    for cn in range(len(cons)-1):
                        conpair = cons[cn:cn+2]
                        conpairs[cid] = conpair
                        for con in conpair:
                            chargecongroups[con][cid] = [intensitymeandiff, massmeandiff, generalizedoverlap]
                        cid += 1
                    chargeconsets.add(tuple(sorted(cons)))
    
    print(time() - t2, 'charge match scoring')
    t3 = time()
    
    #score balancing
    prioritycharges = []
    secondpriorities = []
    for con, congroups in chargecongroups.items():
        if len(congroups) > 1:
            csums = np.array(list(congroups.values())).sum(axis=0)
            #0s here would actually be pretty welcome
            csums[csums == 0] = 1
            n1, n2, n3 = csums
            for congroup, (s1, s2, s3) in congroups.items():
                newconscore = sum((s1/n1, s2/n2, s3/n3))
                prioritycharges.append([congroup, newconscore])
        else:
            congroup, conscore = list(congroups.items())[0]
            conscore = conscore[:-1] #going to exclude the rt overlap here as there's no competition I suppose
            secondpriorities.append([congroup, sum(conscore)])
    
    rankedcharges = sorted(prioritycharges, key=lambda x: x[1]) #competition among matched chargegroups where matching can be done easily
    secondprioritycharges = sorted(secondpriorities, key=lambda x: x[1]) #no competition happening here, the order doesn't actually matter
    rankedcharges.extend(secondprioritycharges)
    
    blocked = set()
    preservedchargeranks = []
    for group, score in rankedcharges:
        if group not in blocked:
            preservedchargeranks.append(conpairs[group])
            blocked.add(group)
    
    print(time() - t3, 'charge priority ranking')
    t4 = time()
    
    sn = 0
    distributionsbychargegroup = defaultdict(set) #analyteid: [dist ids]
    groupsofdists = {} #dist: analyteid
    for dists in preservedchargeranks:
        locs = set()
        for i in dists:
            if i in groupsofdists:
                locs.add(groupsofdists[i])
        if locs:
            joiner = min(locs)
            if len(locs) > 1:
                tempchargegroup = set(dists)
                for l in locs:
                    tempchargegroup.update(distributionsbychargegroup[l])
                if tuple(sorted(tempchargegroup)) in chargeconsets:
                    for oldlocs in locs.difference([joiner]):
                        for ol in distributionsbychargegroup[oldlocs]:
                            groupsofdists[ol] = joiner
                        distributionsbychargegroup[joiner].update(distributionsbychargegroup.pop(oldlocs))
                    distributionsbychargegroup[joiner].update(dists)
                    for i in dists:
                        groupsofdists[i] = joiner
            else:
                if tuple(sorted(distributionsbychargegroup[joiner].union(dists))) in chargeconsets:
                    distributionsbychargegroup[joiner].update(dists)
                    for i in dists:
                        groupsofdists[i] = joiner
        else:
            joiner = sn
            sn += 1
            distributionsbychargegroup[joiner].update(dists)
            for i in dists:
                groupsofdists[i] = joiner
    
    chargedistgroups = defaultdict(dict) #analyteid: charge: distributionid
    chargegroupsbyline = {} #line: analyteid, doubles as blocking list
    for groupid, dists in distributionsbychargegroup.items():
        for dist in dists:
            charge = distributioncharges[dist]
            chargedistgroups[groupid][charge] = dist
            distlines = linesofdistributions[dist]
            #chargedistlines[groupid][charge] = distlines
            for line in distlines:
                chargegroupsbyline[line] = groupid
    
    print(time() - t4, 'charge group assembly')
    t5 = time()
    
    additionaldistributions = set(chargegroupsbyline).intersection(nodists) #nodists that made it into distributions will be given their own distribution
    distributionchangesetup = {k:distributionsoflines[k] for k in additionaldistributions}
    
    #removing all previous nodists from any distribution dicts
    #finaldefiniteind = max(max(i.keys()) for i in solodists.values()))
    indcopy = finaldefiniteind
    while indcopy in distributionmasses:
        del distributionmasses[indcopy]
        del distributiontimelimits[indcopy]
        del distributionintensities[indcopy]
        distributionsbycharge[distributioncharges[indcopy]].remove(indcopy)
        for line in linesofdistributions[indcopy]:
            del distributionsoflines[line]
        del linesofdistributions[indcopy]
        indcopy += 1
    
    distributionchanges = {}
    for line in additionaldistributions:
        dreg = regions[line]
        cgroup = chargegroupsbyline[line]
        potentialcharges = distchargesofnodists[line]
        chargekey = distributionsbychargegroup[cgroup].intersection(potentialcharges)
        if len(chargekey) > 1:
            print('big error -', line, 'has multiple potential nodist keys')
        ckey = list(chargekey)[0]
        charge = distchargesofnodists[line][ckey]
        dmass = dreg[7]
        dintensity = dreg[5]
        rtlimit = dreg[2:4]
        minrt = rtlimit.min()
        maxrt = rtlimit.max()
        distributionmasses[finaldefiniteind] = np.array([dmass])
        distributioncharges[finaldefiniteind] = charge
        distributionsoflines[line] = finaldefiniteind
        linesofdistributions[finaldefiniteind] = [line]
        distributiontimelimits[finaldefiniteind] = [minrt, maxrt]
        distributionintensities[finaldefiniteind] = np.array([dintensity])
        #distributionsbycharge[charge][dist] = sortedlines
        distributionsbycharge[charge].append(finaldefiniteind)
        oldkey = distributionchangesetup[line]
        distributionchanges[oldkey] = finaldefiniteind
        solodists[charge][finaldefiniteind] = [line]
        finaldefiniteind += 1 #continuing from solodists count

    #adding the rest of nodists anyways, i added this separate from the rest of the design, it should be fine
    for line in nodists:
        if line not in additionaldistributions:
            dreg = regions[line]
            charge = 2 #assuming
            dmass = dreg[7]
            dintensity = dreg[5]
            rtlimit = dreg[2:4]
            minrt = rtlimit.min()
            maxrt = rtlimit.max()
            distributionmasses[finaldefiniteind] = np.array([dmass])
            distributioncharges[finaldefiniteind] = charge
            distributionsoflines[line] = finaldefiniteind
            linesofdistributions[finaldefiniteind] = [line]
            distributiontimelimits[finaldefiniteind] = [minrt, maxrt]
            distributionintensities[finaldefiniteind] = np.array([dintensity])
            distributionsbycharge[charge].append(finaldefiniteind)
            solodists[charge][finaldefiniteind] = [line]
            finaldefiniteind += 1 #continuing from solodists count
    
    distributionregions = []
    for k, v in solodists.items():
        for sk, sv in v.items():
            masses = distributionmasses[sk]
            massmax = masses.max()
            massmin = masses.min()
            intensities = distributionintensities[sk]
            mainmass = masses[intensities.argmax()]
            mintime, maxtime = distributiontimelimits[sk]
            signalsum = defaultdict(float) #time: total intensity
            for line in sv:
                data = trackedgroups[line]
                for m, t, i in data.tolist():
                    signalsum[t] += i
            signals = np.array(list(signalsum.items()))
            area = np.trapezoid(signals[:,1], signals[:,0])
            el = [massmin, massmax, mintime, maxtime, len(sv), area, k, mainmass, sk]
            distributionregions.append(el)
    
    distributionregions = np.array(distributionregions)
    distributionregions = distributionregions[distributionregions[:,8].argsort()]
    
    print(time() - t5, 'finalized distribution information')
    t6 = time()
    
    chargestatelines = {} #needs confirmation -> lineid: chargeid
    chargeregions = []
    for k, v in chargedistgroups.items():
        mincharge = min(v)
        maxcharge = max(v)
        times = []
        signalsum = defaultdict(float) #time: total intensity
        changers = {}
        for vc, sv in v.items():
            if sv in distributionchanges:
                changers[vc] = sv
        for vc, sv in changers.items():
            del chargedistgroups[k][vc]
            sv = distributionchanges[sv]
            vc = distributioncharges[sv]
            chargedistgroups[k][vc] = sv
        for vc, sv in v.items():
            times.extend(distributiontimelimits[sv])
            for line in linesofdistributions[sv]:
                data = trackedgroups[line]
                for m, t, i in data.tolist():
                    signalsum[t] += i
                chargestatelines[line] = vc
        signals = np.array(list(signalsum.items()))
        sa = signals[:,0].argsort()
        signals = signals[sa]
        area = np.trapezoid(signals[:,1], signals[:,0])
        distinds = np.array(list(v.values()))
        charges = np.array(list(v))
        maincharge = charges[distributionregions[distinds,5].argmax()]
        mintime = min(times)
        maxtime = max(times)
        el = [mincharge, maxcharge, mintime, maxtime, len(v), area, maincharge, k]
        chargeregions.append(el)
    chargeregions = np.array(chargeregions)
    
    print(time() - t6, 'finalized charge state information')
    print(len(distributionregions), 'distributions')
    print(len(chargeregions), 'ions with multiple charge states')
    t17 = time()
    
    #analyte id == chargegroup id if there is one, continue using sn for this count
    analytekeys = defaultdict(dict) #analyte id: distid: charge, this can be one or multiple distributions based on whether it has multiple charge states
    analytedistributions = defaultdict(dict) #analyte id: ordered masses: AUC of merged isotopomers, weighted means [via intensity] across isotopomers from every charge state, isotopomer datapoints merged across each charge state - if there are any
    analytesbydistribution = {} #distid: analyte id
    linesofanalytes = {} #analyteid: [[lines across charge states at this position], [...]]
    chargesoflines = {} #line: distribution charge
    blocked = set()
    for k, v in chargedistgroups.items():
        massholder = []
        lineholder = []
        for charge, distid in v.items():
            analytekeys[k][distid] = charge
            dmasses = distributionmasses[distid] * charge - proton * charge
            dlines = linesofdistributions[distid]
            for line in dlines:
                chargesoflines[line] = charge
            massholder.append(dmasses)
            lineholder.append(dlines)
            analytesbydistribution[distid] = k
            blocked.add(distid)
        
        sortedlines = sorted(lineholder, key=lambda x: len(x))
        sortedmasses = sorted(massholder, key=lambda x: x.size)
        lineup = sortedmasses[-1][:,None]
        distinds = []
        for sm in sortedmasses:
            sdiff = np.abs(lineup - sm)
            alignmentloc = np.argwhere(sdiff == sdiff.min())[0]
            alignmentloc -= alignmentloc.min()
            luind = alignmentloc[0] - alignmentloc[1]
            outinds = [luind, luind + sm.size]
            distinds.append(outinds)
        
        spaceorganizer = defaultdict(list) #position: [lineuids]
        spacemasses = defaultdict(list) #position: [basemasses]
        for di, sl, sm in zip(distinds, sortedlines, sortedmasses):
            dspaces = range(*di)
            for space, line, mass in zip(dspaces, sl, sm):
                spaceorganizer[space].append(line)
                spacemasses[space].append(mass)
        
        linesofanalytes[k] = []
        massesandintensities = [[], []]
        for sk in sorted(spaceorganizer): #why does this need to be sorted? its putting massesandintensities in order
            slines = tuple(spaceorganizer[sk])
            smasses = np.array(spacemasses[sk])
            areavals = regions[slines,5]
            weightedmass = (smasses * areavals).sum() / areavals.sum()
            signalsum = defaultdict(float) #time: sum intensity
            for line in slines:
                linegroup = trackedgroups[line]
                for m, t, i in linegroup.tolist():
                    signalsum[t] += i
            signals = np.array(list(signalsum.items()))
            sa = signals[:,0].argsort()
            signals = signals[sa]
            area = np.trapezoid(signals[:,1], signals[:,0])
            massesandintensities[0].append(weightedmass)
            massesandintensities[1].append(area)
            linesofanalytes[k].append(slines)
        analytedistributions[k] = np.array(massesandintensities)
    
    #adding any distribution without multiple charge states to analytedistributions
    for distid, masses in distributionmasses.items():
        if distid not in blocked:
            charge = distributioncharges[distid]
            basemasses = masses * charge - proton * charge
            intensities = distributionintensities[distid]
            analytekeys[sn][distid] = charge
            analytesbydistribution[distid] = sn
            massesandintensities = np.array([basemasses, intensities])
            analytedistributions[sn] = massesandintensities
            dlines = linesofdistributions[distid]
            for line in dlines:
                chargesoflines[line] = charge
            linesofanalytes[sn] = [tuple([i]) for i in dlines]
            sn += 1
    
    print(time() - t17, 'summarized analyte distributions')
    
    saverloc = ''.join((processingdirectory, 'distributionregions.pickle'))
    with open(saverloc, 'wb') as pick:
        pickle.dump(distributionregions, pick)
    
    saverloc = ''.join((processingdirectory, 'chargeregions.pickle'))
    with open(saverloc, 'wb') as pick:
        pickle.dump(chargeregions, pick)
    
    saverloc = ''.join((processingdirectory, 'solodists.pickle'))
    with open(saverloc, 'wb') as pick:
        pickle.dump(solodists, pick)
    
    saverloc = ''.join((processingdirectory, 'chargedistgroups.pickle'))
    with open(saverloc, 'wb') as pick:
        pickle.dump(chargedistgroups, pick)
    
    saverloc = ''.join((processingdirectory, 'chargesoflines.pickle'))
    with open(saverloc, 'wb') as pick:
        pickle.dump(chargesoflines, pick)

    saverloc = ''.join((processingdirectory, 'distributioncharges.pickle'))
    with open(saverloc, 'wb') as pick:
        pickle.dump(distributioncharges, pick)

    saverloc = ''.join((processingdirectory, 'analytefactors.pickle'))
    savedbits = [analytekeys, analytedistributions, analytesbydistribution, distributionsoflines, linesofdistributions, linesofanalytes]
    with open(saverloc, 'wb') as pick:
        pickle.dump(savedbits, pick)
