from elementalcomponents import proton
from generalfunctions import boundary_stack
from database import environment

from itertools import chain, combinations
from collections import defaultdict
from time import time
import numpy as np
import pickle

def overlap_counts(rbounds):
    finalscore = 0
    for (lmin, lmax), (rmin, rmax) in combinations(rbounds, 2):
        #first check if they even overlap, if not -> set up for a subtraction
        #encompassment -> flanks subtracted from each other -> this difference will be subtracted from the total overlap -> final value
        #double flank, both has a flank -> both flanks added together -> subtracted from overlap -> final value
        #single flank -> flank is subtracted from overlap -> final value
        #non-overlap -> minus the full range of both
        if rmax > lmin and lmax > rmin:
            #overlap exists
            if rmax <= lmax:
                minmax = rmax
                maxmax = lmax
                maxl = True
            else:
                minmax = lmax
                maxmax = rmax
                maxl = False
            if rmin >= lmin:
                maxmin = rmin
                minmin = lmin
                minl = True
            else:
                maxmin = lmin
                minmin = rmin
                minl = False
            if maxl and minl:
                if rmax == lmax and rmin == lmin:
                    #equal boundaries
                    overlap = rmax - rmin
                    fullrange = lmax - lmin
                    finalscore += overlap
                elif rmax == lmax:
                    #min l flank
                    overlap = lmax - rmin
                    fullrange = rmax - lmin
                    minflank = rmin - lmin
                    score = overlap - minflank
                    finalscore += score
                elif rmin == lmin:
                    #max l flank
                    overlap = rmax - rmin
                    fullrange = lmax - lmin
                    maxflank = lmax - rmax
                    score = overlap - maxflank
                    finalscore += score
                else:
                    #l encompasses r
                    overlap = rmax - rmin
                    fullrange = lmax - lmin
                    minflank = rmin - lmin
                    maxflank = lmax - rmax
                    flankdiff = abs(maxflank - minflank)
                    score = overlap - flankdiff
                    finalscore += score
            elif maxl:
                if rmax == lmax:
                    #r min flank
                    overlap = rmax - lmin
                    fullrange = lmax - rmin
                    minflank = lmin - rmin
                    score = overlap - minflank
                    finalscore += score
                else:
                    #l max flank
                    #r min flank
                    overlap = rmax - lmin
                    fullrange = lmax - rmin
                    minflank = lmin - rmin
                    maxflank = lmax - rmax
                    flankdiff = abs(maxflank - minflank)
                    score = overlap - flankdiff
                    finalscore += score
            elif minl:
                if rmin == lmin:
                    #r max flank
                    overlap = lmax - rmin
                    fullrange = rmax - lmin
                    maxflank = rmax - lmax
                    score = overlap - maxflank
                    finalscore += score
                else:
                    #r max flank
                    #l min flank
                    overlap = lmax - rmin
                    fullrange = rmax - minl
                    minflank = rmin - lmin
                    maxflank = rmax - lmax
                    flankdiff = abs(maxflank - minflank)
                    score = overlap - flankdiff
                    finalscore += score
            else:
                #r encompasses l
                overlap = lmax - lmin
                fullrange = rmax - rmin
                maxflank = rmax - lmax
                minflank = rmin - lmin
                flankdiff = abs(maxflank - minflank)
                score = overlap - flankdiff
                finalscore += score
        else:
            #there's no overlap
            loverlap = lmax - lmin
            roverlap = rmax - rmin
            #subtract both overlaps
            finalscore -= loverlap + roverlap
    return finalscore

def distribution_assembly(minpoints, chargetolerance, librarylocation, processingdirectory, proteome):
    
    regionfile = ''.join((processingdirectory, 'regions.pickle'))
    with open(regionfile, 'rb') as pick:
        regions = pickle.load(pick)
    
    roundcutfile = ''.join((processingdirectory, 'roundcutoff.pickle'))
    with open(roundcutfile, 'rb') as pick:
        roundcutoff = pickle.load(pick)
    
    with environment(librarylocation) as env:
        proteomedb = env.open_db((proteome + '.info').encode())
        with env.begin(write=False) as txn:
            with txn.cursor(proteomedb) as cursor:
                uppermasslimit = float(cursor.get('uppermasslimit'.encode()).decode())
    
    subisomax = 0.01337851739
    newinclimit = 0.1
    steplimit = 0.5
    subisomax = subisomax + subisomax * chargetolerance
    
    t1 = time()
    
    rsorted = regions[regions[:,7].argsort()]
    regiter = rsorted[rsorted[:,4] >= minpoints]
    
    pairkeys = {} #pairkey: pair
    previousdecrease = {} #connectionindex: True, exists if something is increasing
    
    di = 0 #connectionindexs, keeping this incorrect plural makes it searchable
    paircharges = {} #connection: charge
    scoresbypair = {} #pair: [[absacdiff, datapercdiff, rtoffset],]
    pairsbyline = defaultdict(str) #line: 'pair,pair,'
    
    masswidthlimit = roundcutoff * 2
    
    connectionspine = defaultdict(str) #connectionindex: 'pairkey,pairkey,'
    latestconnections = defaultdict(lambda: defaultdict(str)) #masskey: charge: [current connection indices that end in masskey], a mass is always active as long as it's in masspool, this keeps track of all the groups that end in any given mass
    latestmass = {} #connectionindex: latest masskey
    
    pi = 0 #pairkeys
    masspool = []
    for reg in regiter.tolist():
        npoints = reg[4]
        nm = reg[7]
        nmlt = reg[2]
        nmrt = reg[3]
        nmassmin = reg[0]
        nmassmax = reg[1]
        nmwidth = nmassmax - nmassmin
        if nmwidth > masswidthlimit:
            #this is strict, but seems to help in the end, meaning that essentially the means are rather important and useful
            nmassmin = nm - masswidthlimit
            nmassmax = nm + masswidthlimit
        nkey = int(reg[8])
        nintensity = reg[5]
        masspoolremovals = []
        nlrange = nmrt - nmlt
        for okey in masspool:
            oreg = regions[okey]
            om = oreg[7]
            omassmin = oreg[0]
            omassmax = oreg[1]
            omwidth = omassmax - omassmin
            if omwidth > masswidthlimit:
                omassmin = om - masswidthlimit
                omassmax = om + masswidthlimit
            if nmassmin - omassmax <= proton:
                omlt = oreg[2]
                omrt = oreg[3]
                olrange = omrt - omlt
                overpass = False
                if nmlt < omrt and nmrt > omlt: #rt's overlap
                    #if you want to allow non-overlaps to match, and there are good reasons to, you need percentoverlap to become negative so that its later addition is valued worse
                    overlap = min(omrt, nmrt) - max(omlt, nmlt)
                    if nmlt > omlt and nmrt < omrt: #new encompassed
                        combinedrange = nlrange + olrange
                        percentoverlap = (overlap * 2) / combinedrange
                        #using newinclimit here is appropriate because it's directly related to expected differences in adjacent isotopomer quantities, which is what the overlap is also ~somewhat assessing
                        if percentoverlap > newinclimit:
                            overpass = True
                    elif omlt > nmlt and omrt < nmrt: #old encompassed
                        combinedrange = nlrange + olrange
                        percentoverlap = (overlap * 2) / combinedrange
                        if percentoverlap > newinclimit:
                            overpass = True
                    else:
                        fullrange = max(omrt, nmrt) - min(omlt, nmlt)
                        percentoverlap = overlap / fullrange
                        if percentoverlap > 0.5: #this is super lenient I think
                            overpass = True
                #else:
                    #negative percentoverlap values for non-overlaps within some 1.75 range or something goes here
                if overpass:
                    diff = nm - om
                    #sometimes a charge registers as one thing but the value is way closer to +/-1
                    initialcharge = round(proton / diff)
                    if initialcharge > 1: #nothing close enough to care about for 1, plus the first bit would go to 0 and cause annoying zerodivision warnings
                        chargespread = np.linspace(initialcharge - 1, initialcharge + 1, 3)
                        expspread = proton / chargespread
                        minexpind = np.abs(diff - expspread).argmin()
                        charge = int(chargespread[minexpind])
                    else:
                        charge = initialcharge
                    actualmass = om * charge - proton * charge
                    #min/max peptide mass from the isotope database, anything with a basemmass larger than that max shouldn't be considered, it would help keep out ridiculously large charges - instead having them be a part of nodists, the way they probably should
                    if actualmass <= uppermasslimit:
                        expdiff = proton / charge
                        acdiff = expdiff - diff
                        diffcut = expdiff * chargetolerance
                        if acdiff > -1 * (diffcut * chargetolerance + masswidthlimit):
                        #^a lighter acceptance for distances above the charge dist - it does happen from time to time on lower-hanging straggler isotopomers
                            if acdiff <= diffcut + masswidthlimit:
                                absacdiff = abs(acdiff) * charge #normalizing -> distance to proton
                                ointensity = oreg[5]
                                ncons = 0
                                intensitypercdiff = abs(nintensity - ointensity) / (nintensity + ointensity) / 2
                                subcheck = False
                                if okey in latestconnections:
                                    if charge in latestconnections[okey]:
                                        subcheck = True
                                decreasecheck = nintensity < ointensity
                                if subcheck:
                                    csubs = map(int, latestconnections[okey][charge].split(',')[:-1])
                                    for adi in csubs:
                                        ratiocheck = steplimit
                                        if adi in previousdecrease:
                                            if not decreasecheck:
                                                #intensity is increasing
                                                ratiocheck = newinclimit
                                        if intensitypercdiff <= ratiocheck:
                                            if latestmass[adi] == okey:
                                                xdi = adi
                                            else:
                                                xdi = di
                                                opairs = pairsbyline[okey].split(',')[:-1]
                                                spinecopy = connectionspine[adi].split(',')[:-1]
                                                spinds = []
                                                for op in opairs:
                                                    if op in spinecopy:
                                                        spinds.append(spinecopy.index(op))
                                                spi = max(spinds)
                                                spinecopy = spinecopy[:spi]
                                                connectionspine[xdi] = ','.join((spinecopy)) + ','
                                                di += 1
                                            connectionspine[xdi] += str(pi) + ','
                                            latestmass[xdi] = nkey
                                            latestconnections[nkey][charge] += str(xdi) + ','
                                            if decreasecheck:
                                                previousdecrease[xdi] = True
                                            ncons += 1
                                else:
                                    #no previous subgroup 
                                    if intensitypercdiff <= steplimit:
                                        connectionspine[di] += str(pi) + ','
                                        latestmass[di] = nkey
                                        latestconnections[nkey][charge] += str(di) + ','
                                        if decreasecheck:
                                            previousdecrease[di] = True
                                        ncons += 1
                                        di += 1
                                if ncons > 0:
                                    opoints = oreg[4]
                                    #dpercdiff = abs(npoints - opoints) / (npoints + opoints) / 2
                                    dpercdiff = abs(nintensity - ointensity) / (nintensity + ointensity) / 2
                                    lpair = (okey, nkey)
                                    pairkeys[pi] = lpair
                                    paircharges[pi] = charge
                                    scorelist = [absacdiff, dpercdiff, percentoverlap, ~decreasecheck]
                                    for p in lpair:
                                        pairsbyline[p] += str(pi) + ','
                                    scoresbypair[pi] = scorelist
                                    if subcheck:
                                        connectionspine[di] += str(pi) + ','
                                        latestmass[di] = nkey
                                        latestconnections[nkey][charge] += str(di) + ','
                                        if decreasecheck:
                                            previousdecrease[di] = True
                                        di += 1
                                    pi += 1
            else:
                #om is past proton distance, remove om from mass pool
                masspoolremovals.append(okey)
        for mpr in masspoolremovals:
            masspool.remove(mpr)
            if mpr in latestconnections:
                for charge, cons in latestconnections[mpr].items():
                    for con in map(int, cons.split(',')[:-1]):
                        if mpr == latestmass[con]:
                            del latestmass[con]
                            try:
                                del previousdecrease[con]
                            except KeyError: #no previous decrease, idc really
                                pass
                del latestconnections[mpr]
        masspool.append(nkey)
    
    del latestconnections
    del previousdecrease
    del latestmass
    
    print(time() - t1, 'mass processing')
    print('connectionspine -', len(connectionspine))
    t2 = time()
    
    flatdistgroups = set()

    distributionscoresbyline = defaultdict(list) #linekey: [pairkeys] -> without a set there is redundancy in here that disrupts downstream
    distributionscoredict = defaultdict(list) #pairkey: [[scores],]
    for distkey, pairspine in connectionspine.items():
        pairspine = map(int, pairspine.split(',')[:-1])
        activepairlist = []
        scores = []
        for pairs in pairspine:
            activepairlist.append(pairs)
            scores.append(scoresbypair[pairs])
            flatdist = []
            competinglines = set()
            for pair in activepairlist:
                pk = pairkeys[pair]
                flatdist.extend(pk)
                for p in pk:
                    if set(map(int, pairsbyline[p].split(',')[:-1])).difference(activepairlist): #checking that the line is in more distributions than whats generated via this spine
                        competinglines.add(p)
            rbounds = regions[flatdist,2:4]
            rstack = overlap_counts(rbounds.tolist())
            if rstack > 0:
                flatdistgroups.add(tuple(sorted(set(flatdist))))
                if competinglines and len(activepairlist) > 1:
                    #unstackedboundaries = rbounds - rbounds.min(axis=1)[:,None]
                    #ustack = boundary_stack(unstackedboundaries)
                    scorearray = np.array(scores)
                    distmean = scorearray[:,0].mean()
                    #rtmultiplier = rstack / ustack
                    rtmultiplier = rstack
                    decreasingmultiplier = scorearray[:,3].sum() + 1
                    slen = len(scorearray)
                    for pair, score in zip(activepairlist, scores): 
                        if paircharges[pair] > 1: # a lot of bad 1+ matches get high priority from this, I essentially want less 1+ than 3+ and this helps
                            dist, ddiff, rtoffset, decs = score
                            meandiff = abs(distmean - dist) / slen
                            distdiff = meandiff * (2**decreasingmultiplier)
                            datadiff = ddiff / rtmultiplier
                            scorelist = tuple([distdiff, datadiff])
                            distributionscoredict[pair].append(scorelist)
                            for p in pairkeys[pair]:
                                if p in competinglines:
                                    if not pair in distributionscoresbyline[p]: #avoiding set use to save memory
                                        distributionscoresbyline[p].append(pair)
    
    print(time() - t2, 'distribution scoring')
    print('flatdistgroups', len(flatdistgroups))
    t3  = time()
    
    preservedpairs = set()
    #lines being comparison-ranked across all individual distribution scores it participates in, if it has multiple
    rankedpairs = [] #[pair, minval]
    for line, pairs in distributionscoresbyline.items():
        vals = []
        pairexpansion = []
        for pair in pairs:
            pairscores = distributionscoredict[pair]
            vals.extend(pairscores)
            pairexpansion.extend([pair for _ in range(len(pairscores))])
        vals = np.array(vals)
        vsums = vals.sum(axis=0)
        vsums[vsums == 0] = 1 #0-sums make nans but this doesn't change the final answer of 0
        vpercs = vals / vsums
        sumvals = vpercs.sum(axis=1)
        rankedpairs.extend(list(zip(pairexpansion, sumvals.tolist())))
        preservedpairs.update(pairs)

    secondpriorities = [] #[pair, score]
    thirdpriorities = [] #[pair, score]
    for line, pairs in pairsbyline.items():
        pairs = list(map(int, pairs.split(',')[:-1]))
        plen = len(pairs)
        if plen > 1:
            scores = [scoresbypair[i] for i in pairs]
            scorearray = np.array(scores)
            scorearray[:,3] += 1
            scoresums = scorearray.sum(axis=0)
            scoresums[scoresums == 0] = 1 #0-sums make nans but this doesn't change the final answer of 0
            normscores = scorearray / scoresums
            offsetnorms = normscores[:,:2] * scorearray[:,3,None] / scorearray[:,2,None]
            for pair, score in zip(pairs, offsetnorms.sum(axis=1).tolist()):
                if pair not in preservedpairs: #memory saving
                    secondpriorities.append([pair, score])
        else:
            pairs = pairs[0]
            if pairs not in preservedpairs:
                dist, ddiff, rtoffset, dec = scoresbypair[pairs]
                dec += 1 #adds 1 to things that dec'd and makes non-dec's 1 -> no change
                equalizednorm = abs(ddiff - dist)
                outscore = equalizednorm * dec / rtoffset
                thirdpriorities.append([pairs, outscore])

    firstranks = sorted(rankedpairs, key=lambda x: x[1])
    secondranks = sorted(secondpriorities, key=lambda x: x[1])
    thirdranks = sorted(thirdpriorities, key=lambda x: x[1])

    #secondrankindices = {i[0]: i[1] for i in secondranks}

    #the other option is to remove the zeros as they're unreliable and probably not as important due to the fact the rt overlap being 1 probably means they're small lines
    #forcedcompetition = []
    #for p, s in firstranks:
    #    if s == 0:
    #        forcedcompetition.append([p, secondrankindices[p]])
    #    else:
    #        break
    if firstranks[0][1] == 0:
        print('need forced competition')

    #sortedranks = sorted(forcedcompetition, key=lambda x: x[1])
    
    sortedranks = []
    sortedranks.extend(firstranks)
    sortedranks.extend(secondranks)
    sortedranks.extend(thirdranks)

    preservedpairs = set()
    preservedranks = []
    for pairkey, score in sortedranks:
        if pairkey not in preservedpairs:
            preservedranks.append([pairkeys[pairkey], paircharges[pairkey]])
            preservedpairs.add(pairkey)

    print(time() - t3, 'pair ranking setup')
    #print('forcedcompetition', len(forcedcompetition))
    print('firstranks', len(firstranks))
    print('secondranks', len(secondranks))
    print('thirdranks', len(thirdranks))
    print('preservedranks', len(preservedranks))
    t4 = time()
    
    sn = 0
    distsets = defaultdict(set) #distloc: [distributions]
    linelocations = {} #masskey: distloc
    setcharges = defaultdict(set) #index of distsets: [charges]
    for pair, paircharge in preservedranks:
        pairset = set(pair)
        locs = set()
        for line in pairset:
            if line in linelocations:
                locs.add(linelocations[line])
        if locs:
            joiner = min(locs)
            if len(locs) > 1:
                dist = pairset.copy()
                for l in locs:
                    dist.update(distsets[l])
                if tuple(sorted(dist)) in flatdistgroups:
                    for l in locs.difference([joiner]):
                        for line in distsets[l]:
                            linelocations[line] = joiner
                        setcharges[joiner].update(setcharges.pop(l))
                        distsets[joiner].update(distsets.pop(l))
                    distsets[joiner].update(pairset)
                    for line in pairset:
                        linelocations[line] = joiner
            else:
                if tuple(sorted(distsets[joiner].union(pairset))) in flatdistgroups:
                    distsets[joiner].update(pairset)
                    setcharges[joiner].add(paircharge)
                    for line in pairset:
                        linelocations[line] = joiner
        else:
            joiner = sn
            distsets[joiner].update(pairset)
            setcharges[joiner].add(paircharge)
            for line in pairset:
                linelocations[line] = joiner
            sn += 1

    print(time() - t4, 'pair ranking')
    t5 = time()
    
    dr =  0 #easier to work with dr here rather than re-index the finalized distribution just to get the max index for finaldefiniteind afterwards, the order of these aren't important
    solodists = defaultdict(dict) #charge: distid: lines
    for distindex, dist in distsets.items():
        charges = setcharges[distindex]
        if dist:
            charge = max(charges)
            solodists[charge][dr] = list(dist)
            dr += 1
    finaldefiniteind = dr
    
    foundvals = []
    for charge, sgd in solodists.items():
        foundvals.extend(list(chain(*sgd.values())))
    specvals = regiter[:,8].astype(int)
    nodists = np.setdiff1d(specvals, foundvals)
    
    print(time() - t5, 'distribution assembling')
    print(sum(len(i) for i in solodists.values()), 'initial distributions')
    print(len(nodists), 'nodists')
    
    nodistfile = ''.join((processingdirectory, 'nodists.pickle'))
    with open(nodistfile, 'wb') as pick:
        pickle.dump(nodists, pick)
    
    solodistfile = ''.join((processingdirectory, 'solodists.pickle'))
    with open(solodistfile, 'wb') as pick:
        pickle.dump(solodists, pick)
    
    finalindfile = ''.join((processingdirectory, 'finalindex.pickle'))
    with open(finalindfile, 'wb') as pick:
        pickle.dump(finaldefiniteind, pick)
    
    #return finaldefiniteind, processingdirectory
