from generalfunctions import intersection_merge
from database import environment

from collections import defaultdict
from itertools import chain
from pyteomics import mzml
from time import time
import numpy as np
import pickle
import math
import sys
import os

def line_model(mzmlfile, minpoints, minmovinginds, deadsignal, chargetolerance, librarylocation, processingdirectory, proteome):
    
    mslevelfile = ''.join((processingdirectory, 'centroid.ms1.pickle'))
    with open(mslevelfile, 'rb') as pick:
        ms1scans = pickle.load(pick)

    retentiontimesbyscanfile = ''.join((processingdirectory, 'retentiontimesbyscan.pickle'))
    with open(retentiontimesbyscanfile, 'rb') as pick:
        retentiontimesbyscan = pickle.load(pick)
    
    with environment(librarylocation) as env:
        proteomedb = env.open_db((proteome + '.info').encode())
        with env.begin(write=False) as txn:
            with txn.cursor(proteomedb) as cursor:
                uppermasslimit = float(cursor.get('uppermasslimit'.encode()).decode())
    
    subisomax = 0.01337851739
    newinclimit = 0.1
    steplimit = 0.5
    print('subisomax', subisomax)
    print('newinclimit', newinclimit)
    print('steplimit', steplimit)
    print('uppermasslimit', uppermasslimit)
    
    subisomax = subisomax + subisomax * chargetolerance
    
    t1 = time()
    #msrun = mzml.MzML(mzmlfile, dtype=np.float64)
    msiter = iter(ms1scans)
    
    timearray = []
    
    #scan = next(msrun)
    scanindex = next(msiter)
    mza, intensities = ms1scans[scanindex].values()
    #mza = scan['m/z array']
    previousdata = mza.copy()
    scancount = 1
    
    #rt = scan['scanList']['scan'][0]['scan start time'].real
    rt = retentiontimesbyscan[scanindex]
    timearray.append(rt)
    
    mlen = mza.size
    #intensities = scan['intensity array'] #no IT normalization needed, pre-normalized by the instrument software apparently
    retentiontimes = np.repeat(rt, mlen)
    
    coords = np.stack((mza, retentiontimes, intensities), axis=1).reshape(mza.size, 1, 3).tolist()
    
    uids = np.arange(mlen).tolist()
    uidcount = max(uids) + 1
    
    trackedgroups = {} #uniqueid: [[masses], [rt-inds], [intensities], [intensities/injection times], [percent intensities of scans]]
    trackedma = {} #latest moving average mass of trackedgroup: lineuid
    linedeletioncounter = defaultdict(int) #lineuid: notmatched count
    groupmovingaverages = {} #lineuid: latest moving average of line
    groupdifftoma = {} #lineid: moving difference to moving average
    groupranges = {} #uniqueid: [minmass, maxmass]
    modeltracking = {} #scan: number of masses being [added, matched, nonmatched, removed]
    
    modeltracker = [0, 0, 0, 0]
    modeltracker[0] += mlen
    modeltracking[scanindex] = modeltracker
    
    flatmasslist = mza.flatten().tolist()
    trackedma.update(zip(flatmasslist, uids))
    trackedgroups.update(zip(uids, coords))
    groupmovingaverages.update(zip(uids, flatmasslist))
    elen = len(uids)
    groupdifftoma.update(zip(uids, np.zeros(elen).tolist()))
    groupranges.update(zip(uids, np.stack((mza, mza), axis=1).tolist()))
    
    modify = False
    widestmassrange = 0 #a tracked float of the widest mass range
    wides = []
    linecorrections = []
    
    roundcutoff = 0
    #for scan in msrun:
    for scanindex in msiter:
        mza, intensities = ms1scans[scanindex].values()
        #scanlist = scan['scanList']['scan'][0]
        #rt = scanlist['scan start time'].real
        rt = retentiontimesbyscan[scanindex]
        timearray.append(rt)
        
        trackedkeys = {} #latest mass in a trackedgroup: lineid
        
        #intensities = scan['intensity array']
        
        modeltracker = [0, 0, 0, 0]
        #mza = scan['m/z array']
        
        baseind = 0
        catches = []
        massdist = []
        #a k=1 nearest neighbors for signal-processing
        #this is iterating over numpy arrays because its slower to convert them to lists and slower to index a list
        #picking whatevers closer in intensity might start to fail as a concept if the ms1 scans are more spaced out, or perhaps boxcar'd
        for fn, f in enumerate(mza.tolist()):
            mindist = np.inf
            for n, b in enumerate(previousdata[baseind:]):
                dist = abs(b-f)
                if dist < mindist:
                    minind = n + baseind
                    mindist = dist
                elif dist == mindist:
                    #two new masses have symmetrical distances to existing moving average
                    #choose whichever is within the original lines range
                    currentind = trackedma[b]
                    linerange = groupranges[currentind]
                    othermass = previousdata[minind]
                    currentmatch = b > linerange[0] and b < linerange[1]
                    othermatch = othermass > linerange[0] and othermass < linerange[1]
                    if currentmatch and not othermatch:
                        #new match wins
                        minind = n + baseind
                        mindist = dist
                        #no other distances will be closer
                        break
                    elif othermatch and not currentmatch:
                        #old match wins
                        #no other distances will be closer
                        break
                    else:
                        #either both or neither are within the range
                        #switch from comparing masses to comparing intensities
                        currentintensity = trackedgroups[currentind][-1][2]
                        otherintensity = trackedgroups[minind][-1][2]
                        massintensity = intensities[fn]
                        cabs = abs(massintensity - currentintensity)
                        oabs = abs(massintensity - otherintensity)
                        if cabs < oabs:
                            #current mass wins out
                            minind = n + baseind
                            mindist = dist
                            #no other masses will be closer
                            break
                        else:
                            #other mass wins out
                            #no other masses will be closer
                            break
                else:
                    break
            catches.append(minind)
            massdist.append(mindist)
            baseind = minind
        massdist = np.array(massdist)
        catches = np.array(catches)
        
        found = previousdata[catches] #check for duplicates here
        uf, ufc = np.unique(found, return_counts=True)
        ub = ufc > 1
        redundants = np.any(ub)
        #finding redundant matches
        if redundants:
            removals = []
            for umatch in uf[ub].tolist():
                mwhere = np.where(found == umatch)[0]
                mwdists = massdist[mwhere]
                mwdargmin = mwdists.argmin()
                removals.extend(np.delete(mwhere, mwdargmin).tolist())
        
        #removing redundant matches -> the line ended for these ones, or its skipping an index
        retentiontimes = np.repeat(rt, mza.size)
        coords = np.stack((mza, retentiontimes, intensities), axis=1).reshape(mza.size, 1, 3)
        fmassdist = massdist.copy()
        if redundants:
            #things that had a redundant match, and weren't taken, from mza are put up as new lines
            ecoords = coords[removals]
            flatmasslist = ecoords[:,0,0]
            elen = len(ecoords)
            uids = np.arange(elen) + uidcount
            uidcount += elen
            trackedgroups.update(zip(uids, ecoords.tolist()))
            trackedkeys.update(zip(flatmasslist, uids))
            trackedma.update(zip(flatmasslist, uids))
            groupmovingaverages.update(zip(uids, flatmasslist))
            groupdifftoma.update(zip(uids, np.zeros(elen)))
            groupranges.update(zip(uids, np.stack((flatmasslist, flatmasslist), axis=1).tolist()))
            modeltracker[0] += elen #newly added
            #fixing the originals
            found = np.delete(found, removals)
            fmassdist = np.delete(fmassdist, removals)
            coords = np.delete(coords, removals, axis=0)
        found = found.flatten().tolist()
        
        sorteddistances = np.sort(massdist)
        mbool = np.arange(sorteddistances.size)[::-1] + 1
        countsums = mbool / sorteddistances.size
        sumcounts = sorteddistances.cumsum() / sorteddistances.sum()
        mincomboind = (countsums + sumcounts).argmin()
        mincombo = sorteddistances[mincomboind]
        #moving average of average of dists under mincombo
        explicitcutoff = sorteddistances[sorteddistances <= mincombo].mean()
        roundcutoff = (roundcutoff * scancount + explicitcutoff) / (scancount + 1)
        
        #modifying things that are already being tracked
        fmzaremovals = []
        foundremovals = []
        for c, f, d in zip(coords.tolist(), found, fmassdist):
            modify = False
            nf = c[0][0]
            tid = trackedma[f]
            tgroup = trackedgroups[tid]
            tlen = len(tgroup)
            lastmass = tgroup[-1][0]
            rmin, rmax = groupranges[tid]
            grange = rmax - rmin
            #this is for when the moving decision fails for something within the existing range, ain't no thang
            rangepass = nf <= rmax and nf >= rmin
            distancepass = abs(nf - lastmass) < grange / 2
            if rangepass or distancepass:
                oldma = groupmovingaverages[tid]
                nma = (oldma * tlen + nf) / (tlen + 1)
                nmadiff = abs(oldma - nma)
                groupmovingaverages[tid] = nma
                madiff = groupdifftoma[tid]
                groupdifftoma[tid] = (madiff * (tlen - 1) + nmadiff) / tlen
                modify = True
            #generally, this is good for on-the-fly decision making when the moving target is outside the existing mass range. This dominates later on, where it's more robust
            elif tlen >= minmovinginds:
                oldma = groupmovingaverages[tid]
                madiff = groupdifftoma[tid]
                nma = (oldma * tlen + nf) / (tlen + 1)
                nmadiff = abs(oldma - nma)
                #previous methods
                #if nmadiff <= np.mean(madiff): #max(madiff) + (2*np.mean(madiff)):
                #if nmadiff <= np.mean(madiff):
                if nmadiff <= madiff:
                    groupmovingaverages[tid] = nma
                    groupdifftoma[tid] = (madiff * (tlen - 1) + nmadiff) / tlen
                    modify = True
            else:
                if tlen > 1:
                    grange = rmax - rmin
                    if d <= roundcutoff + grange:
                            oldma = groupmovingaverages[tid]
                            nma = (oldma * tlen + nf) / (tlen + 1)
                            nmadiff = abs(oldma - nma)
                            groupmovingaverages[tid] = nma
                            madiff = groupdifftoma[tid]
                            groupdifftoma[tid] = (madiff * (tlen - 1) + nmadiff) / tlen
                            modify = True
                else: #first one's not free, but comes at a discount
                    #if (d - roundcutoff) - d * roundcutoff <= roundcutoff: #1154745
                    if d <= roundcutoff * 2: #1154797 and way less lenient
                        oldma = groupmovingaverages[tid]
                        nma = (oldma * tlen + nf) / (tlen + 1)
                        nmadiff = abs(oldma - nma)
                        groupmovingaverages[tid] = nma
                        groupdifftoma[tid] = nmadiff
                        modify = True
            if modify:
                trackedkeys[nf] = tid
                trackedma[nma] = trackedma.pop(f)
                trackedgroups[tid].append(c[0])
                linedeletioncounter[tid] //= 2
                if nf < rmin:
                    groupranges[tid][0] = nf
                if nf > rmax:
                    groupranges[tid][1] = nf
                gmin, gmax = groupranges[tid]
                grange = gmax - gmin
                if grange > widestmassrange:
                    widestmassrange = grange
                modeltracker[1] += 1 #matched
            else:
                trackedgroups[uidcount] = c
                trackedkeys[nf] = uidcount
                trackedma[nf] = uidcount
                groupmovingaverages[uidcount] = nf
                groupdifftoma[uidcount] = 0 #this zero won't bog down any averages, same principle new mechanics
                groupranges[uidcount] = [nf, nf]
                uidcount += 1
                modeltracker[0] += 1 #newly added
                foundremoval = groupmovingaverages[tid]
                foundremovals.append(foundremoval)
        
        for fr in foundremovals:
            found.remove(fr)
        nonmatched = np.setdiff1d(previousdata, found)
        nmlen = nonmatched.size - 1
        mzlen = mza.size
        newmodelremovals = []
        #things from previousdata not in found gets +1 to linedeletioncounter
        for n, nm in enumerate(nonmatched.tolist()): #flattened above
            linekey = trackedma[nm]
            linedeletioncounter[linekey] += 1
            if linedeletioncounter[linekey] > deadsignal:
                #determine, out of all matched and nonmatched, which fall into a +/- subisomax distance to this movingma, put the lineuids together in a list to be later intersection_merged for line corrections
                newmodelremovals.append(n)
                modeltracker[3] += 1 #removed
                
                trackedgroups[linekey] = np.array(trackedgroups[linekey]) #more efficient memory storage now that it doesn't need to be appended, not much speed compromised but it is a little slower
                
                #collection subisomax-radius lines with overlapping rt's for line corrections
                correctionradius = set()
                nmoverlaps = nonmatched[np.abs(nm - nonmatched) <= subisomax].tolist()
                mzoverlaps = mza[np.abs(nm - mza) <= subisomax].tolist()
                nmkeys = list(map(trackedma.get, nmoverlaps)) #goes to trackedma
                mzkeys = list(map(trackedkeys.get, mzoverlaps)) #goes to trackedkeys
                correctionradius.update(nmkeys)
                correctionradius.update(mzkeys)
                linecorrections.append(tuple(correctionradius))
            else:
                modeltracker[2] += 1 #nonmatched
        
        wides.append(widestmassrange)
        nonmatched = np.delete(nonmatched, newmodelremovals)
        currentmasskeys = list(map(trackedkeys.get, mza.flatten().tolist()))
        currentmasses = np.array(list(map(groupmovingaverages.get, currentmasskeys)))
        #newtrain = np.append(currentmasses, nonmatched, axis=0)
        #model = spatial.KDTree(newtrain)
        previousdata = np.sort(np.append(currentmasses, nonmatched))
        modeltracking[scanindex] = modeltracker
        scancount += 1

    for cmk in currentmasskeys:
        trackedgroups[cmk] = np.array(trackedgroups[cmk])

    nonmatchedkeys = list(map(trackedma.get, nonmatched.flatten().tolist()))
    for nmk in nonmatchedkeys:
        trackedgroups[nmk] = np.array(trackedgroups[nmk])
    
    timearray = np.array(timearray)

    print(time() - t1, 'line model')
    t2 = time()

    regions = [] #t, b, l, r
    for k, a in trackedgroups.items():
        minmass, mintime, mii = a.min(axis=0)
        maxmass, maxtime, mai = a.max(axis=0)
        wmean = (a[:,0] * a[:,2]).sum() / a[:,2].sum()
        regions.append([minmass, maxmass, mintime, maxtime, wmean, k])

    regions = np.array(regions)
    #regions = regions[regions[:,5].argsort()]
    #(np.arange(regions.shape[0]) == regions[:,8]).all() #passes!
    
    print(time() - t2, 'initial regions -', len(regions))
    t3 = time()
    
    correctiongroups = intersection_merge(linecorrections)
    correctiongroups = [list(i) for i in correctiongroups if len(i) > 1]
    
    timeextension = np.diff(timearray).mean() * minpoints #i'm not a huge fan of this because of the potential for it to connect to completely different things, but there's some shit i just need to connect also...
    linecorrections = []
    for cg in correctiongroups:
        torder = regions[cg,4].argsort()
        ncg = np.array(cg)[torder]
        tregs = regions[ncg]
        masstable = tregs[:,:2]
        masswidths = np.diff(masstable)
        moverlaps = np.logical_and(masstable[:,0] - masswidths.flatten() <= masstable[:,1,None] + masswidths, masstable[:,1] + masswidths.flatten() >= masstable[:,0,None] - masswidths)
        timetable = tregs[:,2:4]
        toverlaps = np.logical_and(timetable[:,0] - timeextension <= timetable[:,1,None] + timeextension, timetable[:,1] + timeextension >= timetable[:,0,None] - timeextension)
        overlaps = np.logical_and(moverlaps, toverlaps)
        overwheres = np.argwhere(overlaps).tolist()
        ogroups = intersection_merge(overwheres)
        for ogs in ogroups:
            if len(ogs) > 1:
                og = list(ogs)
                tmatches = tregs[og]
                tmkeys = ncg[og]
                matchtimes = [trackedgroups[i][:,1].tolist() for i in tmkeys]
                flattimes = list(chain(*matchtimes))
                if len(flattimes) == len(set(flattimes)):
                    linecorrections.append(tmkeys.tolist())
                else:
                    #below is making sure only appropriately ordered links are concected and further connections aren't skipping over other lines. it's essentially assuring the intersection merge below works on the basis of this as a directional graph
                    linkedpairs = {} #pair: distance of means
                    uppers = {} #mkeys: upper of the pair
                    uplinks = {} #mkey: closest above line
                    updists = {} #mkey: distance to closest above
                    downers = {} #mkeys: lower of the pair
                    downlinks = {} #mkey: closest below line
                    downdists = {} #mkey: distance to closest below
                    for ow in overwheres:
                        if len(ogs.intersection(ow)) == 2:
                            mregs = tregs[ow]
                            mkeys = tuple(mregs[:,5].astype(int).tolist())
                            l, r = mregs[:,4]
                            lk, rk = mkeys
                            massdiff = abs(l - r)
                            if l > r:
                                upkey = lk
                                downkey = rk
                                uppers[mkeys] = lk
                                downers[mkeys] = rk
                            else:
                                upkey = rk
                                downkey = lk
                                uppers[mkeys] = rk
                                downers[mkeys] = lk
                            if upkey in downlinks:
                                if massdiff < downdists[upkey]:
                                    downlinks[upkey] = downkey
                                    downdists[upkey] = massdiff
                            else:
                                downlinks[upkey] = downkey
                                downdists[upkey] = massdiff
                            if downkey in uplinks:
                                if massdiff < updists[downkey]:
                                    uplinks[downkey] = upkey
                                    updists[downkey] = massdiff
                            else:
                                uplinks[downkey] = upkey
                                updists[downkey] = massdiff
                            #check timepoints here to make sure they're qualified for merging?
                            linkedpairs[mkeys] = massdiff
                    
                    sortedpairs = sorted(linkedpairs.items(), key=lambda x: x[1])
                    #below is making sure no directly adjacent connections have timepoint redundancy
                    passedpairs = []
                    linkedtimes = {} #pair: [mergedtimes]
                    for mkeys, score in sortedpairs:
                        if uplinks[downers[mkeys]] == uppers[mkeys] or downlinks[uppers[mkeys]] == downers[mkeys]:
                            mtimes = [trackedgroups[i][:,1].tolist() for i in mkeys]
                            mergedtimes = mtimes[0] + mtimes[1]
                            if len(mergedtimes) == len(set(mergedtimes)):
                                passedpairs.append(mkeys)
                                linkedtimes[mkeys] = set(mergedtimes)
                    #intersection merge with non-redundancy requirement for timepoints
                    #this doesn't expand mass-ranges as the signals expand, might be a flaw of this whole process I suppose
                    sn = 0
                    itemgroups = defaultdict(set) #groupn: [members]
                    itemtimes = defaultdict(set) #groupn: [covered timepoints]
                    groupsofitems = {} #line: groupn
                    for items in passedpairs: 
                        locs = set()
                        for i in items: 
                            if i in groupsofitems:
                                locs.add(groupsofitems[i])
                            else:
                                otherline = i
                        combine = False
                        if locs:
                            joiner = min(locs)
                            if len(locs) > 1:
                                combine = True
                                for oldlocs in locs.difference([joiner]):
                                    #the timepoint checks here are to check that non-adjacent connections aren't redundant
                                    if oldlocs in itemtimes:
                                        oldtimes = itemtimes[oldlocs]
                                    else:
                                        oldtimes = linkedtimes[items]
                                    if not itemtimes[joiner].intersection(oldtimes):
                                        for ol in itemgroups[oldlocs]:
                                            groupsofitems[ol] = joiner
                                        itemgroups[joiner].update(itemgroups.pop(oldlocs))
                                        #there should only be one oldloc no matter the iteration, so this all operates without the need for complete prior checking of all timepoint redundancy
                                        itemtimes[joiner].update(oldtimes)
                                        if oldlocs in itemtimes:
                                            del itemtimes[oldlocs]
                                    else:
                                        combine = False
                            else:
                                #check that the non-loc'd item isn't redundant for timepoints
                                oldtimes = trackedgroups[otherline][:,1].tolist()
                                if not itemtimes[joiner].intersection(oldtimes):
                                    combine = True
                        else:
                            joiner = sn
                            sn += 1
                            combine = True
                        if combine:
                            itemgroups[joiner].update(items)
                            itemtimes[joiner].update(linkedtimes[items])
                            for i in items:
                                groupsofitems[i] = joiner
                    linecorrections.extend([list(i) for i in itemgroups.values()])
            #else:
                #nopes.append(tmkeys.tolist())
                #add to the later < minpoint + within massrange check
                #I'll forget about anything else for now, this is good enough. Only lone datapoints < minpoints or signals that didn't pass deadsignal would be left here. minimal error if any

    for lines in linecorrections:
        linegrid = []
        for line in lines:
            for c in trackedgroups[line].tolist():
                linegrid.append(list(c))
            del trackedgroups[line]
        linegrid = np.array(sorted(linegrid, key=lambda x: x[1]))
        trackedgroups[uidcount] = np.array(linegrid)
        uidcount += 1

    groupholder = {}

    kl = list(trackedgroups.keys())
    for k in kl:
        groupholder[k] = trackedgroups.pop(k)

    for n, k in enumerate(kl):
        trackedgroups[n] = groupholder.pop(k)

    print(time() - t3, 'line corrections')
    t4 = time()

    startingpoints = defaultdict(list)
    regions = [] #t, b, l, r
    for k, a in trackedgroups.items():
        minmass, mintime, mii = a.min(axis=0)
        maxmass, maxtime, mai = a.max(axis=0)
        wmean = (a[:,0] * a[:,2]).sum() / a[:,2].sum()
        peakarea = np.trapezoid(a[:,2], a[:,1])
        maxintensity = a[:,2].max()
        regions.append([minmass, maxmass, mintime, maxtime, len(a), peakarea, maxintensity, wmean, k])
        startingpoints[mintime].append(k)
    
    regions = np.array(regions)
    
    newwides = []
    maxrange = 0
    for t in sorted(startingpoints):
        for line in startingpoints[t]:
            minmass, maxmass = regions[line,:2]
            massrange = maxmass - minmass
            if massrange > maxrange:
                maxrange = massrange
        newwides.append(maxrange)

    print('old max mass width', max(wides))
    print('new max mass width', max(newwides))

    print(time() - t4, 'new regions -', len(regions))
    
    regionfile = ''.join((processingdirectory, 'regions.pickle'))
    with open(regionfile, 'wb') as pick:
        pickle.dump(regions, pick)
    
    saverloc = ''.join((processingdirectory, 'trackedgroups.pickle'))
    with open(saverloc, 'wb') as pick:
        pickle.dump(trackedgroups, pick)
    
    saverloc = ''.join((processingdirectory, 'modelinfo.pickle'))
    savedbits = [modeltracking, timearray]
    with open(saverloc, 'wb') as pick:
        pickle.dump(savedbits, pick)
    
    saverloc = ''.join((processingdirectory, 'roundcutoff.pickle'))
    with open(saverloc, 'wb') as pick:
        pickle.dump(roundcutoff, pick)
