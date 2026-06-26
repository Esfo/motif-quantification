from elementalcomponents import nfragmentcompositions, cfragmentcompositions, proton
from database import environment

from itertools import chain, product, combinations
from collections import Counter, defaultdict
from decimal import Decimal, getcontext
from contextlib import contextmanager
from scipy.stats import linregress
from multiprocessing import Lock
from functools import partial
from time import time, sleep
import multiprocessing as mp
from bisect import bisect
import numpy as np
import pickle
import heapq
import lmdb
import math
import csv
import os

getcontext().prec = 50

class LimitedProcessPool:
    def __init__(self, max_processes):
        self.max_processes = max_processes
        self.active_processes = []
        self.manager = mp.Manager()
        self.results_list = self.manager.list()  # Shared list for results

    def add_process(self, target, args):
        while len(self.active_processes) >= self.max_processes:
            self._cleanup_finished_processes()
            if len(self.active_processes) >= self.max_processes:
                sleep(0.1)  # Wait a bit before checking again
        
        print(args[1], 'started')
        wrapped_target = self._wrap_target(target)
        p = mp.Process(target=wrapped_target, args=(args, self.results_list))
        p.start()
        self.active_processes.append(p)

    def _wrap_target(self, target):
        def wrapper(args, results_list):
            result = target(*args)
            results_list.append(result)
        return wrapper

    def _cleanup_finished_processes(self):
        self.active_processes = [p for p in self.active_processes if p.is_alive()]

    def join_all(self):
        for p in self.active_processes:
            p.join()

    def get_results(self):
        return list(self.results_list)  # Convert manager list to regular list

@contextmanager
def limited_process_pool(max_processes):
    pool = LimitedProcessPool(max_processes)
    try:
        yield pool
    finally:
        pool.join_all()

class FragmentOrganizer:
    def __init__(self, librarylocation, processingdirectory, proteome, nprocs, ppmtol, ions):
        self.peptidefilename = processingdirectory + 'peptiderankings.csv'
        self.distributionfilename = processingdirectory + 'distributionrankings.csv'
        self.scanfilename = processingdirectory + 'linerankings.csv'
        self.scanfilename = processingdirectory + 'scanrankings.csv'
        self.subformulafilename = processingdirectory + 'subformularankings.csv'
        self.peptidefilelock = Lock()
        self.distributionfilelock = Lock()
        self.scanfilelock = Lock()
        self.scanfilelock = Lock()
        self.subformulafilelock = Lock()
        self.nprocs = nprocs
        self.ppmmod = ppmtol / 1000000
        self.countrange = 3 #number of flanking AAs to cluster by
        
        divisionfile = ''.join((processingdirectory, 'dividedgroups.pickle'))
        with open(divisionfile, 'rb') as pick:
            self.dividedgroups = pickle.load(pick)
        
        chargesoflinesfile = ''.join((processingdirectory, 'chargesoflines.pickle'))
        with open(chargesoflinesfile, 'rb') as pick:
            self.chargesoflines = pickle.load(pick)
        #chargesoflines = line: charge
        
        analytefile = ''.join((processingdirectory, 'analytefactors.pickle'))
        with open(analytefile, 'rb') as pick:
            self.analytesbydistribution, self.distributionsoflines = pickle.load(pick)[2:4]
        #analytesbydistribution = {} #distid: analyte id
        #distributionsoflines = {} #lineid: distid
        
        linesbyscanbysubformulafile = ''.join((processingdirectory, 'linesbyscanbysubformula.pickle'))
        with open(linesbyscanbysubformulafile, 'rb') as pick:
            self.linesbyscanbysubformula = pickle.load(pick)
        #linesbyscanbysubformula = {} #subformula: scan: [lines]
        
        scoredscansfile = ''.join((processingdirectory, 'scored.ms2.pickle'))
        with open(scoredscansfile, 'rb') as pick:
            self.scoredms2scans = pickle.load(pick)
        #scoredms2scans = {} #scan: [[masses], [intensities], [ion scores]]
        
        self.intensityaverages = {} #scan: average intensity
        for scan, (masses, intensities, scores) in self.scoredms2scans.items():
            self.intensityaverages[scan] = np.mean(intensities)
        
        elementsofprobindicesfile = ''.join((processingdirectory, 'elementsofprobabilityindices.pickle'))
        with open(elementsofprobindicesfile, 'rb') as pick:
            self.elementsofprobabilityindices = pickle.load(pick)
        #elementsofprobabilityindices = {} #prob index: e
    
        probabilityorganizerfile = ''.join((processingdirectory, 'probabilityorganizer.pickle'))
        with open(probabilityorganizerfile, 'rb') as pick:
            self.probabilityorganizer = pickle.load(pick)
        #probabilityorganizer = defaultdict(dict) #prob index: iso: prob
        
        matchprobfile = ''.join((processingdirectory, 'matchprobabilities.pickle'))
        with open(matchprobfile, 'rb') as pick:
            self.matchprobabilities = pickle.load(pick)
        #matchprobabilities = defaultdict(list) #subformula: [prob indices]
        
        subformulasubindsfile = ''.join((processingdirectory, 'subformulasubindices.pickle'))
        with open(subformulasubindsfile, 'rb') as pick:
            self.subformulasubindices = pickle.load(pick)
        #subformulasubindices = defaultdict(list) #subformula: [sub match indices]
        
        submatchsequencesfile = ''.join((processingdirectory, 'submatchsequences.pickle'))
        with open(submatchsequencesfile, 'rb') as pick:
            self.submatchsequences = pickle.load(pick)
        #submatchsequences = {} #submatchindex: sequence
        
        subformulapercentfile = ''.join((processingdirectory, 'subformulapercent.pickle'))
        with open(subformulapercentfile, 'rb') as pick:
            self.subformulapercent = pickle.load(pick)
        #subformulapercent = defaultdict(dict) #subformula: sequence: (subiso abundance rank, subiso abundance)
        
        maxintensitylinesofdistsfile = ''.join((processingdirectory, 'maxintensitylinesofdists.pickle'))
        with open(maxintensitylinesofdistsfile, 'rb') as pick:
            self.maxintensitylinesofdists = pickle.load(pick)
        #maxintensitylinesofdists = defaultdict(dict) #distid: (line, scan)
        
        subformulasofsequencedistributionfile = ''.join((processingdirectory, 'subformulasofsequencedistribution.pickle'))
        with open(subformulasofsequencedistributionfile, 'rb') as pick:
            self.subformulasofsequencedistribution = pickle.load(pick)
        #subformulasofsequencedistribution = {} #seq: dist: subformula

        lineintensitiesofscansfile = ''.join((processingdirectory, 'lineintensitiesofscans.pickle'))
        with open(lineintensitiesofscansfile, 'rb') as pick:
            self.lineintensitiesofscans = pickle.load(pick)
        #lineintensitiesofscans = defaultdict(lambda: defaultdict(dict)) #scan: line: raw area of ms1 line, average of 2 flanking points
        
        self.processingdirectory = processingdirectory
        
        with environment(librarylocation) as env:
            proteomedb = env.open_db((proteome + '.info').encode())
            with env.begin(write=False) as txn:
                with txn.cursor(proteomedb) as cursor:
                    self.aminoacidcomposition = eval(cursor.get('aminoacidcomposition'.encode()).decode())
                    self.elementalmasses = eval(cursor.get('elementalmasses'.encode()).decode())
            defaults = env.open_db('defaults'.encode())
            with env.begin(write=False) as txn:
                with txn.cursor(defaults) as cursor:
                    self.dividingthreshold = float(cursor.get('dividingthreshold'.encode()).decode())
        
        ionlist = list(ions)
        self.ndict = {k: nfragmentcompositions[k] for k in ionlist if k in nfragmentcompositions}
        self.cdict = {k: cfragmentcompositions[k] for k in ionlist if k in cfragmentcompositions}
    
    #this only catches nearest 2 at most
    def nearest_neighbors_ppm_tolerance(self, baselist, flylist):
        indices = {} #baseindex: [flyindex] or [flyindex1, flyindex2]
        distances = {} #baseindex: distance
        for bn, rightfn in enumerate(np.searchsorted(flylist, baselist).tolist()): #iter the tolist for profiling
            b = baselist[bn]
            btol = b * self.ppmmod
            bmin = b - btol
            bmax = b + btol
            
            leftfn = rightfn - 1 #worse case scenario this is -1 -> left = False
            left = False
            leftf = flylist[leftfn]
            if leftf > bmin and leftf < bmax:
                left = True
            
            right = False
            try:
                rightf = flylist[rightfn]
                if rightf > bmin and rightf < bmax:
                    right = True
            except IndexError:
                #rightfn == len(flylist), the iteration is over
                if not left:
                    return indices

            if left and right:
                leftdist = b - leftf
                rightdist = rightf - b
                if leftdist < rightdist:
                    indices[bn] = [leftfn]
                    distances[bn] = leftdist
                elif rightdist < leftdist:
                    indices[bn] = [rightfn]
                    distances[bn] = rightdist
                elif leftdist == rightdist:
                    indices[bn] = [leftfn, rightfn]
                    distances[bn] = leftdist
            elif left:
                leftdist = b - leftf
                indices[bn] = [leftfn]
                distances[bn] = leftdist
            elif right:
                rightdist = rightf - b
                indices[bn] = [rightfn]
                distances[bn] = rightdist
        return indices
    
    #full radius allowing more than 2 matches
    def radius_neighbors_ppm_tolerance(self, baselist, flylist):
        #the ppm on this currently goes out of bounds (by a lot), and i probably won't use this function anyways, replacing it for the above
        f = 0
        pool = []
        matches = {} #baselist index: [flylist indices]
        fiter = enumerate(flylist)
        for bn, b in enumerate(baselist):
            btol = b * self.ppmmod
            bmin = b - btol
            bmax = b + btol
            removals = []
            submatches = []
            for fi, pf in pool:
                if pf < bmin:
                    removals.append([fi, pf])
                elif pf <= bmax:
                    submatches.append(fi)
            for r in removals:
                pool.remove(r)
            while f <= bmax:
                try:
                    i, f = next(fiter)
                    if f >= bmin:
                        pool.append([i, b])
                        if f <= bmax:
                            submatches.append(i)
                except StopIteration:
                    break
            if submatches:
                matches[bn] = submatches
        return matches

    def timeshift(self, experimental, theoretical):
        #normalize both arrays to preserve their relationship
        tlen = len(theoretical)
        ratiolen = tlen - 1
        thmean = sum(theoretical) / tlen
        exmean = sum(experimental) / tlen
        scalefactor = thmean / exmean
        exnorm = [e * scalefactor for e in experimental]
        
        #calculate relative differences
        relativedifferences = [abs((e - t) / t) for e, t in zip(exnorm, theoretical)]
        
        #calculate ratios between consecutive points
        thratio = [theoretical[i+1] / theoretical[i] for i in range(ratiolen)]
        exratio = [exnorm[i+1] / exnorm[i] for i in range(ratiolen)]
        
        #compare ratios
        ratiodiffs = [abs(t - e) / t for t, e in zip(thratio, exratio)]
        
        #combine relative differences and ratio differences
        reldiffmean = sum(relativedifferences) / tlen
        ratiodiffmean = sum(ratiodiffs) / ratiolen
        combinationdiffs = (reldiffmean + ratiodiffmean) / 2
        
        return combinationdiffs

    def difference_maximization(self, arr, double):
        #function for refined two-phase greedy difference maximization
        sarr, sdouble = map(list, zip(*sorted(zip(arr, double))))
        
        #initialize sequence with the largest and smallest elements
        sequence = [sarr.pop(0), sarr.pop(-1)]
        sequencedouble = [sdouble.pop(0), sdouble.pop(-1)]
        
        while sarr:
            #compute the difference for adding either the next smallest or next largest element
            min_value, max_value = sarr[0], sarr[-1]
            
            #compare adding to both ends with the smallest and largest elements
            add_to_left_diff_min = abs(min_value - sequence[0])
            add_to_right_diff_min = abs(min_value - sequence[-1])
            
            add_to_left_diff_max = abs(max_value - sequence[0])
            add_to_right_diff_max = abs(max_value - sequence[-1])
            
            #decide to place the minimum or maximum value based on the maximum possible gain
            if add_to_left_diff_min >= add_to_right_diff_min and add_to_left_diff_min >= add_to_left_diff_max and add_to_left_diff_min >= add_to_right_diff_max:
                sequence.insert(0, sarr.pop(0))
                sequencedouble.insert(0, sdouble.pop(0))
            elif add_to_right_diff_min >= add_to_left_diff_min and add_to_right_diff_min >= add_to_left_diff_max and add_to_right_diff_min >= add_to_right_diff_max:
                sequence.append(sarr.pop(0))
                sequencedouble.append(sdouble.pop(0))
            elif add_to_left_diff_max >= add_to_right_diff_max and add_to_left_diff_max >= add_to_left_diff_min and add_to_left_diff_max >= add_to_right_diff_min:
                sequence.insert(0, sarr.pop(-1))
                sequencedouble.insert(0, sdouble.pop(-1))
            else:
                sequence.append(sarr.pop(-1))
                sequencedouble.append(sdouble.pop(-1))
        
        return sequence, sequencedouble

    def sequence_geometry(self, seq, ioncoverage):
        slen = len(seq)
        maxncoverage = 0
        maxccoverage = 0
        dividers = set()
        ntermcoverage = []
        ctermcoverage = []
        for ion in ioncoverage:
            iontype = ion[0]
            ioncount = int(ion[1:])
            if iontype in 'abc': #nterm
                dividers.add(ioncount)
                pseq = seq[:ioncount]
                ntermcoverage.append(ioncount)
                if ioncount > maxccoverage:
                    maxccoverage = ioncount
            elif iontype in 'xyz': #cterm
                dividers.add(slen - ioncount)
                pseq = seq[slen-ioncount:]
                ctermcoverage.append(slen-ioncount)
                if ioncount > maxncoverage:
                    maxncoverage = ioncount
        dividers = sorted(dividers)
        coverageweight = 1 / (maxncoverage + maxccoverage)

        #isolation counts need to be robust against redundant pseqs
        ind = 0
        ddiff = np.diff(dividers, prepend=0).tolist()
        #dividerstring = ''
        partialseqs = defaultdict(int) #index-pseq: count #the index safeguards against multiple isolations of the same partial sequence, the defaultdict rather than a Counter keeps the keys in order so i can view it easier in regards to the order of the sequence
        for d in ddiff:
            pseq = seq[ind:ind+d]
            #dividerstring += pseq + '|'
            ntermcovers = [i for i in ntermcoverage if i > ind]
            ctermcovers = [i for i in ctermcoverage if i <= ind]
            covers = len(ntermcovers) + len(ctermcovers)
            if covers > 0:
                label = str(ind) + '-' + pseq
                partialseqs[label] += covers
            ind += d
        pseq = seq[ind:]
        #dividerstring += pseq
        ntermcovers = [i for i in ntermcoverage if i > ind]
        ctermcovers = [i for i in ctermcoverage if i <= ind]
        covers = len(ntermcovers) + len(ctermcovers)
        if covers > 0:
            label = str(ind) + '-' + pseq
            partialseqs[label] += covers
        pairsum = 0
        matchcounts = len(set(ioncoverage))
        isolationlengthweight = 1
        for indseq, count in partialseqs.items():
            ind, pseq = indseq.split('-')
            ind = int(ind)
            #i could use this index to weight based on distance from the ends i guess?
            #isolationlengthweight *= len(pseq) / len(seq)
            isolationlengthweight *= 1 / len(pseq) / len(seq) #this 1 / provides an additional layer of geometric success to this scheme, grants success where there was previously failure
            plen = len(pseq)
            matchcounts += plen * count
        dividerweight = 1 / len(partialseqs)
        out1 = 1 / (dividerweight + isolationlengthweight + coverageweight)
        out2 = (1 / dividerweight) + (1 / isolationlengthweight) + (1 / coverageweight)
        return out1, out2, matchcounts
    
    def fragmentation_compositions(self, seq):
        fragments = {}

        #calculate the compositions of the n-term fragments
        fragcomp_n = {}
        for n, aa in enumerate(seq[:-1]):  
            aa_composition = self.aminoacidcomposition[aa]
            for k in aa_composition:
                fragcomp_n[k] = fragcomp_n.get(k, 0) + aa_composition.get(k, 0)
            for ion, modcomp in self.ndict.items():
                fragment_composition = fragcomp_n.copy()
                for k in modcomp:
                    fc = fragment_composition.get(k, 0) + modcomp.get(k, 0)
                    if fc > 0:
                        fragment_composition[k] = fc
                    else:
                        del fragment_composition[k]
                fragments[ion + str(n + 1)] = fragment_composition
        
        #calculate the compositions of the c-term fragments
        fragcomp_c = {}
        for n, aa in enumerate(seq[::-1][:-1]): 
            aa_composition = self.aminoacidcomposition[aa]
            for k in aa_composition:
                fragcomp_c[k] = fragcomp_c.get(k, 0) + aa_composition.get(k, 0)
            for ion, modcomp in self.cdict.items():
                fragment_composition = fragcomp_c.copy()
                for k in modcomp:
                    fc = fragment_composition.get(k, 0) + modcomp.get(k, 0)
                    if fc > 0:
                        fragment_composition[k] = fc
                    else:
                        del fragment_composition[k]
                fragments[ion + str(n + 1)] = fragment_composition
        
        #aa = seq[0]
        #aa_composition = self.aminoacidcomposition[aa]
        #for k in aa_composition:
        #    fragcomp_c[k] = fragcomp_c.get(k, 0) + aa_composition.get(k, 0)
        #fragcomp_c['H'] += 2
        #fragcomp_c['O'] += 1
        #fragments['precursor'] = fragcomp_c
        
        return fragments
    
    def fragment_element_binomial_walk(self, e, acount, fragprobabilities):
        nvector = []
        fragmentvectorpositions = {} #iso: position in vector, replacing nvectorpositions
        fragmentelementpositions = {} #position: iso
        maxinitial = 0
        for n, (iso, prob) in enumerate(fragprobabilities.items()):
            nvector.append(0)
            fragmentvectorpositions[iso] = n
            fragmentelementpositions[n] = iso
            if prob > maxinitial:
                maxinitial = prob
                mk = iso
        lesserfragmentisotopes = [i for i in fragprobabilities if i != mk] #replacing nonmonoisotopicgroups
        elementlist = []
        mainheap = []
        vectorsets = defaultdict(set) #element: set of used vectors
        nvector[fragmentvectorpositions[mk]] += acount
        flen = len(fragprobabilities)
        if flen > 2:
            baseprob = fragprobabilities[mk] ** acount
            preheap = []
            preheap.append([baseprob, acount * self.elementalmasses[mk], e, nvector.copy()])
            greater = True
            lastprob = baseprob
            while greater:
                greater = False
                for iso in lesserfragmentisotopes:
                    newelementvector = nvector.copy()
                    newelementvector[fragmentvectorpositions[mk]] -= 1
                    if newelementvector[fragmentvectorpositions[mk]] > -1:
                        newelementvector[fragmentvectorpositions[iso]] += 1
                        vectorsets[e].add(tuple(newelementvector))
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(newelementvector):
                            loopiso = fragmentelementpositions[n]
                            newelementmass += self.elementalmasses[loopiso] * c
                            newelementprob *= fragprobabilities[loopiso]**c
                            if n > 0:
                                newelementprob *= math.comb(acount-pn, c)
                                pn += c
                        preheap.append([newelementprob, newelementmass, e, newelementvector.copy()])
                        if newelementprob > lastprob:
                            lastprob = newelementprob
                            greater = True
            preheap = sorted(preheap)
            maxiso = preheap[-1]
            maxprob, m, e, nv = maxiso
            elementlist.append([-1, maxprob, m, e, nv])
            maxprob *= -1
            preheap = preheap[:-1]
            for h in preheap:
                r = h[0] / maxprob
                h.insert(0, r)
                heapq.heappush(mainheap, h)
            for iso in lesserfragmentisotopes:
                v = nv.copy()
                v[fragmentvectorpositions[mk]] -= 1
                if v[fragmentvectorpositions[mk]] > -1:
                    v[fragmentvectorpositions[iso]] += 1
                    tuplevec = tuple(v)
                    if tuplevec not in vectorsets[e]:
                        vectorsets[e].add(tuplevec)
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(v):
                            loopiso = fragmentelementpositions[n]
                            newelementmass += self.elementalmasses[loopiso] * c
                            newelementprob *= fragprobabilities[loopiso]**c
                            if n > 0:
                                newelementprob *= math.comb(acount-pn, c)
                                pn += c
                        heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
        else:
            preheap = []
            baseprob = fragprobabilities[mk] ** acount
            preheap.append([baseprob, acount * self.elementalmasses[mk], e, nvector.copy()])
            greater = True
            lastprob = baseprob
            iso = lesserfragmentisotopes[0]
            while greater:
                greater = False
                nvector[fragmentvectorpositions[mk]] -= 1
                if nvector[fragmentvectorpositions[mk]] > -1:
                    nvector[fragmentvectorpositions[iso]] += 1
                    vectorsets[e].add(tuple(nvector))
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(nvector):
                        loopiso = fragmentelementpositions[n]
                        newelementmass += self.elementalmasses[loopiso] * c
                        newelementprob *= fragprobabilities[loopiso]**c
                        if n > 0:
                            newelementprob *= math.comb(acount-pn, c)
                            pn += c
                    preheap.append([newelementprob, newelementmass, e, nvector.copy()])
                    if newelementprob > lastprob:
                        lastprob = newelementprob
                        greater = True
            preheap = sorted(preheap)
            maxiso = preheap[-1]
            maxprob, m, e, nv = maxiso
            elementlist.append([-1, maxprob, m, e, nv])
            maxprob *= -1
            preheap = preheap[:-1]
            for h in preheap:
                r = h[0] / maxprob
                h.insert(0, r)
                heapq.heappush(mainheap, h)
            v = nv.copy()
            v[fragmentvectorpositions[mk]] -= 1
            if v[fragmentvectorpositions[mk]] > -1:
                v[fragmentvectorpositions[iso]] += 1
                tuplevec = tuple(v)
                if tuplevec not in vectorsets[e]:
                    vectorsets[e].add(tuplevec)
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(v):
                        loopiso = fragmentelementpositions[n]
                        newelementmass += self.elementalmasses[loopiso] * c
                        newelementprob *= fragprobabilities[loopiso]**c
                        if n > 0:
                            newelementprob *= math.comb(acount-pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
        
        cutoff = -maxprob * self.dividingthreshold

        r, p, m, e, v = heapq.heappop(mainheap)
        elementlist.append([r, p, m, e, v])
        if flen > 2:
            while p > cutoff:
                for iso in lesserfragmentisotopes:
                    newelementvector = v.copy()
                    newelementvector[fragmentvectorpositions[mk]] -= 1
                    if newelementvector[fragmentvectorpositions[mk]] > 0:
                        newelementvector[fragmentvectorpositions[iso]] += 1
                        tuplevec = tuple(newelementvector)
                        if tuplevec not in vectorsets[e]:
                            vectorsets[e].add(tuplevec)
                            pn = 0
                            newelementmass = 0
                            newelementprob = 1
                            for n, c in enumerate(newelementvector):
                                loopiso = fragmentelementpositions[n]
                                newelementmass += self.elementalmasses[loopiso] * c
                                newelementprob *= fragprobabilities[loopiso]**c
                                if n > 0:
                                    newelementprob *= math.comb(acount-pn, c)
                                    pn += c
                            heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, newelementvector.copy()])
                r, p, m, e, v = heapq.heappop(mainheap)
                elementlist.append([r, p, m, e, v])
                try:
                    r, p, m, e, v = heapq.heappop(mainheap)
                    elementlist.append([r, p, m, e, v])
                except IndexError:
                    #mainheap is empty, this can happen when count is low and probabilities are evenly split. When this happened it was in the below loop, but I'll keep this here too just in case
                    break
        else:
            iso = lesserfragmentisotopes[0]
            while p > cutoff:
                nvector = v.copy()
                nvector[fragmentvectorpositions[mk]] -= 1
                if nvector[fragmentvectorpositions[mk]] > 0:
                    nvector[fragmentvectorpositions[iso]] += 1
                    tuplevec = tuple(nvector)
                    if tuplevec not in vectorsets[e]:
                        vectorsets[e].add(tuplevec)
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(nvector):
                            loopiso = fragmentelementpositions[n]
                            newelementmass += self.elementalmasses[loopiso] * c
                            newelementprob *= fragprobabilities[loopiso]**c
                            if n > 0:
                                newelementprob *= math.comb(acount-pn, c)
                                pn += c
                        heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, nvector.copy()])
                try:
                    r, p, m, e, v = heapq.heappop(mainheap)
                    elementlist.append([r, p, m, e, v])
                except IndexError:
                    #mainheap is empty, this can happen when count is low and probabilities are evenly split
                    break
        heapq.heapify(elementlist)
        return elementlist, fragmentelementpositions
    
    def fragment_descending_partial_products(self, elementalorganizer, fragmentpositions):
        mainpool = defaultdict(list) #things already popped from elementalorganizer
        for k in elementalorganizer:
            mainpool[k].append(heapq.heappop(elementalorganizer[k]))

        subformulas = []
        sumabundances = []
        massnumberindices = {} #mass number: index in the 2 above lists
        
        formula = ''
        maxprob = 1
        mainmass = 0
        massnumber = 0
        for b in sorted(mainpool):
            for r, p, m, e, v in mainpool[b]:
                for n, c in enumerate(v):
                    if c > 0:
                        iso = fragmentpositions[e][n]
                        splitval = 0
                        #for handling elements with multiple letters
                        while True:
                            if iso[splitval].isalpha():
                                splitval += 1
                            else:
                                break
                        massnumber += int(iso[splitval:]) * c
                        formula += f'{iso}({c})'
                maxprob *= p
                mainmass += m
        
        massnumberindices[massnumber] = 0
        subformulas.append(formula)
        sumabundances.append([mainmass * maxprob, maxprob])
        
        cutoff = maxprob * self.dividingthreshold
        mainheap = list(chain(*elementalorganizer.values()))
        heapq.heapify(mainheap)
        
        vectorpool = set()
        multinomialpath = [] #sublists not in mainpool
        probabilityranking = [] #representative lists of ratio probability to sort multinomialpath
        while mainheap:
            r, p, m, e, v = heapq.heappop(mainheap)
            baseiter = {k: v for k, v in mainpool.items() if k != e}
            baseiter[e] = [(r, p, m, e, v)]
            
            formula = ''
            prob = 1
            mass = 0
            massnumber = 0
            for b in sorted(baseiter):
                for sr, sp, sm, se, sv in baseiter[b]:
                    for n, c in enumerate(sv):
                        if c > 0:
                            iso = fragmentpositions[se][n]
                            splitval = 0
                            #for handling elements with multiple letters
                            while True:
                                if iso[splitval].isalpha():
                                    splitval += 1
                                else:
                                    break
                            massnumber += int(iso[splitval:]) * c
                            formula += f'{iso}({c})'
                    prob *= sp
                    mass += sm
            
            try:
                index = massnumberindices[massnumber]
                subformulas[index] += '-' + formula
                sumabundances[index][0] += mass * prob
                sumabundances[index][1] += prob
            except KeyError: #not in there
                index = len(massnumberindices)
                massnumberindices[massnumber] = index
                subformulas.append(formula)
                sumabundances.append([mass * prob, prob])
            if prob < cutoff:
                break
            
            tsv = tuple(v)
            if tsv not in vectorpool:
                ind = bisect(probabilityranking, r)
                probabilityranking.insert(ind, r)
                multinomialpath.insert(ind, (r, p, m, e, v))
                vectorpool.add(tsv)
            
            checkedcombos = set()
            for path in multinomialpath.copy():
                multielement = False
                match path[1]:
                    case tuple():
                        multielement = True
                        sepool = set()
                        sepool.add(e)
                        seformulas = []
                        multipath = []
                        nsr = 1
                        for sr, sp, sm, se, sv in path[1:]:
                            if se not in sepool:
                                nsr *= sr
                                sepool.add(se)
                                sef = ''
                                for n, c in enumerate(sv):
                                    if c > 0:
                                        sef += f'{fragmentpositions[se][n]}({c})'
                                seformulas.append(sef)
                                multipath.append((sr, sp, sm, se, sv))
                        checkformula = ''.join((sorted(seformulas)))
                        if checkformula in checkedcombos:
                            continue
                        else:
                            checkedcombos.add(checkformula)
                        if len(multipath) == 0:
                            continue
                    case _:
                        sr, sp, sm, se, sv = path
                        sef = ''.join((f'{se}{str(n)}{(val)}' for n, val in enumerate(sv)))
                        if sef in checkedcombos:
                            continue
                        else:
                            checkedcombos.add(sef)
                        if se == e:
                            continue
                        nsr = sr
                newratio = nsr * r
                if newratio > 0:
                    newratio *= -1
                if -newratio >= self.dividingthreshold:
                    if multielement:
                        seformula = ''
                        newprob = 1
                        newmass = 0
                        newmassnum = 0
                        newiter = {k: v for k, v in baseiter.items() if k not in sepool}
                        newiter[e] = [(r, p, m, e, v)]
                        for ir, ip, im, ie, iv in multipath:
                            newiter[ie] = [(ir, ip, im, ie, iv)]
                        for b in sorted(newiter):
                            for ir, ip, im, ie, iv in newiter[b]:
                                for n, c in enumerate(iv):
                                    if c > 0:
                                        iso = fragmentpositions[ie][n]
                                        splitval = 0
                                        #for handling elements with multiple letters
                                        while True:
                                            if iso[splitval].isalpha():
                                                splitval += 1
                                            else:
                                                break
                                        newmassnum += int(iso[splitval:]) * c
                                        seformula += f'{iso}({c})'
                                newprob *= ip
                                newmass += im
                    else:
                        newiter = {k: v for k, v in baseiter.items() if k != se}
                        newiter[se] = [(sr, sp, sm, se, sv)]
                        seformula = ''
                        newprob = 1
                        newmass = 0
                        newmassnum = 0
                        for b in sorted(newiter):
                            for ir, ip, im, ie, iv in newiter[b]:
                                for n, c in enumerate(iv):
                                    if c > 0:
                                        iso = fragmentpositions[ie][n]
                                        splitval = 0
                                        #for handling elements with multiple letters
                                        while True:
                                            if iso[splitval].isalpha():
                                                splitval += 1
                                            else:
                                                break
                                        newmassnum += int(iso[splitval:]) * c
                                        seformula += f'{iso}({c})'
                                newprob *= ip
                                newmass += im
                    if newprob >= cutoff:
                        try:
                            index = massnumberindices[newmassnum]
                            subformulas[index] += '-' + seformula
                            sumabundances[index][0] += newmass * newprob
                            sumabundances[index][1] += newprob
                        except KeyError: #not in there
                            index = len(massnumberindices)
                            massnumberindices[newmassnum] = index
                            subformulas.append(seformula)
                            sumabundances.append([newmass * newprob, newprob])
                        if multielement:
                            ind = bisect(probabilityranking, newratio)
                            probabilityranking.insert(ind, newratio)
                            multinomialpath.insert(ind, (newratio, *multipath))
                        else: #this is rarely ever needed, but it is needed
                            newmulti = []
                            tsv = tuple(sv)
                            #should this one be first? does it matter? i don't believe it does
                            if tsv not in vectorpool:
                                newmulti.append((sr, sp, sm, se, sv))
                                vectorpool.add(tsv)
                            tvv = tuple(v)
                            if tvv not in vectorpool:
                                newmulti.append((r, p, m, e, v))
                                vectorpool.add(tvv)
                            if newmulti:
                                ind = bisect(probabilityranking, newratio)
                                probabilityranking.insert(ind, newratio)
                                multinomialpath.insert(ind, (newratio, *newmulti))
                else:
                    break

        subformulas = np.array(subformulas, dtype='S')
        massesandabundances = np.array(sumabundances)
        massesandabundances[:,0] /= massesandabundances[:,1]
        #sorting by intensity
        subformulas = subformulas[massesandabundances[:,1].argsort()[::-1]].tolist()
        massesandabundances = massesandabundances[massesandabundances[:,1].argsort()[::-1]]
        return subformulas, massesandabundances
    
    def group_fragmentation(self, group, count):
        searchtime = 0
        filtertime = 0
        fraglens = 0
        scanlens = 0
        chargeiterations = 0
        nt = time()
        positioncache = {}
        elementalcache = {}
        descentcache = {}
        groupseqs = []
        groupsubformulas = []
        #postfragmenttypes = Counter()
        #postfragmentcounts = defaultdict(lambda: Counter())
        initialmatches = 0
        finalmatches = 0
        subformulaoutput = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))))) #seq: distid: line: scan: subformula: ion charge: ion: [metrics]
        for member in group:
            if '(' in member:
                groupsubformulas.append(member)
            else:
                groupseqs.append(member)
        fragments = {}
        for seq in groupseqs:
            fragments[seq] = self.fragmentation_compositions(seq)
        for subformula in groupsubformulas:
            probindices = {self.elementsofprobabilityindices[i]: self.probabilityorganizer[i] for i in self.matchprobabilities[subformula]}
            subindices = self.subformulasubindices[subformula]
            output, fragmasses = [], []
            for submatchindex in subindices:
                seq = self.submatchsequences[submatchindex]
                for ion, fragcomp in fragments[seq].items():
                    elementalorganizer = {} #element: [[iso heaps]]
                    fragmentpositions = {} #element: position: iso
                    fragstrings = ''
                    for e, c in fragcomp.items():
                        fragprobs = probindices[e]
                        fragstring = str(c) + '/' + '/'.join(('/'.join((k, str(v))) for k, v in probindices[e].items()))
                        fragstrings += fragstring
                        if len(fragprobs) > 1:
                            #try/except is faster than an if/else, so i might as well
                            try:
                                elementlist = elementalcache[fragstring]
                                positions = positioncache[fragstring]
                            except KeyError: #not in cache
                                elementlist, positions = self.fragment_element_binomial_walk(e, c, fragprobs)
                                elementalcache[fragstring] = elementlist
                                positioncache[fragstring] = positions
                            elementalorganizer[e] = elementlist.copy()
                            fragmentpositions[e] = positions
                        else: #no need for cache, only 1 iso
                            iso = list(fragprobs)[0]
                            elementalorganizer[e] = [[-1, 1, self.elementalmasses[iso]*c, e, [c]]]
                            fragmentpositions[e] = {0: iso}
                    try:
                        fragformulas, massesandabundances = descentcache[fragstrings]
                    except KeyError: #not done prior
                        fragformulas, massesandabundances = self.fragment_descending_partial_products(elementalorganizer, fragmentpositions)
                        descentcache[fragstrings] = fragformulas, massesandabundances
                    for n, (m, i) in enumerate(massesandabundances.tolist()):
                        out = (seq, ion, fragformulas[n].decode(), n, i)
                        output.append(out)
                        fragmasses.append(m)
            fragmasses, output = zip(*sorted(zip(fragmasses, output)))
            fragmasses = np.array(fragmasses)
            fraglens += fragmasses.size
            st = time()
            for scan, lines in self.linesbyscanbysubformula[subformula].items():
                if len(lines) > 1:
                    #analyteid = '_'.join((str(self.analytesbydistribution[self.distributionsoflines[i]]) for i in lines))
                    chargeset = set(self.chargesoflines[i] for i in lines)
                    ##if len(chargeset) > 1: -> test passes
                    ##    print('PROBLEM, chargeset len > 1')
                    maxcharge = max(chargeset)
                    #CHECK if ^this is ever different, i'm pretty sure its always the same charge, there shouldn't be different ones
                    #^because even if the same subformula is in the same scan more than once, it will never be of a different charge than itself in another distribution..
                    #linestring = '_'.join((str(i) for i in lines))
                    linesofmatchdistributions = defaultdict(list) #distid: [lines]
                    for line in lines:
                        linesofmatchdistributions[self.distributionsoflines[line]].append(line)
                else:
                    line = lines[0]
                    #analyteid = self.analytesbydistribution[self.distributionsoflines[lines]]
                    maxcharge = self.chargesoflines[line]
                    #linestring = str(lines)
                    linesofmatchdistributions = {self.distributionsoflines[line]: [line]}
                #put fragmasses here -> append precursor ion of the line? because i dont want to calculate precursors via the above dists, but i want to match them
                #^i might just make a special precursor search because this gets too retarded
                #i'm removing precursors from the above calculations for now
                ms2masses, ms2intensities, ms2scores = self.scoredms2scans[scan]
                scanlens += ms2masses.size
                #ms2masses = ms2masses.tolist()
                outputorganizer = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))) #sequence: charge: ion: fragrank: [[metrics],]
                for charge in range(1, maxcharge+1):
                    chargeiterations += 1
                    chargedfragments = ((fragmasses + proton * charge) / charge)
                    #i'm going to ditch the radius neighbors for a nearest neighbors concept instead, ought to have a minor speedboost at least
                        #when MS resolution is low -> you won't get anything close enough to have more than 1 thing in a radius
                        #if its high -> you CAN have this, but you should also expect the masses to be accurate
                            #i can see this in the ms1 vs ms2 data for the fr400 file
                        #this will match at most 2 ions if they both have the same distance to a theoretical fragment
                    matches = self.nearest_neighbors_ppm_tolerance(chargedfragments, ms2masses)
                    #matches = self.radius_neighbors_ppm_tolerance(chargedfragments.tolist(), ms2masses)
                    for fragindex, scanindices in matches.items(): #frag index: [mass index] or [mass index 1, mass index 2]
                        #a scanmass can match to multiple generated fragment ions
                        for scanindex in scanindices:
                            experimentalmass = ms2masses[scanindex]
                            ionscore = ms2scores[scanindex]
                            experimentalintensity = ms2intensities[scanindex]
                            theoreticalmass = chargedfragments[fragindex]
                            #ppmerror = ((experimentalmass - theoreticalmass) / experimentalmass) * 1000000 #maybe this could be done somewhere below instead
                            seq, ion, fragformula, fragrank, theoreticalabundance = output[fragindex]
                            metrics = [fragrank, theoreticalmass, experimentalmass, theoreticalabundance, experimentalintensity, scanindex, fragformula, ionscore]
                            outputorganizer[seq][charge][ion][fragrank].append(metrics)
                            initialmatches += 1
                nst = time()
                for seq, ioncharges in outputorganizer.items():
                    for ioncharge, ions in ioncharges.items():
                        for ion, fragranks in ions.items():
                            fcount = 0
                            for fragrank in sorted(fragranks):
                                if fragrank == fcount:
                                    fcount += 1
                            if fcount > 0:
                                combinationoutputs = defaultdict(list) #either 0: [single ion] or 1: [multiple ions]
                                #^multiple ions will always be chosen over the single in ranked order
                                #^if there's multiple single ions, whichever is closer in intensity to the average intensity of that scan will be chosen
                                if fcount > 1:
                                    #process multiple potential fragiso ranks
                                    #iterate products i guess
                                    #assemble a list of the number of total combinations, check if len > 1 or not
                                    isoiterators = {}
                                    for c in range(fcount):
                                        isoiterators[c] = fragranks[c]
                                    #this is determining which out of all possible frag iso combos makes the best distribution for this matched frag dist
                                    for rankcombos in product(*isoiterators.values()):
                                        #assess on charge distance, its consistency, and abundance modeling
                                        rankcombinations = sorted(rankcombos) #not guaranteed to be most-intense ion as first mass
                                        scorearray = np.array([i[:4] for i in rankcombinations])
                                        rankpairs = np.array([(rankcombinations[i][6], rankcombinations[i+1][6]) for i in range(len(rankcombinations)-1)])
                                        theoreticalabundances = scorearray[:,2]
                                        experimentalintensities = scorearray[:,3]
                                        lowerthbounds = theoreticalabundances[:-1] / theoreticalabundances[1:]
                                        #the closer the max / min ratio of theoretical abundance are to 1, the more accurate the ratio between the 2 needs to be ->> square the theoretical ratio, that's the limit for the experimental
                                        #the min would be half the initial ratio as long as its > 1? maybe
                                        #and the max would be the squared value
                                        upperthbounds = lowerthbounds ** 2 #this decreases lower bounds and increases upper ones.. i think this is fine?
                                        #model by mass but rank by ranks i guess?
                                        #start at the 0 to 1 rank linkage
                                        #and if linkages dont always connect next to each other you can always link to a rank previously seen in this ranking
                                        extimeshift = experimentalintensities[:-1] / experimentalintensities[1:]
                                        acceptablepairs = np.sort(rankpairs[np.logical_and(extimeshift > lowerthbounds, extimeshift < upperthbounds)]).tolist()
                                        try:
                                            firstgroup = acceptablepairs[0]
                                        except IndexError:
                                            #useless, no good matches here
                                            continue
                                        isoindices = set()
                                        #the first two ranks MUST have the correct rank order intensities, or else it just won't be taken, take only the top rank instead
                                        if firstgroup == [0, 1]:
                                            #good, start distassembling
                                            isoindices.update(firstgroup)
                                            #from these rankpairs just accumulate anything adjacent to whatevers already in there
                                            for l, r in acceptablepairs[1:]:
                                                if l in isoindices or r in isoindices:
                                                    #connect anything adjacent
                                                    isoindices.add(l)
                                                    isoindices.add(r)
                                                else:
                                                    #finished accumulating
                                                    break
                                        else:
                                            #mismatch, no good
                                            if 0 in firstgroup:
                                                #take 0
                                                combinationoutputs[0].append([abs(self.intensityaverages[scan]-rankcombos[0][3]), rankcombos[0]])
                                                continue
                                            else:
                                                #useless
                                                continue
                                        finalindices = sorted(isoindices)
                                        #score the dist and add it to the final list
                                        massdiffs = np.diff(scorearray[finalindices,1])
                                        avgmdiff = np.abs(massdiffs.mean() - massdiffs).mean() #mass distance consistency measure
                                        theoreticalabundances = scorearray[finalindices,2].tolist()
                                        experimentalintensities = scorearray[finalindices,3].tolist()
                                        shiftdeviance = self.timeshift(experimentalintensities, theoreticalabundances) #time-series comparison
                                        #shiftdeviance = linregress(experimentalintensities, theoreticalabundances).pvalue
                                        combinationoutputs[1].append([avgmdiff * shiftdeviance, [rankcombos[i] for i in finalindices]])
                                else: #fcount == 1
                                    fullmetrics = fragranks[0]
                                    if len(fullmetrics) > 1:
                                        #multiple matches to this fragrank
                                        #pick whichever is closer in intensity to the average intensity of the scan
                                        scanav = self.intensityaverages[scan]
                                        avdiffs = [abs(scanav-i[3]) for i in fullmetrics]
                                        minav = min(avdiffs)
                                        finalmetric = fullmetrics[avdiffs.index(minav)]
                                        combinationoutputs[0].append([minav, finalmetric])
                                    else:
                                        #single, take it
                                        minav = abs(self.intensityaverages[scan] - fullmetrics[0][3])
                                        combinationoutputs[0].append([minav, fullmetrics[0]])
                                #with either of these results below i'm assuming an equal score would only be given to matches that are exactly the same
                                if 1 in combinationoutputs:
                                    #sort and pick best
                                    selection = min(combinationoutputs[1])
                                elif 0 in combinationoutputs:
                                    #take whichever is nearest to the mean intensity of the scan
                                    selection = min(combinationoutputs[0])
                                else:
                                    #got nothing
                                    continue
                                for distid, alines in linesofmatchdistributions.items():
                                    for line in alines:
                                        #just a quick test -> test passes -> the order is different now but it should still be fine
                                        #if seq in subformulaoutput:
                                        #    if distid in subformulaoutput[seq]:
                                        #        if scan in subformulaoutput[seq][distid]:
                                        #            if line in subformulaoutput[seq][distid][scan]:
                                        #                if subformula in subformulaoutput[seq][distid][scan][line]:
                                        #                    if ioncharge in subformulaoutput[seq][distid][scan][line][subformula]:
                                        #                        if ion in subformulaoutput[seq][distid][scan][line][subformula][ioncharge]:
                                        #                            if selection == subformulaoutput[seq][distid][scan][subformula][ioncharge][ion]:
                                        #                                print('selection present')
                                        #                            else:
                                        #                                print('different selection present')
                                        subformulaoutput[seq][self.analytesbydistribution[distid]][distid][line][scan][subformula][ioncharge][ion] = selection
                            else:
                                #nada
                                continue
                filtertime += time() - nst
            searchtime += time() - st
        fragtime = time() - nt - searchtime
        searchtime -= filtertime
        sct = time()
        peptideleveloutput = []
        distributionleveloutput = []
        scanleveloutput = []
        for seq, analyteids in subformulaoutput.items():
            for analyteid, distributions in analyteids.items():
                analytescore = 1
                analyteionscore = 0
                analyteppm = 0
                analyteintensity = 0
                ioncoverage = set()
                scanindexstring = ''
                intensityratios = []
                abundanceratios = []
                analytefragmentindices = defaultdict(set) #scan: [scan indices]
                for distid, lines in distributions.items():
                    distioncoverage = set()
                    #this ion superset samples the most intense ion of a distribution and determines a top-down superset of all the fragmenting ions to be imposed on every other subformula and MS2 sampling taken at lesser intensities, if an ion doesn't show up here then it's not allowed to contribute to the rest of the scoring/ID process as its inconsistent and probably not real
                    #so ie this assumes all subformulas fragment similar enough for it to matter despite slight isotopic differences - which i think should be ok
                    #starting with the line and scan where this distribution sampled the largest MS1 intensity
                    #line, scan = self.maxintensitylinesofdists[distid]
                    ##if scan in scans:
                    #if line in lines:
                    #    #lines = scans[scan]
                    #    scans = lines[line]
                    #    #if line in lines:
                    #    if scan in scans:
                    #        #this is what should be the most abundant subformula at that position
                    #        subformula = self.subformulasofsequencedistribution[distid][seq]
                    #        subformulas = scans[scan]
                    #        #if subformula in subformulas:
                    #        ioncharges = subformulas[subformula]
                    #        ionsuperset = defaultdict(set) #charge: [ions]
                    #        for ioncharge, ions in ioncharges.items():
                    #            #add to these ions to the superset of this identification instance
                    #            ionsuperset[ioncharge].update(ions)
                    #                #for ion, metrics in ions.items():
                    #        else:
                    #            #no superset to be made, the supposed best match isn't there
                    #            continue
                    #    else:
                    #        continue
                    #else:
                    #    continue
                    scanorder = []
                    ms1intensities = []
                    ms2intensities = []
                    scanlineintensitiesbyion = defaultdict(lambda: Counter()) #line-scan: ion: intensity
                    distppm = 0
                    distscore = 1
                    distionscore = 0
                    distintensity = 0
                    #for scan, lines in scans.items():
                    for line, scans in lines.items():
                        #for line, subformulas in lines.items():
                        for scan, subformulas in scans.items():
                            fragmentindices = set() #all fragmass indices in a scan
                            linescan = str(line) + '-' + str(scan)
                            #sort subformulas using subformulapercent, take only adjacent matches, if something has no superset matches -> break
                            subformulalist = sorted((self.subformulapercent[i][seq], i) for i in subformulas)
                            #subformulalist = sorted((*self.subformulapercent[i][seq], subformula) for i in subformulas)
                            #^which is faster?
                            #main scoring mechanisms:
                                #fragdist multiple -> here
                                #cross-scan consistency -> here as a time series across scans
                                #sequence geometry -> here
                                #cross-subformula entropy -> here
                                #intensity pair entropy -> here -> nah im not implementing this
                                #MS1/MS2 intensity entropy by scan % -> next script -> groups into dists i suppose
                            #aiming to make lower scores better in every case to multiply them all together
                            subformulamassindices = defaultdict(lambda: defaultdict(dict)) #ioncharge: ion: subformula: [masses]
                            subformulaintensities = defaultdict(lambda: defaultdict(dict)) #ioncharge: ion: subformula: intensity sum
                            abundanceofsubformulas = {} #subformula: abundance
                            #fragdistmultiple = 1
                            #fragioncount = 0
                            fragmentintensities = set() #assuming each intensity is unique which is actually false apparently, but within scans maybe more probable
                            #ionsubformulastring = '' #subformula^subformularank%ioncharge&ion&scanindex&ppmerror_ioncharge&ion...-subformula&... in order of subformula abundance
                            #gonna capture intensities across ions and check a timeshift
                            qcount = 0
                            sumppm = 0
                            sumintensity = 0
                            scanionscore = 0
                            scanions = set()
                            for (qualrank, abundance), subformula in subformulalist:
                                #this is multiplying the fragdist scores and assembling info to be used for everything else
                                #iterating and applying multiples in order of decreasing subformula abundance
                                if qualrank == qcount:
                                    #substring = '-' + subformula + '^' + str(qualrank) + '%'
                                    ionmatches = False
                                    ioncharges = subformulas[subformula]
                                    abundanceofsubformulas[subformula] = abundance
                                    for ioncharge, ions in ioncharges.items():
                                        #if ioncharge in ionsuperset:
                                        if True:
                                            sortedions = sorted(ions.items(), key=lambda x: x[1][0]) #sorting by fragdist score, so the fragment indices go to better ion first in case they overlap
                                            #for ion, metrictuple in ions.items():
                                            for ion, metrictuple in sortedions:
                                                #if ion in ionsuperset[ioncharge]:
                                                if True:
                                                    fragdistscore, metrics = metrictuple
                                                    fragintensitysum = 0
                                                    fragmassindices = []
                                                    #i only want to take these ions if i haven't seen that scanindex in this scan yet? or i need a way to determine which i prefer it to be labeled as
                                                    match metrics[0]:
                                                        case list():
                                                            metrics = sorted(metrics) #sorting by fragrank
                                                            for fragrank, theoreticalmass, experimentalmass, theoreticalabundance, experimentalintensity, scanindex, fragformula, ionscore in metrics:
                                                                if scanindex not in analytefragmentindices[scan]:
                                                                    analytefragmentindices[scan].add(scanindex)
                                                                    ppmerror = ((experimentalmass - theoreticalmass) / experimentalmass) * 1000000
                                                                    #substring += f'_{ioncharge}&{ion}&{fragrank}&{scanindex}&{round(ppmerror,8)}&{round(theoreticalmass,8)}&{round(experimentalmass,8)}&{round(theoreticalabundance,8)}&{round(experimentalintensity,8)}&{fragformula}'
                                                                    fragintensitysum += experimentalintensity
                                                                    fragmassindices.append(scanindex)
                                                                    fragmentintensities.add(experimentalintensity)
                                                                    #if qualrank == 0:
                                                                    scanlineintensitiesbyion[linescan][ion] += experimentalintensity
                                                                    abserror = abs(ppmerror)
                                                                    sumppm += abserror
                                                                    distppm += abserror
                                                                    analyteppm += abserror
                                                                    sumintensity += experimentalintensity
                                                                    distintensity += experimentalintensity
                                                                    analyteintensity += experimentalintensity
                                                                    scanionscore += ionscore
                                                                    distionscore += ionscore
                                                                    analyteionscore += ionscore
                                                                    scanions.add(ion)
                                                                    ioncoverage.add(ion)
                                                                    distioncoverage.add(ion)
                                                                    finalmatches += 1
                                                                else:
                                                                    #fragranks finished
                                                                    break
                                                        case float:
                                                            if scanindex not in analytefragmentindices[scan]:
                                                                analytefragmentindices[scan].add(scanindex)
                                                                fragrank, theoreticalmass, experimentalmass, theoreticalabundance, experimentalintensity, scanindex, fragformula, ionscore = metrics
                                                                ppmerror = ((experimentalmass - theoreticalmass) / experimentalmass) * 1000000
                                                                #substring += f'_{ioncharge}&{ion}&{fragrank}&{scanindex}&{round(ppmerror,8)}&{round(theoreticalmass,8)}&{round(experimentalmass,8)}&{round(theoreticalabundance,8)}&{round(experimentalintensity,8)}&{fragformula}'
                                                                fragintensitysum += experimentalintensity
                                                                fragmassindices.append(scanindex)
                                                                fragmentintensities.add(experimentalintensity)
                                                                #if qualrank == 0:
                                                                scanlineintensitiesbyion[linescan][ion] += experimentalintensity
                                                                abserror = abs(ppmerror)
                                                                sumppm += abserror
                                                                distppm += abserror
                                                                analyteppm += abserror
                                                                sumintensity += experimentalintensity
                                                                distintensity += experimentalintensity
                                                                analyteintensity += experimentalintensity
                                                                scanionscore += ionscore
                                                                distionscore += ionscore
                                                                analyteionscore += ionscore
                                                                scanions.add(ion)
                                                                ioncoverage.add(ion)
                                                                distioncoverage.add(ion)
                                                                finalmatches += 1
                                                    #fragdistmultiple *= fragdistscore
                                                    #if ioncharge in subformulamassindices: -> test passes
                                                    #    if ion in subformulamassindices[ioncharge]:
                                                    #        if subformula in subformulamassindices[ioncharge][ion]:
                                                    #            print('problem')
                                                    if fragmassindices:
                                                        subformulamassindices[ioncharge][ion][subformula] = fragmassindices
                                                        subformulaintensities[ioncharge][ion][subformula] = fragintensitysum
                                                        fragmentindices.update(fragmassindices)
                                                    #ionmatches = True
                                                    #fragioncount += 1
                                    #if ionmatches:
                                    #    ionsubformulastring += substring
                                    qcount += 1
                                else:
                                    #descending order of subformulas is finished
                                    break
                            #for ion in ioncoverage:
                            #    #tracking fragment ions to infer the likelihood of their surroundings
                            #    iontype = ion[0]
                            #    postfragmenttypes[iontype] += 1
                            #    if iontype in 'abc': #nfrags
                            #        ioncount = int(ion[1:])
                            #        npartialseq = seq[ioncount-self.countrange:ioncount]
                            #        if len(npartialseq) == self.countrange: #partialseq isn't cut off by the end of the sequence
                            #            postfragmentcounts['n'][npartialseq] = 1
                            #        cpartialseq = seq[ioncount:ioncount+self.countrange+1]
                            #        if len(cpartialseq) == self.countrange: #partialseq isn't cut off by the end of the sequence
                            #            postfragmentcounts['n'][cpartialseq] = 1
                            #    elif iontype in 'xyz': #cfrags
                            #        ioncount = len(seq) - int(ion[1:])
                            #        npartialseq = seq[ioncount-self.countrange:ioncount]
                            #        if len(npartialseq) == self.countrange: #partialseq isn't cut off by the end of the sequence
                            #            postfragmentcounts['c'][npartialseq] = 1
                            #        cpartialseq = seq[ioncount:ioncount+self.countrange+1]
                            #        if len(cpartialseq) == self.countrange: #partialseq isn't cut off by the end of the sequence
                            #            postfragmentcounts['c'][cpartialseq] = 1
                            
                            #if len(fragmentintensities) > 1:
                            if len(scanions) > 1:
                                #distionscore += scanionscore
                                #avgppm = sumppm / sumintensity
                                #distppm += sumppm
                                #distintensity += sumintensity
                                
                                indexstring = '/'.join((map(str, sorted(fragmentindices))))
                                scanindexstring += f'{scan}[{indexstring}]'
                                
                                ms1intensities.append(self.lineintensitiesofscans[scan][line])
                                matchedintensitysum = sum(fragmentintensities)
                                ms2intensities.append(matchedintensitysum)
                                scanorder.append(scan)
                                
                                #intensity scoring based on whether the entirety of the intensity is kept by just 1 ion or if its well-dispersed (which is presumably better)
                                #case: 1 ion takes up 20% of the total intensity of 9 ions -> 1/(.8*8)
                                #case: 1 ion takes up 80% of the total intensity of 4 ions -> 1/(.2*3)
                                #etc
                                #matchesintensitymax = max(fragmentintensities)
                                #intensitydispersion = 1 / (((matchedintensitysum - matchesintensitymax) / matchedintensitysum) * (len(fragmentintensities) - 1)) #intensity dispersion
                                #invertedsum = 1 / matchedintensitysum
                                #analytescore *= intensitydispersion * invertedsum
                                #analytescore *= avgppm
                                g1, g2, g3 = self.sequence_geometry(seq, scanions)
                                #scanoutput = f'{seq},{analyteid},{distid},{line},{scan},{scanionscore},{avgppm},{intensitydispersion},{invertedsum}'
                                scanoutput = f'{seq},{analyteid},{distid},{line},{scan},{indexstring},{scanionscore},{sumppm},{sumintensity},{g1},{g2},{g3}'
                                scanleveloutput.append(scanoutput)
                            #else:
                            #    #intensitycoverage = 1
                            #    if sumintensity:
                            #        scanoutput = f'{seq},{analyteid},{distid},{line},{scan},{scanionscore},{avgppm},1,1'
                            #    else:
                            #        scanoutput = f'{seq},{analyteid},{distid},{line},{scan},{scanionscore},0,1,1'
                            #if fragmentindices:
                            
                                #assess cross-subformula entropy if its present
                                #maybe these don't need to be subformulas that are only at the same location, this could be done across all lines? you'd need to normalize by line intensities
                                for ioncharge, ions in subformulamassindices.items():
                                    #this is taking all mass-shifted ion matches and puttig them in a time series as identifying the same ion from different subformulas at different masses is something i consider to be good evidence of a match
                                    #scanions.update(ions)
                                    for ion, ionsubformulas in ions.items():
                                        if len(ionsubformulas) > 1:
                                            for l, r in combinations(ionsubformulas, 2):
                                                #it could possibly be slower to make these into sets while adding them to subformulamassindices but i think this wouldn't happen enough to make that faster, test it later if you need
                                                if not set(ionsubformulas[l]).intersection(ionsubformulas[r]):
                                                    #all indices are different -> compare total intensities
                                                    intensityratio = subformulaintensities[ioncharge][ion][l] / subformulaintensities[ioncharge][ion][r]
                                                    abundanceratio = abundanceofsubformulas[l] / abundanceofsubformulas[r]
                                                    #if abundanceratio > 1:
                                                    #    #maintain the scaling abundanceratio to be between 0 and 1
                                                    #    intensityratio = 1 / intensityratio
                                                    #    abundanceratio = 1 / abundanceratio
                                                    intensityratios.append(intensityratio)
                                                    abundanceratios.append(abundanceratio)
                                                #else:
                                                #    print(len(set(ionsubformulas[l]).intersection(ionsubformulas[r])), 'overlapped shifts')
                                                    #if you can't do either of those because the indices overlap -> don't add the comparison to the score
                                                    #i'm not going to do max-max comparisons because not all distributions of the same ion across subformulas are similar enough in nature for this to be wise
                                
                                #if len(intensityratios) > 1:
                                #    #i guess i could make this > 2 and make a 2 option but i'll be lazy for now
                                #    #this maximizes the adjacent differences of each matched ion reflected across the abundance of the subformula it belongs to
                                #    #it turns it into a hard-to-fake time series, then i use a timeshift concept to score the ion matches based on their experimental intensity ratios compared to the theoretical ones that ought to be present based on theoretical subformula abundance
                                #    #print(ilen, 'length intensity / abundance ratios') #not seen this happen yet
                                #    maximalorderedabundances, maximalorderedintensities = self.difference_maximization(abundanceratios, intensityratios)
                                #    maximaldiffsubformularatiotimeshift = self.timeshift(maximalorderedintensities, maximalorderedabundances)
                                #    #maximaldiffsubformularatiotimeshift = linregress(maximalorderedintensities, maximalorderedabundances).pvalue
                                #    analytescore *= maximaldiffsubformularatiotimeshift
                                #    
                                #    #rawsubformularatiotimeshift = self.timeshift(intensityratios, abundanceratios)
                                #    
                                #    #sortedabundanceratios, sortedintensityratios = zip(*sorted(zip(abundanceratios, intensityratios)))
                                #    #sortedsubformularatiotimeshift = self.timeshift(sortedintensityratios, sortedabundanceratios)
                                #    #sortedsubformularatiotimeshift = linregress(sortedintensityratios, sortedabundanceratios).pvalue
                                #    
                                #    scanoutput += f',{maximaldiffsubformularatiotimeshift}'
                                #    scanleveloutput.append(scanoutput)
                                #else:
                                #    #if not scanoutput.endswith('0,1,1'):
                                #    #    scanoutput += ',1,1,1,1'
                                #    #    scanleveloutput.append(scanoutput)
                                #    scanoutput += ',1'
                                #    scanleveloutput.append(scanoutput)
                                #elif ilen == 1:
                                #    #idk if this one is a good idea?
                                #    analytescore *= abs(intensityratios[0] - abundanceratios[0]) / abundanceratios[0]
                                #else:
                                #    #length 0
                                #    subformulashiftdeviance = 1
                                #subformulaidentificationmultiple = 1 / qcount #idk about this one, can't necessarily verify the subformula presence with this
                                #finalscore = geometry * subformulashiftdeviance * intensitycoverage
                                #finaldistoutput = f'{seq},{analyteid},{distid},{scan},{line},{ionsubformulastring},{qcount+1},{ilen},{fragioncount},{len(ioncoverage)},{len(fragmentindices)},{matchedintensitysum},{geometry},{subformulashiftdeviance},{intensitycoverage},{finalscore}' + '\n'
                                #subformuladistoutput.append(finaldistoutput)
                    if len(distioncoverage) > 1:
                        #analyteionscore += distionscore
                        #analyteintensity += distintensity
                        #avgdistppm = distppm / distintensity
                        #scanlineintensitiesbyion = defaultdict(lambda: defaultdict(lambda: defaultdict(float))) #line-scan: ion: intensity
                        ctshifts = 1
                        for l, r in combinations(scanlineintensitiesbyion, 2):
                            llen = len(scanlineintensitiesbyion[l])
                            rlen = len(scanlineintensitiesbyion[r])
                            keys = set(scanlineintensitiesbyion[l]).intersection(scanlineintensitiesbyion[r])
                            if len(keys) > 1:
                                llist, rlist = [], []
                                for k in keys:
                                    llist.append(scanlineintensitiesbyion[l][k])
                                    rlist.append(scanlineintensitiesbyion[r][k])
                                iontimeshift = self.timeshift(llist, rlist)
                                #iontimeshift = linregress(llist, rlist).pvalue
                                if iontimeshift > 0:
                                    ctshifts *= iontimeshift
                        #this is a cross-scan consistency of matched fragion intensity that doubles as MS1 time series entropy
                        #scanorder, ms1intensities, ms2intensities = zip(*sorted(zip(scanorder, ms1intensities, ms2intensities)))
                        #intensitytimeshift = self.timeshift(ms2intensities, ms1intensities)
                        #intensitytimeshift = linregress(ms2intensities, ms1intensities).pvalue
                        dg1, dg2, dg3 = self.sequence_geometry(seq, distioncoverage)
                        #analytescore *= intensitytimeshift * ctshifts
                        #analytescore *= ctshifts
                        #analyteppm += distppm
                        #should probably factor in the length of each series for this timeshift i guess
                        distoutput = f'{seq},{analyteid},{distid},{distionscore},{distppm},{distintensity},{dg1},{dg2},{dg3},{ctshifts}'
                        distributionleveloutput.append(distoutput)
                if len(ioncoverage) > 1:
                    ilen = len(intensityratios)
                    if ilen > 1:
                        maximalorderedabundances, maximalorderedintensities = self.difference_maximization(abundanceratios, intensityratios)
                        subformularatiotimeshift = self.timeshift(maximalorderedintensities, maximalorderedabundances)
                        analytescore *= subformularatiotimeshift
                    elif ilen == 1:
                        subformularatiotimeshift = abs(intensityratios[0] - abundanceratios[0]) / abundanceratios[0]
                    else:
                        subformularatiotimeshift = 1
                    #applying geometric logic to the coverage of the sequence based on the matched ions
                    #ioncoverage = list(ioncoverage)
                    ag1, ag2, ag3 =  self.sequence_geometry(seq, ioncoverage)
                    #analytescore *= analytegeometry
                    #avganalyteppm = analyteppm / analyteintensity
                    #analytescore *= avganalyteppm
                    ioncoveragestring = '/'.join(map(str, sorted(ioncoverage)))
                    #finalscore = Decimal(1) / Decimal(analytescore)
                    #analytescore *= analyteionscore
                    #finalscore = Decimal(analytescore)
                    outputstring = f'{seq},{analyteid},{ioncoveragestring},{scanindexstring},{analyteionscore},{analyteppm},{analyteintensity},{ag1},{ag2},{ag3},{subformularatiotimeshift}'
                    peptideleveloutput.append(outputstring)
                #else:
                #    geometry = 1
        #explore:
            #how many subformula ion subsets can you actually find in the data
                #regarding the basis for subformula ion differences:
                    #there does tend to be areas for difference, usually 1-2 dalton in some direction, but not that often, usually only 1-2 times in a single set of subisos
            #can this subformula ion filtering present more b/y ions in an abcxyz search?
        #next script:
            #multi-seq entropy
                #minimized scanmass index overlap
                #intensity pair sum division -> intensity pair entropy
                #MS1 entropy across sequences
                #(intra-subformula entropy across mass-shifted ions is done in this file)
            #split lines and analyteids multiples into individual pieces
        #confidence:
            #i'm thinking of looking at all distributions of each kind of score individually
            #then taking the area under that curve
            #and ranking each individual score along what percentile of area that is away from the best score
        #ranking
            #peptide -> accumulates lines/analytes in scans -> different lines can be labeled as the same subformula, and fixed later if its even listed
        scoretime = time() - sct
        with self.peptidefilelock:
            with open(self.peptidefilename, 'a', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                for piece in peptideleveloutput:
                    writer.writerow(piece.split(','))
        
        with self.distributionfilelock:
            with open(self.distributionfilename, 'a', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                for piece in distributionleveloutput:
                    writer.writerow(piece.split(','))
        
        with self.scanfilelock:
            with open(self.scanfilename, 'a', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                for piece in scanleveloutput:
                    writer.writerow(piece.split(','))
        
        #for term, pseqs in postfragmentcounts.items():
        #    postfragmentcounts[term] = dict(pseqs)
        #postfragmentcounts = dict(postfragmentcounts)
        
        #postfragmenttypes = dict(postfragmenttypes)
        
        print(f'{round(time() - nt, 2)} - group {count} -> fragtime: {round(fragtime, 2)}, searchtime: {round(searchtime, 2)}, filtertime: {round(filtertime, 2)}, scoretime {round(scoretime, 2)}, generated fragments: {fraglens}, scan masses: {scanlens}, charge-iterations: {chargeiterations}, initial fragment matches: {initialmatches}, finalized fragment matches: {finalmatches}')
        #return postfragmentcounts, postfragmenttypes
    
    def group_processing(self):
        
        peptideheaders = ['sequence', 'analyteid', 'ion_coverage', 'scan_indices', 'analyte_ion_score', 'analyte_ppm_error', 'analyte_match_intensity', 'ag1', 'ag2', 'ag3', 'subformula_ratio_shift']
        with open(self.peptidefilename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(peptideheaders)
        
        distributionheaders = ['sequence', 'analyteid', 'distribution', 'dist_ion_score', 'dist_ppm_error', 'dist_match_intensity', 'dg1', 'dg2', 'dg3', 'adjacent_intensity_timeshift']
        with open(self.distributionfilename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(distributionheaders)
        
        scanheaders = ['sequence', 'analyteid', 'distribution', 'line', 'scan', 'scan_indices', 'scan_ion_score', 'ppm_error', 'matched_intensity_sum', 'sg1', 'sg2', 'sg3']
        with open(self.scanfilename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(scanheaders)
        
        print(len(self.dividedgroups), 'dividedgroups total')
        
        with limited_process_pool(self.nprocs) as pool:
            for count, divgroup in enumerate(self.dividedgroups):
                #pool.add_process(target=self.group_fragmentation, args=(divgroup, count))
                self.group_fragmentation(divgroup, count)
        #precounts, postcounts, pretypes, posttypes = zip(*pool.get_results())
        #postcounts, posttypes = zip(*pool.get_results())
        
        #prefragmentcounts = defaultdict(lambda: Counter())
        #for pc in precounts:
        #    for term, pseqs in pc.items():
        #        for pseq, count in pseqs.items():
        #            prefragmentcounts[term][pseq] += count
        #for term in prefragmentcounts:
        #    prefragmentcounts[term] = dict(prefragmentcounts[term])
        #prefragmentcounts = dict(prefragmentcounts)
        
        #postfragmentcounts = defaultdict(lambda: Counter())
        #for pc in postcounts:
        #    for term, pseqs in pc.items():
        #        for pseq, count in pseqs.items():
        #            postfragmentcounts[term][pseq] += count
        #for term in postfragmentcounts:
        #    postfragmentcounts[term] = dict(postfragmentcounts[term])
        #postfragmentcounts = dict(postfragmentcounts)
        
        #prefragmenttypes = Counter()
        #for pt in pretypes:
        #    for ion, count in pt.items():
        #        prefragmenttypes[ion] += count
        #prefragmenttypes = dict(prefragmenttypes)

        #postfragmenttypes = Counter()
        #for pt in posttypes:
        #    for ion, count in pt.items():
        #        postfragmenttypes[ion] += count
        #postfragmenttypes = dict(postfragmenttypes)
        
        #prefragmentcountsfile = ''.join((self.processingdirectory, 'prefragmentcounts.pickle'))
        #with open(prefragmentcountsfile, 'wb') as pick:
        #    pickle.dump(prefragmentcounts, pick)
        ##prefragmentcounts = defaultdict(lambda: Counter()) #n/c-term: ion: count
        
        #postfragmentcountsfile = ''.join((self.processingdirectory, 'postfragmentcounts.pickle'))
        #with open(postfragmentcountsfile, 'wb') as pick:
        #    pickle.dump(postfragmentcounts, pick)
        ##postfragmentcounts = defaultdict(lambda: Counter()) #n/c-term: ion: count
        
        #prefragmenttypesfile = ''.join((self.processingdirectory, 'prefragmenttypes.pickle'))
        #with open(prefragmenttypesfile, 'wb') as pick:
        #    pickle.dump(prefragmenttypes, pick)
        ##prefragmenttypes = Counter() #ion type: count
        
        #postfragmenttypesfile = ''.join((self.processingdirectory, 'postfragmenttypes.pickle'))
        #with open(postfragmenttypesfile, 'wb') as pick:
        #    pickle.dump(postfragmenttypes, pick)
        ##postfragmenttypes = Counter() #ion type: count

def fragment_writer(librarylocation, processingdirectory, proteome, nprocs, ppmtol, ions):
    nt = time()
    
    model = FragmentOrganizer(librarylocation, processingdirectory, proteome, nprocs, ppmtol, ions)
    model.group_processing()
    
    print(time() - nt, 'fragment writing complete')


librarylocation = '/home/sfo/data/proteomics/fastas/search-db/'
proteome = 'Human_Homo_sapien-NoTremb'
processingdirectory = '/home/sfo/store/flowcharacterizations/round3/fileprocessing/200901_fR_400/'
ppmtol = 25
ions = 'by'
nprocs = os.cpu_count()
self = FragmentOrganizer

fragment_writer(librarylocation, processingdirectory, proteome, nprocs, ppmtol, ions)
