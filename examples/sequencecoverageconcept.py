import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
from colour import Color
import psutil
import asyncio
import aiofiles
from pyteomics import mzml
import csv
import bisect
import heapq
from time import time
import pandas as pd
import gc
import concurrent.futures
import multiprocessing as mp
from collections import Counter, defaultdict
from scipy import sparse, integrate, spatial, stats, special
from pandas.api.types import CategoricalDtype
from sklearn.neighbors import NearestNeighbors
from more_itertools import sort_together
from distinctipy import distinctipy as dp
from functools import partial
from pickleshare import PickleShareDB
import math
import zlib
import lmdb
import random
import itertools
import string
import pickle
import sys
import os
import warnings
#warnings.filterwarnings("error")
np.set_printoptions(suppress=True)
gc.enable()

#for more https://matplotlib.org/stable/tutorials/introductory/customizing.html
if os.uname()[1] == 'toaster':
    plt.rcParams['figure.dpi'] = 180
elif os.uname()[1] == 'box':
    plt.rcParams['figure.dpi'] = 300
plt.rcParams['axes.facecolor'] = 'gray'
plt.rcParams['figure.facecolor'] = 'gray'
plt.rcParams['axes.labelcolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'white'
plt.rcParams['ytick.labelcolor'] = 'white'
plt.rcParams['xtick.labelcolor'] = 'white'
plt.rcParams['ytick.color'] = 'white'
plt.rcParams['xtick.color'] = 'white'
chexes = ['#ffffff',
        '#8ff6ff',
        '#ff9f9c',
        '#2ded8d',
        '#fbffb3',
        '#ea68f2',
        '#7d26ff',
        ]
plt.rcParams['axes.prop_cycle'] = plt.cycler('color', chexes)
#n = 0
#for c in chexes:
#    plt.plot([n, n+1], [n, n+1], color=c, label=c)
#    n += 1
#plt.legend()
#plt.show()

subisotopomericdepth = 0.8
proton = 1.007276554940804
dividingthreshold = 0.1
ions = 'by'

peptidelength = 7
nseqs = 1000

npicks = 2
ppmtol = 25 #ppm

minintensity = 1e2
maxintensity = 1e5


proton = 1.007276554940804

elementalprobabilities = { #isotope: abundance
        'H1': 0.999885,
        'H2': 0.000115,
        'C12': 0.9893,
        'C13': 0.0107,
        'N14': 0.99636,
        'N15': 0.00364,
        'O16': 0.99757,
        'O17': 0.00038,
        'O18': 0.00205,
        'S32': 0.9499,
        'S33': 0.0075,
        'S34': 0.0425,
        'S36': 0.0001}

elementalmasses = { #isotope: mass
            'H1': 1.00782503223,
            'H2': 2.01410177812,
            'C12': 12.0000000, 
            'C13': 13.00335483507,
            'N14': 14.00307400443,
            'N15': 15.00010889888,
            'O16': 15.99491461957,
            'O17': 16.99913175650,
            'O18': 17.99915961286,
            'S32': 31.9720711744,
            'S33': 32.9714589098,
            'S34': 33.967867004,
            'S36': 35.96708071}

vectorrangesbyelement = {'H': range(0,2),
                         'C': range(2,4),
                         'N': range(4,6),
                         'O': range(6,9),
                         'S': range(9,13)}

nonmonoisotopicelements = {'H2', 'C13', 'N15', 'O17', 'O18', 'S33', 'S34', 'S36'}

isotopesbyelement = { #element: isotopes
        'H': ('H1', 'H2'),
        'C': ('C12', 'C13'),
        'N': ('N14', 'N15'),
        'O': ('O16', 'O17', 'O18'),
        'S': ('S32', 'S33', 'S34', 'S36')}

monoisotopickeys = { #element: monoisotopic element
        'H': 'H1',
        'C': 'C12',
        'N': 'N14',
        'O': 'O16',
        'S': 'S32'}

nonmonoisotopicgroups = { #element: nonmonoisotopic elements
        'H': ('H2',),
        'C': ('C13',),
        'N': ('N15',),
        'O': ('O17', 'O18'),
        'S': ('S33', 'S34', 'S36')}

elementvectors = {}
vectorpositions = {}
elementpositions = {}
for e, isos in isotopesbyelement.items():
    elementvectors[e] = [0 for _ in range(len(isos))]
    vectorpositions[e] = {k: n for n, k in enumerate(isos)}
    elementpositions[e] = {n: k for n, k in enumerate(isos)}

aminoacidcomposition = {
        'A': {'C': 3, 'H': 5, 'N': 1, 'O': 1},
        'R': {'C': 6, 'H': 12, 'N': 4, 'O': 1},
        'N': {'C': 4, 'H': 6, 'N': 2, 'O': 2},
        'D': {'C': 4, 'H': 5, 'N': 1, 'O': 3},
        'C': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1},
        'Q': {'C': 5, 'H': 8, 'N': 2, 'O': 2},
        'E': {'C': 5, 'H': 7, 'N': 1, 'O': 3},
        'G': {'C': 2, 'H': 3, 'N': 1, 'O': 1},
        'H': {'C': 6, 'H': 7, 'N':3, 'O': 1},
        'I': {'C': 6, 'H': 11, 'N':1, 'O': 1},
        'L': {'C': 6, 'H': 11, 'N':1, 'O': 1},
        'K': {'C': 6, 'H': 12, 'N': 2, 'O': 1},
        'M': {'C': 5, 'H': 9, 'N':1, 'O': 1, 'S': 1},
        'F': {'C': 9, 'H': 9, 'N':1, 'O': 1},
        'P': {'C': 5, 'H': 7, 'N':1, 'O': 1},
        'S': {'C': 3, 'H': 5, 'N':1, 'O': 2},
        'T': {'C': 4, 'H': 7, 'N':1, 'O': 2},
        'W': {'C': 11, 'H': 10, 'N': 2, 'O': 1},
        'Y': {'C': 9, 'H': 9, 'N': 1, 'O': 2},
        'V': {'C': 5, 'H': 9, 'N': 1, 'O': 1}
        }

nfragmentcompositions = {}
nfragmentcompositions['a'] = Counter({'C': -1, 'O': -1})
nfragmentcompositions['b'] = Counter()
nfragmentcompositions['c'] = Counter({'N': 1, 'H': 3})

cfragmentcompositions = {}
cfragmentcompositions['x'] = Counter({'C': 1, 'O': 2})
cfragmentcompositions['y'] = Counter({'H': 2, 'O': 1})
cfragmentcompositions['z'] = Counter({'N': -1, 'H': -1, 'O': 1})
#cfragmentcompositions['z•'] = Counter({'N': -1, 'O': 1})

def individual_element_binomial_walk(dividingthreshold, e, acount):
    #elementalorganizer = defaultdict(list)
    elementlist = []
    mainheap = []
    vectorsets = defaultdict(set) #element: set of used vectors
    #for e, acount in atomiccomposition.items():
    mk = monoisotopickeys[e]
    nvector = elementvectors[e].copy()
    nvector[vectorpositions[e][mk]] += acount
    if len(isotopesbyelement[e]) > 2:
        baseprob = elementalprobabilities[mk] ** acount
        preheap = []
        preheap.append([baseprob, acount * elementalmasses[mk], e, nvector.copy()])
        greater = True
        lastprob = baseprob
        while greater:
            greater = False
            for iso in nonmonoisotopicgroups[e]:
                newelementvector = nvector.copy()
                newelementvector[vectorpositions[e][mk]] -= 1
                if newelementvector[vectorpositions[e][mk]] > -1:
                    newelementvector[vectorpositions[e][iso]] += 1
                    vectorsets[e].add(tuple(newelementvector))
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(newelementvector):
                        loopiso = elementpositions[e][n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= elementalprobabilities[loopiso]**c
                        if loopiso in nonmonoisotopicelements:
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
        for iso in nonmonoisotopicgroups[e]:
            v = nv.copy()
            v[vectorpositions[e][mk]] -= 1
            if v[vectorpositions[e][mk]] > -1:
                v[vectorpositions[e][iso]] += 1
                tuplevec = tuple(v)
                if tuplevec not in vectorsets[e]:
                    vectorsets[e].add(tuplevec)
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(v):
                        loopiso = elementpositions[e][n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= elementalprobabilities[loopiso]**c
                        if loopiso in nonmonoisotopicelements:
                            newelementprob *= math.comb(acount-pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
    else:
        preheap = []
        baseprob = elementalprobabilities[mk] ** acount
        preheap.append([baseprob, acount * elementalmasses[mk], e, nvector.copy()])
        greater = True
        lastprob = baseprob
        iso = nonmonoisotopicgroups[e][0]
        while greater:
            greater = False
            nvector[vectorpositions[e][mk]] -= 1
            if nvector[vectorpositions[e][mk]] > -1:
                nvector[vectorpositions[e][iso]] += 1
                vectorsets[e].add(tuple(nvector))
                pn = 0
                newelementmass = 0
                newelementprob = 1
                for n, c in enumerate(nvector):
                    loopiso = elementpositions[e][n]
                    newelementmass += elementalmasses[loopiso] * c
                    newelementprob *= elementalprobabilities[loopiso]**c
                    if loopiso in nonmonoisotopicelements:
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
        v[vectorpositions[e][mk]] -= 1
        if v[vectorpositions[e][mk]] > -1:
            v[vectorpositions[e][iso]] += 1
            tuplevec = tuple(v)
            if tuplevec not in vectorsets[e]:
                vectorsets[e].add(tuplevec)
                pn = 0
                newelementmass = 0
                newelementprob = 1
                for n, c in enumerate(v):
                    loopiso = elementpositions[e][n]
                    newelementmass += elementalmasses[loopiso] * c
                    newelementprob *= elementalprobabilities[loopiso]**c
                    if loopiso in nonmonoisotopicelements:
                        newelementprob *= math.comb(acount-pn, c)
                        pn += c
                heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
    
    cutoff = -maxprob * dividingthreshold

    r, p, m, e, v = heapq.heappop(mainheap)
    elementlist.append([r, p, m, e, v])
    if len(isotopesbyelement[e]) > 2:
        while p > cutoff:
            for iso in nonmonoisotopicgroups[e]:
                newelementvector = v.copy()
                newelementvector[vectorpositions[e][mk]] -= 1
                if newelementvector[vectorpositions[e][mk]] > 0:
                    newelementvector[vectorpositions[e][iso]] += 1
                    tuplevec = tuple(newelementvector)
                    if tuplevec not in vectorsets[e]:
                        vectorsets[e].add(tuplevec)
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(newelementvector):
                            loopiso = elementpositions[e][n]
                            newelementmass += elementalmasses[loopiso] * c
                            newelementprob *= elementalprobabilities[loopiso]**c
                            if loopiso in nonmonoisotopicelements:
                                newelementprob *= math.comb(acount-pn, c)
                                pn += c
                        heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, newelementvector.copy()])
            r, p, m, e, v = heapq.heappop(mainheap)
            elementlist.append([r, p, m, e, v])
    else:
        iso = nonmonoisotopicgroups[e][0]
        while p > cutoff:
            nvector = v.copy()
            nvector[vectorpositions[e][mk]] -= 1
            if nvector[vectorpositions[e][mk]] > 0:
                nvector[vectorpositions[e][iso]] += 1
                tuplevec = tuple(nvector)
                if tuplevec not in vectorsets[e]:
                    vectorsets[e].add(tuplevec)
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(nvector):
                        loopiso = elementpositions[e][n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= elementalprobabilities[loopiso]**c
                        if loopiso in nonmonoisotopicelements:
                            newelementprob *= math.comb(acount-pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, nvector.copy()])
            r, p, m, e, v = heapq.heappop(mainheap)
            elementlist.append([r, p, m, e, v])
    heapq.heapify(elementlist)
    return elementlist

def descending_partial_products(dividingthreshold, elementalorganizer):
    for k in elementalorganizer:
        heapq.heapify(elementalorganizer[k])
    
    mainpool = defaultdict(list) #things already popped from elementalorganizer
    for k in elementalorganizer:
        mainpool[k].append(heapq.heappop(elementalorganizer[k]))
    
    formula = ''
    maxprob = 1
    mainmass = 0
    finalabundances = {} #subformula: prob
    for b in sorted(mainpool):
        for r, p, m, e, v in mainpool[b]:
            for n, c in enumerate(v):
                if c > 0:
                    formula += f'{elementpositions[e][n]}({c})'
            maxprob *= p
            mainmass += m

    finalabundances[formula] = [mainmass, maxprob]

    cutoff = maxprob * dividingthreshold
    mainheap = list(itertools.chain(*elementalorganizer.values()))
    heapq.heapify(mainheap)

    multinomialpath = [] #sublists not in mainpool
    probabilityranking = [] #representative lists of ratio probability to sort multinomialpath
    while mainheap:
        r, p, m, e, v = heapq.heappop(mainheap)
        baseiter = {k: v for k, v in mainpool.items() if k != e}
        baseiter[e] = [[r, p, m, e, v]]
        
        formula = ''
        prob = 1
        mass = 0
        for b in sorted(baseiter):
            for sr, sp, sm, se, sv in baseiter[b]:
                for n, c in enumerate(sv):
                    if c > 0:
                        formula += f'{elementpositions[se][n]}({c})'
                prob *= sp
                mass += sm
        
        finalabundances[formula] = [mass, prob]
        if prob < cutoff:
            break
        
        ind = bisect.bisect(probabilityranking, r)
        probabilityranking.insert(ind, r)
        multinomialpath.insert(ind, [r, p, m, e, v])
        
        checkedcombos = set()
        for path in multinomialpath.copy():
            multielement = False
            match path[1]:
                case list():
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
                                    sef += f'{elementpositions[se][n]}({c})'
                            seformulas.append(sef)
                            multipath.append([sr, sp, sm, se, sv])
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
            if -newratio >= dividingthreshold:
                if multielement:
                    seformula = ''
                    newprob = 1
                    newmass = 0
                    newiter = {k: v for k, v in baseiter.items() if k not in sepool}
                    newiter[e] = [[r, p, m, e, v]]
                    for ir, ip, im, ie, iv in multipath:
                        newiter[ie] = [[ir, ip, im, ie, iv]]
                    for b in sorted(newiter):
                        for ir, ip, im, ie, iv in newiter[b]:
                            for n, c in enumerate(iv):
                                if c > 0:
                                    seformula += f'{elementpositions[ie][n]}({c})'
                            newprob *= ip
                            newmass += im
                else:
                    newiter = {k: v for k, v in baseiter.items() if k != se}
                    newiter[se] = [[sr, sp, sm, se, sv]]
                    seformula = ''
                    newprob = 1
                    newmass = 0
                    for b in sorted(newiter):
                        for ir, ip, im, ie, iv in newiter[b]:
                            for n, c in enumerate(iv):
                                if c > 0:
                                    seformula += f'{elementpositions[ie][n]}({c})'
                            newprob *= ip
                            newmass += im
                if newprob >= cutoff:
                    finalabundances[seformula] = [newmass, newprob]
                    if multielement:
                        ind = bisect.bisect(probabilityranking, newratio)
                        probabilityranking.insert(ind, newratio)
                        multinomialpath.insert(ind, [newratio, *multipath])
                    else:
                        ind = bisect.bisect(probabilityranking, newratio)
                        probabilityranking.insert(ind, newratio)
                        multinomialpath.insert(ind, [newratio, [sr, sp, sm, se, sv], [r, p, m, e, v]])
            else:
                break

    subformulas, massesandabundances = list(zip(*finalabundances.items()))
    subformulas = np.array(subformulas, dtype='S')
    massesandabundances = np.array(massesandabundances).transpose()
    #sorting by mass
    subformulas = subformulas[massesandabundances[0].argsort()].tolist()
    massesandabundances = massesandabundances[:,massesandabundances[0].argsort()]
    return subformulas, massesandabundances

def distribution_generation(dividingthreshold, atomiccomposition):
    elementalorganizer = {} #element: [[preheaps]]
    for e, acount in atomiccomposition.items():
        elementstring = e + str(acount)
        elementlist = individual_element_binomial_walk(dividingthreshold, e, acount)
        elementalorganizer[e] = elementlist.copy() #don't need to copy the insides
    subformulas, massesandabundances = descending_partial_products(dividingthreshold, elementalorganizer)
    return subformulas, massesandabundances

def fragment_element_binomial_walk(dividingthreshold, e, acount, fragprobabilities):
    #elementalorganizer = defaultdict(list)
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
    #for e, acount in atomiccomposition.items():
    #mk = monoisotopickeys[e]
    #nvector = elementvectors[e].copy()
    #nvector[nvectorpositions[e][mk]] += acount
    nvector[fragmentvectorpositions[mk]] += acount
    #if len(isotopesbyelement[e]) > 2:
    flen = len(fragprobabilities)
    if flen > 2:
        baseprob = fragprobabilities[mk] ** acount
        preheap = []
        preheap.append([baseprob, acount * elementalmasses[mk], e, nvector.copy()])
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
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= fragprobabilities[loopiso]**c
                        #if loopiso in nonmonoisotopicelements:
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
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= fragprobabilities[loopiso]**c
                        #if loopiso in nonmonoisotopicelements:
                        if n > 0:
                            newelementprob *= math.comb(acount-pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
    else:
        preheap = []
        baseprob = fragprobabilities[mk] ** acount
        preheap.append([baseprob, acount * elementalmasses[mk], e, nvector.copy()])
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
                    newelementmass += elementalmasses[loopiso] * c
                    newelementprob *= fragprobabilities[loopiso]**c
                    #if loopiso in nonmonoisotopicelements:
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
                    newelementmass += elementalmasses[loopiso] * c
                    newelementprob *= fragprobabilities[loopiso]**c
                    #if loopiso in nonmonoisotopicelements:
                    if n > 0:
                        newelementprob *= math.comb(acount-pn, c)
                        pn += c
                heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
    
    cutoff = -maxprob * dividingthreshold

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
                            newelementmass += elementalmasses[loopiso] * c
                            newelementprob *= fragprobabilities[loopiso]**c
                            #if loopiso in nonmonoisotopicelements:
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
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= fragprobabilities[loopiso]**c
                        #if loopiso in nonmonoisotopicelements:
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

def fragment_descending_partial_products(dividingthreshold, elementalorganizer, fragmentpositions):
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
                    massnumber += int(iso[1:]) * c
                    formula += f'{iso}({c})'
            maxprob *= p
            mainmass += m

    massnumberindices[massnumber] = 0
    subformulas.append(formula)
    sumabundances.append([mainmass * maxprob, maxprob])
    
    cutoff = maxprob * dividingthreshold
    mainheap = list(itertools.chain(*elementalorganizer.values()))
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
                        massnumber += int(iso[1:]) * c
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
            ind = bisect.bisect(probabilityranking, r)
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
            if -newratio >= dividingthreshold:
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
                                    newmassnum += int(iso[1:]) * c
                                    seformula += f'{iso}({c})'
                                    #seformula += f'{fragmentpositions[ie][n]}({c})'
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
                                    newmassnum += int(iso[1:]) * c
                                    seformula += f'{iso}({c})'
                                    #seformula += f'{fragmentpositions[ie][n]}({c})'
                            newprob *= ip
                            newmass += im
                if newprob >= cutoff:
                    #finalabundances[seformula] = [newmass, newprob]
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
                        ind = bisect.bisect(probabilityranking, newratio)
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
                            ind = bisect.bisect(probabilityranking, newratio)
                            probabilityranking.insert(ind, newratio)
                            multinomialpath.insert(ind, (newratio, *newmulti))
            else:
                break

    #subformulas, massesandabundances = list(zip(*finalabundances.items()))
    subformulas = np.array(subformulas, dtype='S')
    #massesandabundances = np.array(massesandabundances).transpose()
    #massesandabundances = np.array(massesandabundances) #this is smaller in memory footprint somehow when its [[m, a], [m, a]] rather than [[m, m], [a, a]] for larger arrays
    massesandabundances = np.array(sumabundances)
    massesandabundances[:,0] /= massesandabundances[:,1]
    #sorting by intensity
    subformulas = subformulas[massesandabundances[:,1].argsort()[::-1]].tolist()
    massesandabundances = massesandabundances[massesandabundances[:,1].argsort()[::-1]]
    return subformulas, massesandabundances

def fragmentation_compositions(aminocomps, nfrags, cfrags, seq):
    fragments = {}

    #calculate the compositions of the n-term fragments
    fragcomp_n = {}
    for n, aa in enumerate(seq[:-1]):  
        aa_composition = aminocomps[aa]
        for k in aa_composition:
            fragcomp_n[k] = fragcomp_n.get(k, 0) + aa_composition.get(k, 0)
        for ion, modcomp in nfrags.items():
            fragment_composition = fragcomp_n.copy()
            for k in modcomp:
                fragment_composition[k] = fragment_composition.get(k, 0) + modcomp.get(k, 0)
            fragments[ion + str(n + 1)] = fragment_composition

    #calculate the compositions of the c-term fragments
    fragcomp_c = {}
    for n, aa in enumerate(seq[::-1][:-1]): 
        aa_composition = aminocomps[aa]
        for k in aa_composition:
            fragcomp_c[k] = fragcomp_c.get(k, 0) + aa_composition.get(k, 0)
        for ion, modcomp in cfrags.items():
            fragment_composition = fragcomp_c.copy()
            for k in modcomp:
                fragment_composition[k] = fragment_composition.get(k, 0) + modcomp.get(k, 0)
            fragments[ion + str(n + 1)] = fragment_composition
    
    #aa = seq[0]
    #aa_composition = aminocomps[aa]
    #for k in aa_composition:
    #    fragcomp_c[k] = fragcomp_c.get(k, 0) + aa_composition.get(k, 0)
    #fragcomp_c['H'] += 2
    #fragcomp_c['O'] += 1
    #fragments['precursor'] = fragcomp_c
    
    return fragments

def intersection_merge(mergable_items):
    sn = 0
    groupsofitems = {} #member: group
    itemgroups = defaultdict(set) #group: [members]
    for items in mergable_items:
        locs = set()
        for i in items:
            if i in groupsofitems:
                locs.add(groupsofitems[i])
        if locs:
            joiner = min(locs)
            if len(locs) > 1:
                for oldloc in locs.difference([joiner]):
                    for ol in itemgroups[oldloc]:
                        groupsofitems[ol] = joiner
                    itemgroups[joiner].update(itemgroups.pop(oldloc))
        else:
            joiner = sn
            sn += 1
        itemgroups[joiner].update(items)
        for i in items:
            groupsofitems[i] = joiner
    return list(itemgroups.values())

def shuffle_string(s):
    char_list = list(s)
    random.shuffle(char_list)
    return ''.join(char_list)

originalseq = ''.join((random.choices(list(aminoacidcomposition.keys()), k=peptidelength)))

atomiccomposition = Counter()
for aa in originalseq:
    atomiccomposition += aminoacidcomposition[aa]
#no OH loss on last residue, no H lost on first residue
atomiccomposition['H'] += 2
atomiccomposition['O'] += 1
formulastring = ''.join((''.join((k, str(v))) for k, v in atomiccomposition.items()))

subformulas, massesandabundances = distribution_generation(dividingthreshold, atomiccomposition)
subformulas = [i.decode() for i in subformulas]

massgroups = defaultdict(list) #massnumber: [masses]
intensitygroups = defaultdict(list) #massnumber: [abundances]
masses, abundances = massesandabundances
for n, s in enumerate(subformulas):
    massnumber = 0
    for ss in s.split(')')[:-1]:
        i1, i2 = map(int, ss[1:].split('('))
        massnumber += i1 * i2
    massgroups[massnumber].append(masses[n])
    intensitygroups[massnumber].append(abundances[n])

meansofmasses = []
sumsofabundances = []
subisodepthindices = [] #coordinates that reset for each proton location
condensationindices = []
for mn, m in massgroups.items():
    condensationindices.append(len(m))
    a = intensitygroups[mn]
    totalabundance = sum(a)
    weightedmass = 0
    cumulativeabundance = 0
    subinds = []
    #cumulatively adding intensities to determine which subisos can be added for ms2 searches
    subisocontinuation = True
    for n, (sm, sa) in sorted(enumerate(zip(m, a)), key=lambda x: -x[1][1]):
        weightedmass += sm * sa
        cumulativeabundance += sa
        cumpercent = cumulativeabundance / totalabundance
        if cumpercent <= subisotopomericdepth:
            subinds.append(n)
        elif subisocontinuation and cumpercent >= subisotopomericdepth:
            #subisodepth threshold breached, end it with this one
            subinds.append(n)
            subisocontinuation = False
        #    break
        #NO, stupid, you stopped the weighted mass means from being right by breaking here
    meansofmasses.append(weightedmass / totalabundance)
    sumsofabundances.append(totalabundance)
    subisodepthindices.append(tuple(subinds))
sumabundancedist = np.array([meansofmasses, sumsofabundances])

conlengths = np.array(condensationindices)
conends = conlengths.cumsum()
constarts = conends - conlengths

ionlist = list(ions)
ndict = {k: nfragmentcompositions[k] for k in ionlist if k in nfragmentcompositions}
cdict = {k: cfragmentcompositions[k] for k in ionlist if k in cfragmentcompositions}

seq = originalseq
seqfragments = {} #seq: [[subformula, mass, intensity]]
for _ in range(nseqs):
    fragments = fragmentation_compositions(aminoacidcomposition, ndict, cdict, seq)
    
    fragmasses = defaultdict(list) #subformula: [masses]
    #fragintensities = []
    #fragmentsubformulas = defaultdict(list) #subformula: [ions]
    for n, bi in enumerate(constarts.tolist()):
        for sq in subisodepthindices[n]:
            subindex = bi + sq
            subformula = subformulas[subindex]
            isocounts = set()
            competing = set()
            competitors = {}
            isosums = {}
            for ss in subformula.split(')')[:-1]:
                iso, c = ss.split('(')
                c = int(c)
                e = iso[0]
                if e in isocounts:
                    competing.add(e)
                    competitors[e][iso] = c
                    isosums[e] += c
                else:
                    isocounts.add(e)
                    competitors[e] = {iso: c}
                    isosums[e] = c
            isoprobs = defaultdict(dict) #e: probs of that element
            for e, v in competitors.items():
                if e in competing:
                    for iso, c in v.items():
                        prob = c / isosums[e]
                        isoprobs[e][iso] = prob
                else:
                    for iso in v:
                        isoprobs[e][iso] = 1
            for ion, fragcomp in fragments.items():
                elementalorganizer = {} #element: [[iso heaps]]
                fragmentpositions = {} #element: position: iso
                for e, c in fragcomp.items():
                    if e in isoprobs: #element is in subformula
                        fragprobs = isoprobs[e]
                        if len(fragprobs) > 1:
                            elementlist, positions = fragment_element_binomial_walk(dividingthreshold, e, c, fragprobs)
                            elementalorganizer[e] = elementlist.copy()
                            fragmentpositions[e] = positions
                        else: #no need for cache, only 1 iso
                            iso = list(fragprobs)[0]
                            elementalorganizer[e] = [[-1, 1, elementalmasses[iso]*c, e, [c]]]
                            fragmentpositions[e] = {0: iso}
                fragformulas, fragmassesandabundances = fragment_descending_partial_products(dividingthreshold, elementalorganizer, fragmentpositions)
                #simplifying this for now bc idc about the entire dists for this exploration
                mainmass, mainintensity = fragmassesandabundances[0]
                #fragmasses[subformula].append([mainmass, ion, fragformulas[0]])
                fragmasses[subformula].append([mainmass, ion])
                #fragintensities.append(mainintensity)
                #fragmentsubformulas[subformula].append(fragformulas[0])
    #seqfragments[seq] = list(zip(fragmentsubformulas, fragmasses))
    seqfragments[seq] = fragmasses
    seq = shuffle_string(originalseq)

#from the original seq pick 3-8 fragment ions and work off of those, make the masses
#then from all the matching ms1 dists, find every potential seq
#match all potential seq fragments to the masses you found, stick to just the most abundant isotopomer like 90% of the time or something

def coverage_print(seq, ionlist):
    olen = len(seq)
    maxncoverage = 0
    maxccoverage = 0
    dividers = []
    coverage = []
    ntermcoverage = []
    ctermcoverage = []
    for ion in ionlist:
        iontype = ion[0]
        ioncount = int(ion[1:])
        if iontype in 'abc': #nterm
            dividers.append(ioncount)
            pseq = seq[:ioncount]
            coverage.append(pseq + '_' * (olen - ioncount))
            ntermcoverage.append(ioncount)
            if ioncount > maxccoverage:
                maxccoverage = ioncount
        elif iontype in 'xyz': #cterm
            dividers.append(olen - ioncount)
            pseq = seq[olen-ioncount:]
            coverage.append((olen - ioncount) * '_' + pseq)
            ctermcoverage.append(olen-ioncount)
            if ioncount > maxncoverage:
                maxncoverage = ioncount
    dividers = sorted(set(dividers))

    #isolation counts need to be robust against redundant pseqs
    ind = 0
    ddiff = np.diff(dividers, prepend=0).tolist()
    dividerstring = ''
    partialseqs = defaultdict(int) #index-pseq: count #the index safeguards against multiple isolations of the same partial sequence, the defaultdict rather than a Counter keeps the keys in order so i can view it easier in regards to the order of the sequence
    for d in ddiff:
        pseq = seq[ind:ind+d]
        dividerstring += pseq + '|'
        ntermcovers = sum(1 for i in ntermcoverage if i > ind)
        ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
        covers = ntermcovers + ctermcovers
        if covers > 0:
            partialseqs[str(ind) + '-' + pseq] += covers
        ind += d
    pseq = seq[ind:]
    dividerstring += pseq
    ntermcovers = sum(1 for i in ntermcoverage if i > ind)
    ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
    covers = ntermcovers + ctermcovers
    if covers > 0:
        partialseqs[str(ind) + '-' + pseq] += covers
    for c in coverage:
        print(c)
    print(dividerstring)


maxloc = sumabundancedist[1].argmax()
maxsubformulas = [subformulas[i+constarts[maxloc]] for i in subisodepthindices[maxloc]]

mainsubformula = maxsubformulas[0]
mainpicks = seqfragments[originalseq][mainsubformula]

pickedindices = np.unique(np.random.randint(0, len(mainpicks)-1, size=npicks)) #excluding precursor
picks = [mainpicks[i] for i in pickedindices]
intensityvals = np.random.uniform(minintensity, maxintensity, size=len(picks))

#omasses, oions, oformulas = zip(*picks)
omasses, oions = zip(*picks)
omasses = np.array(omasses)[:,None]
model = spatial.KDTree(omasses)

nomatches = 0
finalmetrics = {} #seq: [metrics]
secondfinalmetrics = {} #seq: score
thirdfinalmetrics = {} #seq: score
finalions = {}
seqcoverage = {} #seq: [ion coverage]
for seq, isubformulas in seqfragments.items():
    #masses, fragions, formulas = zip(*seqfragments[seq][mainsubformula])
    masses, fragions = zip(*seqfragments[seq][mainsubformula])
    masses = np.array(masses)[:-1,None]
    radius = (masses / 1000000 * ppmtol).flatten()
    matches = model.query_ball_point(masses, radius).tolist()
    coverage = []
    #matchedintensities = []
    for n, m in enumerate(matches):
        if m:
            coverage.append(fragions[n])
            #matchedintensities.append(intensityvals[m[0]])
    if coverage:
        print(seq)
        coverage_print(seq, coverage)
        seqcoverage[seq] = coverage
        print('~~~~~~~~~~')
        #matchedsum = sum(matchedintensities)
        #maxintensity = max(matchedintensities)
        #intensitycoverage = 1 / (((matchedsum - maxintensity) / matchedsum) * (len(matchedintensities) - 1))
        slen = len(seq)
        maxncoverage = 0
        maxccoverage = 0
        dividers = []
        ntermcoverage = []
        ctermcoverage = []
        for ion in coverage:
            iontype = ion[0]
            ioncount = int(ion[1:])
            if iontype in 'abc': #nterm
                dividers.append(ioncount)
                pseq = seq[:ioncount]
                ntermcoverage.append(ioncount)
                if ioncount > maxccoverage:
                    maxccoverage = ioncount
            elif iontype in 'xyz': #cterm
                dividers.append(slen - ioncount)
                pseq = seq[slen-ioncount:]
                ctermcoverage.append(slen-ioncount)
                if ioncount > maxncoverage:
                    maxncoverage = ioncount
        dividers = sorted(set(dividers))
        coverageweight = 1 / (maxncoverage + maxccoverage)

        #isolation counts need to be robust against redundant pseqs
        ind = 0
        ddiff = np.diff(dividers, prepend=0).tolist()
        #dividerstring = ''
        partialseqs = defaultdict(int) #index-pseq: count #the index safeguards against multiple isolations of the same partial sequence, the defaultdict rather than a Counter keeps the keys in order so i can view it easier in regards to the order of the sequence
        #isolationintensities = defaultdict(list) #index-pseq: [individual ion intensities]
        for d in ddiff:
            pseq = seq[ind:ind+d]
            #dividerstring += pseq + '|'
            ntermcovers = [i for i in ntermcoverage if i > ind]
            ctermcovers = [i for i in ctermcoverage if i <= ind]
            covers = len(ntermcovers) + len(ctermcovers)
            if covers > 0:
                label = str(ind) + '-' + pseq
                partialseqs[label] += covers
                #isolationintensities[label] = [i[1] for i in ntermcovers] + [i[1] for i in ctermcovers]
            ind += d
        pseq = seq[ind:]
        #dividerstring += pseq
        ntermcovers = [i for i in ntermcoverage if i > ind]
        ctermcovers = [i for i in ctermcoverage if i <= ind]
        covers = len(ntermcovers) + len(ctermcovers)
        if covers > 0:
            label = str(ind) + '-' + pseq
            partialseqs[label] += covers
            #isolationintensities[label] = [i[1] for i in ntermcovers] + [i[1] for i in ctermcovers]
        pairsum = 0
        #uniquesum = sum(set((itertools.chain.from_iterable(isolationintensities.values()))))
        matchcounts = len(set(coverage))
        isolationlengthweight = 1
        for indseq, count in partialseqs.items():
            ind, pseq = indseq.split('-')
            ind = int(ind)
            #i could use this index to weight based on distance from the ends i guess?
            #isolationlengthweight *= len(pseq) / len(seq)
            isolationlengthweight *= 1 / len(pseq) / len(seq) #this 1 / provides an additional layer of geometric success to this scheme, grants success where there was previously failure
            #if len(isolationintensities[indseq]) > 1:
            #    for il, ir in itertools.combinations(isolationintensities[indseq], 2):
            #        #add pairs
            #        pairsum += il + ir
            #else:
            #    #add single i guess
            #    pairsum += isolationintensities[indseq][0]
            plen = len(pseq)
            matchcounts += plen * count
        dividerweight = 1 / len(partialseqs)
        finalmetrics[seq] = [dividerweight * isolationlengthweight * coverageweight, dividerweight, isolationlengthweight, coverageweight, len(partialseqs), len(seq)]
        finalions[seq] = [i for i in coverage]
        secondfinalmetrics[seq] = matchcounts
        thirdfinalmetrics[seq] = 1 / (dividerweight + isolationlengthweight + coverageweight)
print(nomatches, 'nonmatches')

print(oions)

olen = len(originalseq)
maxncoverage = 0
maxccoverage = 0
dividers = []
coverage = []
ntermcoverage = []
ctermcoverage = []
for ion in oions:
    iontype = ion[0]
    ioncount = int(ion[1:])
    if iontype in 'abc': #nterm
        dividers.append(ioncount)
        pseq = originalseq[:ioncount]
        coverage.append(pseq + '_' * (olen - ioncount))
        ntermcoverage.append(ioncount)
        if ioncount > maxccoverage:
            maxccoverage = ioncount
    elif iontype in 'xyz': #cterm
        dividers.append(olen - ioncount)
        pseq = originalseq[olen-ioncount:]
        coverage.append((olen - ioncount) * '_' + pseq)
        ctermcoverage.append(olen-ioncount)
        if ioncount > maxncoverage:
            maxncoverage = ioncount
dividers = sorted(set(dividers))

#isolation counts need to be robust against redundant pseqs
ind = 0
ddiff = np.diff(dividers, prepend=0).tolist()
dividerstring = ''
partialseqs = defaultdict(int) #index-pseq: count #the index safeguards against multiple isolations of the same partial sequence, the defaultdict rather than a Counter keeps the keys in order so i can view it easier in regards to the order of the sequence
for d in ddiff:
    pseq = originalseq[ind:ind+d]
    dividerstring += pseq + '|'
    ntermcovers = sum(1 for i in ntermcoverage if i > ind)
    ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
    covers = ntermcovers + ctermcovers
    if covers > 0:
        partialseqs[str(ind) + '-' + pseq] += covers
    ind += d
pseq = originalseq[ind:]
dividerstring += pseq
ntermcovers = sum(1 for i in ntermcoverage if i > ind)
ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
covers = ntermcovers + ctermcovers
if covers > 0:
    partialseqs[str(ind) + '-' + pseq] += covers

print('\n')
for c in coverage:
    print(c)
print(dividerstring)

print(partialseqs)
print('\n')

sortedmetrics = sorted(finalmetrics.items(), key=lambda x: x[1][0])
bestseq = sortedmetrics[0][0]
if bestseq == originalseq:
    if len(sortedmetrics) > 1:
        if sortedmetrics[0][1][0] == sortedmetrics[1][1][0]:
            print('success but tied!')
            print(sortedmetrics[1][0])
            coverage_print(sortedmetrics[1][0], seqcoverage[sortedmetrics[1][0]])
            psum = sum(1 for i in sortedmetrics if i[1] == sortedmetrics[0][1])
            print(f'{psum} / {len(sortedmetrics)} tied')
        else:
            print('success!')
        print(sortedmetrics[0][1][0], len(finalions[bestseq]))
        if len(sortedmetrics) > 1:
            print(sortedmetrics[1][1][0], len(finalions[sortedmetrics[1][0]]))
    else:
        print('success!')
        print(sortedmetrics[0][1][0], len(finalions[bestseq]))
else:
    print('FAILURE :(')
    for n, s in enumerate(sortedmetrics):
        if s[0] == originalseq:
            n += 1
            break
    if n > 10:
        print('big failure')
    else:
        for s in sortedmetrics[:n]:
            print(s[0], s[1][0], ' '.join((map(str, map(lambda x: round(x, 3), s[1:-2])))), len(finalions[s[0]]))
        seq = sortedmetrics[0][0]
        masses, fragions = zip(*seqfragments[seq][mainsubformula])
        masses = np.array(masses)[:-1,None]
        radius = (masses / 1000000 * ppmtol).flatten()
        matches = model.query_ball_point(masses, radius).tolist()
        nseqcoverage = []
        #matchedintensities = []
        for n, m in enumerate(matches):
            if m:
                nseqcoverage.append(fragions[n])
                #matchedintensities.append(intensityvals[m[0]])
        slen = len(seq)
        maxncoverage = 0
        maxccoverage = 0
        dividers = []
        coverage = []
        ntermcoverage = []
        ctermcoverage = []
        for ion in nseqcoverage:
            iontype = ion[0]
            ioncount = int(ion[1:])
            if iontype in 'abc': #nterm
                dividers.append(ioncount)
                pseq = seq[:ioncount]
                coverage.append(pseq + '_' * (slen - ioncount))
                ntermcoverage.append(ioncount)
                if ioncount > maxccoverage:
                    maxccoverage = ioncount
            elif iontype in 'xyz': #cterm
                dividers.append(slen - ioncount)
                pseq = seq[slen-ioncount:]
                coverage.append((slen - ioncount) * '_' + pseq)
                ctermcoverage.append(slen-ioncount)
                if ioncount > maxncoverage:
                    maxncoverage = ioncount
        dividers = sorted(set(dividers))

        #isolation counts need to be robust against redundant pseqs
        ind = 0
        ddiff = np.diff(dividers, prepend=0).tolist()
        dividerstring = ''
        partialseqs = defaultdict(int) #index-pseq: count #the index safeguards against multiple isolations of the same partial sequence, the defaultdict rather than a Counter keeps the keys in order so i can view it easier in regards to the order of the sequence
        for d in ddiff:
            pseq = seq[ind:ind+d]
            dividerstring += pseq + '|'
            ntermcovers = sum(1 for i in ntermcoverage if i > ind)
            ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
            covers = ntermcovers + ctermcovers
            if covers > 0:
                partialseqs[str(ind) + '-' + pseq] += covers
            ind += d
        pseq = seq[ind:]
        dividerstring += pseq
        ntermcovers = sum(1 for i in ntermcoverage if i > ind)
        ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
        covers = ntermcovers + ctermcovers
        if covers > 0:
            partialseqs[str(ind) + '-' + pseq] += covers
        
        print('\n')
        for c in coverage:
            print(c)
        print('\n')
        print(dividerstring)
        
        print(partialseqs)

print('~~~~~~~~~~')
print('second try!')
sortedmetrics = sorted(secondfinalmetrics.items(), key=lambda x: x[1], reverse=True)
bestseq = sortedmetrics[0][0]
if bestseq == originalseq:
    if len(sortedmetrics) > 1:
        if sortedmetrics[0][1] == sortedmetrics[1][1]:
            print('success but tied!')
            print(sortedmetrics[1][0])
            coverage_print(sortedmetrics[1][0], seqcoverage[sortedmetrics[1][0]])
            psum = sum(1 for i in sortedmetrics if i[1] == sortedmetrics[0][1])
            print(f'{psum} / {len(sortedmetrics)} tied')
        else:
            print('success!')
        print(sortedmetrics[0][1], len(finalions[bestseq]))
        if len(sortedmetrics) > 1:
            print(sortedmetrics[1][1], len(finalions[sortedmetrics[1][0]]))
    else:
        print('success!')
        print(sortedmetrics[0][1], len(finalions[bestseq]))
else:
    print('FAILURE :(')
    for n, s in enumerate(sortedmetrics):
        if s[0] == originalseq:
            n += 1
            break
    if n > 10:
        print('big failure')
    else:
        for s in sortedmetrics[:n]:
            print(s[0], s[1], len(finalions[s[0]]))
        seq = sortedmetrics[0][0]
        masses, fragions = zip(*seqfragments[seq][mainsubformula])
        masses = np.array(masses)[:-1,None]
        radius = (masses / 1000000 * ppmtol).flatten()
        matches = model.query_ball_point(masses, radius).tolist()
        nseqcoverage = []
        #matchedintensities = []
        for n, m in enumerate(matches):
            if m:
                nseqcoverage.append(fragions[n])
                #matchedintensities.append(intensityvals[m[0]])
        slen = len(seq)
        maxncoverage = 0
        maxccoverage = 0
        dividers = []
        coverage = []
        ntermcoverage = []
        ctermcoverage = []
        for ion in nseqcoverage:
            iontype = ion[0]
            ioncount = int(ion[1:])
            if iontype in 'abc': #nterm
                dividers.append(ioncount)
                pseq = seq[:ioncount]
                coverage.append(pseq + '_' * (slen - ioncount))
                ntermcoverage.append(ioncount)
                if ioncount > maxccoverage:
                    maxccoverage = ioncount
            elif iontype in 'xyz': #cterm
                dividers.append(slen - ioncount)
                pseq = seq[slen-ioncount:]
                coverage.append((slen - ioncount) * '_' + pseq)
                ctermcoverage.append(slen-ioncount)
                if ioncount > maxncoverage:
                    maxncoverage = ioncount
        dividers = sorted(set(dividers))

        #isolation counts need to be robust against redundant pseqs
        ind = 0
        ddiff = np.diff(dividers, prepend=0).tolist()
        dividerstring = ''
        partialseqs = defaultdict(int) #index-pseq: count #the index safeguards against multiple isolations of the same partial sequence, the defaultdict rather than a Counter keeps the keys in order so i can view it easier in regards to the order of the sequence
        for d in ddiff:
            pseq = seq[ind:ind+d]
            dividerstring += pseq + '|'
            ntermcovers = sum(1 for i in ntermcoverage if i > ind)
            ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
            covers = ntermcovers + ctermcovers
            if covers > 0:
                partialseqs[str(ind) + '-' + pseq] += covers
            ind += d
        pseq = seq[ind:]
        dividerstring += pseq
        ntermcovers = sum(1 for i in ntermcoverage if i > ind)
        ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
        covers = ntermcovers + ctermcovers
        if covers > 0:
            partialseqs[str(ind) + '-' + pseq] += covers
        
        print('\n')
        for c in coverage:
            print(c)
        print('\n')
        print(dividerstring)
        
        print(partialseqs)

print('~~~~~~~~~~')
print('third try!')
sortedmetrics = sorted(thirdfinalmetrics.items(), key=lambda x: x[1], reverse=True)
bestseq = sortedmetrics[0][0]
if bestseq == originalseq:
    if len(sortedmetrics) > 1:
        if sortedmetrics[0][1] == sortedmetrics[1][1]:
            print('success but tied!')
            print(sortedmetrics[1][0])
            coverage_print(sortedmetrics[1][0], seqcoverage[sortedmetrics[1][0]])
            psum = sum(1 for i in sortedmetrics if i[1] == sortedmetrics[0][1])
            print(f'{psum} / {len(sortedmetrics)} tied')
        else:
            print('success!')
        print(sortedmetrics[0][1], len(finalions[bestseq]))
        if len(sortedmetrics) > 1:
            print(sortedmetrics[1][1], len(finalions[sortedmetrics[1][0]]))
    else:
        print('success!')
        print(sortedmetrics[0][1], len(finalions[bestseq]))
else:
    print('FAILURE :(')
    for n, s in enumerate(sortedmetrics):
        if s[0] == originalseq:
            n += 1
            break
    if n > 10:
        print('big failure')
    else:
        for s in sortedmetrics[:n]:
            print(s[0], s[1], len(finalions[s[0]]))
        seq = sortedmetrics[0][0]
        masses, fragions = zip(*seqfragments[seq][mainsubformula])
        masses = np.array(masses)[:-1,None]
        radius = (masses / 1000000 * ppmtol).flatten()
        matches = model.query_ball_point(masses, radius).tolist()
        nseqcoverage = []
        #matchedintensities = []
        for n, m in enumerate(matches):
            if m:
                nseqcoverage.append(fragions[n])
                #matchedintensities.append(intensityvals[m[0]])
        slen = len(seq)
        maxncoverage = 0
        maxccoverage = 0
        dividers = []
        coverage = []
        ntermcoverage = []
        ctermcoverage = []
        for ion in nseqcoverage:
            iontype = ion[0]
            ioncount = int(ion[1:])
            if iontype in 'abc': #nterm
                dividers.append(ioncount)
                pseq = seq[:ioncount]
                coverage.append(pseq + '_' * (slen - ioncount))
                ntermcoverage.append(ioncount)
                if ioncount > maxccoverage:
                    maxccoverage = ioncount
            elif iontype in 'xyz': #cterm
                dividers.append(slen - ioncount)
                pseq = seq[slen-ioncount:]
                coverage.append((slen - ioncount) * '_' + pseq)
                ctermcoverage.append(slen-ioncount)
                if ioncount > maxncoverage:
                    maxncoverage = ioncount
        dividers = sorted(set(dividers))

        #isolation counts need to be robust against redundant pseqs
        ind = 0
        ddiff = np.diff(dividers, prepend=0).tolist()
        dividerstring = ''
        partialseqs = defaultdict(int) #index-pseq: count #the index safeguards against multiple isolations of the same partial sequence, the defaultdict rather than a Counter keeps the keys in order so i can view it easier in regards to the order of the sequence
        for d in ddiff:
            pseq = seq[ind:ind+d]
            dividerstring += pseq + '|'
            ntermcovers = sum(1 for i in ntermcoverage if i > ind)
            ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
            covers = ntermcovers + ctermcovers
            if covers > 0:
                partialseqs[str(ind) + '-' + pseq] += covers
            ind += d
        pseq = seq[ind:]
        dividerstring += pseq
        ntermcovers = sum(1 for i in ntermcoverage if i > ind)
        ctermcovers = sum(1 for i in ctermcoverage if i <= ind)
        covers = ntermcovers + ctermcovers
        if covers > 0:
            partialseqs[str(ind) + '-' + pseq] += covers
        
        print('\n')
        for c in coverage:
            print(c)
        print('\n')
        print(dividerstring)
        
        print(partialseqs)

#things must be covered to be isolated
#more things covered -> better
#smaller isolation -> better
#isolated groups further away from either end are more valuable -> im against this now, it actually depends on the context of the other seqs in the library, some positions are more valuable relative to what it has to compete with
#more things covered -> better
#repetetive coverage -> good
#more ions covered during the isolation process -> better
    #ie is this a y1-y2 or a y1-b12?
    #im thinking to do combinatorics of all routes towards coverage
#confidence weighted by intensity
    #for coverage of an isolation -> simply adding the intensity of each ion that covered it is fine
    #% of total scan intensity, as well as % of total intensity used for the calculation as the summed intensities will be used across groups, possibly even the sum of all the intensity pairs used across every group RELATED to that scan, or maybe even related to that peptide
    #and this is all the percentage to be related to MS1 percentages

#the order/basis of execution:
#find covered isolation groups -> their length is a weight for that group
#determine the number of times they're covered -> that's an inverse weight for that group
#determine the number of AAs in their coverage -> redundantly? -> that's a weight for those ions?
#sum unique intensities and combinatoric pair intensity -> assign combinatoric intensity sums via individual mass intensities
    #this operates on the levels of:
        #sum of unique intensities / total scan intensity
        #np.sqrt(sum of redundant intensity pairs / total pair intensity) used for IDing things in that scan
            #should the total pair sum here be from non-redundant sets of IDs for specific lines?
            #or should this be across every potential ID in a scan?
        #the sqrt here is because this is like a squared upper bound due to the redundancy? i think it might match better
        #^no i think you should count the number of total individual intensities used (not pairs) across the whole combinatoric, then divide each group by that
    #how to process them:
        #intensitycoverage
        #case: 1 ion takes up 20% of the total intensity of 9 ions -> 1/(.8*8)
        #case: 1 ion takes up 80% of the total intensity of 4 ions -> 1/(.2*3)
        #case: 1 ion match -> output 1 i guess?
            #1 isn't the highest possible score, but this should only happen when there's no other matches with more ions, this is basically going to be a shit match anyways
    #these are weights for that group


#for this simulation:
#do complete permutations of individual peptide
    #i dont feel the need to do this rn, i think this works optimall and i get the coverage thing:
    #basically, if you have P|E|PT|I|DE
    #then if you have P|E|TP|I|ED or other combinatoric variations, then you'll never be able to distinguish them
    #geometically, i'm maximized the ability of distinction using purely the sequence itself
#introduce background noise -> weight metrics by intensity
    #this will require an intensity concept
    #not too worried about this rn, moving on back to peptidefragmentscoring instead
#assess multiple different lengths and sequence compositions
    #basically what i found, which is just a duh moment, shorter sequences and less fragment ions pose less opportunity for specificity of identification
    #better ID chance even for small peptides when you have way more ions, its basically the only change you have at getting those ones
#maybe guarantee K/R endings too i guess later
    #idc about this
#eventually use linepositionsbyformula to match a wider variety of MS1 dists
    #nah

#for decoy generation:
#i think a neighboring swap/bridge shuffle, like in cards, is better than random shuffling, although id like to compare the 2
#maybe to explore this: do a FUCKTON of shuffling with a massively disproportionate amount of decoys, and get the levenshtein distance compared to the original seq -> which distances give the furthest separation? or any other measurement?
#PEPTIDE becomes EPTPDIE? there should be 2 different shuffles that are the closest and can be distinguished from what you see, if you can beat both of them its probably real
#decoy generation needs to take into consideration groups of peptides in the database that have the same AA composition -> and the decoys should look to distinguish both of these from something false
#^if there's only one of any composition, the decoy should aim to differentiate that one sequence from something false
#decoy generation should take matching distributions (from the database, so across size but matching mass/shape of the dists >>> BUT ONLY for peptide that are in the search in question aka the data)
    #i think this makes a specific problem of: then SPECIFIC ions are needed to prove something is correct
    #and i think id like to mitigate this problem by producing as many DIFFERENT ions of the same mass? it almost doesn't make sense
        #i think as the fragments are small (in the generated decoy) they should have the MOST similarity to the genuine sequence
        #and as the fragments grow towards the middle, and beyond -> they become MORE distinguishable from the actual sequence
        #so as fragments cover more ground, its more likely to brush off any decoys
        #this acts as a strategic mechanism towards coverage and identification itself
        #YES

#i think you want to generate decoys from a birds-eye view on just AAs and their masses alone -> as working off of fragment masses and specific isotopes would be too slow
#the idea is that you want multiple fragments generated from a decoy to be similar enough to fragments generated from the real sequences
#but what will vary is their coverage -> as this is a great avenue for differentiation
#a simple AA comp combination table will lead towards good good selection of how to distinguish this
#also i dont need to do this by mass as i can just do it for individual AAs using shuffling -> dont need MS1 dists PLUS it will always apply by mass and towards every mass match by default

#an exploration:
#heavy combinatorics in both sequence shuffling as well as ion selection of a single sequence in order to determine what actually makes a good decoy

#maybe you can simplify decoys down to decoy fragents?
#as a pseudo de-novo?

#big picture of this file:
#take the biggest group of matchers
#make all theirs dists and fragments
#pick "one true" thing and work from there

#the aim of this concept should be to show:
#how does an increasing/decreasing number of ions affect the specificity of the identification across different-lengthed peptides?
#i also want to explore how often different subisos of the same proton location offer different masses for the same ion
#does including all abc/xyz make this harder to distinguish?
