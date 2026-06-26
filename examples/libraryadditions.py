from elementalcomponents import elementinfo, aminoacidcomposition
from database import environment

from itertools import product, chain, combinations
from collections import Counter, defaultdict
import multiprocessing as mp
from bisect import bisect
from Bio import SeqIO
from math import comb
from time import time
import numpy as np
import random
import string
import heapq
import lmdb
import os

#make proteinseqs from fasta, decoys included
#digest proteinseqs -> return digestsofproteins
#apply switchers to digestsofproteins
# - make and save proteinsofseqs in proteome-proteinsofseqs
# - make flatseqs -> apply variable mods
# > return flatseqs
#put flatseqs into formula_consolidation
# - returns formulas, atomiccomps
# - make seqsbyformula and atomiccompositions -> add seqsbyformula to the database
# - add the list of formulas to the database as being for this proteome
#    > proteome-name: [list of formulas]
# > return atomiccompositions
#check if formulas can be pulled from the database
# - don't generate the ones that do, but retain the formulas
# - generate the isotopic distributions of the formulas
# - add the newly generated fullabundances and abundanceformulas to the db
# > return abundances and subformulas (list version of subformulas is ok)
#generate sumabundances, condensationcoordinates, and subisodepthqualifiers
#make upper/lower mass limits
#^and add these to the db

#so this is a bit slower than the last database, the multiprocessing gets weirder when you try to put everything in a class, but in terms of organizational steps i think this is acceptable
#unfortunately i can no longer import things like elementvectors to make them accessible across all functions more easily and have to refer to them via a self which might slow them down? idk, but whatever
#we've got more consistent modification terminology plus customizable enzymes now so that's a plus

class LibraryOrganizer:
    def __init__(self, librarylocation, proteomefile, minlength, maxlength, missedcleavages, nprocs, variablemodifications=[], staticmodifications=[], enzyme=-1, maxvmods=0):
        if enzyme < 0:
            print('you need to pick an enzyme!')
            self.passing = False
            return
        if set(variablemodifications).intersection(staticmodifications):
            #this doesn't cover for if multiple staticmods are used on the same AA
            print('you have overlapping static and variable modifications, straighten that out')
            self.passing = False
            return
        self.proteome = proteomefile.split('/')[-1].split('.')[0]
        badwords = set(['modifications', 'enzyme', 'defaults'])
        if self.proteome in badwords:
            print('your proteomes name interferes with the internal database, make it realistic')
            self.passing = False
            return
        
        self.passing = True
        self.librarylocation = librarylocation
        self.proteomefile = proteomefile
        self.minlength = minlength
        self.maxlength = maxlength
        self.missedcleavages = missedcleavages
        self.maxvmods = maxvmods
        self.nprocs = nprocs
        
        with environment(self.librarylocation) as env:
            while True:
                try:
                    staticmods = defaultdict(dict)
                    variablemods = defaultdict(dict)
                    enzymedb = env.open_db('enzymes'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(enzymedb) as cursor:
                            enzymename, (self.cutters, self.noncutters) = eval(cursor.get(str(enzyme).encode()).decode())
                    modificationdb = env.open_db('modifications'.encode())
                    encodedvarmods = [str(i).encode() for i in variablemodifications]
                    encodedstaticmods = [str(i).encode() for i in staticmodifications]
                    with env.begin(write=False) as txn:
                        with txn.cursor(modificationdb) as cursor:
                            for k, v in cursor.getmulti(encodedvarmods):
                                mod, aa, modtype = eval(v.decode())
                                variablemods[aa][mod] = modtype
                            for k, v in cursor.getmulti(encodedstaticmods):
                                mod, aa, modtype = eval(v.decode())
                                staticmods[aa][mod] = modtype
                    defaults = env.open_db('defaults'.encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(defaults) as cursor:
                            self.dividingthreshold = float(cursor.get('dividingthreshold'.encode()).decode())
                            self.subisotopomericdepth = float(cursor.get('subisotopomericdepth'.encode()).decode())
                            pcount = int(cursor.get('proteomes'.encode()).decode())
                            cursor.put('proteomes'.encode(), str(pcount + 1).encode())
                    allproteomesdb = env.open_db('proteomes'.encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(allproteomesdb) as cursor:
                            cursor.put(str(pcount + 1).encode(), self.proteome.encode())
                    proteomedb = env.open_db((self.proteome + '.info').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(proteomedb) as cursor:
                            cursor.put('variablemods'.encode(), str(variablemods).encode())
                            cursor.put('staticmods'.encode(), str(staticmods).encode())
                            cursor.put('enzyme'.encode(), str(enzymename).encode())
                            
                            self.aminoacidcomposition = aminoacidcomposition
                            newcharacters = string.ascii_lowercase
                            self.modifiers = defaultdict(set) #existing AA: [new AA letters]
                            modoriginals = {} #new AA letter: existing AA
                            newmods = {} #new AA letter: its atomic composition
                            vn = 0
                            #i need to print all mod/enzyme parameters here for log purposes
                            for aa, mods in variablemods.items():
                                for mod, modtype in mods.items():
                                    match modtype:
                                        case dict():
                                            #comp -> add the comp to old aacomp -> make a new AA
                                            compstring = ''.join((i[0] + f'({i[1]})' for i in sorted(modtype.items())))
                                            print(aa, mod, compstring)
                                            newaa = newcharacters[vn]
                                            self.modifiers[aa].add(newaa)
                                            newcomp = dict(Counter(modtype) + Counter(self.aminoacidcomposition[aa]))
                                            self.aminoacidcomposition[newaa] = newcomp
                                            modoriginals[newaa] = aa
                                            newmods[newaa] = newcomp
                                            vn += 1
                                        case str():
                                            #mass -> make a new atom -> make a new AA
                                            mass = float(modtype)
                                            print(aa, mod, mass)
                                            newaa = newcharacters[vn]
                                            vn += 1
                                            newlement = newcharacters[vn]
                                            vn += 1
                                            elementinfo[newelement] = {newelement: (mass, 1)}
                                            self.modifiers[aa].add(newaa)
                                            newcomp = dict(Counter(modtype) + Counter(self.aminoacidcomposition[aa]))
                                            newcomp = self.aminoacidcomposition[aa].copy()
                                            newcomp[newelement] = 1
                                            self.aminoacidcomposition[newaa] = newcomp
                                            modoriginals[newaa] = aa
                                            newmods[newaa] = newcomp
                            
                            for aa, mods in staticmods.items():
                                for mod, modtype in mods.items():
                                    match modtype:
                                        case dict():
                                            #comp -> simply add to the old AA
                                            compstring = ''.join((i[0] + f'({i[1]})' for i in sorted(modtype.items())))
                                            print(aa, mod, compstring)
                                            newcomp = dict(Counter(self.aminoacidcomposition[aa]) + Counter(modtype))
                                            self.aminoacidcomposition[aa] = newcomp
                                        case str():
                                            #mass -> add new atom and add to the old AA
                                            mass = float(modtype)
                                            print(aa, mod, mass)
                                            newelement = newcharacters[vn]
                                            elementinfo[newelement] = {newelement: (mass, 1)}
                                            self.aminoacidcomposition[aa][newelement] = 1
                                            vn += 1
                            
                            cursor.put('modifiers'.encode(), str(self.modifiers).encode())
                            cursor.put('modoriginals'.encode(), str(modoriginals).encode())
                            cursor.put('newmods'.encode(), str(newmods).encode())
                            cursor.put('elementinfo'.encode(), str(elementinfo).encode())
                            cursor.put('aminoacidcomposition'.encode(), str(self.aminoacidcomposition).encode())
                            
                            self.nonmonoisotopicelements = set()
                            self.isotopesbyelement = {} #element: [isotopes]
                            self.monoisotopickeys = {} #element: this is actually going to be the most abundant mass, but its monoisotopic for most of them
                            self.nonmonoisotopicgroups = {} #element: [non-most abundant masses which are usually nonmonoisotopic]
                            self.elementalmasses = {} #iso: mass
                            self.elementalprobabilities = {} #iso: abundance
                            for e, isos in elementinfo.items():
                                for iso, (mass, prob) in isos.items():
                                    self.elementalprobabilities[iso] = prob
                                    self.elementalmasses[iso] = mass
                                isolist, massandprobs = zip(*isos.items())
                                masses, probs = zip(*massandprobs)
                                maxind = np.argmax(probs)
                                monokey = isolist[maxind]
                                self.monoisotopickeys[e] = monokey
                                self.isotopesbyelement[e] = isolist
                                if len(isolist) > 1:
                                    self.nonmonoisotopicgroups[e] = list(isolist)
                                    self.nonmonoisotopicgroups[e].remove(monokey)
                                    self.nonmonoisotopicgroups[e] = tuple(self.nonmonoisotopicgroups[e])
                                    self.nonmonoisotopicelements.update(self.nonmonoisotopicgroups[e])
                            
                            self.elementvectors = {}
                            self.vectorpositions = {}
                            self.elementpositions = {}
                            for e, isos in self.isotopesbyelement.items():
                                self.elementvectors[e] = [0 for _ in range(len(isos))]
                                self.vectorpositions[e] = {k: n for n, k in enumerate(isos)}
                                self.elementpositions[e] = {n: k for n, k in enumerate(isos)}
                             
                            cursor.put('nonmonoisotopicelements'.encode(), str(self.nonmonoisotopicelements).encode())
                            cursor.put('isotopesbyelement'.encode(), str(self.isotopesbyelement).encode())
                            cursor.put('monoisotopickeys'.encode(), str(self.monoisotopickeys).encode())
                            cursor.put('nonmonoisotopicgroups'.encode(), str(self.nonmonoisotopicgroups).encode())
                            cursor.put('elementalmasses'.encode(), str(self.elementalmasses).encode())
                            cursor.put('elementalprobabilities'.encode(), str(self.elementalprobabilities).encode())
                            cursor.put('elementvectors'.encode(), str(self.elementvectors).encode())
                            cursor.put('vectorpositions'.encode(), str(self.vectorpositions).encode())
                            cursor.put('elementpositions'.encode(), str(self.elementpositions).encode())
                    break
                except lmdb.MapFullError:
                    defaults = env.open_db('defaults'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(defaults) as cursor:
                            mapaddition = int(cursor.get('mapsize'.encode()).decode())
                    newmapsize = env.info()['map_size'] + mapaddition
                    env.set_mapsize(newmapsize)
    
    def shuffle_string(self, seq):
        startchar = seq[0]
        endchar = seq[-1]
        charlist = list(seq[1:-1])
        random.shuffle(charlist)
        return ''.join((startchar, *charlist, endchar))
    
    def fasta_handling(self):
        t1 = time()
        
        fasta = SeqIO.parse(open(self.proteomefile), 'fasta')
        proteinseqs = {} #proteinid: seq
        for f in fasta:
            sequence, idn = str(f.seq), f.id
            proteinseqs[idn] = sequence
            #i'm not going to put too much thought into decoy shuffling, ie if it generates the same peptides as whats in the DB i'm just going to remove them as decoys and not care too much
            #^actually i won't remove them, they'll still count as target peptides but can be used for decoy proteins maybe?
            #i can loop this X amount of times if I want more decoys
            decoyidn = 'decoy_' + idn
            decoysequence = self.shuffle_string(sequence)
            proteinseqs[decoyidn] = decoysequence
        
        print(time() - t1, 'fasta processed')
        t2 = time()
        
        with environment(self.librarylocation) as env:
            while True:
                try:
                    proteomedb = env.open_db((self.proteome + '.info').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(proteomedb) as cursor:
                            cursor.put('nproteins'.encode(), str(len(proteinseqs)).encode())
                    break
                except lmdb.MapFullError:
                    defaults = env.open_db('defaults'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(defaults) as cursor:
                            mapaddition = int(cursor.get('mapsize'.encode()).decode())
                    newmapsize = env.info()['map_size'] + mapaddition
                    env.set_mapsize(newmapsize)
        
        print(time() - t2, 'database updated')
        return proteinseqs

    #its slower than the old function, 1.8ms vs 4.8s, but idc really, i forget how to do regex and its not like this is a bottleneck, this covers all the edgecases i can foresee
    def enzymatic_cleavage(self, proteinseqs):
        t1 = time()
        #appropriately complex example of the below:
        #cutters = {'P': 0, 'T': 1}
        #noncutters = {'P': ('EGP', 'PG'), 'T': ('TK', 'TTTA')}
        #noncutmaps would be: {'P': ((-2, 1), (0, 2)), 'T': ((0, 2), (0, 4))}
        #this approach would start to fail if there's overlap between what something that cuts and something isn't supposed to cut -> user error
        
        noncutmaps = defaultdict(dict) #AA: noncutAAs: (left coord, right coord) for what range to search for a non-cut
        for c in self.cutters:
            if c in self.noncutters:
                for section in self.noncutters[c]:
                    start = -section.index(c)
                    end = len(section) - abs(start)
                    noncutmaps[c][section] = (start, end)
        
        digestsofproteins = {} #proteinid: [peptides]
        for proteinid, seq in proteinseqs.items():
            lastcut = 0
            cutseq = []
            for n, aa in enumerate(seq):
                if aa in self.cutters:
                    if aa in self.noncutters:
                        cutpass = True
                        for section, (start, end) in noncutmaps[aa].items():
                            if seq[n+start:n+end] == section:
                                #print(section, seq[n+start:n+end]) #test passes
                                cutpass = False
                                break
                        if cutpass:
                            #sequence will cut
                            cutseq.append(seq[lastcut:n+1])
                            lastcut = n + 1
                    else:
                        #sequence will cut
                        cutseq.append(seq[lastcut:n+1])
                        lastcut = n + 1
            if lastcut < n:
                cutseq.append(seq[lastcut:n+1])
            extramissedcleavages = []
            if cutseq[0].startswith('M'):
                #adding in the cleaved version of the initial peptide, keeping both for good measure i suppose
                #this needs to be handled uniquely for missed cleavages
                premissedcleavages = cutseq[:self.missedcleavages+1]
                premissedcleavages[0] = premissedcleavages[0][1:]
                joinedseqs = ''
                for mseq in premissedcleavages:
                    joinedseqs += mseq
                    extramissedcleavages.append(joinedseqs)
            missedcleavagelist = []
            for length in range(1, self.missedcleavages+2):
                for position in range(len(cutseq)-length+1):
                    missedcleavagelist.append(''.join((cutseq[position:position+length+1])))
            cutseq.extend(extramissedcleavages)
            cutseq.extend(missedcleavagelist)
            digestsofproteins[proteinid] = list(filter(lambda x: (len(x) >= self.minlength and len(x) <= self.maxlength), cutseq))
        
        print(time() - t1, 'sequences digested')
        return digestsofproteins
    
    def sequence_modifications(self, digestsofproteins):
        t1 = time()

        switchmods = {'B': ('D', 'N'),
                      'J': ('L', 'I'),
                      'Z': ('E', 'Q')}
        
        #doing these before modifications for obvious reasons
        for proteinid, cseqs in digestsofproteins.items():
            removals = []
            for seq in cseqs:
                if 'X' in seq:
                    #not bothering with this one, although its possible
                    removals.append(seq)
                    continue
                switchers = set()
                if 'B' in seq: #D or N
                    switchers.add('B')
                if 'J' in seq: #L or I
                    switchers.add('J')
                if 'Z' in seq: #E or Q
                    switchers.add('Z')
                if switchers:
                    switchinds = {} #index: 
                    for n, aa in enumerate(seq):
                        if aa in switchers:
                            switchinds[n] = switchmods[aa]
                    seqindices, aminoswitches = zip(*switchinds.items())
                    ilen = len(seqindices)
                    seqslices = [slice(0, seqindices[0])]
                    for n, si in enumerate(seqindices):
                        if n < ilen-1:
                            seqslices.append(slice(si+1, seqindices[n+1]))
                    slen = len(seq)
                    modseqs = []
                    for seqmods in product(*aminoswitches):
                        modseq = ''
                        for sl, sm in zip(seqslices, seqmods):
                            modseq += seq[sl] + sm
                        if sl.stop < slen:
                            modseq += seq[sl.stop+1:]
                        modseqs.append(modseq)
                    cseqs.extend(modseqs)
                    removals.append(seq)
                #else:
                    #no problems
            for r in removals:
                while True:
                    try:
                        cseqs.remove(r)
                    except ValueError:
                        #all instances removed
                        break
            digestsofproteins[proteinid] = cseqs
        
        print(time() - t1, 'sequence amino acids adjusted')
        t2 = time()
        
        proteinsofseqs = defaultdict(list) #seq: [proteinids]
        for proteinid, cseqs in digestsofproteins.items():
            for seq in cseqs:
                proteinsofseqs[seq].append(proteinid)
        
        encodedproteinsofseqs = {} #encoded version
        for seq, proteinids in proteinsofseqs.items():
            proteinids = list(set(proteinids))
            encodedproteinsofseqs[seq.encode()] = str(proteinids).encode()
        
        with environment(self.librarylocation) as env:
            while True:
                try:
                    proteomedb = env.open_db((self.proteome + '.proteinsofseqs').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(proteomedb) as cursor:
                            cursor.putmulti(encodedproteinsofseqs.items(), dupdata=False, overwrite=True, append=False)
                    break
                except lmdb.MapFullError:
                    defaults = env.open_db('defaults'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(defaults) as cursor:
                            mapaddition = int(cursor.get('mapsize'.encode()).decode())
                    newmapsize = env.info()['map_size'] + mapaddition
                    env.set_mapsize(newmapsize)
        
        print(time() - t2, 'proteins indexed and saved')
        t3 = time()
        
        #i can just reverse-derive non-modified sequences later to save storage
        flatseqs = list(set(chain.from_iterable(digestsofproteins.values())))
        
        if self.modifiers:
            modseqs = set()
            for seq in flatseqs:
                if any(i in seq for i in self.modifiers):
                    #setup for allowing multiple variable mods of the same AA
                    modlocs = {}
                    modsubs = {}
                    idn = 0
                    for n, s in enumerate(seq):
                        if s in self.modifiers:
                            for mod in self.modifiers[s]:
                                modlocs[idn] = n
                                modsubs[idn] = mod
                                idn += 1
                    
                    for cn in range(1, self.maxvmods+1):
                        for modcombo in combinations(modlocs, cn):
                            if len(set(modlocs[i] for i in modcombo)) == cn:
                                newseq = seq
                                for subid in modcombo:
                                    newseq = ''.join((newseq[:modlocs[subid]], modsubs[subid], newseq[modlocs[subid]+1:]))
                                modseqs.add(newseq)
            flatseqs.extend(modseqs)
            print(time() - t3, 'modifications applied')
        else:
            print(time() - t3, 'sequences organized')
        return flatseqs
    
    def formula_counts(self, seq):
        atomiccomposition = Counter()
        for aa in seq:
            atomiccomposition += self.aminoacidcomposition[aa]
        #no OH loss on last residue, no H lost on first residue
        atomiccomposition['H'] += 2
        atomiccomposition['O'] += 1
        formulastring = ''.join((''.join((k, str(v))) for k, v in atomiccomposition.items()))
        return formulastring, atomiccomposition, seq
    
    def formula_consolidation(self, flatseqs):
        t1 = time()
        
        atomiccompositions = {} #formula string: atomic composition dict
        seqsbyformula = defaultdict(list) #formula string: [seqs]
        with mp.Pool(self.nprocs) as pool:
            for formulastring, atomiccomposition, seq in pool.map(self.formula_counts, flatseqs):
                seqsbyformula[formulastring].append(seq)
                atomiccompositions[formulastring] = atomiccomposition

        print(time() - t1, 'sequences aggregated by formula')
        t2 = time()
        
        encodedseqsbyformula = {}
        for formula, seqs in seqsbyformula.items():
            encodedseqsbyformula[formula.encode()] = str(seqs).encode()
        
        with environment(self.librarylocation) as env:
            while True:
                try:
                    proteomedb = env.open_db((self.proteome + '.seqsbyformula').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(proteomedb) as cursor:
                            cursor.putmulti(encodedseqsbyformula.items(), dupdata=False, overwrite=True, append=False)
                    break
                except lmdb.MapFullError:
                    defaults = env.open_db('defaults'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(defaults) as cursor:
                            mapaddition = int(cursor.get('mapsize'.encode()).decode())
                    newmapsize = env.info()['map_size'] + mapaddition
                    env.set_mapsize(newmapsize)
        
        print(time() - t2, 'aggregated formulas added to database')
        return atomiccompositions
    
    def distribution_handler(self, mainformula, atomiccomposition):
        elementalorganizer = {} #element: [[preheaps]]
        for e, acount in atomiccomposition.items():
            elementstring = e + str(acount)
            #try/except if faster than an if/else, so i might as well
            #try:
            #    #elementlist = fast_nested_copy(elementalcache[elementstring])
            #    elementlist = self.elementalcache[elementstring]
            #except KeyError: #not in cache
            elementlist = self.individual_element_binomial_walk(e, acount)
            #elementalcache[elementstring] = elementlist
            elementalorganizer[e] = elementlist #don't need to copy the insides
        mainformula, subformulas, massesandabundances = self.descending_partial_products(mainformula, elementalorganizer)
        return mainformula, subformulas, massesandabundances
    
    def individual_element_binomial_walk(self, e, acount):
        #elementalorganizer = defaultdict(list)
        elementlist = []
        mainheap = []
        vectorsets = defaultdict(set) #element: set of used vectors
        #for e, acount in atomiccomposition.items():
        mk = self.monoisotopickeys[e]
        nvector = self.elementvectors[e].copy()
        nvector[self.vectorpositions[e][mk]] += acount
        if len(self.isotopesbyelement[e]) > 2:
            baseprob = self.elementalprobabilities[mk] ** acount
            preheap = []
            preheap.append([baseprob, acount * self.elementalmasses[mk], e, nvector.copy()])
            greater = True
            lastprob = baseprob
            while greater:
                greater = False
                for iso in self.nonmonoisotopicgroups[e]:
                    newelementvector = nvector.copy()
                    newelementvector[self.vectorpositions[e][mk]] -= 1
                    if newelementvector[self.vectorpositions[e][mk]] > -1:
                        newelementvector[self.vectorpositions[e][iso]] += 1
                        vectorsets[e].add(tuple(newelementvector))
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(newelementvector):
                            loopiso = self.elementpositions[e][n]
                            newelementmass += self.elementalmasses[loopiso] * c
                            newelementprob *= self.elementalprobabilities[loopiso]**c
                            if loopiso in self.nonmonoisotopicelements:
                                newelementprob *= comb(acount-pn, c)
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
            for iso in self.nonmonoisotopicgroups[e]:
                v = nv.copy()
                v[self.vectorpositions[e][mk]] -= 1
                if v[self.vectorpositions[e][mk]] > -1:
                    v[self.vectorpositions[e][iso]] += 1
                    tuplevec = tuple(v)
                    if tuplevec not in vectorsets[e]:
                        vectorsets[e].add(tuplevec)
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(v):
                            loopiso = self.elementpositions[e][n]
                            newelementmass += self.elementalmasses[loopiso] * c
                            newelementprob *= self.elementalprobabilities[loopiso]**c
                            if loopiso in self.nonmonoisotopicelements:
                                newelementprob *= comb(acount-pn, c)
                                pn += c
                        heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
        else:
            preheap = []
            baseprob = self.elementalprobabilities[mk] ** acount
            preheap.append([baseprob, acount * self.elementalmasses[mk], e, nvector.copy()])
            greater = True
            lastprob = baseprob
            iso = self.nonmonoisotopicgroups[e][0]
            while greater:
                greater = False
                nvector[self.vectorpositions[e][mk]] -= 1
                if nvector[self.vectorpositions[e][mk]] > -1:
                    nvector[self.vectorpositions[e][iso]] += 1
                    vectorsets[e].add(tuple(nvector))
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(nvector):
                        loopiso = self.elementpositions[e][n]
                        newelementmass += self.elementalmasses[loopiso] * c
                        newelementprob *= self.elementalprobabilities[loopiso]**c
                        if loopiso in self.nonmonoisotopicelements:
                            newelementprob *= comb(acount-pn, c)
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
            v[self.vectorpositions[e][mk]] -= 1
            if v[self.vectorpositions[e][mk]] > -1:
                v[self.vectorpositions[e][iso]] += 1
                tuplevec = tuple(v)
                if tuplevec not in vectorsets[e]:
                    vectorsets[e].add(tuplevec)
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(v):
                        loopiso = self.elementpositions[e][n]
                        newelementmass += self.elementalmasses[loopiso] * c
                        newelementprob *= self.elementalprobabilities[loopiso]**c
                        if loopiso in self.nonmonoisotopicelements:
                            newelementprob *= comb(acount-pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
        
        cutoff = -maxprob * self.dividingthreshold
    
        r, p, m, e, v = heapq.heappop(mainheap)
        elementlist.append([r, p, m, e, v])
        if len(self.isotopesbyelement[e]) > 2:
            while p > cutoff:
                for iso in self.nonmonoisotopicgroups[e]:
                    newelementvector = v.copy()
                    newelementvector[self.vectorpositions[e][mk]] -= 1
                    if newelementvector[self.vectorpositions[e][mk]] > 0:
                        newelementvector[self.vectorpositions[e][iso]] += 1
                        tuplevec = tuple(newelementvector)
                        if tuplevec not in vectorsets[e]:
                            vectorsets[e].add(tuplevec)
                            pn = 0
                            newelementmass = 0
                            newelementprob = 1
                            for n, c in enumerate(newelementvector):
                                loopiso = self.elementpositions[e][n]
                                newelementmass += self.elementalmasses[loopiso] * c
                                newelementprob *= self.elementalprobabilities[loopiso]**c
                                if loopiso in self.nonmonoisotopicelements:
                                    newelementprob *= comb(acount-pn, c)
                                    pn += c
                            heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, newelementvector.copy()])
                r, p, m, e, v = heapq.heappop(mainheap)
                elementlist.append([r, p, m, e, v])
        else:
            iso = self.nonmonoisotopicgroups[e][0]
            while p > cutoff:
                nvector = v.copy()
                nvector[self.vectorpositions[e][mk]] -= 1
                if nvector[self.vectorpositions[e][mk]] > 0:
                    nvector[self.vectorpositions[e][iso]] += 1
                    tuplevec = tuple(nvector)
                    if tuplevec not in vectorsets[e]:
                        vectorsets[e].add(tuplevec)
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(nvector):
                            loopiso = self.elementpositions[e][n]
                            newelementmass += self.elementalmasses[loopiso] * c
                            newelementprob *= self.elementalprobabilities[loopiso]**c
                            if loopiso in self.nonmonoisotopicelements:
                                newelementprob *= comb(acount-pn, c)
                                pn += c
                        heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, nvector.copy()])
                r, p, m, e, v = heapq.heappop(mainheap)
                elementlist.append([r, p, m, e, v])
        heapq.heapify(elementlist)
        return elementlist
    
    def descending_partial_products(self, mainformula, elementalorganizer):
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
                        formula += f'{self.elementpositions[e][n]}({c})'
                maxprob *= p
                mainmass += m
    
        finalabundances[formula] = [mainmass, maxprob]
    
        cutoff = maxprob * self.dividingthreshold
        mainheap = list(chain(*elementalorganizer.values()))
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
                            formula += f'{self.elementpositions[se][n]}({c})'
                    prob *= sp
                    mass += sm
            
            finalabundances[formula] = [mass, prob]
            if prob < cutoff:
                break
            
            ind = bisect(probabilityranking, r)
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
                                        sef += f'{self.elementpositions[se][n]}({c})'
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
                if -newratio >= self.dividingthreshold:
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
                                        seformula += f'{self.elementpositions[ie][n]}({c})'
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
                                        seformula += f'{self.elementpositions[ie][n]}({c})'
                                newprob *= ip
                                newmass += im
                    if newprob >= cutoff:
                        finalabundances[seformula] = [newmass, newprob]
                        if multielement:
                            ind = bisect(probabilityranking, newratio)
                            probabilityranking.insert(ind, newratio)
                            multinomialpath.insert(ind, [newratio, *multipath])
                        else:
                            ind = bisect(probabilityranking, newratio)
                            probabilityranking.insert(ind, newratio)
                            multinomialpath.insert(ind, [newratio, [sr, sp, sm, se, sv], [r, p, m, e, v]])
                else:
                    break
    
        subformulas, massesandabundances = list(zip(*finalabundances.items()))
        subformulas = np.array(subformulas, dtype='S')
        massesandabundances = np.array(massesandabundances).transpose()
        subformulas = subformulas[massesandabundances[0].argsort()].tolist()
        massesandabundances = massesandabundances[:,massesandabundances[0].argsort()]
        return mainformula, subformulas, massesandabundances
    
    def isotope_library_generation(self, atomiccompositions):
        t1 = time()
        
        with environment(self.librarylocation) as env:
            while True:
                try:
                    pulledformulas = []
                    encodedformulas = tuple(i.encode() for i in atomiccompositions)
                    fulls = env.open_db('distributions.full'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(fulls) as cursor:
                            for k, v in cursor.getmulti(encodedformulas):
                                pulledformulas.append(k.decode())
                    break
                except lmdb.MapFullError:
                    defaults = env.open_db('defaults'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(defaults) as cursor:
                            mapaddition = int(cursor.get('mapsize'.encode()).decode())
                    newmapsize = env.info()['map_size'] + mapaddition
                    env.set_mapsize(newmapsize)
        
        for k in pulledformulas:
            del atomiccompositions[k]
        
        alen = len(pulledformulas)
        if alen > 0:
            print(time() - t1, '- pulled', alen, 'existing formulas from library')
        t2 = time()
        
        #elementalcache = {} #this worked as a cache with the old functional approach and i have no idea why
        with mp.Pool(self.nprocs) as pool:
            formulastrings, subformulas, massesandabundances = zip(*pool.starmap(self.distribution_handler, atomiccompositions.items()))
        
        #you set this above for db pulling and then delete it here, why?
        abundances = dict(zip(formulastrings, massesandabundances))
        abundanceformulas = dict(zip(formulastrings, subformulas))
        
        encodedabundanceformulas = {}
        for formula, subformulagroup in abundanceformulas.items():
            encodedabundanceformulas[formula.encode()] = str(subformulagroup).encode()

        encodedabundances = {}
        for formula, massesandabundances in abundances.items():
            encodedabundances[formula.encode()] = massesandabundances.tobytes()
        
        print(time() - t2, len(abundances), 'isotopic distributions total')
        t3 = time()
        
        with environment(self.librarylocation) as env:
            while True:
                try:
                    formuladb = env.open_db(('distributions.formulas').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(formuladb) as cursor:
                            cursor.putmulti(encodedabundanceformulas.items(), dupdata=False, overwrite=True, append=False)
                    fulldistributiondb = env.open_db(('distributions.full').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(fulldistributiondb) as cursor:
                            cursor.putmulti(encodedabundances.items(), dupdata=False, overwrite=True, append=False)
                    break
                except lmdb.MapFullError:
                    defaults = env.open_db('defaults'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(defaults) as cursor:
                            mapaddition = int(cursor.get('mapsize'.encode()).decode())
                    newmapsize = env.info()['map_size'] + mapaddition
                    env.set_mapsize(newmapsize)
        
        print(time() - t3, 'distributions and subformulas processed')
        return abundances, subformulas, pulledformulas
    
    def subisotopomer_handler(self, formulawithabundancedist, subformulas):
        formula, abundancedist = formulawithabundancedist
        massgroups = defaultdict(list) #massnumber: [masses]
        intensitygroups = defaultdict(list) #massnumber: [abundances]
        masses, abundances = abundancedist
        for n, s in enumerate(subformulas):
            s = s.decode()
            massnumber = 0
            for ss in s.split(')')[:-1]:
                splitval = 0
                #for handling elements with multiple letters
                while True:
                    if ss[splitval].isalpha():
                        splitval += 1
                    else:
                        break
                i1, i2 = list(map(int, ss[splitval:].split('(')))
                massnumber += i1 * i2
            massgroups[massnumber].append(masses[n])
            intensitygroups[massnumber].append(abundances[n])
        meansofmasses = []
        maxabundances = []
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
            maxab = 0
            for n, (sm, sa) in sorted(enumerate(zip(m, a)), key=lambda x: -x[1][1]): #sorting the enumeration is done purposefully here to preserve the original order i think
                if sa > maxab:
                    maxab = sa
                weightedmass += sm * sa
                cumulativeabundance += sa
                cumpercent = cumulativeabundance / totalabundance
                if cumpercent <= self.subisotopomericdepth:
                    subinds.append(n)
                elif subisocontinuation and cumpercent >= self.subisotopomericdepth:
                    #subisodepth threshold breached, end it with this one, but the rest continues for the weighted average to be correct
                    subisocontinuation = False
                    subinds.append(n)
            maxabundances.append(maxab)
            meansofmasses.append(weightedmass / totalabundance)
            sumsofabundances.append(totalabundance)
            subisodepthindices.append(tuple(subinds))
        sumabundancedist = np.array([meansofmasses, sumsofabundances])
        return formula, np.array(maxabundances), sumabundancedist, np.array(condensationindices), subisodepthindices
    
    def library_processing(self):
        abundances, subformulas, pulledformulas = self.isotope_library_generation(self.formula_consolidation(self.sequence_modifications(self.enzymatic_cleavage(self.fasta_handling()))))
        t1 = time()
        
        encodedmaxabundancedistributions = {}
        encodedsumabundancedistributions = {}
        encodedcondensationcoordinates = {} #formula: [# isotopomers per proton-step]
        encodedsubisodepthqualifiers = {} #formula: [[top n subisos of proton locations in descending order of abundance]]
        with mp.Pool(self.nprocs) as pool:
            for k, maxes, sums, coords, dsubisos in pool.starmap(self.subisotopomer_handler, zip(abundances.items(), subformulas)):
                ek = k.encode()
                encodedmaxabundancedistributions[ek] = maxes.tobytes()
                encodedsumabundancedistributions[ek] = sums.tobytes()
                encodedcondensationcoordinates[ek] = coords.tobytes()
                encodedsubisodepthqualifiers[ek] = str(dsubisos).encode()
        
        print(time() - t1, 'subisotopic differences processed')
        
        uppermass = max(v[0].max() for v in abundances.values())
        rum = round(uppermass)
        rlen = len(str(rum))
        #ceiling of the largest place
        uppermasslimit = (int(str(rum)[0]) + 1) * 10**(rlen-1)
        
        print('upper mass limit', uppermasslimit)
        t2 = time()

        pulledformulas.extend(abundances) #full formula list
    
        with environment(self.librarylocation) as env:
            while True:
                try:
                    formuladb = env.open_db(('proteomes.formulalist').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(formuladb) as cursor:
                            cursor.put(self.proteome.encode(), str(pulledformulas).encode())
                    proteomedb = env.open_db((self.proteome + '.info').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(proteomedb) as cursor:
                            cursor.put('uppermasslimit'.encode(), str(uppermasslimit).encode())
                            cursor.put('nformulas'.encode(), str(len(pulledformulas)).encode())
                    sumdistributiondb = env.open_db(('distributions.sum').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(sumdistributiondb) as cursor:
                            cursor.putmulti(encodedsumabundancedistributions.items(), dupdata=False, overwrite=True, append=False)
                    sumdistributiondb = env.open_db(('distributions.max').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(sumdistributiondb) as cursor:
                            cursor.putmulti(encodedmaxabundancedistributions.items(), dupdata=False, overwrite=True, append=False)
                    condensationdb = env.open_db(('distributions.condensationcoordinates').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(condensationdb) as cursor:
                            cursor.putmulti(encodedcondensationcoordinates.items(), dupdata=False, overwrite=True, append=False)
                    subisodb = env.open_db(('distributions.subisodepthqualifiers').encode())
                    with env.begin(write=True) as txn:
                        with txn.cursor(subisodb) as cursor:
                            cursor.putmulti(encodedsubisodepthqualifiers.items(), dupdata=False, overwrite=True, append=False)
                    break
                except lmdb.MapFullError:
                    defaults = env.open_db('defaults'.encode())
                    with env.begin(write=False) as txn:
                        with txn.cursor(defaults) as cursor:
                            mapaddition = int(cursor.get('mapsize'.encode()).decode())
                    newmapsize = env.info()['map_size'] + mapaddition
                    env.set_mapsize(newmapsize)
        print(time() - t2, 'proteome saved and processed')

def library_additions(librarylocation, proteomefile, minlength, maxlength, missedcleavages, nprocs, variablemodifications=[], staticmodifications=[], enzyme=-1, maxvmods=0):
    t1 = time()
    
    if nprocs == 0:
        nprocs = os.cpu_count()
    
    lib = LibraryOrganizer(librarylocation, proteomefile, minlength, maxlength, missedcleavages, nprocs, enzyme=enzyme, variablemodifications=variablemodifications, staticmodifications=staticmodifications, maxvmods=maxvmods)
    
    if lib.passing:
        lib.library_processing()
    
    print(time() - t1, 'total')

#librarylocation = '/home/sfo/data/proteomics/fastas/search-db'
#minlength = 6
#maxlength = 30
#maxvmods = 3
#staticmodifications = []
#variablemodifications = [0]
#enzyme = 0
#missedcleavages = 1
#proteomefile = '/home/sfo/data/proteomics/fastas/proteomes/Human_Homo_sapien-NoTremb.fasta'
