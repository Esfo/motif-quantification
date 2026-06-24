from collections import defaultdict, Counter

#source:
#https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
#https://physics.nist.gov/cgi-bin/cuu/Value?mp

proton = 1.007276554940804

elementinfo = {'H': {'H1': (1.00782503223, 0.999885),
                     'H2': (2.01410177812, 0.000115)},
               'C': {'C12': (12.0000000, 0.9893),
                     'C13': (13.00335483507, 0.0107)},
               'N': {'N14': (14.00307400443, 0.99636),
                     'N15': (15.00010889888, 0.00364)},
               'O': {'O16': (15.99491461957, 0.99757),
                     'O17': (16.99913175650, 0.00038),
                     'O18': (17.99915961286, 0.00205)},
               'S': {'S32': (31.9720711744, 0.9499),
                     'S33': (32.9714589098, 0.0075),
                     'S34': (33.967867004, 0.0425),
                     'S36': (35.96708071, 0.0001)},
               'P': {'P31': (30.97376199842, 1)},
               'Se': {'Se74': (73.922475934, 0.0089),
                      'Se76': (75.919213704, 0.0937),
                      'Se77': (76.919914154, 0.0763),
                      'Se78': (77.91730928, 0.2377),
                      'Se80': (79.9165218, 0.4961),
                      'Se82': (81.9166995, 0.0873)}}

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
        'V': {'C': 5, 'H': 9, 'N': 1, 'O': 1},
        'U': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'Se': 1}, #selenocysteine
        'O': {'C': 12, 'H': 19, 'N': 3, 'O': 2} #pyrrolysine
        }

#to selectively pick ions below, you can iterate these dicts in order and cumulatively combine until you hit an ion you want to generate, then use those cumulative +/-s as just a single dict entry each in fragmentation_compositions. so you would generate the dicts you plan on using in this file
nfragmentcompositions = {}
nfragmentcompositions['a'] = Counter({'C': -1, 'O': -1})
nfragmentcompositions['b'] = Counter()
nfragmentcompositions['c'] = Counter({'N': 1, 'H': 3})

cfragmentcompositions = {}
cfragmentcompositions['x'] = Counter({'C': 1, 'O': 2})
cfragmentcompositions['y'] = Counter({'H': 2, 'O': 1})
cfragmentcompositions['z'] = Counter({'N': -1, 'H': -1, 'O': 1})
#cfragmentcompositions['z•'] = Counter({'N': -1, 'O': 1})
