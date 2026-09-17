import math
import re
import random
import numpy as np
from itertools import permutations, product

a_lattice = 2.867  # Lattice constant for Fe

# define shells with bcc symmetry, up to a certain shell
def shells_from_model(model_path, nshells):
    """Read the first nshells neighbour shells from a Yell model.txt.

    In:  model_path, nshells
    Out: keys    - representative vector per shell, in doubled integer coords
                   (x2 so BCC half-integer sites stay integers), components
                   sorted large->small so equivalent vectors share a key
         shell_r - matching shell radii in Angstroms

    Use: pass keys to shell_vectors() to get each shell's neighbour offsets.
    """

    # shell order taken from model.txt block order
    text = open(model_path).read()
    blocks = re.findall(
        r'\[\s*\(([^)]+)\)(.*?)SubstitutionalCorrelation\(\s*FeAl\s*,\s*FeAl\s*,\s*([^)]+?)\)',
        text, re.DOTALL)

    # convert string to doubled integer (to avoid fractions)
    def dbl(s):
        s = s.strip()
        if '/' in s:
            n, d = s.split('/'); return int(round(2 * int(n) / int(d)))
        return 2 * int(s)

    # get tuple representative of each shell
    keys = []
    for vec_str, _body, arg in blocks:
        if not re.search(r'(?:sub_|dp)\d+', arg):
            continue # skip the (0,0,0) self-correlation block
        keys.append(tuple(sorted((abs(dbl(c)) for c in vec_str.split(',')), reverse=True)))

    # return the first nshells of model.txt
    keys = keys[:nshells]
    shell_r = [math.sqrt(sum(c*c for c in k) / 4) * a_lattice for k in keys] # /4 undoes the coordinate doubling
    return keys, shell_r

def shell_vectors(key):
    """All symmetry-equivalent vectors of one shell representative.

    In:  key - doubled-coordinate representative, e.g. (2,0,0)
    Out: sorted list of every distinct permutation/sign image; length is the
         shell multiplicity (6 for (2,0,0), 8 for (1,1,1))

    Both v and -v are included, so pair sums over these need a factor of 0.5.
    """
    
    # all sign/permutation images of the doubled representative
    return sorted({tuple(s * c for s, c in zip(signs, perm))
                   for perm in permutations(key)
                   for signs in product((1, -1), repeat=3)})