import math
from typing import List
from .nn_params import (
    stack37, bulge37,
    dangle5_37, dangle3_37,
    lxc37, TerminalAU37, mismatchH37, mismatchM37,
    Triloops, Tetraloops, Hexaloops,
    hairpin37, Tetraloop37, Hexaloop37, Triloop37,
    ML_intern37, ML_closing37,
    int11_37, int21_37, int22_37, internal_loop37,
    mismatchI37, mismatch1nI37, mismatch23I37, mismatchExt37,
    ninio37, MAX_NINIO,
)

NOTON = 5
NOTOND = 25
NOTONT = 125
EXPLICIT_MAX_LEN = 4
SINGLE_MIN_LEN = 0
SINGLE_MAX_LEN = 30
MULTI_MAX_LEN = 30
HAIRPIN_MAX_LEN = 30
BULGE_MAX_LEN = SINGLE_MAX_LEN
INTERNAL_MAX_LEN = SINGLE_MAX_LEN
SYMMETRIC_MAX_LEN = 15
ASYMMETRY_MAX_LEN = 28

def GET_ACGU_NUM(x):
    return 0 if x=='A' else (1 if x=='C' else (2 if x=='G' else (3 if x=='U' else 4)))

def num_to_nuc(x):
    return -1 if x == -1 else 0 if x == 4 else x + 1

def num_to_pair(x, y):
    return {
        0: 5 if y == 3 else 0,
        1: 1 if y == 2 else 0,
        2: 2 if y == 1 else 3 if y == 3 else 0,
        3: 4 if y == 2 else 6 if y == 0 else 0,
    }.get(x, 0)

# def nuc_to_pair(x, y):
#     return {
#         1: 5 if y == 4 else 0,
#         2: 1 if y == 3 else 0,
#         3: 2 if y == 2 else 3 if y == 4 else 0,
#         4: 4 if y == 3 else 6 if y == 1 else 0,
#     }.get(x, 0)

_allowed_pairs = [[False]*NOTON for _ in range(NOTON)]
_helix_stacking = [[[[False]*NOTON for _ in range(NOTON)] for _ in range(NOTON)] for _ in range(NOTON)]
cache_single = [[0.0]*(SINGLE_MAX_LEN+1) for _ in range(SINGLE_MAX_LEN+1)]

def HELIX_STACKING_OLD(x, y, z, w, val):
    global _helix_stacking
    _helix_stacking[GET_ACGU_NUM(x)][GET_ACGU_NUM(y)][GET_ACGU_NUM(z)][GET_ACGU_NUM(w)] = val

def init_stacking_and_pairs():
    global _helix_stacking, _allowed_pairs
    # for i, a in enumerate('ACGU'):
    #     for j, b in enumerate('ACGU'):
    #         _allowed_pairs[i][j] = allowed_pairs.get(a) == b
    _allowed_pairs[GET_ACGU_NUM('A')][GET_ACGU_NUM('U')] = True
    _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('A')] = True
    _allowed_pairs[GET_ACGU_NUM('C')][GET_ACGU_NUM('G')] = True
    _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('C')] = True
    _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('U')] = True
    _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('G')] = True

    # helix_stacking = {'AU': ['AU', 'CG', 'GC', 'GU', 'UA', 'UG'],
    #                   'CG': ['AU', 'CG', 'GC', 'GU', 'UA', 'UG'],
    #                   'GC': ['AU', 'CG', 'GU', 'UA'],
    #                   'GU': ['AU', 'UG', 'GC']}
    # for s1 in helix_stacking:
    #     for c1 in s1:
    #         for s2 in helix_stacking[s1]:
    #             for c2 in s2:
    #                 _helix_stacking[GET_ACGU_NUM(c1)][GET_ACGU_NUM(c2)][GET_ACGU_NUM(s1[0])][GET_ACGU_NUM(s2[1])] = True

    HELIX_STACKING_OLD('A', 'U', 'A', 'U', True)
    HELIX_STACKING_OLD('A', 'U', 'C', 'G', True)
    HELIX_STACKING_OLD('A', 'U', 'G', 'C', True)
    HELIX_STACKING_OLD('A', 'U', 'G', 'U', True)
    HELIX_STACKING_OLD('A', 'U', 'U', 'A', True)
    HELIX_STACKING_OLD('A', 'U', 'U', 'G', True)
    HELIX_STACKING_OLD('C', 'G', 'A', 'U', True)
    HELIX_STACKING_OLD('C', 'G', 'C', 'G', True)
    HELIX_STACKING_OLD('C', 'G', 'G', 'C', True)
    HELIX_STACKING_OLD('C', 'G', 'G', 'U', True)
    HELIX_STACKING_OLD('C', 'G', 'U', 'G', True)
    HELIX_STACKING_OLD('G', 'C', 'A', 'U', True)
    HELIX_STACKING_OLD('G', 'C', 'C', 'G', True)
    HELIX_STACKING_OLD('G', 'C', 'G', 'U', True)
    HELIX_STACKING_OLD('G', 'C', 'U', 'G', True)
    HELIX_STACKING_OLD('G', 'U', 'A', 'U', True)
    HELIX_STACKING_OLD('G', 'U', 'G', 'U', True)
    HELIX_STACKING_OLD('G', 'U', 'U', 'G', True)
    HELIX_STACKING_OLD('U', 'A', 'A', 'U', True)
    HELIX_STACKING_OLD('U', 'A', 'G', 'U', True)
    HELIX_STACKING_OLD('U', 'G', 'G', 'U', True)

def v_init_tetra_hex_tri(seq: str, seq_length: int, if_tetraloops: List[int], if_hexaloops: List[int], if_triloops: List[int]) -> None:
    # Triloops
    if seq_length-4 > 0:
        del if_triloops[seq_length-4:]
    for i in range(seq_length-4):
        if not ((seq[i] == 'C' and seq[i+4] == 'G') or (seq[i] == 'G' and seq[i+4] == 'C')):
            continue
        ts = Triloops.find(seq[i:i+5])
        if ts != -1:
            if_triloops[i] = ts // 6

    # TetraLoops
    if seq_length-5 > 0:
        del if_tetraloops[seq_length-5:]
    for i in range(seq_length-5):
        if seq[i] != 'C' or seq[i+5] != 'G':
            continue
        ts = Tetraloops.find(seq[i:i+6])
        if ts != -1:
            if_tetraloops[i] = ts // 7

    # Hexaloops
    if seq_length-7 > 0:
        del if_hexaloops[seq_length-7:]
    for i in range(seq_length-7):
        if seq[i] != 'A' or seq[i+7] != 'U':
            continue
        ts = Hexaloops.find(seq[i:i+8])
        if ts != -1:
            if_hexaloops[i] = ts // 9

def v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri_index=-1):
    size = j - i - 1
    type_ = num_to_pair(nuci, nucj)
    si1 = num_to_nuc(nuci1)
    sj1 = num_to_nuc(nucj_1)

    energy = 0

    if size <= 30:
        energy = hairpin37[size]
    else:
        energy = hairpin37[30] + int(lxc37 * math.log(size / 30))

    if size < 3:
        return energy

    if size == 4 and tetra_hex_tri_index > -1:
        return Tetraloop37[tetra_hex_tri_index]
    elif size == 6 and tetra_hex_tri_index > -1:
        return Hexaloop37[tetra_hex_tri_index]
    elif size == 3:
        if tetra_hex_tri_index > -1:
            return Triloop37[tetra_hex_tri_index]
        return energy + (TerminalAU37 if type_ > 2 else 0)

    energy += mismatchH37[type_][si1][sj1]
    return energy

# multi_loop
def E_MLstem(type_, si1, sj1):
    energy = 0

    if si1 >= 0 and sj1 >= 0:
        energy += mismatchM37[type_][si1][sj1]
    elif si1 >= 0:
        energy += dangle5_37[type_][si1]
    elif sj1 >= 0:
        energy += dangle3_37[type_][sj1]

    if type_ > 2:
        energy += TerminalAU37

    energy += ML_intern37

    return energy

def v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, len_):
    # Note: i, j, len_ unused
    tt = num_to_pair(nucj, nuci)
    si1 = num_to_nuc(nuci1)
    sj1 = num_to_nuc(nucj_1)

    return E_MLstem(tt, sj1, si1) + ML_closing37

def MIN2(a, b):
    return a if a <= b else b

MAXLOOP = 30
def v_score_single(i, j, p, q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1):
    si1 = num_to_nuc(nuci1)
    sj1 = num_to_nuc(nucj_1)
    sp1 = num_to_nuc(nucp_1)
    sq1 = num_to_nuc(nucq1)
    type_ = num_to_pair(nuci, nucj)
    type_2 = num_to_pair(nucq, nucp)
    n1,n2 = p-i-1, j-q-1
    if n1 > n2:
        nl,ns = n1,n2
    else:
        nl,ns = n2,n1
    u = 0
    energy = 0

    # stack
    if nl == 0:
        return stack37[type_][type_2]

    # bulge
    if ns == 0:
        energy = bulge37[nl] if nl <= MAXLOOP else (bulge37[30] + int(lxc37*math.log(nl/30.)))
        if nl == 1:
            energy += stack37[type_][type_2]
        else:
            if type_ > 2:
                energy += TerminalAU37
            if type_2 > 2:
                energy += TerminalAU37
        return energy

    # interior loops
    else:
        if ns == 1:
            # 1x1 loop
            if nl == 1:
                return int11_37[type_][type_2][si1][sj1]
            # 2x1 loop
            if nl == 2:
                if n1 == 1:
                    return int21_37[type_][type_2][si1][sq1][sj1]
                else:
                    return int21_37[type_2][type_][sq1][si1][sp1]
            # 1xn loop
            else:
                energy = internal_loop37[nl+1] if nl + 1 <= MAXLOOP else internal_loop37[30] + int(lxc37*math.log((nl+1)/30.))
                energy += MIN2(MAX_NINIO, (nl-ns)*ninio37)
                energy += mismatch1nI37[type_][si1][sj1] + mismatch1nI37[type_2][sq1][sp1];
                return energy

        if ns == 2:
            # 2x2 loop
            if nl == 2:
                return int22_37[type_][type_2][si1][sp1][sq1][sj1]
            # 2x3 loop
            if nl == 3:
                energy = internal_loop37[5] + ninio37
                energy += mismatch23I37[type_][si1][sj1] + mismatch23I37[type_2][sq1][sp1]
                return energy

        # generic interior loop (no else here!)
        u = nl + ns
        energy = internal_loop37[u] if u <= MAXLOOP else internal_loop37[30]+int(lxc37*math.log(u/30.))
        energy += MIN2(MAX_NINIO, (nl-ns)*ninio37)
        energy += mismatchI37[type_][si1][sj1] + mismatchI37[type_2][sq1][sp1]

    return energy

def v_score_M1(i, j, k, nuci_1, nuci, nuck, nuck1, len_):
    p, q = i, k
    tt = num_to_pair(nuci, nuck)
    sp1 = num_to_nuc(nuci_1)
    sq1 = num_to_nuc(nuck1)

    return E_MLstem(tt, sp1, sq1)

def v_score_external_paired(i, j, nuci_1, nuci, nucj, nucj1, len_):
    type_ = num_to_pair(nuci, nucj)
    si1 = num_to_nuc(nuci_1)
    sj1 = num_to_nuc(nucj1)
    energy = 0

    if si1 >= 0 and sj1 >= 0:
        energy += mismatchExt37[type_][si1][sj1]
    elif si1 >= 0:
        energy += dangle5_37[type_][si1]
    elif sj1 >= 0:
        energy += dangle3_37[type_][sj1]

    if type_ > 2:
        energy += TerminalAU37

    return energy