import math
from collections import defaultdict

from .helpers import quickselect, fast_log_plus_equals, fast_exp
from .constants import VALUE_MIN, NEG_INF
from .vienna_scoring import (
    _allowed_pairs,
    SINGLE_MAX_LEN, NOTON, GET_ACGU_NUM,
    init_stacking_and_pairs,
)

from .contrafold_scoring import (
    initialize_cachesingle,
    score_multi_unpaired, score_multi, score_helix, score_hairpin,
    score_M1, score_junction_B, score_single_without_junctionB,
    score_external_paired, score_external_unpaired,
)

class State:
    def __init__(self):
        self.alpha = VALUE_MIN
        self.beta = VALUE_MIN

class ContrafoldLinearPartition:
    def __init__(self, seq, beam_size, bpp_cutoff, no_sharp_turn):
        self.seq = seq
        self.seq_length = len(seq)
        self.nucs = [GET_ACGU_NUM(seq[i]) for i in range(self.seq_length)]
        seq_length = self.seq_length

        self.no_sharp_turn = no_sharp_turn

        # Pair probs
        self.Pij = {}
        self.bpp_cutoff = bpp_cutoff

        # Structures for Contrafold scoring model
        initialize_cachesingle()
        init_stacking_and_pairs()

        # Pair traversal structure (Uses _allowed_pairs from init_stacking_and_pairs)
        self.next_pair = [[] for _ in range(NOTON)]
        for nuci in range(NOTON):
            self.next_pair[nuci] = [-1] * seq_length
            next_ = -1
            for j in range(seq_length - 1, -1, -1):
                self.next_pair[nuci][j] = next_
                if _allowed_pairs[nuci][self.nucs[j]]:
                    next_ = j

        # Filled out in recursions
        self.beam = beam_size
        self.bestH = [defaultdict(State) for _ in range(seq_length + 1)]
        self.bestP = [defaultdict(State) for _ in range(seq_length + 1)]
        self.bestM2 = [defaultdict(State) for _ in range(seq_length + 1)]
        self.bestMulti = [defaultdict(State) for _ in range(seq_length + 1)]
        self.bestM = [defaultdict(State) for _ in range(seq_length + 1)]
        self.bestC = [State() for _ in range(seq_length + 1)]

    def beam_prune(self, beamstep):
        scores = []
        for i, cand in beamstep.items():
            k = i - 1
            newalpha = (self.bestC[k].alpha if k >= 0 else 0.0) + cand.alpha
            scores.append((newalpha, i))

        if len(scores) <= self.beam:
            return NEG_INF
        threshold = quickselect(scores, 0, len(scores) - 1, len(scores) - self.beam)
        for p in scores:
            if p[0] < threshold:
                beamstep.pop(p[1])

        return threshold

    def inside(self):
        seq_length = self.seq_length
        nucs = self.nucs
        next_pair = self.next_pair

        if seq_length > 0:
            self.bestC[0].alpha = fast_log_plus_equals(self.bestC[0].alpha, score_external_unpaired(0, 0))
        if seq_length > 1:
            self.bestC[1].alpha = fast_log_plus_equals(self.bestC[1].alpha, score_external_unpaired(0, 1))

        for j in range(seq_length):
            nucj = nucs[j]
            nucj1 = nucs[j + 1] if j + 1 < seq_length else -1

            beamstepH = self.bestH[j]
            beamstepMulti = self.bestMulti[j]
            beamstepP = self.bestP[j]
            beamstepM2 = self.bestM2[j]
            beamstepM = self.bestM[j]
            beamstepC = self.bestC[j]

            # beam of H
            if self.beam > 0 and len(beamstepH) > self.beam:
                self.beam_prune(beamstepH)
            # for nucj put H(j, j_next) into H[j_next]
            jnext = next_pair[nucj][j]
            if self.no_sharp_turn:
                while jnext - j < 4 and jnext != -1:
                    jnext = next_pair[nucj][jnext]
            if jnext != -1:
                nucjnext = nucs[jnext]
                nucjnext_1 = nucs[jnext - 1] if jnext - 1 > -1 else -1

                newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext)
                self.bestH[jnext][j].alpha = fast_log_plus_equals(self.bestH[jnext][j].alpha, newscore)
                #print(f'1. Setting: bestH[{jnext}][{j}] Alpha after: {self.bestH[jnext][j].alpha}')

            # for every state h in H[j]
            #   1. extend h(i, j) to h(i, jnext)
            #   2. generate p(i, j)
            for i, state in reversed(beamstepH.items()): # note: reversed to match C++ unordered_map iteration pattern
                nuci = nucs[i]
                jnext = next_pair[nuci][j]

                if jnext != -1:
                    nuci1 = nucs[i + 1] if i + 1 < seq_length else -1
                    nucjnext = nucs[jnext]
                    nucjnext_1 = nucs[jnext - 1] if (jnext - 1) > -1 else -1

                    # 1. extend h(i, j) to h(i, jnext)=
                    newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext)
                    self.bestH[jnext][i].alpha = fast_log_plus_equals(self.bestH[jnext][i].alpha, newscore)
                    #print(f'Alpha after: {self.bestH[jnext][i].alpha}')
                    #print(f'2. Setting: bestH[{jnext}][{i}] Alpha after: {self.bestH[jnext][i].alpha}')

                # 2. generate p(i, j)
                # print(f'beamstepP Alpha before: {beamstepP[i].alpha}')
                beamstepP[i].alpha = fast_log_plus_equals(beamstepP[i].alpha, state.alpha)
                # print(f'beamstepP Alpha after: {beamstepP[i].alpha}')
            if j == 0:
                continue

            # beam of Multi
            if self.beam > 0 and len(beamstepMulti) > self.beam:
                self.beam_prune(beamstepMulti)
            for i, state in reversed(beamstepMulti.items()): # note: reversed to match C++ unordered_map iteration pattern
                nuci = nucs[i]
                nuci1 = nucs[i+1]
                jnext = next_pair[nuci][j]

                # 1. extend (i, j) to (i, jnext)
                if jnext != -1:
                    newscore = score_multi_unpaired(j, jnext - 1)
                    self.bestMulti[jnext][i].alpha = fast_log_plus_equals(self.bestMulti[jnext][i].alpha, state.alpha + newscore)
                    #print(f'1. Setting: bestMulti[{jnext}][{i}] Alpha after: {self.bestMulti[jnext][i].alpha}')

                # 2. generate P (i, j)
                newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length)
                beamstepP[i].alpha = fast_log_plus_equals(beamstepP[i].alpha, state.alpha + newscore)

            # beam of P
            if self.beam > 0 and len(beamstepP) > self.beam:
                self.beam_prune(beamstepP)
            # for every state in P[j]
            # 1. generate new helix/bulge
            # 2. M = P
            # 3. M2 = M + P
            # 4. C = C + P
            for i, state in reversed(beamstepP.items()): # note: reversed to match C++ unordered_map iteration pattern
                nuci = nucs[i]
                nuci_1 = nucs[i-1] if i > 0 else -1

                # 1. generate new helix / single_branch
                # new state is of shape p..i..j..q
                if i > 0 and j < seq_length-1:
                    precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci)
                    for p in range(i-1, max(i-SINGLE_MAX_LEN-1, -1), -1):
                        nucp = nucs[p]
                        nucp1 = nucs[p+1]
                        q = next_pair[nucp][j]
                        while q != -1 and (i-p)+(q-j)-2 <= SINGLE_MAX_LEN:
                            nucq = nucs[q]
                            nucq_1 = nucs[q-1]

                            if p == i-1 and q == j+1:
                                # helix
                                newscore = score_helix(nucp, nucp1, nucq_1, nucq)
                                self.bestP[q][p].alpha = fast_log_plus_equals(self.bestP[q][p].alpha, state.alpha + newscore)
                            else:
                                # single branch
                                newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) + \
                                           precomputed + \
                                           score_single_without_junctionB(p, q, i, j,
                                                               nuci_1, nuci, nucj, nucj1)
                                self.bestP[q][p].alpha = fast_log_plus_equals(self.bestP[q][p].alpha, state.alpha + newscore)
                            q = next_pair[nucp][q]

                # 2. M = P
                if i > 0 and j < seq_length - 1:
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length)
                    beamstepM[i].alpha = fast_log_plus_equals(beamstepM[i].alpha, state.alpha + newscore)

                # 3. M2 = M + P
                k = i - 1
                if k > 0 and len(self.bestM[k]) > 0:
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length)
                    m1_alpha = state.alpha + newscore
                    for newi, m_state in reversed(self.bestM[k].items()):
                        beamstepM2[newi].alpha = fast_log_plus_equals(beamstepM2[newi].alpha, m_state.alpha + m1_alpha)

                # 4. C = C + P
                k = i - 1
                if k >= 0:
                    prefix_C = self.bestC[k]
                    nuck = nuci_1
                    nuck1 = nuci
                    newscore = score_external_paired(k+1, j, nuck, nuck1,
                                                         nucj, nucj1, seq_length)
                    beamstepC.alpha = fast_log_plus_equals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore)
                else:
                    newscore =  score_external_paired(0, j, -1, nucs[0],
                                                         nucj, nucj1, seq_length)
                    beamstepC.alpha = fast_log_plus_equals(beamstepC.alpha, state.alpha + newscore)

            # beam of M2
            if self.beam > 0 and len(beamstepM2) > self.beam:
                self.beam_prune(beamstepM2)
            for i, state in reversed(beamstepM2.items()): # note: reversed to match C++ unordered_map iteration pattern
                # 1. multi-loop
                for p in range(i-1, max(i-SINGLE_MAX_LEN, 0)-1, -1):
                    nucp = nucs[p]
                    q = next_pair[nucp][j]
                    if q != -1 and (i - p - 1) <= SINGLE_MAX_LEN:
                        newscore = score_multi_unpaired(p+1, i-1) + \
                                   score_multi_unpaired(j+1, q-1)
                        self.bestMulti[q][p].alpha = fast_log_plus_equals(self.bestMulti[q][p].alpha, state.alpha + newscore)

                # 2. M = M2
                beamstepM[i].alpha = fast_log_plus_equals(beamstepM[i].alpha, state.alpha)

            # beam of M
            if self.beam > 0 and len(beamstepM) > self.beam:
                self.beam_prune(beamstepM)
            for i, state in reversed(beamstepM.items()): # note: reversed to match C++ unordered_map iteration pattern
                if j < seq_length-1:
                    newscore = score_multi_unpaired(j + 1, j + 1)
                    self.bestM[j+1][i].alpha = fast_log_plus_equals(self.bestM[j+1][i].alpha, state.alpha + newscore)

            # beam of C
            # C = C + U
            if j < seq_length-1:
                newscore = score_external_unpaired(j+1, j+1)
                self.bestC[j+1].alpha = fast_log_plus_equals(self.bestC[j+1].alpha, beamstepC.alpha + newscore)

    def outside(self):
        next_pair = self.next_pair

        seq_length = self.seq_length
        nucs = self.nucs

        self.bestC[seq_length-1].beta = 0.0

        # from right to left
        for j in range(seq_length-1, 0, -1):
            nucj = nucs[j]
            nucj1 = nucs[j + 1] if j + 1 < seq_length else -1

            # beamstepH = self.bestH[j]
            beamstepMulti = self.bestMulti[j]
            beamstepP = self.bestP[j]
            beamstepM2 = self.bestM2[j]
            beamstepM = self.bestM[j]
            beamstepC = self.bestC[j]

            # beam of C
            if j < seq_length-1:
                newscore = score_external_unpaired(j+1, j+1)
                beamstepC.beta = fast_log_plus_equals(beamstepC.beta, self.bestC[j+1].beta + newscore)

            # beam of M
            for i, state in reversed(beamstepM.items()): # note: reversed to match C++ unordered_map iteration pattern
                if j < seq_length-1:
                    newscore = score_multi_unpaired(j + 1, j + 1)
                    state.beta = fast_log_plus_equals(state.beta, self.bestM[j+1][i].beta + newscore)

            # beam of M2
            for i, state in reversed(beamstepM2.items()): # note: reversed to match C++ unordered_map iteration pattern
                # 1. multi-loop
                for p in range(i-1, max(i-SINGLE_MAX_LEN, 0)-1, -1):
                    nucp = nucs[p]
                    q = next_pair[nucp][j]
                    if q != -1 and (i - p - 1) <= SINGLE_MAX_LEN:
                        newscore = score_multi_unpaired(p+1, i-1) + \
                                   score_multi_unpaired(j+1, q-1)
                        state.beta = fast_log_plus_equals(state.beta, self.bestMulti[q][p].beta + newscore)

                # 2. M = M2
                state.beta = fast_log_plus_equals(state.beta, beamstepM[i].beta)

            # beam of P
            for i, state in reversed(beamstepP.items()): # note: reversed to match C++ unordered_map iteration pattern
                nuci = nucs[i]
                nuci_1 = nucs[i-1] if i > 0 else -1

                # 1. generate new helix / single_branch
                # new state is of shape p..i..j..q
                if i > 0 and j < seq_length-1:
                    precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci)
                    for p in range(i-1, max(i-SINGLE_MAX_LEN-1, -1), -1):
                        nucp = nucs[p]
                        nucp1 = nucs[p+1]
                        q = next_pair[nucp][j]
                        while q != -1 and (i-p)+(q-j)-2 <= SINGLE_MAX_LEN:
                            nucq = nucs[q]
                            nucq_1 = nucs[q-1]

                            if p == i-1 and q == j+1:
                                # helix
                                newscore = score_helix(nucp, nucp1, nucq_1, nucq)
                                state.beta = fast_log_plus_equals(state.beta, self.bestP[q][p].beta + newscore)
                            else:
                                # single branch
                                newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) + \
                                           precomputed + \
                                           score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1)
                                state.beta = fast_log_plus_equals(state.beta, self.bestP[q][p].beta + newscore)
                            q = next_pair[nucp][q]

                # 2. M = P
                if i > 0 and j < seq_length - 1:
                    newscore =  score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length)
                    state.beta = fast_log_plus_equals(state.beta, beamstepM[i].beta + newscore)

                # 3. M2 = M + P
                k = i - 1
                if k > 0 and len(self.bestM[k]) > 0:
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length)
                    m1_alpha = newscore
                    m1_plus_P_alpha = state.alpha + m1_alpha
                    for newi, m_state in reversed(self.bestM[k].items()):
                        state.beta = fast_log_plus_equals(state.beta, beamstepM2[newi].beta + m_state.alpha + m1_alpha)
                        m_state.beta = fast_log_plus_equals(m_state.beta, beamstepM2[newi].beta + m1_plus_P_alpha)

                # 4. C = C + P
                k = i - 1
                if k >= 0:
                    nuck = nuci_1
                    nuck1 = nuci
                    newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length)
                    external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore
                    self.bestC[k].beta = fast_log_plus_equals(self.bestC[k].beta, state.alpha + external_paired_alpha_plus_beamstepC_beta)
                    state.beta = fast_log_plus_equals(state.beta, self.bestC[k].alpha + external_paired_alpha_plus_beamstepC_beta)
                else:
                    newscore = score_external_paired(0, j, -1, nucs[0],
                                                         nucj, nucj1, seq_length)
                    state.beta = fast_log_plus_equals(state.beta, beamstepC.beta + newscore)

            # beam of Multi
            for i, state in reversed(beamstepMulti.items()): # note: reversed to match C++ unordered_map iteration pattern
                nuci = nucs[i]
                nuci1 = nucs[i+1]
                jnext = next_pair[nuci][j]

                # 1. extend (i, j) to (i, jnext)
                if jnext != -1:
                    newscore = score_multi_unpaired(j, jnext - 1)
                    state.beta = fast_log_plus_equals(state.beta, self.bestMulti[jnext][i].beta + newscore)
                    #print(f'1. Setting: bestMulti[{jnext}][{i}] Alpha after: {self.bestMulti[jnext][i].alpha}')

                # 2. generate P (i, j)
                newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                state.beta = fast_log_plus_equals(state.beta, beamstepP[i].beta + newscore)

    def cal_PairProb(self):
        viterbi = self.bestC[self.seq_length-1]
        for j in range(self.seq_length):
            for i, state in reversed(self.bestP[j].items()): # note: reversed to match C++ unordered_map iteration pattern
                temp_prob_inside = state.alpha + state.beta - viterbi.alpha
                if temp_prob_inside > -9.91152:
                    prob = fast_exp(temp_prob_inside)
                    if prob > 1.0:
                        prob = 1.0
                    if prob < self.bpp_cutoff:
                        continue
                    self.Pij[(i+1, j+1)] = prob

    def log_partition_coeff(self):
        viterbi = self.bestC[self.seq_length-1]
        return viterbi.alpha