def back_trace(i, j, back_pointer, seq_length):
    if i > j:
        return ""
    elif back_pointer[i][j] == -1:
        if i == j:
            return "."
        else:
            return "." + back_trace(i + 1, j, back_pointer, seq_length)
    elif back_pointer[i][j] != 0:
        k = back_pointer[i][j]
        assert 0 < k + 1 <= seq_length
        if k == j:
            temp = ""
        else:
            temp = back_trace(k + 1, j, back_pointer, seq_length)
        return "(" + back_trace(i + 1, k - 1, back_pointer, seq_length) + ")" + temp
    assert False

def PairProb_MEA(model, gamma):
    seq_length = model.seq_length
    OPT = [[0] * seq_length for _ in range(seq_length)]
    P = [[0] * seq_length for _ in range(seq_length)]
    back_pointer = [[0] * seq_length for _ in range(seq_length)]
    paired = [[] for _ in range(seq_length)]
    Q = [1.0] * seq_length

    for (i,j), score in model.Pij.items():
        i, j = i-1, j-1
        P[i][j] = score
        paired[i].append(j)
        Q[i] -= score
        Q[j] -= score

    for i in range(seq_length):
        paired[i].sort()
    for l in range(seq_length):
        for i in range(seq_length - l):
            j = i + l
            if i == j:
                OPT[i][j] = Q[i]
                back_pointer[i][j] = -1
                continue
            OPT[i][j] = OPT[i][i] + OPT[i + 1][j]
            back_pointer[i][j] = -1
            for k in paired[i]:
                if k > j:
                    break
                if k < j:
                    temp_OPT_k1_j = OPT[k + 1][j]
                else:
                    temp_OPT_k1_j = 0.0
                temp_score = 2 * gamma * P[i][k] + OPT[i + 1][k - 1] + temp_OPT_k1_j
                if OPT[i][j] < temp_score:
                    OPT[i][j] = temp_score
                    back_pointer[i][j] = k
    structure = back_trace(0, seq_length - 1, back_pointer, seq_length)
    return structure