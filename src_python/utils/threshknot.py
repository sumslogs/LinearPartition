def dotbracket_to_pairs(structure):
    pairs = {}
    s = []
    index = 1
    for elem in structure:
        if elem == '(':
            s.append(index)
        elif elem == ')':
            pre_index = s.pop()
            pairs[pre_index] = index
            pairs[index] = pre_index
        index += 1
    return pairs

def output_bpseq(file_name, pairs, seq, mode='w'):
    if file_name:
        print(f"Outputting base pairs in bpseq format to {file_name}...")
        with open(file_name, mode) as fptr:
            for i in range(1, len(seq)+1):
                j = pairs.get(i, 0)
                nuc = seq[i-1]
                print(f"{i} {nuc} {j}", file=fptr)
            print("", file=fptr)
    else:
        for i in range(1, len(seq)+1):
            j = pairs.get(i, 0)
            nuc = seq[i-1]
            print(f"{i} {nuc} {j}")

def threshknot_pairs(model, threshknot_threshold):
    rowprob = {}
    prob_list = []

    pairs = {}
    visited = set()

    for (i,j), score in model.Pij.items():
        if score < threshknot_threshold:
            continue
        prob_list.append((i, j, score))
        rowprob[i] = max(rowprob.get(i, float('-inf')), score)
        rowprob[j] = max(rowprob.get(j, float('-inf')), score)

    for elem in prob_list:
        i, j, score = elem
        if score == rowprob[i] and score == rowprob[j]:
            if i in visited or j in visited:
                continue
            visited.add(i)
            visited.add(j)

            pairs[i] = j
            pairs[j] = i

    return pairs