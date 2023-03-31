def print_states(fptr, states, j, label, inside_only, threshold):
    for i, state in reversed(states.items()):
        if inside_only:
            print(f"{label} {i+1} {j+1} {state.alpha:.5f}", file=fptr)
        elif state.alpha + state.beta > threshold: # lhuang : alpha + beta - totalZ < ...
            print(f"{label} {i+1} {j+1} {state.alpha:.5f} {state.beta:.5f}", file=fptr)

def dump_forest(model, file_name, inside_only=False):
    print(f"Dumping ({'Inside-Only' if inside_only else 'Inside-Outside'}) Forest to {file_name}...")
    with open(file_name, "w") as fptr:
        seq_length = len(model.seq)

        print(model.seq, file=fptr)
        for j in range(seq_length):
            if inside_only:
                print(f"E {j+1} {model.bestC[j].alpha:.5f}", file=fptr)
            else:
                print(f"E {j+1} {model.bestC[j].alpha:.5f} {model.bestC[j].beta:.5f}", file=fptr)

        threshold = model.bestC[seq_length-1].alpha - 9.91152 # lhuang -9.xxx or ?
        for j in range(seq_length):
            print_states(fptr, model.bestP[j], j, "P", inside_only, threshold)
            print_states(fptr, model.bestM[j], j, "M", inside_only, threshold)
            print_states(fptr, model.bestM2[j], j, "M2", inside_only, threshold)
            print_states(fptr, model.bestMulti[j], j, "Multi", inside_only, threshold)