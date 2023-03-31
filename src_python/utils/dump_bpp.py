def output_bpp_to_file(model, file_name):
    seq_length = model.seq_length

    with open(file_name, 'w') as f:
        turn = 3 if model.no_sharp_turn else 0
        for i in range(1, seq_length+1):
            for j in range(i + turn + 1, seq_length+1):
                key = (i, j)
                if key in model.Pij:
                    f.write(f"{i} {j} {model.Pij[key]:.4e}\n")
        f.write("\n")
        print(f"Wrote bpp to {file_name}")