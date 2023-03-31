def load_shape_data(shape_file_path):
    # Load shape data if provided
    with open(shape_file_path, 'r') as f:
        SHAPE_data = [l.split("\t") for l in f.read().split("\n")]
        SHAPE_data = [float(l[1]) if l[1] != 'NA' else -1.0 for l in SHAPE_data if len(l) == 2]
    return SHAPE_data