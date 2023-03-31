import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
from collections import defaultdict

def draw_heatmap(matrix, figsize=(10,10)):
    seq_length = len(matrix)
    matrix_ticks = pd.DataFrame(data=matrix,
                                columns=range(1,seq_length+1),
                                index=range(1,seq_length+1))

    # Generate a mask for the upper triangle
    mask = np.tril(np.ones_like(matrix, dtype=np.bool))

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=figsize)

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    sns.heatmap(matrix_ticks, mask=mask, cmap=cmap, vmax=1., center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5},
                xticklabels="auto", yticklabels="auto")

    #plt.savefig(bpp_file+"_heatmap", format="pdf")
    #plt.savefig("heatmap", format="pdf")

    return f, ax

def drawarc_clockwise(a, b, deg, style, length, lengthfix, rotate): ## counterclock wise
    angle_a = 360./(length+lengthfix)*a
    angle_b = 360./(length+lengthfix)*b
    angle_a = angle_a-70 + rotate
    angle_b = angle_b-70 + rotate

    alpha = (angle_b-angle_a)
    pa = angle_a * math.pi / 180.
    pb = angle_b * math.pi / 180.
    palpha = alpha * math.pi / 180.
    arc_txt = ''
    if alpha < 170:
        arc_txt += f"\\draw[line width = 0.0mm] {style} ([shift=({angle_b}:10cm)]0,0) arc ({angle_b + 90}:{angle_a + 270}:{10 * math.tan(palpha/2.)}cm); %% {a} {b} {length}\n"
    elif alpha > 190:
        beta = 360. - alpha
        pbeta = beta * math.pi / 180.
        arc_txt += f"\\draw[line width = 0.0mm] {style} ([shift=({angle_a}:10cm)]0,0) arc ({angle_a + 90}:{angle_b - 90}:{10 * math.tan(pbeta/2.)}cm); %% {a} {b} {length}\n"
    else:
        arc_txt += f"\\draw[line width = 0.0mm] {style} ({a}) to [bend left={2*deg:.1f}] ({b});\n"

    return arc_txt

def drawarc_counterclockwise(a, b, deg, style, length, lengthfix, rotate):
    a,b = b,a
    angle_a = 360./(length+lengthfix)*a
    angle_b = 360./(length+lengthfix)*b
    angle_a = 450 - angle_a -20
    angle_b = 450 - angle_b -20

    alpha = (angle_b-angle_a)
    pa = angle_a * math.pi / 180.
    pb = angle_b * math.pi / 180.
    palpha = alpha * math.pi / 180.
    arc_txt = ''
    if alpha < 170:
        arc_txt += f"\\draw[line width = 0.0mm] {style} ([shift=({angle_b}:10cm)]0,0) arc ({angle_b + 90}:{angle_a + 270}:{10*math.tan(palpha/2.)}cm); %% {a} {b} {length}\n"
    elif alpha > 190:
        beta = 360. - alpha
        pbeta = beta * math.pi / 180.
        arc_txt += f"\\draw[line width = 0.0mm] {style} ([shift=({angle_a}:10cm)]0,0) arc ({angle_a + 90}:{angle_b - 90}:{10*math.tan(pbeta/2.)}cm); %% {a} {b} {length}\n"
    else:
        arc_txt += f"\\draw[line width = 0.0mm] {style} ({a}) to [bend left={2*deg:.1f}] ({b});\n"

def agree(pres, pref, a, b): ## pres[a] = b
    if pref[a] == b:
        return True
    elif pref.get(a-1,-1) == b or pref.get(a+1,-1) == b:
        return True
    elif pref.get(b-1,-1) == a or pref.get(b+1,-1) == a:
        return True
    else:
        return False

def pair_agree(pres, pref, a, b):
    if pres[a] == b or pres.get(a-1,-1) == b or pres.get(a+1,-1) == b or pres.get(b-1,-1) == a or pres.get(b+1,-1) == a:
        return True
    else:
        return False

doc_preamble = r'''
\documentclass{standalone}
%\usepackage{fullpage}
\pagestyle{empty}
\usepackage{tikz}
\usepackage{tkz-euclide}
\usepackage{siunitx}
\usetikzlibrary{shapes, shapes.multipart}
\usepackage{xcolor}

\usepackage{verbatim}
\usepackage{lipsum}
\begin{document}
'''
doc_ending = "\\end{document}\n"

pic_preamble = r'''
%\hspace{-3cm}
%\resizebox{1.2\textwidth}{!}{
\begin{tikzpicture}[darkstyle/.style={}, scale=2] %circle,draw,fill=gray!10}]
'''
def generate_arc_plot_tex(seq, struct, bpp_dict,
                  counter_clockwise=False,
                  rotate=180-5):
    circular = True # Circular = False wasn't actually implemented
    seq_length = len(seq)

    MINLEN,MAXLEN = 0, 5650
    if seq_length > MAXLEN:
        raise Exception("Seq too long (%d)" % (seq_length))
    if seq_length < MINLEN:
        raise Exception("Seq too short (%d)" % (seq_length))

    ref = '.'*seq_length

    txt = ''
    txt += doc_preamble

    notes = ""

    lbs = ['(', '[', '{', '<']
    rbs = [')', ']', '}', '>']

    pairs = []
    goldpairs = []
    pairset = set()
    respair = defaultdict(lambda: -1)
    goldpair = defaultdict(lambda: -1)

    stacks = []
    for _ in range(len(lbs)):
        stacks.append([])
    for i, item in enumerate(struct):
        if item in lbs:
            stackindex = lbs.index(item)
            stacks[stackindex].append(i)
        elif item in rbs:
            stackindex = rbs.index(item)
            left = stacks[stackindex][-1]
            stacks[stackindex] = stacks[stackindex][:-1]
            pairs.append((left+1,i+1, stackindex))
            respair[left+1] = i+1
            respair[i+1] = left+1
            pairset.add((left+1,i+1))
    notes += ";pair=%d" % (len(respair)//2)

    stacks = []
    for _ in range(len(lbs)):
        stacks.append([])
    for i, item in enumerate(ref):
        if item in lbs:
            stackindex = lbs.index(item)
            stacks[stackindex].append(i)
        elif item in rbs:
            stackindex = rbs.index(item)
            left = stacks[stackindex][-1]
            stacks[stackindex] = stacks[stackindex][:-1]
            goldpairs.append((left+1,i+1, stackindex))
            goldpair[left+1] = i+1
            goldpair[i+1] = left+1

    bases = [''] + list(seq)

    txt += pic_preamble + "\n"

    lengthfix = int(seq_length/9.0)
    for i, base in enumerate(bases[1:], 1):
        if circular:
            angle = 360./(seq_length+lengthfix)*i
            if counter_clockwise: # hzhang
                angle = 450-angle -20
            else:
                angle = angle-70 + rotate

            txt += "\\node [darkstyle]  (%d) at (%f:10cm) {};\n" % (i, angle)
            if seq_length <= 100:
                gap = 5
            elif seq_length <= 200:
                gap = 10
            elif seq_length <= 300:
                gap = 20
            elif seq_length <= 700:
                gap = 50
            elif seq_length <= 2000:
                gap = 100
            elif seq_length <= 3000:
                gap = 200
            elif seq_length <= 5000:
                gap = 300
            else:
                gap = 400

            if i % gap == 0 and i < len(bases)-10:
                txt += "\\node [scale=2]           (%d,1) at (%f:10.8cm) {\Huge %d};\n" % (i, angle, i)
            if i > 1:
                txt += "\\draw (%d.center) -- (%d.center);\n" % (i, i-1)

        else:
            txt += "\\node [darkstyle]  (%d) at (%d,0) {};\n" % (i, i)
            if i % 5 == 0:
                txt += "\\node []           (%d,1) at (%d,-1) {%d};\n" % (i, i, i)

    if circular:
        for j in range(4):
            i += 1
            angle = 360./(seq_length+lengthfix)*(i-4)

            if counter_clockwise: # hzhang
                angle = 450-angle -20
            else:
                angle = angle-70 + rotate

            txt += "\\node [darkstyle]  (%d) at (%f:10cm) {};\n" % (i, angle)

        angle = 360./(seq_length+lengthfix)
        angle = 450-angle
        if seq_length <= 50:
            angle5 = angle - 30
            angle3 = angle + 30
        elif seq_length <= 100:
            angle5 = angle - 20
            angle3 = angle + 20
        elif seq_length <= 200:
            angle5 = angle - 17
            angle3 = angle + 17
        else:
            angle5 = angle - 13
            angle3 = angle + 13
        txt += "\\node [scale=2](3prime) at (%f:10cm) {\LARGE \\textbf{%s}};\n" % (angle3, "5'")
        txt += "\\node [right=9.5cm of 3prime, scale=2] {\LARGE \\textbf{%s}};\n" % "3'"


    # add for bpp
    #####################################
    for a, b, prob in bpp_dict:
        color = "red" # not in gold
        if pair_agree(goldpair, respair, a, b): # in gold
            color = "blue"
        lw = ",thick"
        prob *= 100
        color += "!" + str(prob)
        style = "[" + color + lw + "]"

        if circular:
            dist = b - a
            revdist = seq_length+lengthfix - dist
            deg = 90 * (0.5-(dist+.0) / (seq_length+lengthfix+.0))
        else:
            deg = 20

        if counter_clockwise:
            txt += drawarc_counterclockwise(a, b, deg, style, seq_length, lengthfix, rotate)
        else:
            txt += drawarc_clockwise(a, b, deg, style, seq_length, lengthfix, rotate)

    txt += "\\end{tikzpicture}\n"

    txt += doc_ending
    return txt