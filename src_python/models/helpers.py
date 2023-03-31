import math

def quickselect_partition(scores, lower, upper):
    pivot = scores[upper][0]
    while lower < upper:
        while scores[lower][0] < pivot:
            lower += 1
        while scores[upper][0] > pivot:
            upper -= 1
        if scores[lower][0] == scores[upper][0]:
            lower += 1
        elif lower < upper:
            scores[lower], scores[upper] = scores[upper], scores[lower]
    return upper

def quickselect(scores, lower, upper, k):
    if lower == upper:
        return scores[lower][0]
    split = quickselect_partition(scores, lower, upper)
    length = split - lower + 1
    if length == k:
        return scores[split][0]
    elif k < length:
        return quickselect(scores, lower, split - 1, k)
    else:
        return quickselect(scores, split + 1, upper, k - length)

def fast_log_exp_plus_one(x):
    # Bounds for tolerance of 7.05e-06: (0, 11.8625)
    # Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    # Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    # Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    # Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    # Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    # Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    # Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    # Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    # 8 polynomials needed.
    assert 0.0 <= x <= 11.8624794162, "Argument out-of-range."
    if x < 3.3792499610:
        if x < 1.6320158198:
            if x < 0.6615367791:
                return ((-0.0065591595*x+0.1276442762)*x+0.4996554598)*x+0.6931542306
            return ((-0.0155157557*x+0.1446775699)*x+0.4882939746)*x+0.6958092989
        if x < 2.4912588184:
            return ((-0.0128909247*x+0.1301028251)*x+0.5150398748)*x+0.6795585882
        return ((-0.0072142647*x+0.0877540853)*x+0.6208708362)*x+0.5909675829

    if x < 5.7890710412:
        if x < 4.4261691294:
            return ((-0.0031455354*x+0.0467229449)*x+0.7592532310)*x+0.4348794399
        return ((-0.0010110698*x+0.0185943421)*x+0.8831730747)*x+0.2523695427
    if x < 7.8162726752:
        return ((-0.0001962780*x+0.0046084408)*x+0.9634431978)*x+0.0983148903
    return ((-0.0000113994*x+0.0003734731)*x+0.9959107193)*x+0.0149855051

def fast_log_plus_equals(x, y, verb=False):
    if verb:
        print(f"Before: x={x}, y={y}")
    if x < y:
        x, y = y, x
    if y > -math.inf/2 and x-y < 11.8624794162:
        x = fast_log_exp_plus_one(x-y) + y
    if verb:
        print(f"After: x={x}, y={y}")
    return x

def fast_exp(x):
    #return math.exp(x)

    # Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    # Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    # Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    # Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    # Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    # Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    # Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    # 6 polynomials needed.
    if x < -2.4915033807:
        if x < -5.8622823336:
            if x < -9.91152:
                return 0.0
            return ((0.0000803850*x+0.0021627428)*x+0.0194708555)*x+0.0588080014
        if x < -3.8396630909:
            return ((0.0013889414*x+0.0244676474)*x+0.1471290604)*x+0.3042757740
        return ((0.0072335607*x+0.0906002677)*x+0.3983111356)*x+0.6245959221
    if x < -0.6725053211:
        if x < -1.4805375919:
            return ((0.0232410351*x+0.2085645908)*x+0.6906367911)*x+0.8682322329
        return ((0.0573782771*x+0.3580258429)*x+0.9121133217)*x+0.9793091728
    if x < 0:
        return ((0.1199175927*x+0.4815668234)*x+0.9975991939)*x+0.9999505077

    return 1e20 if x > 46.052 else math.exp(x)

# import numpy as np
# r = np.linspace(-10, 10, 5000)
# %timeit for v in r: fast_exp(v)
# %timeit for v in r: np.exp(v)
# %timeit for v in r: math.exp(v)

# r = np.linspace(0.0, 11.8624794162, 5000)
# %timeit for v in r: fast_log_exp_plus_one(v)
# %timeit for v in r: np.log(np.exp(v) + 1)
# %timeit for v in r: math.log(math.exp(v) + 1)