from math import e
import sympy as sp
import math


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[90m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    # Background colors:
    GREYBG = '\033[100m'
    REDBG = '\033[101m'
    GREENBG = '\033[102m'
    YELLOWBG = '\033[103m'
    BLUEBG = '\033[104m'
    PINKBG = '\033[105m'
    CYANBG = '\033[106m'

def EvaluateError(startPoint, endPoint):
    """
    This function Helps us find out if we can find a root in a limited amount of tries in a specific range.
    :param startPoint: start of range.
    :param endPoint: end of range.
    :return:
    """
    exp = pow(10, -10)
    if endPoint - startPoint == 0:
        return 100
    return ((-1)*math.log(exp/(endPoint-startPoint), e))/math.log(2, e)


def NewtonsMethod(f, x0, tries=100, epsylon=0.0001, symbole=sp.symbols('x')):
    """
    This function find a root to a function by using the newton raphson method by a given first guess.
    :param f: The function with sympy symbols.
    :param x0: The first guess.
    :param tries: Number of tries to find the root.
    :param symbole: The symbol you entered in the function (Default is lower case x)
    :param epsylon: The tolerance of the deviation of the solution ;
    How precise you want the solution (the smaller the better).
    :return:Returns the local root by raphson method ,
    if it fails to find a solutions in the given tries limit it will return None .
    """
    if f.subs(symbole, x0) == 0:
        return 0
    for i in range(tries):
        print(bcolors.OKBLUE, "Attempt #", i + 1, ":", bcolors.ENDC)
        print("f({0}) = {1} = {2}".format(x0, f, round(f.subs(symbole, x0), 2)))
        print("f'({0}) = {1} = {2}".format(x0, sp.diff(f, symbole),
                                           round(sp.diff(f, symbole).subs(symbole, x0), 2)))
        if sp.diff(f, symbole).subs(symbole, x0) == 0.0:
            continue
        next_x = (x0 - f.subs(symbole, x0) / sp.diff(f, symbole).subs(symbole, x0))

        print("next_X = ", round(next_x, 2))
        t = abs(next_x - x0)
        if t < epsylon:
            print(bcolors.OKGREEN, "Found a Root Solution ; X =", round(next_x, 8), bcolors.ENDC)
            return next_x
        x0 = next_x
    print(bcolors.FAIL, "Haven't Found a Root Solution ; (returning None)", bcolors.ENDC)
    return None


