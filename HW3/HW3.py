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


def NewtonsMethodInRangeIterations(f, check_range, tries=10, epsilon=0.0001, symbol=sp.symbols('x')):
    """
    This function find a root to a function by using the newton raphson method by a given list of guesses.
    :param f: The function with sympy symbols.
    :param check_range: List of guesses.
    :param tries: Number of tries to find the root.
    :param symbole: The symbol you entered in the function (Default is lower case x)
    :param epsylon: The tolerance of the deviation of the solution ;
    How precise you want the solution (the smaller the better).
    :return:Returns a list roots by raphson method ,
    if it fails to find a solutions in the given tries limit it will return an empty list .
    """
    roots = []
    for i in check_range:
        if i == check_range[-1]:
            break
        for sep in range(0, 10):
            check_number = round(i + (sep * 0.1), 2)
            print(bcolors.HEADER, "First guess:", check_number , bcolors.ENDC)
            local_root = NewtonsMethod(f, check_number, tries, epsilon, symbol)
            if local_root is not None and round(local_root,6) not in roots:
                roots+=[round(local_root,6)]
            else:
                print(bcolors.FAIL, "Already found that root.", bcolors.ENDC)
    return roots


def SecantMethodInRangeIterations(f, check_range, epsilon=0.0001):
    """
    This function find a root to a function by using the secant method by a given list of values to check beetween.
    :param f: The function (as a python function).
    :param check_range: List of values to check between ; e.g (1,2,3,4,5) it will check between 1-2,2-3,....
    :param epsylon: The tolerance of the deviation of the solution ;
    How precise you want the solution (the smaller the better).
    :return:Returns a list roots by secant method ,
    if it fails to find a solutions in the given tries limit it will return an empty list .
    """
    roots = []
    iterCounter = 0
    for i in check_range:
        if i == check_range[-1]:
            break
        for sep in range(0, 10):
            startPoint = round(i + (sep * 0.1), 2)
            endPoint = round(i + ((sep+1) * 0.1), 2)
            print(bcolors.HEADER, "Checked range:", startPoint, "-",endPoint, bcolors.ENDC)
            local_root = SecantMethod(f, startPoint, endPoint, epsilon, iterCounter)
            if local_root is not None and round(local_root,6) not in roots:
                roots += [round(local_root,6)]
            else:
                print(bcolors.FAIL, "Already found that root.",bcolors.ENDC)
    return roots


def SecantMethod(polynomial,firstGuess, secondGuess,epsilon, iterCounter):
    """
     This function find a root to a function by using the SecantMethod method by a given tow guess.
    :param polynomial: The function on which the method is run
    :param firstGuess: The first guess
    :param secondGuess: The second guess
    :param epsilon: The tolerance of the deviation of the solution
    :param iterCounter: number of tries until the function found the root.
    :return:Returns the local root by Secant method ,
    if it fails to find a solutions in the given tries limit it will return None .
    """
    if iterCounter > 100:
        return

    if abs(secondGuess - firstGuess) < epsilon: #Stop condition
        print("after ", iterCounter, "iterations The root found is: ", bcolors.OKBLUE, round(secondGuess, 6), bcolors.ENDC)
        return round(secondGuess, 6) # Returns the root found

    next_guess = (firstGuess * polynomial(secondGuess) - secondGuess * polynomial(firstGuess))/(polynomial(secondGuess) - polynomial(firstGuess))

    return SecantMethod(polynomial, secondGuess, next_guess, epsilon, iterCounter+1)

