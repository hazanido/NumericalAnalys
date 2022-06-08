import numpy as np


# Limits
A = 0
B = np.pi

# Function
F = np.sin

# Maximum of absolute value of function in [A, B]
M = 1

# Number of partitions
N = 4

# ===== END OF INPUT =====


def main():
    print(f'FUNCTION: {F.__name__}\n')
    # Calculate integrals
    step = (B - A) / N
    for method in (trapezoid_method, simpson_method):
        total_sum = 0
        for i in np.arange(A, B, step):
            ai, bi = i, i + step
            total_sum += method(ai, bi, F)

        print(f'{method.__name__}')
        print(f'\tValue: {total_sum:.7f}')
        print(f'\tTheoretical estimated error: '
              f'{theoretical_error_estimation(M, step, A, B, method.__name__):.7f}')




def rename(newname):
    def decorator(f):
        f.__name__ = newname
        return f
    return decorator


@rename('Trapezoid method')
def trapezoid_method(a, b, func):
    return abs((b-a) * (func(a) + func(b)) / 2)


@rename('Simpson method')
def simpson_method(a, b, func):
    return abs((b - a) * (func(a) + 4 * func((a + b) / 2) + func(b)) / 6)


def theoretical_error_estimation(m, h, a, b, type_):
    if type_ == 'Trapezoid method':
        return m * (h ** 2) * (b - a) / 12
    elif type_ == 'Simpson method':
        return m * (h ** 2) * (b - a) / 2880


if __name__ == '__main__':
    main()