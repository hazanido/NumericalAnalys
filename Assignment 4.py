"""
 * Authors: Ofek Elgozi (ID: 318432085), Ido Hazan (ID:316613769), Stav Sharabi (ID:208163840)
"""

from sympy import tan
from math import e
import sympy as sp
import math


def PrintMatrix(matrix):
    """
    Matrix Printing Function
    :param matrix: Matrix nxn
    """
    for line in matrix:
        line.append('|')
        line.insert(0, '|')
        print('  '.join(map(str, line)))
        line.remove('|', )
        line.remove('|', )


def Determinant(matrix, mul):
    """
    Recursive function for determinant calculation
    :param matrix: Matrix nxn
    :param mul: The double number
    :return: determinant of matrix
    """
    width = len(matrix)
    # Stop Conditions
    if width == 1:
        return mul * matrix[0][0]
    else:
        sign = -1
        det = 0
        for i in range(width):
            m = []
            for j in range(1, width):
                buff = []
                for k in range(width):
                    if k != i:
                        buff.append(matrix[j][k])
                m.append(buff)
            # Change the sign of the multiply number
            sign *= -1
            #  Recursive call for determinant calculation
            det = det + mul * Determinant(m, sign * matrix[0][i])
    return det


def MaxNorm(matrix):
    """
    Function for calculating the max-norm of a matrix
    :param matrix: Matrix nxn
    :return:max-norm of a matrix
    """
    max_norm = 0
    for i in range(len(matrix)):
        norm = 0
        for j in range(len(matrix)):
            # Sum of organs per line with absolute value
            norm += abs(matrix[i][j])
        # Maximum row amount
        if norm > max_norm:
            max_norm = norm

    return max_norm


def MultiplyMatrix(matrixA, matrixB):
    """
    Function for multiplying 2 matrices
    :param matrixA: Matrix nxn
    :param matrixB: Matrix nxn
    :return: Multiplication between 2 matrices
    """
    # result matrix initialized as singularity matrix
    result = [[0 for y in range(len(matrixB[0]))] for x in range(len(matrixA))]
    for i in range(len(matrixA)):
        # iterate through columns of Y
        for j in range(len(matrixB[0])):
            # iterate through rows of Y
            for k in range(len(matrixB)):
                result[i][j] += matrixA[i][k] * matrixB[k][j]
    return result


def MakeIMatrix(cols, rows):
    # Initialize a identity matrix
    return [[1 if x == y else 0 for y in range(cols)] for x in range(rows)]


def InverseMatrix(matrix, vector):
    """
    Function for calculating an inverse matrix
    :param matrix:  Matrix nxn
    :return: Inverse matrix
    """
    # Unveri reversible matrix
    if Determinant(matrix, 1) == 0:
        print("Error,Singular Matrix\n")
        return
    # result matrix initialized as singularity matrix
    result = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # turn the pivot into 1 (make elementary matrix and multiply with the result matrix )
        # pivoting process
        matrix, vector = RowXchange(matrix, vector)
        elementary = MakeIMatrix(len(matrix[0]), len(matrix))
        elementary[i][i] = 1 / matrix[i][i]
        result = MultiplyMatrix(elementary, result)
        matrix = MultiplyMatrix(elementary, matrix)
        # make elementary loop to iterate for each row and subtracrt the number below (specific) pivot to zero  (make
        # elementary matrix and multiply with the result matrix )
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            elementary[j][i] = -(matrix[j][i])
            matrix = MultiplyMatrix(elementary, matrix)
            result = MultiplyMatrix(elementary, result)

    # after finishing with the lower part of the matrix subtract the numbers above the pivot with elementary for loop
    # (make elementary matrix and multiply with the result matrix )
    for i in range(len(matrix[0]) - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            elementary[j][i] = -(matrix[j][i])
            matrix = MultiplyMatrix(elementary, matrix)
            result = MultiplyMatrix(elementary, result)

    return result


def RowXchange(matrix, vector):
    """
    Function for replacing rows with both a matrix and a vector
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Replace rows after a pivoting process
    """

    for i in range(len(matrix)):
        max = abs(matrix[i][i])
        for j in range(i, len(matrix)):
            # The pivot member is the maximum in each column
            if abs(matrix[j][i]) > max:
                temp = matrix[j]
                temp_b = vector[j]
                matrix[j] = matrix[i]
                vector[j] = vector[i]
                matrix[i] = temp
                vector[i] = temp_b
                max = abs(matrix[i][i])

    return [matrix, vector]


def GaussJordanElimination(matrix, vector):
    """
    Function for moding a linear equation using gauss's elimination method
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=A(-1)b
    """
    # Pivoting process
    matrix, vector = RowXchange(matrix, vector)
    # Inverse matrix calculation
    invert = InverseMatrix(matrix, vector)
    return MulMatrixVector(invert, vector)


def MulMatrixVector(InversedMat, b_vector):
    """
    Function for multiplying a vector matrix
    :param InversedMat: Matrix nxn
    :param b_vector: Vector n
    :return: Result vector
    """
    result = []
    # Initialize the x vector
    for i in range(len(b_vector)):
        result.append([])
        result[i].append(0)
    # Multiplication of inverse matrix in the result vector
    for i in range(len(InversedMat)):
        for k in range(len(b_vector)):
            result[i][0] += InversedMat[i][k] * b_vector[k][0]
    return result


def UMatrix(matrix, vector):
    """
    :param matrix: Matrix nxn
    :return:Disassembly into a  U matrix
    """
    # result matrix initialized as singularity matrix
    U = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # pivoting process
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            # Finding the M(ij) to reset the organs under the pivot
            elementary[j][i] = -(matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)
    # U matrix is a doubling of elementary matrices that we used to reset organs under the pivot
    U = MultiplyMatrix(U, matrix)
    return U


def LMatrix(matrix, vector):
    """
       :param matrix: Matrix nxn
       :return:Disassembly into a  L matrix
       """
    # Initialize the result matrix
    L = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # pivoting process
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            # Finding the M(ij) to reset the organs under the pivot
            elementary[j][i] = -(matrix[j][i]) / matrix[i][i]
            # L matrix is a doubling of inverse elementary matrices
            L[j][i] = (matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)

    return L


def RowXchageZero(matrix, vector):
    """
      Function for replacing rows with both a matrix and a vector
      :param matrix: Matrix nxn
      :param vector: Vector n
      :return: Replace rows after a pivoting process
      """

    for i in range(len(matrix)):
        for j in range(i, len(matrix)):
            # The pivot member is not zero
            if matrix[i][i] == 0:
                temp = matrix[j]
                temp_b = vector[j]
                matrix[j] = matrix[i]
                vector[j] = vector[i]
                matrix[i] = temp
                vector[i] = temp_b

    return [matrix, vector]


def SolveLU(matrix, vector):
    """
    Function for deconstructing a linear equation by ungrouping LU
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=U(-1)L(-1)b
    """
    matrixU = UMatrix(matrix)
    matrixL = LMatrix(matrix)
    return MultiplyMatrix(InverseMatrix(matrixU), MultiplyMatrix(InverseMatrix(matrixL), vector))


def Cond(matrix, invert):
    """
    :param matrix: Matrix nxn
    :param invert: Inverted matrix
    :return: CondA = ||A|| * ||A(-1)||
    """
    print("\n|| A ||max = ", MaxNorm(matrix))
    print("\n|| A(-1) ||max = ", MaxNorm(invert))
    return MaxNorm(matrix) * MaxNorm(invert)


def solveMatrix(matrixA, vectorb):
    detA = Determinant(matrixA, 1)
    print("\nMatrix A: \n")
    PrintMatrix(matrixA)
    print("\nVector b: \n")
    PrintMatrix(vectorb)
    print("\nDET(A) = ", detA)
    if detA != 0:
        print("\nCondA = ", Cond(matrixA, InverseMatrix(matrixA, vectorb)))
        print("\nGaussJordanElimination")
        print("Solve Ax = b \n")
        result = GaussJordanElimination(matrixA, vectorb)
        PrintMatrix(result)
        return result
    else:
        print("Singular Matrix\n")
        print("\nLU Decomposition \n")
        print("Matrix U: \n")
        PrintMatrix(UMatrix(matrixA, vectorb))
        print("\nMatrix L: \n")
        PrintMatrix(LMatrix(matrixA, vectorb))
        print("\nMatrix A=LU: \n")
        result = MultiplyMatrix(LMatrix(matrixA, vectorb), UMatrix(matrixA, vectorb))
        PrintMatrix(result)
        return result


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


def LinearInterpolation(others_table_points, point):
    p = []
    result = 0
    flag = 1
    for i in range(len(others_table_points)):
        p.append(others_table_points[i][0])
    for i in range(len(p) - 1):
        if i <= point <= i + 1:
            x1 = others_table_points[i][0]
            x2 = others_table_points[i + 1][0]
            y1 = others_table_points[i][1]
            y2 = others_table_points[i + 1][1]
            result = (((y1 - y2) / (x1 - x2)) * point) + ((y2 * x1) - (y1 * x2)) / (x1 - x2)
            print(bcolors.OKBLUE, "(LinearInterpolation) The approximation (interpolation) of the point ", point,
                  " is: ", round(result, 4),
                  bcolors.ENDC)
            flag = 0
    if flag:
        x1 = others_table_points[0][0]
        x2 = others_table_points[1][0]
        y1 = others_table_points[0][1]
        y2 = others_table_points[1][1]
        m = (y1 - y2) / (x1 - x2)
        result = y1 + m * (point - x1)
        print(bcolors.OKBLUE, "(LinearInterpolation) The approximation (extrapolation) of the point ", point, " is: ",
              round(result, 4),
              bcolors.ENDC)


def PolynomialMethod(others_table_points, x):
    matrix = [[point[0] ** i for i in range(len(others_table_points))] for point in
              others_table_points]  # Makes the initial matrix
    b = [[point[1]] for point in others_table_points]

    print(bcolors.OKBLUE, "The matrix is- ", bcolors.ENDC)
    PrintMatrix(matrix)
    print(bcolors.OKBLUE, "b vector is- ", bcolors.ENDC)
    print(b)
    print(bcolors.OKBLUE, "The Solotion:", bcolors.ENDC)
    matrixSol = solveMatrix(matrix, b)

    result = sum([matrixSol[i][0] * (x ** i) for i in range(len(matrixSol))])
    print("\n")
    print('+'.join([' (' + str(matrixSol[i][0]) + ') * x ^ ' + str(i) for i in range(len(matrixSol))]), '= ')
    print(bcolors.OKBLUE, "(PolynomialMethod) The approximation (interpolation) of the point ", x, " is: ",
          round(result, 4), bcolors.ENDC)
    return result


def NevAlgorithm(others_table_points, x):
    tot = 0
    size = len(others_table_points)
    temp_table = [[point[0], point[1]] for point in others_table_points]
    for j in range(1, size):
        for i in range(size - 1, j - 1, -1):
            temp_table[i][1] = ((x - temp_table[i - j][0]) * temp_table[i][1] - (x - temp_table[i][0]) *
                                temp_table[i - 1][1]) / (temp_table[i][0] - temp_table[i - j][0])
    tot = temp_table[size - 1][1]
    print(bcolors.OKBLUE, "(NevAlgorithm) The approximation (interpolation) of the point ", x, " is: ", round(tot, 4),
          bcolors.ENDC)
    return tot


def LagrangeMethod(others_table_points, point):
    xp = []
    yp = []
    result = 0
    for i in range(len(others_table_points)):
        xp.append(others_table_points[i][0])
        yp.append(others_table_points[i][1])
    for i in range(len(others_table_points)):
        lagrange_i = 1
        for j in range(len(others_table_points)):
            if i != j:
                lagrange_i = lagrange_i * (point - xp[j]) / (xp[i] - xp[j])
        result += lagrange_i * yp[i]

    print(bcolors.OKBLUE, "(LagrangeMethod) The approximation (interpolation) of the point ", point, " is: ",
          round(result, 4), bcolors.ENDC)


def cubic_spilne_solveMatrix(matrix, size):
    for i in range(size):
        # preprocess the matrix
        currentPivot = abs(matrix[i][i])
        maxRow = i
        row = i + 1
        while row < size:
            if abs(matrix[row][i]) > currentPivot:
                currentPivot = abs(matrix[row][i])
                maxRow = row
            row += 1
        matrix = cubic_spilne_swapRows(matrix, i, maxRow, size)

        matrix = cubic_spilne_PivotToOne(matrix, i, i, size)
        row = i + 1
        while row < size:
            matrix = cubic_spilne_nullify(matrix, row, i, size, matrix[i][i])
            row += 1
    for i in range(1, size):
        row = i - 1
        while row >= 0:
            matrix = cubic_spilne_nullify(matrix, row, i, size, matrix[i][i])
            row -= 1

    return matrix


def cubic_spilne_identityMatrix(size):
    matrix = []
    for i in range(size):
        matrix.append([])
        for j in range(size):
            matrix[i].append(0)

    for i in range(size):
        matrix[i][i] = 1
    return matrix


def cubic_spilne_swapRows(matrix, row1, row2, size):
    identity = cubic_spilne_identityMatrix(size)
    identity[row1], identity[row2] = identity[row2], identity[row1]
    return cubic_spilne_multiplyMatrices(identity, matrix, size)


def cubic_spilne_editMatrix(matrix, x, y, val):
    matrix[x][y] = val
    return matrix


# assuming nxn+1 matrix
def cubic_spilne_multiplyMatrices(matrix1, matrix2, size):
    mat = []
    for i in range(size + 1):
        mat.append([])

    for i in range(size):
        for j in range(size + 1):
            sum = 0
            for k in range(size):
                sum += matrix1[i][k] * matrix2[k][j]
            mat[i].append(sum)
    return mat


def cubic_spilne_nullify(matrix, x, y, size, pivot):
    identity = cubic_spilne_identityMatrix(size)
    return cubic_spilne_multiplyMatrices(cubic_spilne_editMatrix(identity, x, y, -1 * matrix[x][y] / pivot), matrix,
                                         size)


def cubic_spilne_PivotToOne(matrix, x, y, size):
    identity = cubic_spilne_identityMatrix(size)
    return cubic_spilne_multiplyMatrices(cubic_spilne_editMatrix(identity, x, y, 1 / matrix[x][y]), matrix, size)


class cubic_spilne_Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def cubic_spline_interpolation(table, target):
    gamma = []
    mu = []
    d = []
    h = []

    for i in range(0, len(table) - 1):
        h.append(table[i + 1].x - table[i].x)

    for i in range(1, len(table) - 1):
        gamma.append(h[i] / (h[i] + h[i - 1]))
        mu.append(1 - h[i] / (h[i] + h[i - 1]))
        d.append(
            (6 / (h[i] + h[i - 1]) * ((table[i + 1].y - table[i].y) / h[i] - (table[i].y - table[i - 1].y) / h[i - 1])))

    # build matrix
    mat = cubic_spilne_identityMatrix(len(d))
    for i in range(len(d)):
        mat[i][i] = 2
        if i != 0:
            mat[i][i - 1] = mu[i]
        if i != len(d) - 1:
            mat[i][i + 1] = gamma[i]
        mat[i].append(d[i])

    # extract result
    m = [0]
    res = cubic_spilne_solveMatrix(mat, len(d))
    for x in range(len(res) - 1):
        m.append(res[x][-1])
    m.append(0)

    for y in range(len(table) - 1):
        if target > table[y].x:
            if target < table[y + 1].x:
                point = ((table[y + 1].x - target) ** 3 * m[y] + (target - table[y].x) ** 3 * m[y + 1]) / (6 * h[y]) \
                        + ((table[y + 1].x - target) * table[y].y + (target - table[y].x) * table[y + 1].y) / (h[y]) \
                        - h[y] * ((table[y + 1].x - target) * m[y] + (target - table[y].x) * m[y + 1]) / 6
                return print(bcolors.OKBLUE, "(Natural Cubic Spline) The approximation (interpolation) of the point ",
                             target, " is: ", round(point, 4))
    print("Target out of bounds")


def MainFunction():
    others_table_points = [(0, 0), (math.pi / 6, 0.5), (math.pi / 4, 0.7072), (math.pi / 2, 1)]
    p0 = cubic_spilne_Point(0, 0)
    p1 = cubic_spilne_Point(math.pi / 6, 0.5)
    p2 = cubic_spilne_Point(math.pi / 4, 0.7072)
    p3 = cubic_spilne_Point(math.pi / 2, 1)
    cubic_spline_table = [p0, p1, p2, p3]
    x = math.pi / 3

    print(bcolors.OKBLUE, "Interpolation & Extrapolation Methods\n", bcolors.ENDC)
    print(bcolors.OKGREEN, "Table Points", others_table_points, bcolors.ENDC)
    print("\n")
    print(bcolors.OKGREEN, "Finding an approximation to the point ", x, bcolors.ENDC)
    print("\n")
    PolynomialMethod(others_table_points, x)
    print("\n")
    LinearInterpolation(others_table_points, x)
    print("\n")
    LagrangeMethod(others_table_points, x)
    print("\n")
    NevAlgorithm(others_table_points, x)
    print("\n")
    cubic_spline_interpolation(cubic_spline_table, x)


MainFunction()




