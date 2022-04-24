from Bohan_1 import Gaussian_Elimination

def ZeroMatrix(n):
    mat = [([0] * n) for i in range(n)]  # initialize the matrix with zeros
    return mat


def printMatrix(matrix):
    for line in matrix:
        print('  '.join(map(str, line)))
    print("\n")


def strongU(matrix):
    mat = ZeroMatrix(len(matrix))
    for r in range(len(matrix)):
        for c in range(len(matrix)):
            if c > r:
                mat[r][c] = matrix[r][c]
    return mat


def strongD(matrix):
    mat = ZeroMatrix(len(matrix))
    for r in range(len(matrix)):
        mat[r][r] = matrix[r][r]
    return mat


def strongL(matrix):
    mat = ZeroMatrix(len(matrix))
    for r in range(len(matrix)):
        for c in range(len(matrix)):
            if c < r:
                mat[r][c] = matrix[r][c]
    return mat


def sumMatrix(mat1, mat2):
    assert(len(mat1) == len(mat2)), "matrix not at the same size"
    n = len(mat1)
    mat = Gaussian_Elimination.Identity(n)
    for r in range(n):
        for c in range(n):
            mat[r][c] = mat1[r][c] + mat2[r][c]
    return mat

def isDomDiagonal(mat):
    n = len(mat)
    for r in range(n):
        sumline = 0
        pivot = abs(mat[r][r])
        for c in range(n):
            if c != r:
                sumline += abs(mat[r][c])
        if sumline >= pivot:
            return False
    return True

def negativeMatrix(mat):
    n = len(mat)
    for r in range(n):
        for c in range(n):
            mat[r][c] = -mat[r][c]
    return mat

def getG(matrix):
    return Gaussian_Elimination.Matrix_multiplication(negativeMatrix(Gaussian_Elimination.getInverseMatrix(sumMatrix(strongL(matrix), strongD(matrix)))), strongU(matrix))

def getH(matrix):
    return Gaussian_Elimination.getInverseMatrix(sumMatrix(strongL(matrix), strongD(matrix)))

def DominantDiagonalFix(matrix):
    n = len(matrix)
    dom = [0] * n
    domDiagonalMatrix = []
    for r in range(len(matrix)):
        for c in range(n):
            if (matrix[r][c] > sum(map(abs,matrix[r]))-matrix[r][c]):
                dom[r] = c
                break
    for r in range(n):
        domDiagonalMatrix.append([])
        if r not in dom:
            return matrix
    for r in range(n):
        domDiagonalMatrix[dom[r]] = matrix[r]
    return domDiagonalMatrix

#   #   #   #   #   #   #   #   #   #   #   #   #   #

#matrix = [[1,2,3],[0,1,4],[5,6,0]]
matrix = [[4,1,9],[9,2,4],[3,9,1]]

print("H: ")
matrixH = getH(matrix)
printMatrix(matrixH)

print("G: ")
matrixG = getG(matrix)
printMatrix(matrixG)

if isDomDiagonal(matrix):
    print("Inverse: ")
    printMatrix(Gaussian_Elimination.getInverseMatrix(matrix))
else:
    if isDomDiagonal(DominantDiagonalFix(matrix)):
        mat = DominantDiagonalFix(matrix)
        print("mat: ")
        printMatrix(mat)
        print("Inverse: ")
        printMatrix(Gaussian_Elimination.getInverseMatrix(mat))
    else:
        print("diagonal not dominant")