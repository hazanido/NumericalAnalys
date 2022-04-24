def printMatrix(matrix):
    for line in matrix:
        print('  '.join(map(str, line)))
    print("\n")

def print_mul_format(lmat, rmat, sol):
    assert len(lmat) == len(rmat) and len(lmat) == len(sol), 'Error'
    l = []
    r = []
    s = []
    for row in range(len(lmat)):
        rl = []
        rr = []
        sr = []
        for col in range(len(lmat)):
            rl.append(float(round(lmat[row][col], 2)))
            rr.append(float(round(rmat[row][col], 2)))
            sr.append(float(round(sol[row][col], 2)))
        l.append(rl)
        r.append(rr)
        s.append(sr)
    for line in range(len(lmat)):
        if line == len(lmat)//2:
            print(f'{l[line]} * {r[line]} = {s[line]}')
        else:
            print(f'{l[line]}   {r[line]}   {s[line]}')


def copy_matrix(mat):
    new_mat = []
    for r in range(len(mat)):
        row = []
        for c in range(len(mat)):
            row.append(mat[r][c])
        new_mat.append(row)
    return new_mat


def Identity(n):
    mat = [([0] * n) for i in range(n)]  # initialize the matrix with zeros
    for i in range(0, n):
        mat[i][i] = 1  # the identity matrix includes 1 all over its diagonal, starts at [0][0]
    return mat


def Matrix_multiplication(mat1, mat2):
    if len(mat1[0]) != len(mat2):
        raise Exception("Illegal multiplication between matrix's ")
    result_mat = [([0] * len(mat2[0])) for i in range(len(mat1))]  # initialize the result matrix with zeros
    # iterate through the first matrix rows
    for row1 in range(0, len(mat1)):
        # iterate through the second matrix columns
        for col2 in range(0, len(mat2[0])):
            # iterate through the second matrix rows
            for row2 in range(0, len(mat2)):
                result_mat[row1][col2] += mat1[row1][row2] * mat2[row2][col2]
    return result_mat


def ResetMember(row, col, n, pivot, a):
    elementary_matrix = Identity(n)
    elementary_matrix[row][col] = -(a / pivot)
    return elementary_matrix


def MultiplyRow(row, a, n):
    elementary_matrix = Identity(n)
    elementary_matrix[row][row] = a
    return elementary_matrix


def ExchangeRows(row1, row2, n):
    elementary_matrix = Identity(n)
    elementary_matrix[row1][row1] = 0
    elementary_matrix[row1][row2] = 1
    elementary_matrix[row2][row2] = 0
    elementary_matrix[row2][row1] = 1
    return elementary_matrix


def FinalVector(matrix, b):
    n = len(matrix)
    for j in range(0, n):
        for i in range(0, n):
            if i == j:
                pivot = matrix[i][j]
                for k in range(i + 1, n):
                    if abs(matrix[k][j]) > abs(pivot):  # pivoting
                        originalmatrix = copy_matrix(matrix)
                        elementary_matrix = ExchangeRows(k, i, n)
                        matrix = Matrix_multiplication(elementary_matrix, matrix)
                        pivot = matrix[i][j]
                        b = Matrix_multiplication(elementary_matrix, b)
                        print_mul_format(elementary_matrix, originalmatrix, matrix)
                        print('\n\n')


        for i in range(0, n):
            if i > j:
                if matrix[i][j] != 0:
                    originalmatrix=copy_matrix(matrix)
                    elementary_matrix = ResetMember(i, j, n, pivot, matrix[i][j])
                    matrix = Matrix_multiplication(elementary_matrix, matrix)
                    b = Matrix_multiplication(elementary_matrix, b)
                    print_mul_format(elementary_matrix, originalmatrix, matrix)
                    print('\n\n')


        for i in range(0, n):
            if i < j:
                if matrix[i][j] != 0:
                    originalmatrix = copy_matrix(matrix)
                    elementary_matrix = ResetMember(i, j, n, pivot, matrix[i][j])
                    matrix = Matrix_multiplication(elementary_matrix, matrix)
                    b = Matrix_multiplication(elementary_matrix, b)
                    print_mul_format(elementary_matrix, originalmatrix, matrix)
                    print('\n\n')

    for i in range(0, n):
        if matrix[i][i] != 1:
            if matrix[i][i] < 0:
                originalmatrix = copy_matrix(matrix)
                elementary_matrix = MultiplyRow(i, -1, n)
                matrix = Matrix_multiplication(elementary_matrix, matrix)
                b = Matrix_multiplication(elementary_matrix, b)
                print_mul_format(elementary_matrix, originalmatrix, matrix)
                print('\n\n')

            originalmatrix = copy_matrix(matrix)
            elementary_matrix = MultiplyRow(i, 1 / matrix[i][i], n)
            matrix = Matrix_multiplication(elementary_matrix, matrix)
            b=Matrix_multiplication(elementary_matrix, b)
            print_mul_format(elementary_matrix, originalmatrix, matrix)
            print('\n')
    return b


def getInverseMatrix(matrix):
    n = len(matrix)
    inverseMatrix = Identity(n)

    for j in range(0, n):
        for i in range(0, n):
            if i == j:
                pivot = matrix[i][j]
                for k in range(i + 1, n):
                    if abs(matrix[k][j]) > abs(pivot):  # pivoting
                        originalmatrix = copy_matrix(matrix)
                        elementary_matrix = ExchangeRows(k, i, n)
                        matrix = Matrix_multiplication(elementary_matrix, matrix)
                        pivot = matrix[i][j]
                        inverseMatrix = Matrix_multiplication(elementary_matrix, inverseMatrix)

        for i in range(0, n):
            if i > j:
                if matrix[i][j] != 0:
                    originalmatrix = copy_matrix(matrix)
                    elementary_matrix = ResetMember(i, j, n, pivot, matrix[i][j])
                    matrix = Matrix_multiplication(elementary_matrix, matrix)
                    inverseMatrix = Matrix_multiplication(elementary_matrix, inverseMatrix)

        for i in range(0, n):
            if i < j:
                if matrix[i][j] != 0:
                    originalmatrix = copy_matrix(matrix)
                    elementary_matrix = ResetMember(i, j, n, pivot, matrix[i][j])
                    matrix = Matrix_multiplication(elementary_matrix, matrix)
                    inverseMatrix = Matrix_multiplication(elementary_matrix, inverseMatrix)

    for i in range(0, n):
        if matrix[i][i] != 1:
            if matrix[i][i] < 0:
                originalmatrix = copy_matrix(matrix)
                elementary_matrix = MultiplyRow(i, -1, n)
                matrix = Matrix_multiplication(elementary_matrix, matrix)
                inverseMatrix = Matrix_multiplication(elementary_matrix, inverseMatrix)


            originalmatrix = copy_matrix(matrix)
            elementary_matrix = MultiplyRow(i, 1 / matrix[i][i], n)
            matrix = Matrix_multiplication(elementary_matrix, matrix)
            inverseMatrix = Matrix_multiplication(elementary_matrix, inverseMatrix)

    for r in range(n):
        for c in range(n):
            inverseMatrix[r][c] = (float(round(inverseMatrix[r][c], 2)))
    return inverseMatrix



#Main
'''
R = int(input("Enter the number of rows:"))
C = int(input("Enter the number of columns:"))


#matrix
matrix = []
print("Enter all the matrix:")
for i in range(R):
    a = []
    for j in range(C):
        a.append(float(input()))
    matrix.append(a)


#vector b
b = []
print('Enter vector b: ')
for i in range(R):
    br = []
    for j in range(1):
        br.append(float(input()))
    b.append(br)

b = FinalVector(matrix, b)


#answer
ans=[]
for i in range(R):
    ans.append(float(round(b[i][0],2)))
print(f'x: {ans}')
'''