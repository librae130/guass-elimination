from sympy import symbols, nan, Rational


def getMatrixConsole():
    while True:
        try:
            rows = int(input("Enter the number of rows: "))
            cols = int(input("Enter the number of columns: "))
            if rows > 0 and cols > 0:
                break
            else:
                print("Please enter positive integers for dimensions.")
        except ValueError:
            print("Please enter valid integers.")
    print(f"\nEnter the matrix elements ({rows}x{cols}):")
    print("You can enter decimals, fractions (a/b) or integers.")
    print("Enter elements row by row, separated by spaces.\n")
    matrix = []
    for i in range(rows):
        while True:
            try:
                print(f"Row {i+1}: ", end="")
                rowIn = input().strip().split()
                if len(rowIn) != cols:
                    print(f"Please enter exactly {cols} elements.")
                    continue
                row = []
                for e in rowIn:
                    row.append(Rational(e))
                matrix.append(row)
                break
            except (ValueError, ZeroDivisionError):
                print(" Please enter valid numbers.")
    return matrix

# Counts the number of leading zero(s) for each row into a list, returns that list
def countRowLeadingZero(mat):
    leadingZeroCntRow = [0 for x in range(len(mat))]
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if mat[i][j] != 0:
                break
            leadingZeroCntRow[i] += 1
    return leadingZeroCntRow

# Performs guass elimination on a matrix, returns that matrix
def guassElimination(extMat):
    # Counts and stores the number of leading zero(s) for each row
    # The number of leading zero(s) corresponds with the first non-zero index
    leadingZeroCntRow = countRowLeadingZero(extMat)
    # Chooses pivot row 'i'
    for i in range(len(extMat)):
        # Checks whether pivot has the least leading zeros compare to row(s) underneath
        swapIdx = i
        for j in range(i + 1, len(leadingZeroCntRow)):
            if leadingZeroCntRow[i] > leadingZeroCntRow[j]:
                swapIdx = j
        # Swaps pivot and the row with the least leading zeros, also swaps for counting array
        swap = extMat[i]
        extMat[i] = extMat[swapIdx]
        extMat[swapIdx] = swap
        swap = leadingZeroCntRow[i]
        leadingZeroCntRow[i] = leadingZeroCntRow[swapIdx]
        leadingZeroCntRow[swapIdx] = swap
        # Gets value '1' at first non-zero index of pivot
        divisor = 0
        if leadingZeroCntRow[i] < len(extMat[0]):
            divisor = extMat[i][leadingZeroCntRow[i]]
        for j in range(leadingZeroCntRow[i], len(extMat[0])):
            if divisor != 0:
                extMat[i][j] /= divisor
        # Updates for row(s) underneath, also update counting array
        for j in range(i+1, len(extMat)):
            # Only update for row(s) with equal number of leading zero(s) as pivot
            if leadingZeroCntRow[j] == leadingZeroCntRow[i] and leadingZeroCntRow[j] < len(extMat[0]):
                leadingZeroCurCnt = leadingZeroCntRow[j]
                divisor = extMat[j][leadingZeroCntRow[j]] / extMat[i][leadingZeroCntRow[i]]
                for k in range(leadingZeroCntRow[j], len(extMat[0])):
                    extMat[j][k] -= (extMat[i][k]*divisor)
                    # Store and update counting array later if leading zero's count changes
                    if extMat[j][k] == 0 and leadingZeroCurCnt == k:
                        leadingZeroCurCnt += 1
                leadingZeroCntRow[j] = leadingZeroCurCnt
    return extMat

# Solves linear equations from an extended matrix (A|B) with Ax=B, returns resulting list
def backSubstitution(extMat):
    guassElimination(extMat)
    # Counts and stores the number of leading zero(s) for each row
    # The number of leading zero(s) corresponds with the first non-zero index
    leadingZeroCntRow = countRowLeadingZero(extMat)
    result = [nan for x in range(len(extMat[0])-1)]
    # No solution (example: 0,0,0|a!=0)
    for i in range(len(leadingZeroCntRow)):
        if leadingZeroCntRow[i] == len(extMat[0]) - 1:
            return result
    # Solution exists
    for i in range(len(extMat)-1, -1, -1):
        if leadingZeroCntRow[i] >= len(extMat[0]):
            continue
        calCache = Rational(0)
        for j in range(len(extMat[0])-2, leadingZeroCntRow[i]-1, -1):
            # Already solved for variable at that column
            if result[j] != nan:
                calCache += (extMat[i][j]*result[j])
                continue
            # Solves for unknown variable(s) at column(s) not leading '1'
            if j != leadingZeroCntRow[i]:
                result[j] = symbols('t'+str(j+1))
                calCache += (extMat[i][j]*result[j])
            # Solves for unknown variable(s) at column(s) with leading '1'
            else:
                result[j] = extMat[i][len(extMat[0])-1] - calCache
    return result


matIn = getMatrixConsole()
print(matIn)
print(backSubstitution(matIn))
