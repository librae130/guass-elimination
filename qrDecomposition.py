import numpy
from sympy import Rational, sqrt


def inputMatrixConsole():
  while True:
    try:
      rowCount = input("Enter the number of rows, press 'E' to exit: ")
      if rowCount.upper() == 'E':
        return []
      colCount = input("Enter the number of columns, press 'E' to exit: ")
      if colCount.upper() == 'E':
        return []
      if int(rowCount) > 0 and int(colCount) > 0:
        break
      else:
        print("Please enter positive integers for dimensions.")
    except ValueError:
      print("Please enter valid integers.")
    print("You can enter decimals, fractions (a/b) or integers.")
    print("Enter elements row by row, separated by spaces.\n")
  mat = []
  for i in range(int(rowCount)):
    while True:
      try:
        line = input(f"Row {i+1}, press 'E' to exit: ")
        if line.upper() == 'E':
          return []
        rowIn = line.strip().split()
        if len(rowIn) != int(colCount):
          print(f"Please enter exactly {colCount} elements.")
          continue
        row = []
        for e in rowIn:
          row.append(Rational(e))
        mat.append(row)
        break
      except (ValueError, ZeroDivisionError):
        print(" Please enter valid numbers.")
  return mat


def qrDecomposition(mat):
  if numpy.size(mat) == 0:
    return numpy.array([]), numpy.array([])
  
  row, col = mat.shape
  sqr = min(row, col)

  Q = numpy.zeros((row, sqr), dtype=Rational)

  for i in range(0,sqr):
    Q[:,i] = mat[:,i]
    for j in range(i-1,-1,-1):
      Q[:,i] = Q[:,i] - (mat[:,i].dot(Q[:,j]) / Q[:,j].dot(Q[:,j]) * Q[:,j])

    if(numpy.all(Q[:,i] == 0)):
      raise ValueError("Input matrix has linearly dependent columns")
    
    magnitude = sqrt(numpy.sum(Q[:,i]**2))
    Q[:,i] = Q[:,i] / magnitude

  R = numpy.zeros((sqr,col), dtype=Rational)

  for i in range(0,sqr):
    for j in range(i,col):
         R[i][j] = mat[:, j].dot(Q[:, i])

  return Q, R


matIn = numpy.array(inputMatrixConsole(), dtype=Rational)
try:
  Q, R = qrDecomposition(matIn)
  print("Matrix =\n", matIn)
  print("Q =\n", Q)
  print("R =\n", R)
except ValueError as e:
  print(f"QR decomposition failed: {e}")
