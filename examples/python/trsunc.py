import numpy as np
import trustregion

print("got trustregion package, functions are:")
print(dir(trustregion))
print("testing...")

g = np.array([5.0, 0.0, 4.0])
H = np.array([[1.0, 0.0, 4.0], [0.0, 2.0, 0.0], [4.0, 0.0, 3.0]])
delta = 1.0

s, lda = trustregion.py_trsunc(delta, g, H)

print("TRSUNC easy case")
print("s =", s)
print("lambda = %g" % lda)
print("Note: global min is s = (-1, 0, 0) with lambda=4")

g = np.array([0.0, 2.0, 0.0])
s, lda = trustregion.py_trsunc(delta, g, H)

print("TRSUNC hard case")
print("s =", s)
print("lambda = %g" % lda)
print("Note: global min is s = (+/-0.68926566, -0.48507125, -/+0.53816237) with lambda=2.1231056")

g = np.array([0.0, 2.0, 0.0001])
s, lda = trustregion.py_trsunc(delta, g, H)

print("TRSUNC nearly hard case")
print("s =", s)
print("lambda = %g" % lda)
print("Note: global min is s = (0.6892634, -0.48506297, -0.53817273) with lambda=2.1231760")

