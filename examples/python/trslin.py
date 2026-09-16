import numpy as np
import trustregion

g = np.array([1.0, 0.0, 1.0])
H = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
delta = 2.0

n = len(g)
m = 2 * n

xbase = np.array([1.0, 1.0, 1.0])
sl = xbase + np.array([-0.5, -10.0, -10.0])
su = xbase + np.array([11.0, 11.0, 11.0])

# Fill bounds as linear constraints
A = np.zeros((m, n))
b = np.zeros((m,))

for i in range(n):
    A[i, i] = 1.0
    b[i] = su[i]
    A[n+i, i] = -1.0
    b[n+i] = -sl[i]

s = trustregion.trslin(delta, g, H, A, b, xbase, tol=1e-2)

print("TRSLIN calculated original step")
print("s =", s)

# Replace final constraint xbase[3] + s[3] >= -10 (unused) with a constraint forcing the 
# bound constraint solution s=[-0.5, 0, -0.5] to not be valid
# Use -e * (xopt + s) <= -3, so since e * xopt = 3 this means s[1] + s[2] + s[3] >= 0
A[-1, :] = -1.0
b[-1] = -3.0

s = trustregion.trslin(delta, g, H, A, b, xbase, tol=1e-2)

print("TRSLIN calculated modified step")
print("s =", s)
