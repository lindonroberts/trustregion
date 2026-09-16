import numpy as np
import trustregion

g = np.array([5.0, 0.0, 4.0])
H = np.array([[1.0, 0.0, 4.0], [0.0, 2.0, 0.0], [4.0, 0.0, 3.0]])
delta = 1.0

xbase = np.array([0.0, 0.0, 0.0])
sl = np.array([-0.1, -1.0, -1.0])
su = np.array([1.0, 1.0, 1.0])

s, crvmin = trustregion.trsbox(delta, g, H, sl, su, xbase, tol=1e-2)

print("TRSBOX calculated step")
print("s =", s)
print("crvmin = %g" % crvmin)
print("Note: global min (no box) is s = (-1, 0, 0)")

