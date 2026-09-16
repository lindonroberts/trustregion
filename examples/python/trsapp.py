import numpy as np
import trustregion

g = np.array([5.0, 0.0, 4.0])
H = np.array([[1.0, 0.0, 4.0], [0.0, 2.0, 0.0], [4.0, 0.0, 3.0]])
delta = 1.0

s, crvmin, info = trustregion.trsapp(delta, g, H, tol=1e-2)

print("TRSAPP calculated step")
print("s =", s)
print("crvmin = %g" % crvmin)
print("info = %g" % info)
print("Note: global min is s = (-1, 0, 0) with lambda=4")

