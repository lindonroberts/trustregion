import numpy as np
import trustregion

g = np.array([5.0, 0.0, 4.0])
H = np.array([[1.0, 0.0, 4.0], [0.0, 2.0, 0.0], [4.0, 0.0, 3.0]])
delta = 1.0

s = trustregion.arcunc(delta, g, H)

print("ARCUNC easy case")
print("s =", s)
print("Note: global min is s = (-2.4827526, 0, 1.0418972)")  # lambda=2.6925100

g = np.array([0.0, 2.0, 0.0])
s = trustregion.arcunc(delta, g, H)

print("ARCUNC hard case")
print("s =", s)
print("Note: global min is s = (+/-1.62917555, -0.48507125, -/+1.27203396)")  # lambda=2.1231056

g = np.array([0.0, 2.0, 0.0001])
s = trustregion.arcunc(delta, g, H)

print("ARCUNC nearly hard case")
print("s =", s)
print("Note: global min is s = (1.62920031, -0.48506775, -1.27205329)")  # lambda=2.1231354

