"""
An implementation of the series solution to the quintic from
https://www.tandfonline.com/doi/epdf/10.1080/00029890.2025.2460966?needAccess=true& p.12
The series appears to converge to a solution in around ~50% of cases in the current configuration,
diverging quite badly in other cases.
"""

import random
import mpmath as mp

# Configuration --------------------------------------------------
MAX_M = 15      # truncate the indexes m2..m5 at this value. Runtime is O(MAX_M^4).

PRECISION = 40    # decimal digits of precision for mpmath


MU_C0 = 0.0
MU_C1 = 2.0  # Scale to c_1 = 2.0 for numerical stability.
MU_C2 = 0.0
MU_C3 = 0.0
MU_C4 = 0.0
MU_C5 = 0.0

SIGMA_C0 = 1.0
SIGMA_C1 = 0.0 # Scale to c_1 = 2.0 for numerical stability.
SIGMA_C2 = 1.0
SIGMA_C3 = 1.0
SIGMA_C4 = 1.0
SIGMA_C5 = 1.0

# ---------------------------------------------------------------

mp.mp.dps = PRECISION  # set working precision


def sample_coeffs():
    C0 = random.gauss(MU_C0, SIGMA_C0)
    C1 = random.gauss(MU_C1, SIGMA_C1)
    C2 = random.gauss(MU_C2, SIGMA_C2)
    C3 = random.gauss(MU_C3, SIGMA_C3)
    C4 = random.gauss(MU_C4, SIGMA_C4)
    C5 = random.gauss(MU_C5, SIGMA_C5)
    return [C0, C1, C2, C3, C4, C5]
    


def compute_series_root(c, max_m=MAX_M):
    """
    Evaluate the truncated formal series for x given coefficients c=(c0..c5)
    summing 0 ≤ m2,m3,m4,m5 ≤ max_m.
    """
    c0, c1, c2, c3, c4, c5 = c
    x = mp.mpf('0')
    
    for m2 in range(max_m + 1):
        for m3 in range(max_m + 1):
            for m4 in range(max_m + 1):
                for m5 in range(max_m + 1):
                    top_fact = mp.factorial(2*m2 + 3*m3 + 4*m4 + 5*m5)
                    bot_fact = mp.factorial(1 + m2 + 2*m3 + 3*m4 + 4*m5)
                    
                    num = (
                        top_fact *
                        c0**(1 + m2 + 2*m3 + 3*m4 + 4*m5) *
                        c2**m2 * c3**m3 * c4**m4 * c5**m5
                    )
                    den = (
                        bot_fact *
                        mp.factorial(m2) * mp.factorial(m3) *
                        mp.factorial(m4) * mp.factorial(m5) *
                        c1**(1 + 2*m2 + 3*m3 + 4*m4 + 5*m5)
                    )
                    x += num / den
    return x


def quintic_value(c, x):
    """Evaluate c0 - c1 x + c2 x^2 + c3 x^3 + c4 x^4 + c5 x^5."""
    c0, c1, c2, c3, c4, c5 = c
    return c0 - c1*x + c2*x**2 + c3*x**3 + c4*x**4 + c5*x**5


# ----------- Run a single random test ---------------------------
coeffs = sample_coeffs()
x_est  = compute_series_root(coeffs)
residual = quintic_value(coeffs, x_est)

print("Random coefficients c0..c5:")
for i, c in enumerate(coeffs):
    print(f"  c{i} = {c}")

print("Solving the quintic equation c0 - c1 x + c2 x^2 + c3 x^3 + c4 x^4 + c5 x^5 = 0")

print(f"\nSeries root x (truncated at m≤{MAX_M}):")
print(f"  x ≈ {x_est}")

print("\nQuintic evaluated at x:")
print(f"  c0 - c1 x + ... + c5 x^5 ≈ {residual}")

#### SAMPLE OUTPUT ####
# Seems like in around ~50% of cases in the current configuration the series converges to a solution
# > python test.py
# Random coefficients c0..c5:
#   c0 = -0.5858071533615916
#   c1 = 2.0
#   c2 = 0.41749501779200293
#   c3 = 0.06392557638909224
#   c4 = -0.7570851793744919
#   c5 = 0.29726408382062586
# Solving the quintic equation c0 - c1 x + c2 x^2 + c3 x^3 + c4 x^4 + c5 x^5 = 0

# Series root x (truncated at m≤15):
#   x ≈ -0.2798338480149204031582971929500635272132

# Quintic evaluated at x:
#   c0 - c1 x + ... + c5 x^5 ≈ -0.000000000009992385063606844158896871234617964026211