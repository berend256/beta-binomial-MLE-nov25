import sympy as sp

# ======================================================================
# VARIABLES
# ======================================================================

alpha, beta = sp.symbols('alpha beta')
l0, l1, l2, l3 = sp.symbols('l0 l1 l2 l3')

k = l0 + l1 + l2 + l3

# ======================================================================
# SCORE EQUATIONS (n = 3)
# ======================================================================

eq1 = (
    (l1 + l2 + l3)/alpha
    + (l2 + l3)/(alpha + 1)
    + l3/(alpha + 2)
    - k * (
        1/(alpha + beta)
        + 1/(alpha + beta + 1)
        + 1/(alpha + beta + 2)
    )
)

eq2 = (
    (l0 + l1 + l2)/beta
    + (l0 + l1)/(beta + 1)
    + l0/(beta + 2)
    - k * (
        1/(alpha + beta)
        + 1/(alpha + beta + 1)
        + 1/(alpha + beta + 2)
    )
)

# ======================================================================
# CONVERT TO POLYNOMIALS
# ======================================================================

P1 = sp.expand(sp.together(eq1).as_numer_denom()[0])
P2 = sp.expand(sp.together(eq2).as_numer_denom()[0])

print("=" * 70)
print("P1")
print("=" * 70)
print(P1)

print()
print("=" * 70)
print("P2")
print("=" * 70)
print(P2)

# ======================================================================
# SANITY CHECKS
# ======================================================================

print("\n" + "=" * 70)
print("SANITY CHECKS")
print("=" * 70)

tests = [
    (
        "(1,1,1,1), alpha=1, beta=1",
        {l0: 1, l1: 1, l2: 1, l3: 1, alpha: 1, beta: 1}
    ),
    (
        "(4,3,2,1), alpha=1, beta=2",
        {l0: 4, l1: 3, l2: 2, l3: 1, alpha: 1, beta: 2}
    ),
    (
        "(1,2,3,4), alpha=2, beta=1",
        {l0: 1, l1: 2, l2: 3, l3: 4, alpha: 2, beta: 1}
    )
]

for name, subs in tests:
    print(f"\n{name}")
    print("P1 =", sp.expand(P1.subs(subs)))
    print("P2 =", sp.expand(P2.subs(subs)))

# ======================================================================
# SYMMETRY TEST
# ======================================================================

print("\n" + "=" * 70)
print("SYMMETRY TEST")
print("=" * 70)

a, b = sp.symbols('a b')
m0, m1, m2, m3 = sp.symbols('m0 m1 m2 m3')

P2_swap = (
    P2.subs({
        alpha: a,
        beta: b,
        l0: m0,
        l1: m1,
        l2: m2,
        l3: m3
    })
    .subs({
        a: beta,
        b: alpha,
        m0: l3,
        m1: l2,
        m2: l1,
        m3: l0
    })
)

print("\nP1 - swapped(P2):")
print(sp.expand(P1 - P2_swap))

# ======================================================================
# RESULTANTS
# ======================================================================

print("\n" + "=" * 70)
print("RESULTANTS")
print("=" * 70)

print("\nComputing resultant eliminating beta...")
R_alpha = sp.expand(sp.resultant(P1, P2, beta))

print("\nComputing resultant eliminating alpha...")
R_beta = sp.expand(sp.resultant(P1, P2, alpha))

print("\nDegree of R_alpha in alpha:",
      sp.Poly(R_alpha, alpha).degree())

print("Degree of R_beta in beta:",
      sp.Poly(R_beta, beta).degree())

print("\n" + "=" * 70)
print("FACTORIZATION OF R_alpha")
print("=" * 70)
print(sp.factor(R_alpha))

print("\n" + "=" * 70)
print("FACTORIZATION OF R_beta")
print("=" * 70)
print(sp.factor(R_beta))

# ======================================================================
# EXTRACT THE ESSENTIAL ELIMINATION FACTORS
# ======================================================================
#
# The resultant contains many obvious factors arising from denominator
# clearing in the original likelihood equations.  For n = 3 we expect
# factors such as
#
#     alpha^4 (alpha+1)^3 (alpha+2)^2
#
# and the analogous beta factors.
#
# The mathematically interesting object is the remaining factor.
# We call these
#
#     F_alpha(alpha)
#     F_beta(beta)
#
# and regard them as the genuine elimination polynomials.
#

coeff_alpha, factors_alpha = sp.factor_list(R_alpha)
coeff_beta, factors_beta = sp.factor_list(R_beta)

F_alpha = None
F_beta = None

for f, e in factors_alpha:
    if sp.Poly(f, alpha).degree() > 0:
        deg = sp.Poly(f, alpha).degree()
        if deg == 9:
            F_alpha = sp.expand(f)

# =====================================================================
# EXPORT F_alpha
# =====================================================================
#
# The factor F_alpha is the essential elimination polynomial.
# We store it in a text file so that subsequent experimental programs
# can use it without recomputing the resultant.
#

with open("F_alpha_n_is_3.txt", "w") as f:
    f.write(str(F_alpha))

print("\nSaved F_alpha to F_alpha_n_is_3.txt")

for f, e in factors_beta:
    if sp.Poly(f, beta).degree() > 0:
        deg = sp.Poly(f, beta).degree()
        if deg == 9:
            F_beta = sp.expand(f)

print("\n" + "=" * 70)
print("ESSENTIAL ELIMINATION POLYNOMIALS")
print("=" * 70)

print("\nDegree of F_alpha in alpha:",
      sp.Poly(F_alpha, alpha).degree())

print("Degree of F_beta in beta:",
      sp.Poly(F_beta, beta).degree())

# ======================================================================
# TEST 1: MONOMIAL COUNTS
# ======================================================================
#
# The degree of the elimination polynomial gives only a rough measure
# of complexity.  A second useful statistic is the number of monomials.
#
# This tells us how complicated the symbolic expression really is.
#

print("\n" + "=" * 70)
print("TEST 1: MONOMIAL COUNTS")
print("=" * 70)

num_monomials_alpha = len(sp.Poly(F_alpha, alpha).terms())
num_monomials_beta  = len(sp.Poly(F_beta, beta).terms())

print("Number of monomials in F_alpha:", num_monomials_alpha)
print("Number of monomials in F_beta :", num_monomials_beta)

# ======================================================================
# TEST 2: ELIMINATION-POLYNOMIAL SANITY CHECKS
# ======================================================================
#
# The following datasets are derived from exact beta-binomial models:
#
#    (l0,l1,l2,l3) = (1,1,1,1)  -> (alpha,beta) = (1,1)
#    (l0,l1,l2,l3) = (4,3,2,1)  -> (alpha,beta) = (1,2)
#    (l0,l1,l2,l3) = (1,2,3,4)  -> (alpha,beta) = (2,1)
#
# Since these are genuine solutions of the score equations, the
# essential elimination polynomials F_alpha and F_beta must vanish
# at the corresponding alpha and beta values.
#

print("\n" + "=" * 70)
print("TEST 2: EXACT ROOT CHECKS")
print("=" * 70)

known_cases = [
    (
        "(1,1,1,1)",
        {l0:1, l1:1, l2:1, l3:1},
        1,
        1
    ),
    (
        "(4,3,2,1)",
        {l0:4, l1:3, l2:2, l3:1},
        1,
        2
    ),
    (
        "(1,2,3,4)",
        {l0:1, l1:2, l2:3, l3:4},
        2,
        1
    )
]

for name, counts, alpha_root, beta_root in known_cases:

    Fa_val = sp.expand(
        F_alpha.subs(counts).subs(alpha, alpha_root)
    )

    Fb_val = sp.expand(
        F_beta.subs(counts).subs(beta, beta_root)
    )

    print("\nDataset", name)
    print("F_alpha(alpha =", alpha_root, ") =", Fa_val)
    print("F_beta(beta  =", beta_root, ") =", Fb_val)

# ======================================================================
# TEST 3: FACTORIZATION AFTER SPECIALIZATION
# ======================================================================
#
# We now substitute the exact datasets into the elimination
# polynomials themselves and factor the resulting one-variable
# polynomials.
#
# This reveals:
#
#   1. whether the known exact root appears explicitly;
#   2. whether it appears with multiplicity greater than one;
#   3. how the generic degree-9 elimination polynomial behaves in
#      special situations.
#
# In particular, we expect to see factors such as (alpha-1)
# and (beta-2) in the appropriate examples.
#

print("\n" + "=" * 70)
print("TEST 3: SPECIALIZED FACTORIZATIONS")
print("=" * 70)

for name, counts, alpha_root, beta_root in known_cases:

    print("\nDataset", name)

    Fa_spec = sp.factor(F_alpha.subs(counts))
    Fb_spec = sp.factor(F_beta.subs(counts))

    print("\nF_alpha specialized:")
    print(Fa_spec)

    print("\nF_beta specialized:")
    print(Fb_spec)

# ======================================================================
# TEST 4: SYMMETRY OF THE ESSENTIAL FACTORS
# ======================================================================
#
# The original score equations satisfy the symmetry
#
#     alpha <-> beta
#     l0    <-> l3
#     l1    <-> l2
#
# The essential elimination polynomials should inherit the same
# symmetry.  This test ignores the denominator-generated factors and
# checks only the mathematically meaningful factors F_alpha and F_beta.
#

print("\n" + "=" * 70)
print("TEST 4: SYMMETRY OF ESSENTIAL FACTORS")
print("=" * 70)

a, b = sp.symbols('a b')
m0, m1, m2, m3 = sp.symbols('m0 m1 m2 m3')

F_beta_swap = (
    F_beta.subs({
        beta: b,
        l0: m0,
        l1: m1,
        l2: m2,
        l3: m3
    })
    .subs({
        b: alpha,
        m0: l3,
        m1: l2,
        m2: l1,
        m3: l0
    })
)

print("\nDifference:")
print(sp.expand(F_alpha - F_beta_swap))

print("\nSum:")
print(sp.expand(F_alpha + F_beta_swap))

from sympy import Poly

print("Monomials in F_alpha:",
      len(Poly(F_alpha, alpha, l0, l1, l2, l3).terms()))
