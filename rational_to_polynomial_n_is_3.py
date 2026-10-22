import sympy as sp

# Variables
alpha, beta = sp.symbols('alpha beta')
l0, l1, l2, l3 = sp.symbols('l0 l1 l2 l3')

k = l0 + l1 + l2 + l3

# Score equation 1
eq1 = (
    (l1 + l2 + l3)/alpha
    + (l2 + l3)/(alpha + 1)
    + l3/(alpha + 2)
    - k*(
        1/(alpha + beta)
        + 1/(alpha + beta + 1)
        + 1/(alpha + beta + 2)
    )
)

# Score equation 2
eq2 = (
    (l0 + l1 + l2)/beta
    + (l0 + l1)/(beta + 1)
    + l0/(beta + 2)
    - k*(
        1/(alpha + beta)
        + 1/(alpha + beta + 1)
        + 1/(alpha + beta + 2)
    )
)

# Polynomial numerators
P1 = sp.expand(sp.together(eq1).as_numer_denom()[0])
P2 = sp.expand(sp.together(eq2).as_numer_denom()[0])

print("="*70)
print("P1")
print("="*70)
print(P1)

print()
print("="*70)
print("P2")
print("="*70)
print(P2)

# ======================================================================
# SANITY CHECKS
# ======================================================================

print("\n" + "="*70)
print("SANITY CHECKS")
print("="*70)

tests = [
    (
        "(1,1,1,1), alpha=1, beta=1",
        {l0:1, l1:1, l2:1, l3:1, alpha:1, beta:1}
    ),
    (
        "(4,3,2,1), alpha=1, beta=2",
        {l0:4, l1:3, l2:2, l3:1, alpha:1, beta:2}
    ),
    (
        "(1,2,3,4), alpha=2, beta=1",
        {l0:1, l1:2, l2:3, l3:4, alpha:2, beta:1}
    )
]

for name, subs in tests:
    print(f"\n{name}")
    print("P1 =", sp.expand(P1.subs(subs)))
    print("P2 =", sp.expand(P2.subs(subs)))

# ======================================================================
# SYMMETRY TEST
# ======================================================================

print("\n" + "="*70)
print("SYMMETRY TEST")
print("="*70)

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

print("\nP1 + swapped(P2):")
print(sp.expand(P1 + P2_swap))
