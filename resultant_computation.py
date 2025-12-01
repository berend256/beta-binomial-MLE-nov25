from sympy import symbols, factor, resultant, simplify, pprint

# Declare symbolic variables
alpha, beta = symbols('alpha beta', positive=True)
l0, l1, l2 = symbols('l0 l1 l2', integer=True, positive=True)

# Sum of l_i's
L = l0 + l1 + l2

# Step 0: Show the initial rational equations
print("\n=== Step 0: Starting from Rational Equations ===")
print("For n = 2, the likelihood score equations are:")
print(" (ℓ₁ + ℓ₂)/α + ℓ₂/(α+1) = (ℓ₀ + ℓ₁ + ℓ₂)/(α+β) + (ℓ₀ + ℓ₁ + ℓ₂)/(α+β+1)")
print(" (ℓ₀ + ℓ₁)/β + ℓ₀/(β+1) = (ℓ₀ + ℓ₁ + ℓ₂)/(α+β) + (ℓ₀ + ℓ₁ + ℓ₂)/(α+β+1)")

# Define the simplified polynomial equations (after clearing denominators)
eq1 = - alpha*(alpha + 1)*(2*alpha + 2*beta + 1)*L \
      + l1*(alpha + 1)*(alpha + beta)*(alpha + beta + 1) \
      + l2*(alpha + beta)*(2*alpha + 1)*(alpha + beta + 1)

eq2 = - beta*(beta + 1)*(2*alpha + 2*beta + 1)*L \
      + l0*(alpha + beta)*(2*beta + 1)*(alpha + beta + 1) \
      + l1*(alpha + beta)*(beta + 1)*(alpha + beta + 1)

print("\n=== Step 1: Polynomial Equations ===")
pprint(eq1)
print(" = 0")
pprint(eq2)
print(" = 0")


# Compute resultants
res_beta = resultant(eq1, eq2, beta)
res_alpha = resultant(eq1, eq2, alpha)

# Simplify and factor
res_beta_factored = factor(simplify(res_beta))
res_alpha_factored = factor(simplify(res_alpha))

# Extract essential factor (last parenthesis)
essential_beta = res_beta_factored.as_ordered_factors()[-1]
essential_alpha = res_alpha_factored.as_ordered_factors()[-1]

print("\n=== Step 2: Resultant Eliminating β (full factored) ===")
pprint(res_beta_factored)
print("\nEssential factor for β:")
pprint(essential_beta)

print("\n=== Step 3: Resultant Eliminating α (full factored) ===")
pprint(res_alpha_factored)
print("\nEssential factor for α:")
pprint(essential_alpha)

