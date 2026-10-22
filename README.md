# beta-binomial-MLE-nov25

# Beta-Binomial MLE: Symbolic Computation with SymPy

This repository contains Python code used to compute resultants for the Beta-Binomial Maximum Likelihood Estimation (MLE) problem. The original code supports the derivation of the closed-form solution presented in Theorem 2 of our paper for the case $n=2$, where each item receives exactly two ratings. The repository also contains scripts for the case $n=3$, which derive the polynomial score equations and compute the corresponding resultant and its factorization.

The paper is:
> *On the MLE for the Beta-Binomial Distribution*  
> Daniel Berend, Yuri Chernyavsky, Luba Sapir

---

## **Overview**
The likelihood equations for the Beta-Binomial model can be transformed into a system of two polynomial equations in two variables (α and β). Using elimination theory and resultants, we reduce this system to a single equation in one variable. For n=2, this leads to the closed-form solution presented in Theorem 2.

The repository contains:
- resultant_computation.py - resultant computations for the case n=2.
- rational_to_polynomial_n_is_3.py - derivation of the polynomial score equations for the case n=3.
- resultant_computation_n_is_3.py - computation and factorization of the resultant for the case n=3.


---

## **Requirements**
- Python 3.8+
- https://www.sympy.org/  
Install with:
```bash
pip install sympy
