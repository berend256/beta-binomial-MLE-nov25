# beta-binomial-MLE-nov25

# Beta-Binomial MLE: Symbolic Computation with SymPy

This repository contains Python code used to compute **resultants** for the Beta-Binomial Maximum Likelihood Estimation (MLE) problem in the special case where each item receives exactly two ratings (n = 2). The code supports the derivation of the closed-form solution presented in **Theorem 5** of our paper:

> *On the MLE for the Beta-Binomial Distribution*  
> Daniel Berend, Yuri Chernyavsky, Luba Sapir  
> [Add DOI or arXiv link here]

---

## **Overview**
The likelihood equations for the Beta-Binomial model can be transformed into a system of two polynomial equations in two variables (α and β). Using **elimination theory** and **resultants**, we reduce this system to a single equation in one variable, enabling a closed-form solution.

This script:
- Defines the polynomial system for n = 2.
- Computes the **resultant** with respect to α and β using SymPy.
- Simplifies and factors the resultants for clarity.

---

## **Requirements**
- Python 3.8+
- https://www.sympy.org/  
Install with:
```bash
pip install sympy
