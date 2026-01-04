# Sparse Linear Systems Approach to Sudoku Solving

This repository contains an implementation and analysis of the sparse linear systems approach to solving Sudoku puzzles, based on the paper:

> P. Babu, K. Pelckmans, P. Stoica, and J. Li,  
> *“Linear Systems, Sparse Solutions, and Sudoku”*,  
> IEEE Signal Processing Letters, vol. 17, no. 1, pp. 40–42, January 2010.

The project formulates Sudoku as an underdetermined linear system and recovers the solution via sparse optimization using ℓ₁ minimization.

---

## 📌 Project Overview

Sudoku constraints can be expressed as a linear system  
\[
A x = b
\]
where the solution vector \(x\) is inherently sparse under an indicator-variable representation.

The key idea is that the correct Sudoku solution corresponds to the **sparsest feasible solution**, which can be recovered via convex ℓ₁ minimization (Linpro). For difficult puzzles where standard ℓ₁ fails, an **iterative reweighted ℓ₁** strategy is applied.

---

## 📂 Repository Structure

### Core Solver
- `sudoko_solver.m`  
  Constructs the constraint matrix and solves the ℓ₁ minimization problem using CVX (Linpro).

- `iterative_reweighted_L1.m`  
  Implements iterative reweighted ℓ₁ minimization to improve sparsity recovery for hard puzzles.

---

### Puzzle Generation
- `generate9x9Sudoku.m`  
  Generates valid 9×9 Sudoku puzzles with configurable difficulty.

- `generate4x4Sudoku.m`  
  Generates 4×4 Sudoku puzzles for quick testing and debugging.

---

## ▶️ How to Run

### Requirements
- MATLAB
- **CVX** (convex optimization toolbox for MATLAB)

Install CVX from:  
https://cvxr.com/cvx/

After installation, run:
```matlab
cvx_setup
