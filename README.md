# Some "ForTran" for Saad's Book

Some "Formula Translations" for professor Yousef Saad's book "[Iterative Methods for Sparse Linear Systems (2nd Edition)](http://www-users.cs.umn.edu/~saad/IterMethBook\_2ndEd.pdf)". All codes are written in MATLAB and followed the algorithms given in Saad's book or the listed references. I also referred to Prof. Saad's reference: [SPARSKIT v2](http://www-users.cs.umn.edu/~saad/software/SPARSKIT/). The codes are for my learning and practice, only checked for a few test cases.



## Index for the Algorithms

Chapter 6 & 7: Krylov Subspace Methods, in [this repo](https://github.com/EnigmaHuang/Krylov_Subspace_Methods)

*   (Restarted) Full Orthogonalization Method (FOM);
*   (Restarted) GMRES (with Arnold / Householder orthonormalization);
*   Conjugate Gradient (CG);
*   Conjugate Residual (CR);
*   Biconjugate Gradient without/with Stabilized (BiCG / BiCGStab).



Chapter 9: Preconditioned Iterations, in [this repo](https://github.com/EnigmaHuang/Krylov_Subspace_Methods)

*   Left preconditioned restarted GMRES;
*   Left preconditioned Conjugate Gradient.



Chapter 10: Preconditioning Techniques (in directory `Chapter10`)

*   (Modified) Incomplete LU factorization with zero fill-in (MILU(0) / ILU(0)), using Compressed Sparse Row (CSR) format matrix as input and output;
*   Incomplete LU factorization with dual truncation mechanism (ILUT), using CSR format matrix as input and Modified Sparse Row (MSR) format as output;



Chapter 11: Multigrid Methods

*   Geometric Multigrid with full-weighted restriction operator and Galerkin coarse grid, in [this repo](https://github.com/EnigmaHuang/Poisson_FDM_Multigrid);
*   Classic Algebraic Multigrid, in [this repo](https://github.com/EnigmaHuang/Classic_AMG_Demo).



