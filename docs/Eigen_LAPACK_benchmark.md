LAPACK requires linking to BLAS and libgfortran which can be troublesome to manage. 
LAPACKE is not viable because its routines do not take work arrays but instead dynamically allocate memory, which is too expensive for 2x2 or 3x3 matrices. 
Eigen promises optimized code via lazy evaluation, has statically-sized matrices, is header only and well integrated to CMake. 
We document here the differences observed with LAPACK for future reference, as we intend to make LAPACK optional (and possibly remove it altogether eventually). 


# Eigen decomposition

We begin by testing symmetric matrix eigendecomposition using:
  - The Joachim Kopp implementation of dsyevq, see `src/linalg/dsyevq.(c|h)xx`  and `src/linalg/dsytrd.(c|h)xx`. We originally included this to plug in SANS::SurrealS (automatic differentiation). 
  - LAPACK dsyevq, using the LAPACK libraries installed by homebrew on MacOS. 
  - Eigen with `SelfAdjointEigenSolver`  

The unit test and benchmark is `bunit/test_eigen.cxx`. 
We construct metrics of increased anisotropy ratio defined as the ratio of mesh sizes. 
Recall the eigenvalues are `\lambda_i = 1/h_i^2` so the matrix conditioning is in fact increasing with the square of the anisotropy ratio. 
The largest size is in the order of 1, it is the smallest size (i.e. largest eigenvalue) which is decreased to enforce the ratio. 
We consider ratios from 2 to 2e7 and sample 1e5 metrics per ratio, with a fixed seed. 

