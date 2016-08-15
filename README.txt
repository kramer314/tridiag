Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
See the LICENSE.txt file at the top-level of this distribution.

Description
===========
Fortran library to help handle linear algebra problems involving positive
definite tridiagonal matrices. There are two separate implementations:

[TODO]
* LAPACK interfaces for various classes of positive definite tridiagonal
  matrices. This also includes a set of interfaces to LAPACK's
  eigendecompostion routines for symmetric positive definite tridiagonal
  matrices.

[IMPLEMENTED]
* A "reference" implementation of the Thomas algorithm (modified Gaussian
  elimination) for various classes of positive definite tridiagonal matrices.
  Especially for tridiagonal matrices with symmetry and / or constant elements,
  some of these routines may be preferable to the LAPACK interfaces.

At some point, a ScaLAPACK interface may also be added.

There is also a module which performs optimized matrix-vector multiplication
for tridiagonal matrices.

Required dependencies
=====================
* Fortran compiler (most recently tested on GNU Fortran 6.1)
* (Optional) LAPACK (most recently tested on LAPACK 3.6.0)

Basic setup
===========
The following will download and build all the tridiagonal libraries:

    https://github.com/kramer314/tridiag.git
    cd ./tridiag/
    scons
