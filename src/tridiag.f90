! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods for solving tridiagonal systems using the Thomas algorithm (modified
! Gaussian elimination)
!
! Note that this algorithm is only guaranteed to be stable for diagonally
! dominant tridiagonal systems; that is, when the norm of each diagonal element
! is greater than or equal to the sum of the norms of the off diagonal
! elements.
!
! This module provides the following subroutines:
!
! tridiag_gen(...) -> general tridiagonal system
! tridiag_cnst(...) -> all elements constant, band-wise
! tridiag_cnst_diag -> diagonal elements constant
! tridiag_cnst_odiag -> off-diagonal elements constant, band-wise
!
! tridiag_sym(...) -> symmetric tridiagonal system
! tridiag_sym_cnst(...) -> symmetric, all elements constant, band-wise
! tridiag_sym_cnst_diag(...) -> symmetric, diagonal elements constant
! tridiag_sym_cnst_odiag(...) -> symmetric, off-diagonal elements
!   constant, band-wise
!
! Note that the methods for symmetric tridiagonal systems are simply wrappers
! around the methods for the more general asymmetric case. There are currently
! no symmetry-specific numerical optimizations.
module tridiag


  use iso_fortran_env, only: real32, real64, int32

  implicit none
  private


  public :: tridiag_gen
  interface tridiag_gen
     ! General tridiagonal systems
     module procedure tridiag_gen_real_dp
     module procedure tridiag_gen_cmplx_dp
     module procedure tridiag_gen_real_sp
     module procedure tridiag_gen_cmplx_sp
  end interface tridiag_gen


  public :: tridiag_cnst
  interface tridiag_cnst
     ! Tridiagonal systems where all elements are constant, band-wise
     module procedure tridiag_cnst_real_dp
     module procedure tridiag_cnst_cmplx_dp
     module procedure tridiag_cnst_real_sp
     module procedure tridiag_cnst_cmplx_sp
  end interface tridiag_cnst


  public :: tridiag_cnst_diag
  interface tridiag_cnst_diag
     ! Tridiagonal systems where diagonal elements are constant
     module procedure tridiag_cnst_diag_real_dp
     module procedure tridiag_cnst_diag_cmplx_dp
     module procedure tridiag_cnst_diag_real_sp
     module procedure tridiag_cnst_diag_cmplx_sp
  end interface tridiag_cnst_diag

  public :: tridiag_cnst_odiag
  interface tridiag_cnst_odiag
     ! Tridiagonal systems where off-diagonal elements are constant, band-wise
     module procedure tridiag_cnst_odiag_real_dp
     module procedure tridiag_cnst_odiag_cmplx_dp
     module procedure tridiag_cnst_odiag_real_sp
     module procedure tridiag_cnst_odiag_cmplx_sp
  end interface tridiag_cnst_odiag


  public :: tridiag_sym
  interface tridiag_sym
     ! Symmetric tridiagonal systems where off-band elements are equal,
     ! element-wise
     module procedure tridiag_sym_real_dp
     module procedure tridiag_sym_cmplx_dp
     module procedure tridiag_sym_real_sp
     module procedure tridiag_sym_cmplx_sp
  end interface tridiag_sym


  public :: tridiag_sym_cnst
  interface tridiag_sym_cnst
     ! Symmetric tridiagonal systems where all elements are equal, band-wise
     module procedure tridiag_sym_cnst_real_dp
     module procedure tridiag_sym_cnst_cmplx_dp
     module procedure tridiag_sym_cnst_real_sp
     module procedure tridiag_sym_cnst_cmplx_sp
  end interface tridiag_sym_cnst

  public :: tridiag_sym_cnst_diag
  interface tridiag_sym_cnst_diag
     ! Symmetric tridiagonal systems where diagonal elements are equal
     module procedure tridiag_sym_cnst_diag_real_dp
     module procedure tridiag_sym_cnst_diag_cmplx_dp
     module procedure tridiag_sym_cnst_diag_real_sp
     module procedure tridiag_sym_cnst_diag_cmplx_sp
  end interface tridiag_sym_cnst_diag

  public :: tridiag_sym_cnst_odiag
  interface tridiag_sym_cnst_odiag
     ! Symmetric tridiagonal systems where off-diagonal elements are equal,
     ! band-wise
     module procedure tridiag_sym_cnst_odiag_real_dp
     module procedure tridiag_sym_cnst_odiag_cmplx_dp
     module procedure tridiag_sym_cnst_odiag_real_sp
     module procedure tridiag_sym_cnst_odiag_cmplx_sp
  end interface tridiag_sym_cnst_odiag


  interface tridiag_backsweep
     ! Back-substitution step of Thomas algorithm
     module procedure tridiag_backsweep_real_dp
     module procedure tridiag_backsweep_cmplx_dp
     module procedure tridiag_backsweep_real_sp
     module procedure tridiag_backsweep_cmplx_sp
  end interface tridiag_backsweep


  ! Numeric precision parameters
  integer(int32), parameter :: ip = int32
  integer(ip), parameter :: sp = real32
  integer(ip), parameter :: dp = real64


contains


  ! subroutine tridiag_backsweep
  !
  ! Thomas algorithm back-substitution to obtain solution to A x = b from
  ! Gaussian elimination coefficients
  !
  ! mat_coeff(n - 1) :: Gaussian elimination coefficients for A, where A x = b
  ! vec_coeff(n) :: Gaussian elimination coefficients for b, where A x = b
  ! res(n) :: x (unknown), where A x = b

  ! tridiag_backsweep definition for double precision real values
  subroutine tridiag_backsweep_real_dp(mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: mat_coeff(:), vec_coeff(:)
    real(dp), intent(inout) :: res(:)

    include "./tridiag_src/backsweep.src"
  end subroutine tridiag_backsweep_real_dp

  ! tridiag_backsweep definition for double precision complex values
  subroutine tridiag_backsweep_cmplx_dp(mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: mat_coeff(:), vec_coeff(:)
    complex(dp), intent(inout) :: res(:)

    include "./tridiag_src/backsweep.src"
  end subroutine tridiag_backsweep_cmplx_dp

  ! tridiag_backsweep definition for single precision real values
  subroutine tridiag_backsweep_real_sp(mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: mat_coeff(:), vec_coeff(:)
    real(sp), intent(inout) :: res(:)

    include "./tridiag_src/backsweep.src"
  end subroutine tridiag_backsweep_real_sp

  ! tridiag_backsweep definition for single precision complex values
  subroutine tridiag_backsweep_cmplx_sp(mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: mat_coeff(:), vec_coeff(:)
    complex(sp), intent(inout) :: res(:)

    include "./tridiag_src/backsweep.src"
  end subroutine tridiag_backsweep_cmplx_sp


  ! subroutine tridiag_gen
  !
  ! Solve an n-dimensional tridiagional matrix equation A x = b.
  !
  ! diag(n) :: main diagonal entries
  ! u_diag(n - 1) :: upper diagonal entries
  ! l_diag(n - 1) :: lower diagonal entries
  ! vec(n) :: b (known), where A x = b
  ! mat_coeff(n - 1) :: work array
  ! vec_coeff(n) :: work array
  ! res(:) :: x (unknown), where A x = b

  ! tridiag_gen definition for double precision real values
  subroutine tridiag_gen_real_dp(diag, u_diag, l_diag, vec, mat_coeff, &
    vec_coeff, res)

    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec)), u_diag(size(vec) - 1), &
         l_diag(size(vec) - 1)
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/gen.src"
  end subroutine tridiag_gen_real_dp

  ! tridiag_gen definition for double precision complex values
  subroutine tridiag_gen_cmplx_dp(diag, u_diag, l_diag, vec, mat_coeff, &
    vec_coeff, res)

    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec)), u_diag(size(vec) - 1), &
         l_diag(size(vec) - 1)
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/gen.src"
  end subroutine tridiag_gen_cmplx_dp

  ! tridiag_gen definition for single precision real values
  subroutine tridiag_gen_real_sp(diag, u_diag, l_diag, vec, mat_coeff, &
    vec_coeff, res)

    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec)), u_diag(size(vec) - 1), &
         l_diag(size(vec) - 1)
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/gen.src"
  end subroutine tridiag_gen_real_sp

  ! tridiag_gen definition for single precision complex values
  subroutine tridiag_gen_cmplx_sp(diag, u_diag, l_diag, vec, mat_coeff, &
    vec_coeff, res)

    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec)), u_diag(size(vec) - 1), &
         l_diag(size(vec) - 1)
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/gen.src"
  end subroutine tridiag_gen_cmplx_sp


  ! subroutine tridiag_cnst
  !
  ! Solve an n-dimensional tridiagonal matrix equation A x = b, where all band
  ! elements are equal (band-wise)
  !
  ! diag_cnst :: main diagonal entry
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnst :: lower diagonal entry
  ! vec(n) :: b (known), where A x = b
  ! mat_coeff(n - 1) :: work array
  ! vec_coeff(n) :: work array
  ! res(n) :: x (unknown), where A x = b

  ! tridiag_cnst definition for double-precision real values
  subroutine tridiag_cnst_real_dp(diag_cnst, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst.src"
  end subroutine tridiag_cnst_real_dp

  ! tridiag_cnst definition for double-precision complex values
  subroutine tridiag_cnst_cmplx_dp(diag_cnst, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst.src"
  end subroutine tridiag_cnst_cmplx_dp

  ! tridiag_cnst definition for single-precision real values
  subroutine tridiag_cnst_real_sp(diag_cnst, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst.src"
  end subroutine tridiag_cnst_real_sp

  ! tridiag_cnst definition for single-precision complex values
  subroutine tridiag_cnst_cmplx_sp(diag_cnst, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst.src"
  end subroutine tridiag_cnst_cmplx_sp


  ! subroutine tridiag_cnst_diag
  !
  ! Solve an n-dimensional tridiagonal matrix equation A x = b, where the
  ! diagonal elements are equal.
  !
  ! diag_cnst :: diagonal elements
  ! u_diag(n - 1) :: upper diagonal elements
  ! l_diag(n - 1) :: lower diagonal elements
  ! vec(n) :: b (known), where A x = b
  ! mat_coeff(n - 1) :: work array
  ! vec_coeff(n) :: work array
  ! res(n) :: x (unknown), where A x = b

  ! tridiag_cnst_diag definition for double-precision real values
  subroutine tridiag_cnst_diag_real_dp(diag_cnst, u_diag, l_diag, vec, &
       mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag_cnst
    real(dp), intent(in) :: u_diag(size(vec) - 1), l_diag(size(vec) - 1)
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst_diag.src"
  end subroutine tridiag_cnst_diag_real_dp

  ! tridiag_cnst_diag definition for double-precision complex values
  subroutine tridiag_cnst_diag_cmplx_dp(diag_cnst, u_diag, l_diag, vec, &
       mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag_cnst
    complex(dp), intent(in) :: u_diag(size(vec) - 1), l_diag(size(vec) - 1)
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst_diag.src"
  end subroutine tridiag_cnst_diag_cmplx_dp

  ! tridiag_cnst_diag definition for single-precision real values
  subroutine tridiag_cnst_diag_real_sp(diag_cnst, u_diag, l_diag, vec, &
       mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag_cnst
    real(sp), intent(in) :: u_diag(size(vec) - 1), l_diag(size(vec) - 1)
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst_diag.src"
  end subroutine tridiag_cnst_diag_real_sp

  ! tridiag_cnst_diag definition for single-precision complex values
  subroutine tridiag_cnst_diag_cmplx_sp(diag_cnst, u_diag, l_diag, vec, &
       mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag_cnst
    complex(sp), intent(in) :: u_diag(size(vec) - 1), l_diag(size(vec) - 1)
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst_diag.src"
  end subroutine tridiag_cnst_diag_cmplx_sp


  ! subroutine tridiag_cnst_odiag
  !
  ! Solve an n-dimensional tridiagonal matrix equation A x = b, where the
  ! off-diagonal band elements are equal (band-wise)
  !
  ! diag(n) :: diagonal elements
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnst :: lower diagonal entry
  ! vec(n) :: b (known), where A x = b
  ! mat_coeff(n - 1) :: work array
  ! vec_coeff(n) :: work array
  ! res(n) :: x (unknown), where A x = b

  ! tridiag_cnst_odiag definition for double-precision real values
  subroutine tridiag_cnst_odiag_real_dp(diag, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec)), u_diag_cnst, l_diag_cnst
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst_odiag.src"
  end subroutine tridiag_cnst_odiag_real_dp

  ! tridiag_cnst_odiag definition for double-precision complex values
  subroutine tridiag_cnst_odiag_cmplx_dp(diag, u_diag_cnst, &
       l_diag_cnst, vec, mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec)), u_diag_cnst, l_diag_cnst
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))
;
    include "./tridiag_src/cnst_odiag.src"
  end subroutine tridiag_cnst_odiag_cmplx_dp

  ! tridiag_cnst_odiag definition for single-precision real values
  subroutine tridiag_cnst_odiag_real_sp(diag, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec)), u_diag_cnst, l_diag_cnst
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst_odiag.src"
  end subroutine tridiag_cnst_odiag_real_sp

  ! tridiag_cnst_odiag definition for double-precision complex values
  subroutine tridiag_cnst_odiag_cmplx_sp(diag, u_diag_cnst, &
       l_diag_cnst, vec, mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec)), u_diag_cnst, l_diag_cnst
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/cnst_odiag.src"
  end subroutine tridiag_cnst_odiag_cmplx_sp


  ! subroutine tridiag_sym_cnst
  !
  ! Solve an n-dimensional symmetric tridiagonal matrix equation A x = b, where
  ! all elements are equal band-wise.
  !
  ! diag_cnst :: main diagonal entry
  ! o_diag_cnst :: off diagonal entry
  ! vec(n) :: b (known), where A x = b
  ! mat_coeff(n - 1) :: work array
  ! vec_coeff(n) :: work array
  ! res(n) :: x (unknown), where A x = b

  ! tridiag_syn_cnst definition for double-precision real values
  subroutine tridiag_sym_cnst_real_dp(diag_cnst, o_diag_cnst, vec, &
       mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag_cnst, o_diag_cnst
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst.src"
  end subroutine tridiag_sym_cnst_real_dp

  ! tridiag_syn_cnst definition for double-precision complex values
  subroutine tridiag_sym_cnst_cmplx_dp(diag_cnst, o_diag_cnst, vec, &
       mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag_cnst, o_diag_cnst
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst.src"
  end subroutine tridiag_sym_cnst_cmplx_dp

  ! tridiag_syn_cnst definition for single-precision real values
  subroutine tridiag_sym_cnst_real_sp(diag_cnst, o_diag_cnst, vec, &
       mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag_cnst, o_diag_cnst
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst.src"
  end subroutine tridiag_sym_cnst_real_sp

  ! tridiag_sym_cnst definition for single-precision complex values
  subroutine tridiag_sym_cnst_cmplx_sp(diag_cnst, o_diag_cnst, vec, &
       mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag_cnst, o_diag_cnst
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst.src"
  end subroutine tridiag_sym_cnst_cmplx_sp


  ! subroutine tridiag_sym_cnst_diag
  !
  ! Solve an n-dimensional symmetric tridiagonal matrix equation A x = b, where
  ! all diagonal elements are equal
  !
  ! diag_cnst :: diagonal element
  ! o_diag(n - 1) :: off-diagonal entries
  ! vec(n) :: b (known), where A x = b
  ! mat_coeff(n - 1) :: work array
  ! vec_coeff(n) :: work array
  ! res(n) :: x (unknown), where A x = b

  ! tridiag_sym_cnst_diag for double-precision real values
  subroutine tridiag_sym_cnst_diag_real_dp(diag_cnst, o_diag, vec, &
       mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag_cnst, o_diag(size(vec) - 1)
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst_diag.src"
  end subroutine tridiag_sym_cnst_diag_real_dp

  ! tridiag_sym_cnst_diag for double-precision complex values
  subroutine tridiag_sym_cnst_diag_cmplx_dp(diag_cnst, o_diag, vec, &
       mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag_cnst, o_diag(size(vec) - 1)
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst_diag.src"
  end subroutine tridiag_sym_cnst_diag_cmplx_dp

  ! tridiag_sym_cnst_diag for single-precision real values
  subroutine tridiag_sym_cnst_diag_real_sp(diag_cnst, o_diag, vec, &
       mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag_cnst, o_diag(size(vec) - 1)
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst_diag.src"
  end subroutine tridiag_sym_cnst_diag_real_sp

  ! tridiag_sym_cnst_diag for single-precision complex values
  subroutine tridiag_sym_cnst_diag_cmplx_sp(diag_cnst, o_diag, vec, &
       mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag_cnst, o_diag(size(vec) - 1)
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst_diag.src"
  end subroutine tridiag_sym_cnst_diag_cmplx_sp


  ! subroutine tridiag_sym_cnst_odiag
  !
  ! Solve an n-dimensional symmetric tridiagonal matrix equation A x = b, where
  ! all off-diagonal elements are equal, band-wise.
  !
  ! diag(n) :: diagonal entries
  ! o_diag_cnst :: off diagonal entry
  ! vec(n) :: b (known), where A x = b
  ! mat_coeff(n - 1) :: work array
  ! vec_coeff(n) :: work array
  ! res(n) :: x (unknown), where A x = b

  ! tridiag_sym_cnst_odiag for double-precision real values
  subroutine tridiag_sym_cnst_odiag_real_dp(diag, o_diag_cnst, vec, &
       mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec)), o_diag_cnst
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst_odiag.src"
  end subroutine tridiag_sym_cnst_odiag_real_dp

  ! tridiag_sym_cnst_odiag for double-precision complex values
  subroutine tridiag_sym_cnst_odiag_cmplx_dp(diag, o_diag_cnst, vec, &
       mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec)), o_diag_cnst
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst_odiag.src"
  end subroutine tridiag_sym_cnst_odiag_cmplx_dp

  ! tridiag_sym_cnst_odiag for single-precision real values
  subroutine tridiag_sym_cnst_odiag_real_sp(diag, o_diag_cnst, vec, &
       mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec)), o_diag_cnst
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst_odiag.src"
  end subroutine tridiag_sym_cnst_odiag_real_sp

  ! tridiag_sym_cnst_odiag for single-precision complex values
  subroutine tridiag_sym_cnst_odiag_cmplx_sp(diag, o_diag_cnst, vec, &
       mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec)), o_diag_cnst
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym_cnst_odiag.src"
  end subroutine tridiag_sym_cnst_odiag_cmplx_sp


  ! subroutine tridiag_sym
  !
  ! Solve an n-dimensional symmetric tridiagonal matrix equation A x = b
  !
  ! diag(n) :: diagonal entries
  ! o_diag(n) :: off-diagonal entries
  ! vec(n) :: b (known), where A x = b
  ! mat_coeff(n - 1) :: work array
  ! vec_coeff(n) :: work array
  ! res(n) :: x (unknown), where A x = b

  ! tridiag_sym for double-precision real values
  subroutine tridiag_sym_real_dp(diag, o_diag, vec, mat_coeff, &
       vec_coeff, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec)), o_diag(size(vec))
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym.src"
  end subroutine tridiag_sym_real_dp

  ! tridiag_sym for double-precision complex values
  subroutine tridiag_sym_cmplx_dp(diag, o_diag, vec, mat_coeff, &
       vec_coeff, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec)), o_diag(size(vec))
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym.src"
  end subroutine tridiag_sym_cmplx_dp

  ! tridiag_sym for single-precision real values
  subroutine tridiag_sym_real_sp(diag, o_diag, vec, mat_coeff, &
       vec_coeff, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec)), o_diag(size(vec))
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym.src"
  end subroutine tridiag_sym_real_sp

  ! tridiag_sym for single-precision complex values
  subroutine tridiag_sym_cmplx_sp(diag, o_diag, vec, mat_coeff, &
       vec_coeff, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec)), o_diag(size(vec))
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/sym.src"
  end subroutine tridiag_sym_cmplx_sp


end module tridiag
