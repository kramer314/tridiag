! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods for multiplying tridiagonal matrices by a vector; these are in
! general much more efficient than standard matrix multiplication.
!
! This module provides the following subroutines:
!
! tridiag_matmul_gen(...) -> general tridiagonal matrix
! tridiag_matmul_cnst(...) -> all elements cosntant, band-wise
! tridiag_matmul_cnst_diag(...) -> diagonal elements constant
! tridiag_matmul_cnst_odiag(...) -> off-diagonal elements constant, band-wise
!
! tridiag_sym(...) -> symmetric tridiagonal matrix
! tridiag_sym_cnst(...) -> symmetric, all elements constant, band-wise
! tridiag_sym_cnst_diag(...) -> symmetric, diagonal elements constant
! tridiag_sym_cnst_odiag(...) -> symmetric, off-diagonal elements constant
!   band-wise
!
! Note that the methods for symmetric tridiagonal matrices are simply wrappers
! around the methods for the more general asymmetric case. There are currently
! no symmetry-specific numerical optimizations.
module tridiag_matmul


  use iso_fortran_env, only: real32, real64, int32

  implicit none
  private


  public :: tridiag_matmul_gen
  interface tridiag_matmul_gen
     ! Multiply general tridiagonal matrix by vector
     module procedure tridiag_matmul_gen_real_dp
     module procedure tridiag_matmul_gen_cmplx_dp
     module procedure tridiag_matmul_gen_real_sp
     module procedure tridiag_matmul_gen_cmplx_sp
  end interface tridiag_matmul_gen


  public :: tridiag_matmul_cnst
  interface tridiag_matmul_cnst
     ! Multiply tridiagonal matrix by vector where all elements are equal,
     ! band-wise
     module procedure tridiag_matmul_cnst_real_dp
     module procedure tridiag_matmul_cnst_cmplx_dp
     module procedure tridiag_matmul_cnst_real_sp
     module procedure tridiag_matmul_cnst_cmplx_sp
  end interface tridiag_matmul_cnst


  public :: tridiag_matmul_cnst_diag
  interface tridiag_matmul_cnst_diag
     ! Multiply tridiagonal matrix by vector where diagonal elements are equal
     module procedure tridiag_matmul_cnst_diag_real_dp
     module procedure tridiag_matmul_cnst_diag_cmplx_dp
     module procedure tridiag_matmul_cnst_diag_real_sp
     module procedure tridiag_matmul_cnst_diag_cmplx_sp
  end interface tridiag_matmul_cnst_diag


  public :: tridiag_matmul_cnst_odiag
  interface tridiag_matmul_cnst_odiag
     ! Multiply tridiagonal matrix by vector where off-diagonal elements are
     ! equal, band-wise
     module procedure tridiag_matmul_cnst_odiag_real_dp
     module procedure tridiag_matmul_cnst_odiag_cmplx_dp
     module procedure tridiag_matmul_cnst_odiag_real_sp
     module procedure tridiag_matmul_cnst_odiag_cmplx_sp
  end interface tridiag_matmul_cnst_odiag


  public :: tridiag_matmul_sym
  interface tridiag_matmul_sym
     ! Multiply symmetric tridiagonal matrix by vector
     module procedure tridiag_matmul_sym_real_dp
     module procedure tridiag_matmul_sym_cmplx_dp
     module procedure tridiag_matmul_sym_real_sp
     module procedure tridiag_matmul_sym_cmplx_sp
  end interface tridiag_matmul_sym


  public :: tridiag_matmul_sym_cnst
  interface tridiag_matmul_sym_cnst
     ! Multiply symmetric tridiagonal matrix by vector, where all elements are
     ! equal, band-wise
     module procedure tridiag_matmul_sym_cnst_real_dp
     module procedure tridiag_matmul_sym_cnst_cmplx_dp
     module procedure tridiag_matmul_sym_cnst_real_sp
     module procedure tridiag_matmul_sym_cnst_cmplx_sp
  end interface tridiag_matmul_sym_cnst


  public :: tridiag_matmul_sym_cnst_diag
  interface tridiag_matmul_sym_cnst_diag
     ! Multiply symmetric tridiagonal matrix by vector, where diagonal elements
     ! are constant
     module procedure tridiag_matmul_sym_cnst_diag_real_dp
     module procedure tridiag_matmul_sym_cnst_diag_cmplx_dp
     module procedure tridiag_matmul_sym_cnst_diag_real_sp
     module procedure tridiag_matmul_sym_cnst_diag_cmplx_sp
  end interface tridiag_matmul_sym_cnst_diag


  public :: tridiag_matmul_sym_cnst_odiag
  interface tridiag_matmul_sym_cnst_odiag
     ! Multiply symmetric tridiagonal matrix by vector, where off-diagonal
     ! elements are constant, band-wise
     module procedure tridiag_matmul_sym_cnst_odiag_real_dp
     module procedure tridiag_matmul_sym_cnst_odiag_cmplx_dp
     module procedure tridiag_matmul_sym_cnst_odiag_real_sp
     module procedure tridiag_matmul_sym_cnst_odiag_cmplx_sp
  end interface tridiag_matmul_sym_cnst_odiag


  ! Numeric precision parameters
  integer(int32), parameter :: ip = int32
  integer(ip), parameter :: sp = real32
  integer(ip), parameter :: dp = real64


contains


  ! subroutine tridiag_matmul_gen
  !
  ! Multiply general tridiagonal matrix by vector
  !
  ! diag(n) :: diagonal entries
  ! u_diag(n - 1) :: upper diagonal entries
  ! l_diag(n - 1) :: upper diagonal entries
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_gen definition for double-precision real values
  subroutine tridiag_matmul_gen_real_dp(diag, u_diag, l_diag, vec, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec))
    real(dp), intent(in) :: u_diag(size(diag) - 1)
    real(dp), intent(in) :: l_diag(size(diag) - 1)
    real(dp), intent(inout) :: res(size(diag))

    include "./tridiag_matmul_src/gen.src"
  end subroutine tridiag_matmul_gen_real_dp

  ! tridiag_matmul_gen definition for double-precision complex values
  subroutine tridiag_matmul_gen_cmplx_dp(diag, u_diag, l_diag, vec, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec))
    complex(dp), intent(in) :: u_diag(size(diag) - 1)
    complex(dp), intent(in) :: l_diag(size(diag) - 1)
    complex(dp), intent(inout) :: res(size(diag))

    include "./tridiag_matmul_src/gen.src"
  end subroutine tridiag_matmul_gen_cmplx_dp

  ! tridiag_matmul_gen definition for single-precision real values
  subroutine tridiag_matmul_gen_real_sp(diag, u_diag, l_diag, vec, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec))
    real(sp), intent(in) :: u_diag(size(diag) - 1)
    real(sp), intent(in) :: l_diag(size(diag) - 1)
    real(sp), intent(inout) :: res(size(diag))

    include "./tridiag_matmul_src/gen.src"
  end subroutine tridiag_matmul_gen_real_sp

  ! tridiag_matmul_gen definition for single-precision complex values
  subroutine tridiag_matmul_gen_cmplx_sp(diag, u_diag, l_diag, vec, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec))
    complex(sp), intent(in) :: u_diag(size(diag) - 1)
    complex(sp), intent(in) :: l_diag(size(diag) - 1)
    complex(sp), intent(inout) :: res(size(diag))

    include "./tridiag_matmul_src/gen.src"
  end subroutine tridiag_matmul_gen_cmplx_sp


  ! subroutine tridiag_matmul_cnst
  !
  ! Multiply tridiagonal matrix by a vector, where the band elements are
  ! equal (band-wise)
  !
  ! diag_cnst :: diagonal entry
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnnst :: lower diagonal entry
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_cnst definition for double-precision real values
  subroutine tridiag_matmul_cnst_real_dp(diag_cnst, u_diag_cnst, &
       l_diag_cnst, vec, res)
    real(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst
    real(dp), intent(in) :: vec(:)
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst.src"
  end subroutine tridiag_matmul_cnst_real_dp

  ! tridiag_matmul_cnst definition for double-precision complex values
  subroutine tridiag_matmul_cnst_cmplx_dp(diag_cnst, u_diag_cnst, &
       l_diag_cnst, vec, res)
    complex(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst.src"
  end subroutine tridiag_matmul_cnst_cmplx_dp

  ! tridiag_mmatmul_cnst definition for single-precision real values
  subroutine tridiag_matmul_cnst_real_sp(diag_cnst, u_diag_cnst, &
       l_diag_cnst, vec, res)
    real(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst
    real(sp), intent(in) :: vec(:)
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst.src"
  end subroutine tridiag_matmul_cnst_real_sp

  ! tridiag_matmul_cnst defintion for single-precision complex values
  subroutine tridiag_matmul_cnst_cmplx_sp(diag_cnst, u_diag_cnst, &
       l_diag_cnst, vec, res)
    complex(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst.src"
  end subroutine tridiag_matmul_cnst_cmplx_sp


  ! subroutine tridiag_matmul_cnst_diag
  !
  ! Multiply tridiagonal matrix by a vector, where the diagonal elements are
  ! equal
  !
  ! diag_cnst :: diagonal element
  ! u_diag(n - 1) :: upper diagonal elements
  ! l_diag(n - 1) :: lower diagonal elements
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_cnst_diag definition for double-precision real values
  subroutine tridiag_matmul_cnst_diag_real_dp(diag_cnst, u_diag, l_diag, &
       vec, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag_cnst
    real(dp), intent(in) :: u_diag(size(vec) - 1), l_diag(size(vec) - 1)
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst_diag.src"
  end subroutine tridiag_matmul_cnst_diag_real_dp

  ! tridiag_matmul_cnst_diag definition for double-precision complex values
  subroutine tridiag_matmul_cnst_diag_cmplx_dp(diag_cnst, u_diag, l_diag, &
       vec, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag_cnst
    complex(dp), intent(in) :: u_diag(size(vec) - 1), l_diag(size(vec) - 1)
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst_diag.src"
  end subroutine tridiag_matmul_cnst_diag_cmplx_dp

  ! tridiag_matmul_cnst_diag definition for single-precision real values
  subroutine tridiag_matmul_cnst_diag_real_sp(diag_cnst, u_diag, l_diag, &
       vec, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag_cnst
    real(sp), intent(in) :: u_diag(size(vec) - 1), l_diag(size(vec) - 1)
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst_diag.src"
  end subroutine tridiag_matmul_cnst_diag_real_sp

  ! tridiag_matmul_cnst_diag definition for double-precision complex values
  subroutine tridiag_matmul_cnst_diag_cmplx_sp(diag_cnst, u_diag, l_diag, &
       vec, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag_cnst
    complex(sp), intent(in) :: u_diag(size(vec) - 1), l_diag(size(vec) - 1)
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst_diag.src"
  end subroutine tridiag_matmul_cnst_diag_cmplx_sp


  ! subroutine tridiag_matmul_cnst_offdiag
  !
  ! Multiply tridiagonal matrix by a vector, where the off-diagonal band
  ! elements are equal (band-wise)
  !
  ! diag(n) :: diagonal elements
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnst :: lower diagonal entry
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_cnst_offdiag definition for double-precision real values
  subroutine tridiag_matmul_cnst_odiag_real_dp(diag, u_diag_cnst, &
       l_diag_cnst, vec, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec))
    real(dp), intent(in) :: u_diag_cnst, l_diag_cnst
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst_odiag.src"
  end subroutine tridiag_matmul_cnst_odiag_real_dp

  ! tridiag_matmul_cnst_offdiag definition for double-precision complex
  ! values
  subroutine tridiag_matmul_cnst_odiag_cmplx_dp(diag, u_diag_cnst, &
       l_diag_cnst, vec, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec))
    complex(dp), intent(in) :: u_diag_cnst, l_diag_cnst
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst_odiag.src"
  end subroutine tridiag_matmul_cnst_odiag_cmplx_dp

  ! tridiag_matmul_cnst_offdiag definition for single-precision real values
  subroutine tridiag_matmul_cnst_odiag_real_sp(diag, u_diag_cnst, &
       l_diag_cnst, vec, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec))
    real(sp), intent(in) :: u_diag_cnst, l_diag_cnst
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst_odiag.src"
  end subroutine tridiag_matmul_cnst_odiag_real_sp

  ! tridiag_matmul_cnst_offdiag definition for single-precision complex
  ! values
  subroutine tridiag_matmul_cnst_odiag_cmplx_sp(diag, u_diag_cnst, &
       l_diag_cnst, vec, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec))
    complex(sp), intent(in) :: u_diag_cnst, l_diag_cnst
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/cnst_odiag.src"
  end subroutine tridiag_matmul_cnst_odiag_cmplx_sp


  ! subroutine tridiag_matmul_sym
  !
  ! Multiply symmetric tridiagonal matrix by vector
  !
  ! diag(n) :: diagonal entries
  ! o_diag(n - 1) :: off-diagonal entries
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_sym definition for double-precision real valeus
  subroutine tridiag_matmul_sym_real_dp(diag, o_diag, vec, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec)), o_diag(size(vec) - 1)
    real(dp), intent(inout) :: res(size(diag))

    include "./tridiag_matmul_src/sym.src"
  end subroutine tridiag_matmul_sym_real_dp

  ! tridiag_matmul_sym definition for double-precision complex valeus
  subroutine tridiag_matmul_sym_cmplx_dp(diag, o_diag, vec, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec)), o_diag(size(vec) - 1)
    complex(dp), intent(inout) :: res(size(diag))

    include "./tridiag_matmul_src/sym.src"
  end subroutine tridiag_matmul_sym_cmplx_dp

  ! tridiag_matmul_sym definition for single-precision real valeus
  subroutine tridiag_matmul_sym_real_sp(diag, o_diag, vec, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec)), o_diag(size(vec) - 1)
    real(sp), intent(inout) :: res(size(diag))

    include "./tridiag_matmul_src/sym.src"
  end subroutine tridiag_matmul_sym_real_sp

  ! tridiag_matmul_sym definition for double-precision complex valeus
  subroutine tridiag_matmul_sym_cmplx_sp(diag, o_diag, vec, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec)), o_diag(size(vec) - 1)
    complex(sp), intent(inout) :: res(size(diag))

    include "./tridiag_matmul_src/sym.src"
  end subroutine tridiag_matmul_sym_cmplx_sp


  ! subroutine tridiag_matmul_syn_cnst
  !
  ! Multiply symmetric tridiagonal matrix by vector, where the band elements
  ! equal, band-wise
  !
  ! diag_cnst :: diagonal entry
  ! o_diag_cnst :: off-diagonal entry
  ! vec(n) :: vectory to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_sym_cnst definition for double-precision real values
  subroutine tridiag_matmul_sym_cnst_real_dp(diag_cnst, o_diag_cnst, vec, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag_cnst, o_diag_cnst
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst.src"
  end subroutine tridiag_matmul_sym_cnst_real_dp

  ! tridiagonal_matmul_sym_cnst definition for double-precision complex values
  subroutine tridiag_matmul_sym_cnst_cmplx_dp(diag_cnst, o_diag_cnst, vec, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag_cnst, o_diag_cnst
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst.src"
  end subroutine tridiag_matmul_sym_cnst_cmplx_dp

  ! tridiag_matmul_sym_cnst definition for single-precision real values
  subroutine tridiag_matmul_sym_cnst_real_sp(diag_cnst, o_diag_cnst, vec, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag_cnst, o_diag_cnst
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst.src"
  end subroutine tridiag_matmul_sym_cnst_real_sp

  ! tridiagonal_matmul_sym_cnst definition for single-precision complex values
  subroutine tridiag_matmul_sym_cnst_cmplx_sp(diag_cnst, o_diag_cnst, vec, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag_cnst, o_diag_cnst
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst.src"
  end subroutine tridiag_matmul_sym_cnst_cmplx_sp


  ! subroutine tridiag_matmul_sym_cnst_diag
  !
  ! Multiply symmetric tridiagonal matrix by a vector, where the diagonal
  ! elements are equal
  !
  ! diag_cnst :: diagonal element
  ! o_diag(n - 1) :: off-diagonal elements
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplications

  ! tridiag_matmul_sym_cnst_diag definition for double-precision real values
  subroutine tridiag_matmul_sym_cnst_diag_real_dp(diag_cnst, o_diag, vec, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag_cnst, o_diag(size(vec) - 1)
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst_diag.src"
  end subroutine tridiag_matmul_sym_cnst_diag_real_dp

  ! tridiag_matmul_sym_cnst_diag definition for double-precision complex values
  subroutine tridiag_matmul_sym_cnst_diag_cmplx_dp(diag_cnst, o_diag, vec, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag_cnst, o_diag(size(vec) - 1)
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst_diag.src"
  end subroutine tridiag_matmul_sym_cnst_diag_cmplx_dp

  ! tridiag_matmul_sym_cnst_diag definition for single-precision real values
  subroutine tridiag_matmul_sym_cnst_diag_real_sp(diag_cnst, o_diag, vec, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag_cnst, o_diag(size(vec) - 1)
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst_diag.src"
  end subroutine tridiag_matmul_sym_cnst_diag_real_sp

  ! tridiag_matmul_sym_cnst_diag definition for single-precision complex values
  subroutine tridiag_matmul_sym_cnst_diag_cmplx_sp(diag_cnst, o_diag, vec, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag_cnst, o_diag(size(vec) - 1)
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst_diag.src"
  end subroutine tridiag_matmul_sym_cnst_diag_cmplx_sp


  ! subroutine tridiag_matmul_sym_cnst_odiag
  !
  ! Multiply symmetric tridiagonal matrix by a vector, where the off-diagonal
  ! elements are equal, band-wise
  !
  ! diag(n) :: diagonal elements
  ! o_diag_cnst :: off-diagonal element
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplications

  ! tridiag_matmul_sym_cnst_odiag definition for double-precision real values
  subroutine tridiag_matmul_sym_cnst_odiag_real_dp(diag, o_diag_cnst, vec, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec)), o_diag_cnst
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst_odiag.src"
  end subroutine tridiag_matmul_sym_cnst_odiag_real_dp

  ! tridiag_matmul_sym_cnst_odiag definition for double-precision complex
  ! values
  subroutine tridiag_matmul_sym_cnst_odiag_cmplx_dp(diag, o_diag_cnst, vec, &
       res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec)), o_diag_cnst
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst_odiag.src"
  end subroutine tridiag_matmul_sym_cnst_odiag_cmplx_dp

  ! tridiag_matmul_sym_cnst_odiag definition for single-precision real values
  subroutine tridiag_matmul_sym_cnst_odiag_real_sp(diag, o_diag_cnst, vec, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec)), o_diag_cnst
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst_odiag.src"
  end subroutine tridiag_matmul_sym_cnst_odiag_real_sp

  ! tridiag_matmul_sym_cnst_odiag definition for single-precision complex
  ! values
  subroutine tridiag_matmul_sym_cnst_odiag_cmplx_sp(diag, o_diag_cnst, vec, &
       res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec)), o_diag_cnst
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_matmul_src/sym_cnst_odiag.src"
  end subroutine tridiag_matmul_sym_cnst_odiag_cmplx_sp


end module tridiag_matmul
