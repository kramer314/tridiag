! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods to verify the stability of the Thomas algorithm for solving
! tridiagonal systems. The Thomas algorithm is stable if the tridiagonal
! system to be solved is diagonally dominant (either row- or column-wise), or
! if the system is positive semidefinite. This module does not diagonalize
! the tridiagonal system or compute its Cholesky factorization, etc., so it
! only tests for diagonal dominance.
!
! This module provides the following subroutines:
!
! tridiag_stable_gen(...) -> general tridiagonal matrix
! tridiag_stable_cnst(...) -> all elements constant, band-wise
! tridiag_stable_cnst_diag(...) -> diagonal elements constant
! tridiag_stable_cnst_odiag(...) -> off_diagonal elements constant, band-wise
!
! tridiag_stable_sym(...) -> symmetric tridiagonal matrix
! tridiag_stable_sym_cnst(...) -> symmetric, all elements constant, band-wise
! tridiag_stable_sym_cnst_diag(..) -> symmetric, diagonal elements constant
! tridag_stable_sym_cnst_odiag(...) -> symmetric, off-diagonal elements
!   constant, band-wise
!
! Note that the methods for symmetric tridiagonal matrices are simply wrappers
! around the methods for the more general asymmetric case. There are currently
! no symemtry-specific numerical approximations.
module tridiag_stable


  ! Numeric precision parameters
  use iso_fortran_env, only: sp=>real32, dp=>real64, ip=>int32

  implicit none
  private


  public :: tridiag_stable_gen
  interface tridiag_stable_gen
     ! General tridiagonal matrix
     module procedure tridiag_stable_gen_real_dp
     module procedure tridiag_stable_gen_cmplx_dp
     module procedure tridiag_stable_gen_real_sp
     module procedure tridiag_stable_gen_cmplx_sp
  end interface tridiag_stable_gen


  public :: tridiag_stable_cnst
  interface tridiag_stable_cnst
     ! Tridiagonal matrix where all elements are equal, band-wise
     module procedure tridiag_stable_cnst_real_dp
     module procedure tridiag_stable_cnst_cmplx_dp
     module procedure tridiag_stable_cnst_real_sp
     module procedure tridiag_stable_cnst_cmplx_sp
  end interface tridiag_stable_cnst


  public :: tridiag_stable_cnst_diag
  interface tridiag_stable_cnst_diag
     ! Tridiagonal matrix where diagonal elements are equal
     module procedure tridiag_stable_cnst_diag_real_dp
     module procedure tridiag_stable_cnst_diag_cmplx_dp
     module procedure tridiag_stable_cnst_diag_real_sp
     module procedure tridiag_stable_cnst_diag_cmplx_sp
  end interface tridiag_stable_cnst_diag


  public :: tridiag_stable_cnst_odiag
  interface tridiag_stable_cnst_odiag
     ! Tridiagonal matrix where off-diagonal elements are equal, band-wise
     module procedure tridiag_stable_cnst_odiag_real_dp
     module procedure tridiag_stable_cnst_odiag_cmplx_dp
     module procedure tridiag_stable_cnst_odiag_real_sp
     module procedure tridiag_stable_cnst_odiag_cmplx_sp
  end interface tridiag_stable_cnst_odiag


  public :: tridiag_stable_sym
  interface tridiag_stable_sym
     ! Symmetric tridiagonal matrix
     module procedure tridiag_stable_sym_real_dp
     module procedure tridiag_stable_sym_cmplx_dp
     module procedure tridiag_stable_sym_real_sp
     module procedure tridiag_stable_sym_cmplx_sp
  end interface tridiag_stable_sym


  public :: tridiag_stable_sym_cnst
  interface tridiag_stable_sym_cnst
     ! Symmetric tridiagonal matrix where all elements are equal, band-wise
     module procedure tridiag_stable_sym_cnst_real_dp
     module procedure tridiag_stable_sym_cnst_cmplx_dp
     module procedure tridiag_stable_sym_cnst_real_sp
     module procedure tridiag_stable_sym_cnst_cmplx_sp
  end interface tridiag_stable_sym_cnst


  public :: tridiag_stable_sym_cnst_diag
  interface tridiag_stable_sym_cnst_diag
     ! Symmetric tridiagonal matrix where diagonal elements are equal
     module procedure tridiag_stable_sym_cnst_diag_real_dp
     module procedure tridiag_stable_sym_cnst_diag_cmplx_dp
     module procedure tridiag_stable_sym_cnst_diag_real_sp
     module procedure tridiag_stable_sym_cnst_diag_cmplx_sp
  end interface tridiag_stable_sym_cnst_diag


  public :: tridiag_stable_sym_cnst_odiag
  interface tridiag_stable_sym_cnst_odiag
     ! Symmetric tridiagonal matrix where off-diagonal elements are equal,
     ! band-wise
     module procedure tridiag_stable_sym_cnst_odiag_real_dp
     module procedure tridiag_stable_sym_cnst_odiag_cmplx_dp
     module procedure tridiag_stable_sym_cnst_odiag_real_sp
     module procedure tridiag_stable_sym_cnst_odiag_cmplx_sp
  end interface tridiag_stable_sym_cnst_odiag


contains


  ! subroutine tridiag_stable_gen
  !
  ! General tridiagonal system
  !
  ! diag(n) :: diagonal entries
  ! u_diag(n - 1) :: upper diagonal entries
  ! l_diag(n - 1) :: lower diagonal entries

  ! tridiag_stable_gen definition for double-precision real values
  logical function tridiag_stable_gen_real_dp(diag, u_diag, l_diag) &
       result(val)
    real(dp), intent(in) :: diag(:)
    real(dp), intent(in) :: u_diag(size(diag)), l_diag(size(diag))
    integer(ip), parameter :: fp = dp

    include "./tridiag_stable_src/gen.src"
  end function tridiag_stable_gen_real_dp

  ! tridiag_stable_gen definition for double-precision complex values
  logical function tridiag_stable_gen_cmplx_dp(diag, u_diag, l_diag) &
       result(val)
    complex(dp), intent(in) :: diag(:)
    complex(dp), intent(in) :: u_diag(size(diag)), l_diag(size(diag))
    integer(ip), parameter :: fp = dp

    include "./tridiag_stable_src/gen.src"
  end function tridiag_stable_gen_cmplx_dp

  ! tridiag_stable_gen definition for single-precision real values
  logical function tridiag_stable_gen_real_sp(diag, u_diag, l_diag) &
       result(val)
    real(sp), intent(in) :: diag(:)
    real(sp), intent(in) :: u_diag(size(diag)), l_diag(size(diag))
    integer(ip), parameter :: fp = sp

    include "./tridiag_stable_src/gen.src"
  end function tridiag_stable_gen_real_sp

  ! tridiag_stable_gen definition for single-precision complex values
  logical function tridiag_stable_gen_cmplx_sp(diag, u_diag, l_diag) &
       result(val)
    complex(sp), intent(in) :: diag(:)
    complex(sp), intent(in) :: u_diag(size(diag)), l_diag(size(diag))
    integer(ip), parameter :: fp = sp

    include "./tridiag_stable_src/gen.src"
  end function tridiag_stable_gen_cmplx_sp


  ! subroutine tridiag_stable_cnst
  !
  ! Tridiagonal matrix where all band-elements are equal (band-wise)
  !
  ! diag_cnst :: diagonal entry
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnst :: lower diagonal entry

  ! tridiag_stable_cnst definition for double-precision real values
  logical function tridiag_stable_cnst_real_dp(diag_cnst, u_diag_cnst, &
       l_diag_cnst) result(val)
    real(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst

    include "./tridiag_stable_src/cnst.src"
  end function tridiag_stable_cnst_real_dp

  ! tridiag_stable_cnst definition for double-precision complex values
  logical function tridiag_stable_cnst_cmplx_dp(diag_cnst, u_diag_cnst, &
       l_diag_cnst) result(val)
    complex(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst

    include "./tridiag_stable_src/cnst.src"
  end function tridiag_stable_cnst_cmplx_dp

  ! tridiag_stable_cnst definition for single-precision real values
  logical function tridiag_stable_cnst_real_sp(diag_cnst, u_diag_cnst, &
       l_diag_cnst) result(val)
    real(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst

    include "./tridiag_stable_src/cnst.src"
  end function tridiag_stable_cnst_real_sp

  ! tridiag_stable_cnst definition for single-precision complex values
  logical function tridiag_stable_cnst_cmplx_sp(diag_cnst, u_diag_cnst, &
       l_diag_cnst) result(val)
    complex(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst

    include "./tridiag_stable_src/cnst.src"
  end function tridiag_stable_cnst_cmplx_sp


  ! subroutine tridiag_cnst_diag
  !
  ! Tridiagonal system where diagonal elements are constant
  !
  ! diag_cnst :: diagonal entry
  ! u_diag(n - 1) :: upper diagonal entries
  ! l_diag(n - 1) :: lower diagonal entries

  ! tridiag_stable_cnst_diag definition for double-precision real values
  logical function tridiag_stable_cnst_diag_real_dp(diag_cnst, u_diag, &
       l_diag) result(val)
    real(dp), intent(in) :: u_diag(:)
    real(dp), intent(in) :: l_diag(size(u_diag)), diag_cnst
    integer(ip), parameter :: fp = dp

    include "./tridiag_stable_src/cnst_diag.src"
  end function tridiag_stable_cnst_diag_real_dp

  ! tridiag_stable_cnst_diag definition for double-precision complex values
  logical function tridiag_stable_cnst_diag_cmplx_dp(diag_cnst, u_diag, &
       l_diag) result(val)
    complex(dp), intent(in) :: u_diag(:)
    complex(dp), intent(in) :: l_diag(size(u_diag)), diag_cnst
    integer(ip), parameter :: fp = dp

    include "./tridiag_stable_src/cnst_diag.src"
  end function tridiag_stable_cnst_diag_cmplx_dp

  ! tridiag_stable_cnst_diag definition for single-precision real values
  logical function tridiag_stable_cnst_diag_real_sp(diag_cnst, u_diag, &
       l_diag) result(val)
    real(sp), intent(in) :: u_diag(:)
    real(sp), intent(in) :: l_diag(size(u_diag)), diag_cnst
    integer(ip), parameter :: fp = sp

    include "./tridiag_stable_src/cnst_diag.src"
  end function tridiag_stable_cnst_diag_real_sp

  ! tridiag_stable_cnst_diag definition for single-precision complex values
  logical function tridiag_stable_cnst_diag_cmplx_sp(diag_cnst, u_diag, &
       l_diag) result(val)
    complex(sp), intent(in) :: u_diag(:)
    complex(sp), intent(in) :: l_diag(size(u_diag)), diag_cnst
    integer(ip), parameter :: fp = sp

    include "./tridiag_stable_src/cnst_diag.src"
  end function tridiag_stable_cnst_diag_cmplx_sp


  ! subroutine tridiag_cnst_odiag
  !
  ! Tridiagonal system where off-diagonal elements are constant, band-wise
  !
  ! diag(n) :: diagonal entries
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnst :: lower diagonal entry

  ! tridiag_stable_cnst_odiag definition for double-precision real values
  logical function tridiag_stable_cnst_odiag_real_dp(diag, u_diag_cnst, &
       l_diag_cnst) result(val)
    real(dp), intent(in) :: diag(:), u_diag_cnst, l_diag_cnst
    integer(ip), parameter :: fp = dp

    include "./tridiag_stable_src/cnst_odiag.src"
  end function tridiag_stable_cnst_odiag_real_dp

  ! tridiag_stable_cnst_odiag definition for double-precision complex values
  logical function tridiag_stable_cnst_odiag_cmplx_dp(diag, u_diag_cnst, &
       l_diag_cnst) result(val)
    complex(dp), intent(in) :: diag(:), u_diag_cnst, l_diag_cnst
    integer(ip), parameter :: fp = dp

    include "./tridiag_stable_src/cnst_odiag.src"
  end function tridiag_stable_cnst_odiag_cmplx_dp

  ! tridiag_stable_cnst_odiag definition for single-precision real values
  logical function tridiag_stable_cnst_odiag_real_sp(diag, u_diag_cnst, &
       l_diag_cnst) result(val)
    real(sp), intent(in) :: diag(:), u_diag_cnst, l_diag_cnst
    integer(ip), parameter :: fp = sp

    include "./tridiag_stable_src/cnst_odiag.src"
  end function tridiag_stable_cnst_odiag_real_sp

  ! tridiag_stable_cnst_odiag definition for single-precision complex values
  logical function tridiag_stable_cnst_odiag_cmplx_sp(diag, u_diag_cnst, &
       l_diag_cnst) result(val)
    complex(sp), intent(in) :: diag(:), u_diag_cnst, l_diag_cnst
    integer(ip), parameter :: fp = sp

    include "./tridiag_stable_src/cnst_odiag.src"
  end function tridiag_stable_cnst_odiag_cmplx_sp


  ! subroutine tridiag_stable_sym
  !
  ! Symmetric tridiagonal system
  !
  ! diag(n) :: diagonal entries
  ! o_diag(n - 1) :: off-diagonal entries

  ! tridiag_stable_sym definition for double-precision real values
  logical function tridiag_stable_sym_real_dp(diag, o_diag) result(val)
    real(dp), intent(in) :: diag(:)
    real(dp), intent(in) :: o_diag(size(diag))

    include "./tridiag_stable_src/sym.src"
  end function tridiag_stable_sym_real_dp

  ! tridiag_stable_sym definition for double-precision complex values
  logical function tridiag_stable_sym_cmplx_dp(diag, o_diag) result(val)
    complex(dp), intent(in) :: diag(:)
    complex(dp), intent(in) :: o_diag(size(diag))

    include "./tridiag_stable_src/sym.src"
  end function tridiag_stable_sym_cmplx_dp

  ! tridiag_stable_sym definition for single-precision real values
  logical function tridiag_stable_sym_real_sp(diag, o_diag) result(val)
    real(sp), intent(in) :: diag(:)
    real(sp), intent(in) :: o_diag(size(diag))

    include "./tridiag_stable_src/sym.src"
  end function tridiag_stable_sym_real_sp

  ! tridiag_stable_sym definition for single-precision complex values
  logical function tridiag_stable_sym_cmplx_sp(diag, o_diag) result(val)
    complex(sp), intent(in) :: diag(:)
    complex(sp), intent(in) :: o_diag(size(diag))

    include "./tridiag_stable_src/sym.src"
  end function tridiag_stable_sym_cmplx_sp


  ! subroutine tridiag_stable_sym_cnst
  !
  ! Symmetric tridiagonal system where all elements are equal, band-wise
  !
  ! diag_cnst :: diagonal entry
  ! o_diag_cnst :: off-diagonal entries

  ! tridiag_stable_sym_cnst definition for double-precision real values
  logical function tridiag_stable_sym_cnst_real_dp(diag_cnst, &
       o_diag_cnst) result(val)
    real(dp), intent(in) :: diag_cnst, o_diag_cnst

    include "./tridiag_stable_src/sym_cnst.src"
  end function tridiag_stable_sym_cnst_real_dp

  ! tridiag_stable_sym_cnst definition for double-precision complex values
  logical function tridiag_stable_sym_cnst_cmplx_dp(diag_cnst, &
       o_diag_cnst) result(val)
    complex(dp), intent(in) :: diag_cnst, o_diag_cnst

    include "./tridiag_stable_src/sym_cnst.src"
  end function tridiag_stable_sym_cnst_cmplx_dp

  ! tridiag_stable_sym_cnst definition for single-precision real values
  logical function tridiag_stable_sym_cnst_real_sp(diag_cnst, &
       o_diag_cnst) result(val)
    real(sp), intent(in) :: diag_cnst, o_diag_cnst

    include "./tridiag_stable_src/sym_cnst.src"
  end function tridiag_stable_sym_cnst_real_sp

  ! tridiag_stable_sym_cnst definition for single-precision complex values
  logical function tridiag_stable_sym_cnst_cmplx_sp(diag_cnst, &
       o_diag_cnst) result(val)
    complex(sp), intent(in) :: diag_cnst, o_diag_cnst

    include "./tridiag_stable_src/sym_cnst.src"
  end function tridiag_stable_sym_cnst_cmplx_sp


  ! subroutine tridiag_stable_sym_cnst_diag
  !
  ! Symmetric tridiagonal system with constant diagonal entries
  !
  ! diag_cnst :: diagonal entry
  ! o_diag(n - 1) :: off-diagonal entries

  ! tridiag_stable_sym_cnst_diag definition for double-precision real values
  logical function tridiag_stable_sym_cnst_diag_real_dp(diag_cnst, o_diag) &
       result(val)
    real(dp), intent(in) :: diag_cnst, o_diag(:)

    include "./tridiag_stable_src/sym_cnst_diag.src"
  end function tridiag_stable_sym_cnst_diag_real_dp

  ! tridiag_stable_sym_cnst_diag definition for double-precision complex values
  logical function tridiag_stable_sym_cnst_diag_cmplx_dp(diag_cnst, o_diag) &
       result(val)
    complex(dp), intent(in) :: diag_cnst, o_diag(:)

    include "./tridiag_stable_src/sym_cnst_diag.src"
  end function tridiag_stable_sym_cnst_diag_cmplx_dp

  ! tridiag_stable_sym_cnst_diag definition for single-precision real values
  logical function tridiag_stable_sym_cnst_diag_real_sp(diag_cnst, o_diag) &
       result(val)
    real(sp), intent(in) :: diag_cnst, o_diag(:)

    include "./tridiag_stable_src/sym_cnst_diag.src"
  end function tridiag_stable_sym_cnst_diag_real_sp

  ! tridiag_stable_sym_cnst_diag definition for single-precision complex values
  logical function tridiag_stable_sym_cnst_diag_cmplx_sp(diag_cnst, o_diag) &
       result(val)
    complex(sp), intent(in) :: diag_cnst, o_diag(:)

    include "./tridiag_stable_src/sym_cnst_diag.src"
  end function tridiag_stable_sym_cnst_diag_cmplx_sp


  ! subroutine tridiag_stable_sym_cnst_odiag
  !
  ! Symmetric tridiagonal system with constant off-diagonal entries
  !
  ! diag(n) :: diagonal entries
  ! o_diag_cnst :: off-diagonal entry

  ! tridiag_stable_sym_cnst_odiag definition for double-precision real values
  logical function tridiag_stable_sym_cnst_odiag_real_dp(diag, o_diag_cnst) &
       result(val)
    real(dp), intent(in) :: diag(:), o_diag_cnst

    include "./tridiag_stable_src/sym_cnst_odiag.src"
  end function tridiag_stable_sym_cnst_odiag_real_dp

  ! tridiag_stable_sym_cnst_odiag definition for double-precision complex
  ! values
  logical function tridiag_stable_sym_cnst_odiag_cmplx_dp(diag, o_diag_cnst) &
       result(val)
    complex(dp), intent(in) :: diag(:), o_diag_cnst

    include "./tridiag_stable_src/sym_cnst_odiag.src"
  end function tridiag_stable_sym_cnst_odiag_cmplx_dp

  ! tridiag_stable_sym_cnst_odiag definition for single-precision real values
  logical function tridiag_stable_sym_cnst_odiag_real_sp(diag, o_diag_cnst) &
       result(val)
    real(sp), intent(in) :: diag(:), o_diag_cnst

    include "./tridiag_stable_src/sym_cnst_odiag.src"
  end function tridiag_stable_sym_cnst_odiag_real_sp

  ! tridiag_stable_sym_cnst_odiag definition for single-precision complex
  ! values
  logical function tridiag_stable_sym_cnst_odiag_cmplx_sp(diag, o_diag_cnst) &
       result(val)
    complex(sp), intent(in) :: diag(:), o_diag_cnst

    include "./tridiag_stable_src/sym_cnst_odiag.src"
  end function tridiag_stable_sym_cnst_odiag_cmplx_sp


end module tridiag_stable
