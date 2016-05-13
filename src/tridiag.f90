! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods for solving tridiagonal systems using the Thomas algorithm (modified
! Gaussian elimination)
!
! Note that this algorithm is only guaranteed to be stable for diagonally
! dominant tridiagonal systems; that is, when the norm of each diagonal element
! is greater than or equal to the sum of the norms of the off diagonal
! elements.
module tridiag


  use iso_fortran_env, only: real32, real64, int32

  implicit none
  private


  public :: tridiag_general
  interface tridiag_general
     ! General tridiagonal system
     module procedure tridiag_general_real_dp
     module procedure tridiag_general_cmplx_dp
     module procedure tridiag_general_real_sp
     module procedure tridiag_general_cmplx_sp
  end interface tridiag_general

  
  public :: tridiag_constant
  interface tridiag_constant
     ! Tridiagonal systems where all elements are equal, band-wise
     module procedure tridiag_constant_real_dp
     module procedure tridiag_constant_cmplx_dp
     module procedure tridiag_constant_real_sp
     module procedure tridiag_constant_cmplx_sp
  end interface tridiag_constant


  public :: tridiag_constant_offdiag
  interface tridiag_constant_offdiag
     ! Tridiagonal systems where off-diagonal elements are equal, band-wise
     module procedure tridiag_constant_offdiag_real_dp
     module procedure tridiag_constant_offdiag_cmplx_dp
     module procedure tridiag_constant_offdiag_real_sp
     module procedure tridiag_constant_offdiag_cmplx_sp
  end interface tridiag_constant_offdiag


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

    include "./tridiag_src/tridiag_backsweep.src"
  end subroutine tridiag_backsweep_real_dp

  ! tridiag_backsweep definition for double precision complex values
  subroutine tridiag_backsweep_cmplx_dp(mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: mat_coeff(:), vec_coeff(:)
    complex(dp), intent(inout) :: res(:)

    include "./tridiag_src/tridiag_backsweep.src"
  end subroutine tridiag_backsweep_cmplx_dp

  ! tridiag_backsweep definition for single precision real values
  subroutine tridiag_backsweep_real_sp(mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: mat_coeff(:), vec_coeff(:)
    real(sp), intent(inout) :: res(:)

    include "./tridiag_src/tridiag_backsweep.src"
  end subroutine tridiag_backsweep_real_sp

  ! tridiag_backsweep definition for single precision complex values
  subroutine tridiag_backsweep_cmplx_sp(mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: mat_coeff(:), vec_coeff(:)
    complex(sp), intent(inout) :: res(:)

    include "./tridiag_src/tridiag_backsweep.src"
  end subroutine tridiag_backsweep_cmplx_sp


  ! subroutine tridiag_general
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

  ! tridiag_general definition for double precision real values
  subroutine tridiag_general_real_dp(diag, u_diag, l_diag, vec, mat_coeff, &
    vec_coeff, res)

    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec)), u_diag(size(vec) - 1), &
         l_diag(size(vec) - 1)
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_general.src"
  end subroutine tridiag_general_real_dp

  ! tridiag_general definition for double precision complex values
  subroutine tridiag_general_cmplx_dp(diag, u_diag, l_diag, vec, mat_coeff, &
    vec_coeff, res)

    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec)), u_diag(size(vec) - 1), &
         l_diag(size(vec) - 1)
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_general.src"
  end subroutine tridiag_general_cmplx_dp

  ! tridiag_general definition for single precision real values
  subroutine tridiag_general_real_sp(diag, u_diag, l_diag, vec, mat_coeff, &
    vec_coeff, res)

    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec)), u_diag(size(vec) - 1), &
         l_diag(size(vec) - 1)
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_general.src"
  end subroutine tridiag_general_real_sp

  ! tridiag_general definition for single precision complex values
  subroutine tridiag_general_cmplx_sp(diag, u_diag, l_diag, vec, mat_coeff, &
    vec_coeff, res)

    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec)), u_diag(size(vec) - 1), &
         l_diag(size(vec) - 1)
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_general.src"
  end subroutine tridiag_general_cmplx_sp


  ! subroutine tridiag_constant
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

  ! tridiag_constant definition for double-precision real values
  subroutine tridiag_constant_real_dp(diag_cnst, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_constant.src"
  end subroutine tridiag_constant_real_dp

  ! tridiag_constant definition for double-precision complex values
  subroutine tridiag_constant_cmplx_dp(diag_cnst, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_constant.src"
  end subroutine tridiag_constant_cmplx_dp

  ! tridiag_constant definition for single-precision real values
  subroutine tridiag_constant_real_sp(diag_cnst, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_constant.src"
  end subroutine tridiag_constant_real_sp

  ! tridiag_constant definition for single-precision complex values
  subroutine tridiag_constant_cmplx_sp(diag_cnst, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst, vec(:)
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_constant.src"
  end subroutine tridiag_constant_cmplx_sp


  ! subroutine tridiag_constant_offdiag
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

  ! tridiag_constant_offdiag definition for double-precision real values
  subroutine tridiag_constant_offdiag_real_dp(diag, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec)), u_diag_cnst, l_diag_cnst
    real(dp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(dp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_constant_offdiag.src"
  end subroutine tridiag_constant_offdiag_real_dp

  ! tridiag_constant_offdiag definition for double-precision complex values
  subroutine tridiag_constant_offdiag_cmplx_dp(diag, u_diag_cnst, &
       l_diag_cnst, vec, mat_coeff, vec_coeff, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec)), u_diag_cnst, l_diag_cnst
    complex(dp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(dp), intent(inout) :: res(size(vec))
;
    include "./tridiag_src/tridiag_constant_offdiag.src"
  end subroutine tridiag_constant_offdiag_cmplx_dp

  ! tridiag_constant_offdiag definition for single-precision real values
  subroutine tridiag_constant_offdiag_real_sp(diag, u_diag_cnst, l_diag_cnst, &
       vec, mat_coeff, vec_coeff, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec)), u_diag_cnst, l_diag_cnst
    real(sp), intent(inout) :: mat_coeff(size(vec) - 1), vec_coeff(size(vec))
    real(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_constant_offdiag.src"
  end subroutine tridiag_constant_offdiag_real_sp

  ! tridiag_constant_offdiag definition for double-precision complex values
  subroutine tridiag_constant_offdiag_cmplx_sp(diag, u_diag_cnst, &
       l_diag_cnst, vec, mat_coeff, vec_coeff, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec)), u_diag_cnst, l_diag_cnst
    complex(sp), intent(inout) :: mat_coeff(size(vec) - 1), &
         vec_coeff(size(vec))
    complex(sp), intent(inout) :: res(size(vec))

    include "./tridiag_src/tridiag_constant_offdiag.src"
  end subroutine tridiag_constant_offdiag_cmplx_sp


end module tridiag
