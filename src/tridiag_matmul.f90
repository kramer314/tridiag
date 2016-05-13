! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file in the top-level directory of this distribution

! Methods for multiplying tridiagonal matrices by a vector; these are in
! general much more efficient than standard matrix multiplication.
module tridiag_matmul


  use iso_fortran_env, only: real32, real64, int32

  implicit none
  private


  public :: tridiag_matmul_general
  interface tridiag_matmul_general
     ! Multiply general tridiagonal matrix by vector
     module procedure tridiag_matmul_general_real_dp
     module procedure tridiag_matmul_general_cmplx_dp
     module procedure tridiag_matmul_general_real_sp
     module procedure tridiag_matmul_general_cmplx_sp
  end interface tridiag_matmul_general


  public :: tridiag_matmul_constant
  interface tridiag_matmul_constant
     ! Multiply tridiagonal matrix by vector where all elements are equal,
     ! band-wise
     module procedure tridiag_matmul_constant_real_dp
     module procedure tridiag_matmul_constant_cmplx_dp
     module procedure tridiag_matmul_constant_real_sp
     module procedure tridiag_matmul_constant_cmplx_sp
  end interface tridiag_matmul_constant


  public :: tridiag_matmul_constant_offdiag
  interface tridiag_matmul_constant_offdiag
     ! Multiply tridiagonal matrix by vector where off-diagonal elements are
     ! equal, band-wise
     module procedure tridiag_matmul_constant_real_dp
     module procedure tridiag_matmul_constant_cmplx_dp
     module procedure tridiag_matmul_constant_real_sp
     module procedure tridiag_matmul_constant_cmplx_sp
  end interface tridiag_matmul_constant_offdiag


  ! Numeric precision parameters
  integer(int32), parameter :: ip = int32
  integer(ip), parameter :: sp = real32
  integer(ip), parameter :: dp = real64


contains


  ! subroutine tridiag_matmul
  !
  ! Multiply general tridiagonal matrix by vector
  !
  ! diag(n) :: main diagonal entries
  ! u_diag(n - 1) :: upper diagonal entries
  ! l_diag(n - 1) :: upper diagonal entries
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_general definition for double precision real values
  subroutine tridiag_matmul_general_real_dp(diag, u_diag, l_diag, vec, res)
    real(dp), intent(in) :: vec(:)
    real(dp), intent(in) :: diag(size(vec))
    real(dp), intent(in) :: u_diag(size(diag) - 1)
    real(dp), intent(in) :: l_diag(size(diag) - 1)
    real(dp), intent(out) :: res(size(diag))

    include "./tridiag_matmul_src/tridiag_matmul_general.src"
  end subroutine tridiag_matmul_general_real_dp

  ! tridiag_matmul_general definition for double precision complex values
  subroutine tridiag_matmul_general_cmplx_dp(diag, u_diag, l_diag, vec, res)
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(in) :: diag(size(vec))
    complex(dp), intent(in) :: u_diag(size(diag) - 1)
    complex(dp), intent(in) :: l_diag(size(diag) - 1)
    complex(dp), intent(out) :: res(size(diag))

    include "./tridiag_matmul_src/tridiag_matmul_general.src"
  end subroutine tridiag_matmul_general_cmplx_dp

  ! tridiag_matmul_general definition for single precision real values
  subroutine tridiag_matmul_general_real_sp(diag, u_diag, l_diag, vec, res)
    real(sp), intent(in) :: vec(:)
    real(sp), intent(in) :: diag(size(vec))
    real(sp), intent(in) :: u_diag(size(diag) - 1)
    real(sp), intent(in) :: l_diag(size(diag) - 1)
    real(sp), intent(out) :: res(size(diag))

    include "./tridiag_matmul_src/tridiag_matmul_general.src"
  end subroutine tridiag_matmul_general_real_sp

  ! tridiag_matmul_general definition for single precision complex values
  subroutine tridiag_matmul_general_cmplx_sp(diag, u_diag, l_diag, vec, res)
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(in) :: diag(size(vec))
    complex(sp), intent(in) :: u_diag(size(diag) - 1)
    complex(sp), intent(in) :: l_diag(size(diag) - 1)
    complex(sp), intent(out) :: res(size(diag))

    include "./tridiag_matmul_src/tridiag_matmul_general.src"
  end subroutine tridiag_matmul_general_cmplx_sp


  ! subroutine tridiag_matmul_constant
  !
  ! Multiply tridiagonal matrix by a vector, where the band elements are
  ! equal (band-wise)
  !
  ! diag_cnst :: diagonal entry
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnnst :: lower diagonal entry
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_constant definition for double precision real values
  subroutine tridiag_matmul_constant_real_dp(diag_cnst, u_diag_cnst, &
       l_diag_cnst, vec, res)
    real(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst
    real(dp), intent(in) :: vec(:)
    real(dp), intent(out) :: res(size(vec))

    include "./tridiag_matmul_src/tridiag_matmul_constant.src"
  end subroutine tridiag_matmul_constant_real_dp

  ! tridiag_matmul_constant definition for double precision complex values
  subroutine tridiag_matmul_constant_cmplx_dp(diag_cnst, u_diag_cnst, &
       l_diag_cnst, vec, res)
    complex(dp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(out) :: res(size(vec))

    include "./tridiag_matmul_src/tridiag_matmul_constant.src"
  end subroutine tridiag_matmul_constant_cmplx_dp

  ! tridiag_mmatmul_constant definition for single precision real values
  subroutine tridiag_matmul_constant_real_sp(diag_cnst, u_diag_cnst, &
       l_diag_cnst, vec, res)
    real(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst
    real(sp), intent(in) :: vec(:)
    real(sp), intent(out) :: res(size(vec))

    include "./tridiag_matmul_src/tridiag_matmul_constant.src"
  end subroutine tridiag_matmul_constant_real_sp

  ! tridiag_matmul_constant defintion for single precision complex values
  subroutine tridiag_matmul_constant_cmplx_sp(diag_cnst, u_diag_cnst, &
       l_diag_cnst, vec, res)
    complex(sp), intent(in) :: diag_cnst, u_diag_cnst, l_diag_cnst
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(out) :: res(size(vec))

    include "./tridiag_matmul_src/tridiag_matmul_constant.src"
  end subroutine tridiag_matmul_constant_cmplx_sp


  ! subroutine tridiag_matmul_constant_offdiag
  !
  ! Multiple tridiagonal matrix by a vector, where the off-diagonal band
  ! elements are equal (band-wise)
  !
  ! diag(n) :: diagonal elements
  ! u_diag_cnst :: upper diagonal entry
  ! l_diag_cnst :: lower diagonal entry
  ! vec(n) :: vector to multiply
  ! res(n) :: result of multiplication

  ! tridiag_matmul_constant_offdiag definition for double precision real values
  subroutine tridiag_matmul_constant_offdiag_real_dp(diag, u_diag_cnst, &
       l_diag_cnst, vec, res)
    real(dp), intent(in) :: diag(:)
    real(dp), intent(in) :: u_diag_cnst, l_diag_cnst
    real(dp), intent(in) :: vec(:)
    real(dp), intent(out) :: res(size(vec))

    include "./tridiag_matmul_src/tridiag_matmul_constant_offdiag.src"
  end subroutine tridiag_matmul_constant_offdiag_real_dp

  ! tridiag_matmul_constant_offdiag definition for double precision complex
  ! values
  subroutine tridiag_matmul_constant_offdiag_cmplx_dp(diag, u_diag_cnst, &
       l_diag_cnst, vec, res)
    complex(dp), intent(in) :: diag(:)
    complex(dp), intent(in) :: u_diag_cnst, l_diag_cnst
    complex(dp), intent(in) :: vec(:)
    complex(dp), intent(out) :: res(size(vec))

    include "./tridiag_matmul_src/tridiag_matmul_constant_offdiag.src"
  end subroutine tridiag_matmul_constant_offdiag_cmplx_dp

  ! tridiag_matmul_constant_offdiag definition for single precision real values
  subroutine tridiag_matmul_constant_offdiag_real_sp(diag, u_diag_cnst, &
       l_diag_cnst, vec, res)
    real(sp), intent(in) :: diag(:)
    real(sp), intent(in) :: u_diag_cnst, l_diag_cnst
    real(sp), intent(in) :: vec(:)
    real(sp), intent(out) :: res(size(vec))

    include "./tridiag_matmul_src/tridiag_matmul_constant_offdiag.src"
  end subroutine tridiag_matmul_constant_offdiag_real_sp

  ! tridiag_matmul_constant_offdiag definition for single precision complex
  ! values
  subroutine tridiag_matmul_constant_offdiag_cmplx_sp(diag, u_diag_cnst, &
       l_diag_cnst, vec, res)
    complex(sp), intent(in) :: diag(:)
    complex(sp), intent(in) :: u_diag_cnst, l_diag_cnst
    complex(sp), intent(in) :: vec(:)
    complex(sp), intent(out) :: res(size(vec))

    include "./tridiag_matmul_src/tridiag_matmul_constant_offdiag.src"
  end subroutine tridiag_matmul_constant_offdiag_cmplx_sp


end module tridiag_matmul
