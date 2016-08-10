module cnst_diag


  use precision, only: ip, sp, dp
  use test_vars, only: fmt_test_name, fmt_test_res

  use tridiag, only: tridiag_cnst_diag
  use tridiag_matmul, only: tridiag_matmul_cnst_diag


  implicit none
  private


  public :: cnst_diag_tests

  real(sp) :: dev_sp
  real(dp) :: dev_dp


contains

  subroutine cnst_diag_tests()

    write(*, fmt_test_name, advance="no") "cnst_diag"

    dev_dp = cnst_diag_real_dp()
    write(*, fmt_test_res, advance="no") dev_dp

    dev_sp = cnst_diag_real_sp()
    write(*, fmt_test_res, advance="no") dev_sp

    dev_dp = cnst_diag_cmplx_dp()
    write(*, fmt_test_res, advance="no") dev_dp

    dev_sp = cnst_diag_cmplx_sp()
    write(*, fmt_test_res, advance="no") dev_sp

    write(*,*)

  end subroutine cnst_diag_tests

  real(dp) function cnst_diag_real_dp() result(val)
    integer(ip), parameter :: fp = dp

    include "./cnst_diag_src/real.src"
  end function cnst_diag_real_dp

  real(dp) function cnst_diag_cmplx_dp() result(val)
    integer(ip), parameter :: fp = dp

    include "./cnst_diag_src/cmplx.src"
  end function cnst_diag_cmplx_dp

  real(sp) function cnst_diag_real_sp() result(val)
    integer(ip), parameter :: fp = sp

    include "./cnst_diag_src/real.src"
  end function cnst_diag_real_sp

  real(sp) function cnst_diag_cmplx_sp() result(val)
    integer(ip), parameter :: fp = sp

    include "./cnst_diag_src/cmplx.src"
  end function cnst_diag_cmplx_sp

end module cnst_diag
