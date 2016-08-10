module cnst


  use precision, only: ip, sp, dp
  use test_vars, only: fmt_test_name, fmt_test_res

  use tridiag, only: tridiag_cnst
  use tridiag_matmul, only: tridiag_matmul_cnst

  use tridiag_stable, only: tridiag_stable_cnst


  implicit none
  private


  public :: cnst_tests

  real(sp) :: dev_sp
  real(dp) :: dev_dp


contains

  subroutine cnst_tests()

    write(*, fmt_test_name, advance="no") "cnst"

    dev_dp = cnst_real_dp()
    write(*, fmt_test_res, advance="no") dev_dp

    dev_sp = cnst_real_sp()
    write(*, fmt_test_res, advance="no") dev_sp

    dev_dp = cnst_cmplx_dp()
    write(*, fmt_test_res, advance="no") dev_dp

    dev_sp = cnst_cmplx_sp()
    write(*, fmt_test_res, advance="no") dev_sp

    write(*,*)

  end subroutine cnst_tests


  real(dp) function cnst_real_dp() result(val)
    integer(ip), parameter :: fp = dp

    include "./cnst_src/real.src"
  end function cnst_real_dp

  real(dp) function cnst_cmplx_dp() result(val)
    integer(ip), parameter :: fp = dp

    include "./cnst_src/cmplx.src"
  end function cnst_cmplx_dp

  real(sp) function cnst_real_sp() result(val)
    integer(ip), parameter :: fp = sp

    include "./cnst_src/real.src"
  end function cnst_real_sp

  real(sp) function cnst_cmplx_sp() result(val)
    integer(ip), parameter :: fp = sp

    include "./cnst_src/cmplx.src"
  end function cnst_cmplx_sp

end module cnst
