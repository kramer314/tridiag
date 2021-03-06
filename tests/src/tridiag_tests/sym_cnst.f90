module sym_cnst


  use precision, only: ip, sp, dp
  use test_vars, only: fmt_test_name, fmt_test_res

  use tridiag, only: tridiag_sym_cnst
  use tridiag_matmul, only: tridiag_matmul_sym_cnst

  implicit none

  real(sp) :: dev_sp
  real(dp) :: dev_dp


contains

  subroutine sym_cnst_tests()

    write(*, fmt_test_name, advance="no") "sym_cnst"

    dev_dp = sym_cnst_real_dp()
    write(*, fmt_test_res, advance="no") dev_dp

    dev_sp = sym_cnst_real_sp()
    write(*, fmt_test_res, advance="no") dev_sp

    dev_dp = sym_cnst_cmplx_dp()
    write(*, fmt_test_res, advance="no") dev_dp

    dev_sp = sym_cnst_cmplx_sp()
    write(*, fmt_test_res, advance="no") dev_sp

    write(*,*)

  end subroutine sym_cnst_tests


  real(dp) function sym_cnst_real_dp() result(val)
    integer(ip), parameter :: fp = dp

    include "./sym_cnst_src/real.src"
  end function sym_cnst_real_dp

  real(dp) function sym_cnst_cmplx_dp() result(val)
    integer(ip), parameter :: fp = dp

    include "./sym_cnst_src/cmplx.src"
  end function sym_cnst_cmplx_dp

  real(sp) function sym_cnst_real_sp() result(val)
    integer(ip), parameter :: fp = sp

    include "./sym_cnst_src/real.src"
  end function sym_cnst_real_sp

  real(sp) function sym_cnst_cmplx_sp() result(val)
    integer(ip), parameter :: fp = sp

    include "./sym_cnst_src/cmplx.src"
  end function sym_cnst_cmplx_sp

end module sym_cnst
