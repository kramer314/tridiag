module sym


  use precision, only: ip, sp, dp
  use test_vars, only: fmt_test_name, fmt_test_res

  use tridiag, only: tridiag_sym
  use tridiag_matmul, only: tridiag_matmul_sym

  implicit none

  real(sp) :: dev_sp
  real(dp) :: dev_dp


contains

  subroutine sym_tests()

    write(*, fmt_test_name, advance="no") "sym"

    dev_dp = sym_real_dp()
    write(*, fmt_test_res, advance="no") dev_dp

    dev_sp = sym_real_sp()
    write(*, fmt_test_res, advance="no") dev_sp

    dev_dp = sym_cmplx_dp()
    write(*, fmt_test_res, advance="no") dev_dp

    dev_sp = sym_cmplx_sp()
    write(*, fmt_test_res, advance="no") dev_sp

    write(*,*)

  end subroutine sym_tests


  real(dp) function sym_real_dp() result(val)
    integer(ip), parameter :: fp = dp

    include "./sym_src/real.src"
  end function sym_real_dp

  real(dp) function sym_cmplx_dp() result(val)
    integer(ip), parameter :: fp = dp

    include "./sym_src/cmplx.src"
  end function sym_cmplx_dp

  real(sp) function sym_real_sp() result(val)
    integer(ip), parameter :: fp = sp

    include "./sym_src/real.src"
  end function sym_real_sp

  real(sp) function sym_cmplx_sp() result(val)
    integer(ip), parameter :: fp = sp

    include "./sym_src/cmplx.src"
  end function sym_cmplx_sp

end module sym
