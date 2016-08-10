program test_tridiag


  use test_vars, only: fmt_test_res_header, fmt_test_name

  use cnst, only: cnst_tests
  use cnst_diag, only: cnst_diag_tests
  use cnst_odiag, only: cnst_odiag_tests

  use sym, only: sym_tests
  use sym_cnst, only: sym_cnst_tests
  use sym_cnst_diag, only: sym_cnst_diag_tests
  use sym_cnst_odiag, only: sym_cnst_odiag_tests

  implicit none


  write(*,*) repeat("=", 79)
  write(*,*) "tridiag tests"
  write(*,*) repeat("=", 79)

  write(*,*) "Total deviations from expected values:"

  write(*, fmt_test_name, advance="no") "test name"
  write(*, fmt_test_res_header, advance="no") "real_dp"
  write(*, fmt_test_res_header, advance="no") "real_sp"
  write(*, fmt_test_res_header, advance="no") "cmplx_dp"
  write(*, fmt_test_res_header, advance="no") "cmplx_sp"
  write(*,*)

  call cnst_tests()
  call cnst_diag_tests()
  call cnst_odiag_tests()

  call sym_tests()
  call sym_cnst_tests()
  call sym_cnst_diag_tests()
  call sym_cnst_odiag_tests()

end program test_tridiag
