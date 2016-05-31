PROGRAM test_f_exit
  INTEGER exit_code

  WRITE(*, *) 'Enter desired exit code:'
  READ(*, *) exit_code
  WRITE(*, *) 'Exiting with exit code', exit_code
  call f_exit(exit_code)
END PROGRAM test_f_exit
