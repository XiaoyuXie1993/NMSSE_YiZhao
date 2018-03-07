subroutine diagonal(n_Matrix, Matrix)

  integer, intent(in) :: n_Matrix
  double complex, intent(in) :: Matrix(n_Matrix, n_Matrix)
  double complex :: eigenvalue(n_Matrix)
  double complex :: leigenvector(n_Matrix, n_Matrix), reigenvector(n_Matrix, n_Matrix)
  integer :: info, lwork, lwmax
  parameter(lwmax = 100000)
  double precision :: rwork(2 * n_Matrix)
  double complex :: work(lwmax)
  
  lwork = -1
  
  call zgeev('N', 'N', n_Matrix, Matrix, n_Matrix, eigenvalue, leigenvector, n_Matrix, leigenvector, n_Matrix, work, lwork, rwork, info)
  
  lwork = min(lwmax, int(work(1)))
  
  call zgeev('N', 'N', n_Matrix, Matrix, n_Matrix, eigenvalue, leigenvector, n_Matrix, leigenvector, n_Matrix, work, lwork, rwork, info)
  
end subroutine
