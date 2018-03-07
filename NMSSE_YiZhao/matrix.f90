!! diagonal of real symmetric matrix (output gives range of spectrum, eigenvalue(M) \in [ea - eb, ea + eb])
subroutine diagonal_real(n_Matrix, Matrix, ea, eb)

  integer, intent(in) :: n_Matrix
  double precision, intent(in) :: Matrix(n_Matrix, n_Matrix)
  double precision, intent(out) :: ea, eb
  double precision :: eigenvalue(n_Matrix)
  integer :: info, lwork, liwork, lwmax
  parameter(lwmax = 100000)
  integer :: iwork(lwmax)
  double precision :: work(lwmax)
  
  lwork = -1
  liwork = -1
  
  call dsyevd('N', 'L', n_Matrix, Matrix, n_Matrix, eigenvalue, work, lwork, iwork, liwork, info)
  
  lwork = min(lwmax, int(work(1)))
  liwork = min(lwmax, iwork(1))
  
  call dsyevd('N', 'L', n_Matrix, Matrix, n_Matrix, eigenvalue, work, lwork, iwork, liwork, info)
  
  ea = 0.5d0 * (eigenvalue(n_Matrix) + eigenvalue(1))
  eb = 0.5d0 * (eigenvalue(n_Matrix) - eigenvalue(1))

end subroutine

!! diagonal of complex Hermitian matrix (output gives range of spectrum, eigenvalue(M) \in [ea - eb, ea + eb])
subroutine diagonal_complex(n_Matrix, Matrix, ea, eb)

  integer, intent(in) :: n_Matrix
  double complex, intent(in) :: Matrix(n_Matrix, n_Matrix)
  double precision, intent(out) :: ea, eb
  double precision :: eigenvalue(n_Matrix)
  integer :: info, lwork, lrwork, liwork, lwmax
  parameter(lwmax = 100000)
  double precision :: rwork(lwmax)
  integer :: iwork(lwmax)
  double complex :: work(lwmax)
  
  lwork = -1
  lrwork = -1
  liwork = -1

  call zheevd('N', 'L', n_Matrix, Matrix, n_Matrix, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)
  
  lwork = min(lwmax, int(work(1)))
  lrwork = min(lwmax, int(rwork(1)))
  liwork = min(lwmax, iwork(1))
  
  call zheevd('N', 'L', n_Matrix, Matrix, n_Matrix, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)

  
  ea = 0.5d0 * (maxval(eigenvalue, n_Matrix) + minval(eigenvalue, n_Matrix))
  eb = 0.5d0 * (maxval(eigenvalue, n_Matrix) - minval(eigenvalue, n_Matrix))

end subroutine

!! exp(M) for complex Hermitian matrix
subroutine expMatrix_complex(n_matrix, Matrix, Mx)

  integer, intent(in) :: n_matrix
  double complex, intent(in) :: Matrix(n_matrix, n_matrix)
  double complex, intent(out) :: Mx(n_matrix, n_matrix)
  double complex :: eigenvector(n_matrix, n_matrix), tmp(n_matrix, n_matrix)
  double precision :: eigenvaluex(n_matrix, n_matrix)
  double precision :: eigenvalue(n_matrix)
  integer :: info, lwork, lrwork, liwork, lwmax
  parameter(lwmax = 100000)
  double precision :: rwork(lwmax)
  integer :: iwork(lwmax)
  double complex :: work(lwmax)
  
  eigenvector = Matrix
  
!  write(*, '(4f14.7)') Matrix(1, :)
!  write(*, '(4f14.7)') Matrix(2, :)
!  write(*, *)

  lwork = -1
  lrwork = -1
  liwork = -1

  call zheevd('V', 'L', n_Matrix, eigenvector, n_Matrix, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)
  
  lwork = min(lwmax, int(work(1)))
  lrwork = min(lwmax, int(rwork(1)))
  liwork = min(lwmax, iwork(1))
  
  call zheevd('V', 'L', n_Matrix, eigenvector, n_Matrix, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)

  eigenvaluex = 0.0d0
  do i = 1, n_matrix
    eigenvaluex(i, i) = dexp(eigenvalue(i))
  end do

!  write(*, '(2f14.7)') eigenvaluex(1, :)
!  write(*, '(2f14.7)') eigenvaluex(2, :)
!  write(*, *)
!   write(*, '(4f14.7)') eigenvector(1, :)
!   write(*, '(4f14.7)') eigenvector(2, :)
!   write(*, *)

  call dzgemm('N', 'C', n_matrix, n_matrix, n_matrix, 1.0d0, eigenvaluex, n_matrix, eigenvector, n_matrix, 0.0d0, tmp, n_matrix)
!  write(*, '(4f14.7)') tmp(1, :)
!  write(*, '(4f14.7)') tmp(2, :)
!  write(*, *)
  
  call zgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, eigenvector, n_matrix, tmp, n_matrix, 0.0d0, Mx, n_matrix)
!  write(*, '(4f14.7)') Mx(1, :)
!  write(*, '(4f14.7)') Mx(2, :)
!  write(*, *)
!!  stop

end subroutine