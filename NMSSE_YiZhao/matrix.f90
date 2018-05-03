!! diagonal of real symmetric matrix (output gives range of spectrum, eigenvalue(M) \in [ea - eb, ea + eb])
subroutine diagonal_real(n_Matrix, Matrix, ea, eb)

  integer, intent(in) :: n_Matrix
  double precision, intent(in) :: Matrix(n_Matrix, n_Matrix)
  double precision, intent(out) :: ea, eb
  double precision, allocatable :: eigenvector(:, :)
  double precision, allocatable :: eigenvalue(:)
  integer :: info, lwork, liwork, lwmax
  parameter(lwmax = 100000)
  integer, allocatable :: iwork(:)
  double precision, allocatable :: work(:)
  
  allocate(eigenvalue(n_Matrix))
  allocate(eigenvector(n_Matrix, n_Matrix))

  eigenvector = Matrix
  
  allocate(iwork(lwmax))
  allocate(work(lwmax))
  
  lwork = -1
  liwork = -1
  
  call dsyevd('N', 'L', n_Matrix, eigenvector, n_Matrix, eigenvalue, work, lwork, iwork, liwork, info)
  
  lwork = min(lwmax, int(work(1)))
  liwork = min(lwmax, iwork(1))
  
  call dsyevd('N', 'L', n_Matrix, eigenvector, n_Matrix, eigenvalue, work, lwork, iwork, liwork, info)
  
  deallocate(work)
  deallocate(iwork)
  
  deallocate(eigenvector)
  
  ea = 0.5d0 * (maxval(eigenvalue, n_Matrix) + minval(eigenvalue, n_Matrix))
  eb = 0.5d0 * (maxval(eigenvalue, n_Matrix) - minval(eigenvalue, n_Matrix))

  deallocate(eigenvalue)
  
end subroutine

!! diagonal of complex Hermitian matrix (output gives range of spectrum, eigenvalue(M) \in [ea - eb, ea + eb])
subroutine diagonal_complex(n_Matrix, Matrix, ea, eb)

  integer, intent(in) :: n_Matrix
  double complex, intent(in) :: Matrix(n_Matrix, n_Matrix)
  double precision, intent(out) :: ea, eb
  double complex, allocatable :: eigenvector(:, :)
  double precision, allocatable :: eigenvalue(:)
  integer :: info, lwork, lrwork, liwork, lwmax
  parameter(lwmax = 100000)
  double precision, allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  double complex, allocatable :: work(:)
    
  allocate(eigenvalue(n_Matrix))
  allocate(eigenvector(n_Matrix, n_Matrix))

  eigenvector = Matrix
  
  allocate(rwork(lwmax))
  allocate(iwork(lwmax))
  allocate(work(lwmax))
  
  lwork = -1
  lrwork = -1
  liwork = -1

  call zheevd('N', 'L', n_Matrix, eigenvector, n_Matrix, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)
  
  lwork = min(lwmax, int(work(1)))
  lrwork = min(lwmax, int(rwork(1)))
  liwork = min(lwmax, iwork(1))
  
  call zheevd('N', 'L', n_Matrix, eigenvector, n_Matrix, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)
  
  deallocate(work)
  deallocate(iwork)
  deallocate(rwork)

  deallocate(eigenvector)
  
  ea = 0.5d0 * (maxval(eigenvalue, n_Matrix) + minval(eigenvalue, n_Matrix))
  eb = 0.5d0 * (maxval(eigenvalue, n_Matrix) - minval(eigenvalue, n_Matrix))
  
  deallocate(eigenvalue)
  
end subroutine

!! exp(M) for complex Hermitian matrix
subroutine expMatrix_complex(n_Matrix, Matrix, expM)

  integer, intent(in) :: n_Matrix
  double complex, intent(in) :: Matrix(n_Matrix, n_Matrix)
  double complex, intent(out) :: expM(n_Matrix, n_Matrix)
  double complex, allocatable :: eigenvector(:, :)
  double complex, allocatable :: tmp(:, :)
  double precision, allocatable :: eigenvalue(:)
  double precision, allocatable :: expeigenvalue(:, :)
  double complex :: alpha0, beta0
  integer :: info, lwork, lrwork, liwork, lwmax
  parameter(lwmax = 100000)
  integer, allocatable :: iwork(:)
  double precision, allocatable :: rwork(:)
  double complex, allocatable :: work(:)
  
  allocate(eigenvalue(n_Matrix))
  allocate(eigenvector(n_Matrix, n_Matrix))
  
  eigenvector = Matrix
  
!  write(*, '(4f20.10)') Matrix(1, :)
!  write(*, '(4f20.10)') Matrix(2, :)
!  write(*, *)

!  write(*, '(4f20.10)') eigenvector(1, :)
!  write(*, '(4f20.10)') eigenvector(2, :)
!  write(*, *)
  
  liwork = 5 * n_Matrix + 3
  lrwork =  2 * n_Matrix ** 2 + 5 * n_Matrix + 1
  lwork =  n_Matrix ** 2 + 2 * n_Matrix

  allocate(iwork(liwork))
  allocate(rwork(lrwork))
  allocate(work(lwork))

  call zheevd('V', 'L', n_Matrix, eigenvector, n_Matrix, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)
  
  deallocate(work)
  deallocate(rwork)
  deallocate(iwork)

  allocate(expeigenvalue(n_Matrix, n_Matrix))

  expeigenvalue = 0.0d0
  do i = 1, n_matrix
    expeigenvalue(i, i) = dexp(eigenvalue(i))
!    expeigenvalue(i, i) = eigenvalue(i)
  end do

  deallocate(eigenvalue)
  
!  write(*, '(2f20.10)') expeigenvalue(1, :)
!  write(*, '(2f20.10)') expeigenvalue(2, :)
!  write(*, *)
!  write(*, '(4f20.10)') eigenvector(1, :)
!  write(*, '(4f20.10)') eigenvector(2, :)
!  write(*, *)
!  stop

  allocate(tmp(n_Matrix, n_Matrix))
  alpha0 = dcmplx(1.0d0, 0.0d0)
  beta0 = dcmplx(0.0d0, 0.0d0)
  call dzgemm('N', 'C', n_Matrix, n_Matrix, n_Matrix, alpha0, expeigenvalue, n_Matrix, eigenvector, n_Matrix, beta0, tmp, n_Matrix)
  deallocate(expeigenvalue)
!  write(*, '(4f20.10)') tmp(1, :)
!  write(*, '(4f20.10)') tmp(2, :)
!  write(*, *)

  call zgemm('N', 'N', n_Matrix, n_Matrix, n_Matrix, alpha0, eigenvector, n_Matrix, tmp, n_Matrix, beta0, expM, n_Matrix)
  deallocate(tmp)
  deallocate(eigenvector)
!  write(*, '(4f20.10)') expM(1, :)
!  write(*, '(4f20.10)') expM(2, :)
!  write(*, *)
!  stop

end subroutine