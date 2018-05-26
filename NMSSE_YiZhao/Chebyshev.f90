!! calculate exp(i * M * time) using Chebyshev polynomials (exp(iMt) = sum_n h * i ^ n * Jn(x) * Tn(t)) with -1 <= eigenvalue(M) <= 1
subroutine expansion_real(n_matrix, M, time, expiMt)

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: M(n_matrix, n_matrix)
  double precision, intent(in) :: time
  double complex, intent(out) :: expiMt(n_matrix, n_matrix)
  integer :: truncation0, truncation
  parameter(truncation0 = 1000)
  double precision :: J(truncation0)
  double precision, allocatable :: T(:, :, :)
  double precision :: h
  double complex, allocatable :: pexpansion(:, :)

  call Bessel_function(time, truncation0, truncation, J)
  allocate(T(truncation, n_matrix, n_matrix))
  call Cheyshev_polynomial_real(n_matrix, M, truncation, T)
  
!  do i = 1, truncation
!    write(*, '(i5, 4f10.5)') i, J(i), T(i), (dcmplx(0.0d0, 1.0d0)) ** (i - 1)
!  end do 
  
  allocate(pexpansion(n_matrix, n_matrix))
  
  expiMt = 0.0d0
  do i = 1, truncation
    if(i == 1) then
      h = 1.0d0
    else
      h = 2.0d0
    end if
    pexpansion = h * J(i) * T(i, :, :) * (dcmplx(0.0d0, 1.0d0)) ** (i - 1)
    expiMt = expiMt + pexpansion
  end do
!  write(*, *) expansion
  
  deallocate(pexpansion)
  deallocate(T)

end subroutine

!! Cheyshev polynomials Tn(M) for real symmetric matrix
subroutine Cheyshev_polynomial_real(n_matrix, M, n, results)

  integer, intent(in) :: n_matrix, n
  double precision, intent(in) :: M(n_matrix, n_matrix)
  double precision, intent(out) :: results(n, n_matrix, n_matrix)
  double precision, allocatable :: M_tmp(:, :)
  double precision :: alpha0, beta0

  allocate(M_tmp(n_matrix, n_matrix))
  
  M_tmp = M
  
  results = 0.0d0
  do i = 1, n_matrix
    results(1, i, i) = 1.0d0
  end do

  if(n >= 2) then
  
    results(2, :, :) = M_tmp
  
    alpha0 = 1.0d0
    beta0 = 0.0d0
    do i = 3, n
      call dgemm('N', 'N', n_matrix, n_matrix, n_matrix, alpha0, M_tmp, n_matrix, results(i - 1, :, :), n_matrix, beta0, results(i, :, :), n_matrix)
      results(i, :, :) = 2.0d0 * results(i, :, :) - results(i - 2, :, :)
    end do
  
  end if
  
  deallocate(M_tmp)

end subroutine

!! calculate exp(i * M * time) using Chebyshev polynomials (exp(iMt) = sum_n h * i ^ n * Jn(x) * Tn(t)) with -1 <= eigenvalue(M) <= 1 for complex Hermitian matrix M
subroutine expansion_complex(n_matrix, M, time, expiMt)

  integer, intent(in) :: n_matrix
  double complex, intent(in) :: M(n_matrix, n_matrix)
  double precision, intent(in) :: time
  double complex, intent(out) :: expiMt(n_matrix, n_matrix)
  integer :: truncation0, truncation
  parameter(truncation0 = 100000)
  double precision :: J(truncation0)
  double complex, allocatable :: T(:, :, :)
  double precision :: h
  double complex, allocatable :: pexpansion(:, :)

  call Bessel_function(time, truncation0, truncation, J)
  allocate(T(truncation, n_matrix, n_matrix))
  call Cheyshev_polynomial_complex(n_matrix, M, truncation, T)

!  do i = 1, truncation
!    write(*, '(i5, 4f10.5)') i, J(i), T(i), (dcmplx(0.0d0, 1.0d0)) ** (i - 1)
!  end do 
  
  allocate(pexpansion(n_matrix, n_matrix))
  
  expiMt = 0.0d0
  do i = 1, truncation
    if(i == 1) then
      h = 1.0d0
    else
      h = 2.0d0
    end if
    pexpansion = h * J(i) * T(i, :, :) * (dcmplx(0.0d0, 1.0d0)) ** (i - 1)
    expiMt = expiMt + pexpansion
  end do
!  write(*, *) expansion
  
  deallocate(pexpansion)
  deallocate(T)

end subroutine

!! Cheyshev polynomials Tn(M) for complex Hermitian matrix
subroutine Cheyshev_polynomial_complex(n_matrix, M, n, results)

  integer, intent(in) :: n_matrix, n
  double complex, intent(in) :: M(n_matrix, n_matrix)
  double complex, intent(out) :: results(n, n_matrix, n_matrix)
  double complex, allocatable :: M_tmp(:, :)
  double complex :: alpha0, beta0

  allocate(M_tmp(n_matrix, n_matrix))
  
  M_tmp = M
  
  results = 0.0d0
  do i = 1, n_matrix
    results(1, i, i) = dcmplx(1.0d0, 0.0d0)
  end do
  
  if(n >= 2) then
  
    results(2, :, :) = M_tmp
  
    alpha0 = dcmplx(1.0d0, 0.0d0)
    beta0 = dcmplx(0.0d0, 0.0d0)
    do i = 3, n
      call zgemm('N', 'N', n_matrix, n_matrix, n_matrix, alpha0, M_tmp, n_matrix, results(i - 1, :, :), n_matrix, beta0, results(i, :, :), n_matrix)
      results(i, :, :) = 2.0d0 * results(i, :, :) - results(i - 2, :, :)
    end do
    
  end if

  deallocate(M_tmp)
  
end subroutine
