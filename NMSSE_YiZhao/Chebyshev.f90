!! calculate exp(i * x * time) using Chebyshev polynomials (exp(ixt) = sum_n h * i ^ n * Jn(x) * Tn(t)�� with -1 <= eigenvalue(x) <= 1
subroutine expansion(n_matrix, x, time, expixt)

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: x(n_matrix, n_matrix)
  double precision, intent(in) :: time
  double complex, intent(out) :: expixt(n_matrix, n_matrix)
  integer :: truncation0, truncation
  parameter(truncation0 = 100000)
  double precision :: J(truncation0)
  double precision, allocatable :: T(:, :, :)
  double precision :: h
  double complex :: pexpansion(n_matrix, n_matrix)

  call Bessel_function(time, truncation0, truncation, J)
  allocate(T(truncation, n_matrix, n_matrix))
  call Cheyshev_polynomial_real(n_matrix, x, truncation, T)
  
!  do i = 1, truncation
!    write(*, '(i5, 4f10.5)') i, J(i), T(i), (dcmplx(0.0d0, 1.0d0)) ** (i - 1)
!  end do 
  
  expixt = 0.0d0
  do i = 1, truncation
    if(i == 1) then
      h = 1.0d0
    else
      h = 2.0d0
    end if
    pexpansion = h * J(i) * T(i, :, :) * (dcmplx(0.0d0, 1.0d0)) ** (i - 1)
    expixt = expixt + pexpansion
  end do
!  write(*, *) expansion
  deallocate(T)

end subroutine

!! Cheyshev polynomials Tn(x)
subroutine Cheyshev_polynomial_real(n_matrix, x, n, results)

  integer, intent(in) :: n_matrix, n
  double precision, intent(in) :: x(n_matrix, n_matrix)
  double precision, intent(out) :: results(n, n_matrix, n_matrix)

  results = 0.0d0
  do i = 1, n_matrix
    results(1, i, i) = 1.0d0
  end do
  results(2, :, :) = x
  
  do i = 3, n
    call dgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, x, n_matrix, results(i - 1, :, :), n_matrix, 0.0d0, results(i, :, :), n_matrix)
    results(i, :, :) = 2.0d0 * results(i, :, :) - results(i - 2, :, :)
  end do

end subroutine