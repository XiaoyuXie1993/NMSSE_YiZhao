!! first kind Bessel functions Jn(x) with a truncation (output)
subroutine Bessel_function(x, n, truncation, results)

  integer, intent(in) :: n
  double precision, intent(in) :: x
  integer, intent(out) :: truncation
  double precision, intent(out) :: results(n)
  integer :: truncation0
  parameter(truncation0 = 100000)
  double precision :: threshold, pBessel
  parameter(threshold = 1e-35)
  
  truncation = n
  results = 0.0d0
  do i = 1, n
    do j = 1, truncation0
      if(dabs(pBessel(x, i - 1, j - 1)) < threshold) exit
      results(i) = results(i) + pBessel(x, i - 1, j - 1)
    end do
    if(dabs(results(i)) < threshold) exit
  end do

  truncation = i

end subroutine

double precision function pBessel(x, n, m)

  double precision, intent(in) :: x
  integer, intent(in) :: n, m

  pBessel = 1.0d0

  do i = 1, n + m
    pBessel = pBessel * (x * 0.50d0) / i
  end do
  
  do i = 1, m
    pBessel = pBessel * (-x * 0.50d0) / i
  end do
  
end function