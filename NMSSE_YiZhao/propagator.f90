!! calculate exp(i * H * t / hbar) using Chebyshev polynomials methods for complex H (not Hermitian but Hij = Hji* for i /= j)
!! using exp(i * M * t) = exp(i* (H + i * H_diagonal) * t) = exp(i * H * t)exp(- H_diaganol * t)exp(0.5d0 * i * [H, H_diaganol] * t ^2)
subroutine expansion_Hamiltonian(n_matrix, H, t, expiHt)

  use constants

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: t
  double complex, intent(in) :: H(n_matrix, n_matrix)
  double complex, intent(out) :: expiHt(n_matrix, n_matrix)
  double precision :: real_H(n_matrix, n_matrix), imag_H(n_matrix, n_matrix)
  double precision :: ri_H(n_matrix, n_matrix), ir_H(n_matrix, n_matrix), commutor_H(n_matrix, n_matrix)
  double precision :: expimag(n_matrix, n_matrix)
  double complex :: expreal(n_matrix, n_matrix)!, expinter(n_matrix, n_matrix)
  double complex :: tmp(n_matrix, n_matrix)

  do i1 = 1, n_matrix; do i2 = 1, n_matrix
    real_H(i1, i2) = real(H(i1, i2))
    imag_H(i1, i2) = imag(H(i1, i2))
  end do; end do

  call dgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, real_H, n_matrix, imag_H, n_matrix, 0.0d0, ri_H, n_matrix)
  call dgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, imag_H, n_matrix, real_H, n_matrix, 0.0d0, ir_H, n_matrix)
  commutor_H = ri_H - ir_H
  write(*, '(2f14.5)')  real_H(1, :)
  write(*, '(2f14.5)')  real_H(2, :)
  write(*, *)
  write(*, '(2f14.5)')  imag_H(1, :)
  write(*, '(2f14.5)')  imag_H(2, :)
  write(*, *)
  write(*, '(2f14.5)')  commutor_H(1, :)
  write(*, '(2f14.5)')  commutor_H(2, :)
  write(*, *)
stop
  call expansion_Hamiltonian_real(n_matrix, real_H, t, expreal)
  call expansion_Hamiltonian_real(n_matrix, commutor_H, 0.5d0 * t ** 2, expinter)
  call expMatrix(n_matrix, -imag_H * t, expimag)

  call dzgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, expimag, n_matrix, expinter, n_matrix, 0.0d0, tmp, n_matrix)
  call zgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, expreal, n_matrix, tmp, n_matrix, 0.0d0, expiHt, n_matrix)
!  call dzgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, expimag, n_matrix, expreal, n_matrix, 0.0d0, expiHt, n_matrix)
  
end subroutine

!! calculate exp(i * H * t / hbar) using Chebyshev polynomials methods for real Hamiltonian
subroutine expansion_Hamiltonian_real(n_matrix, H, t, expiHt)

  use constants

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: t
  double precision, intent(in) :: H(n_matrix, n_matrix)
  double complex, intent(out) :: expiHt(n_matrix, n_matrix)
  double precision :: Ha, Hb
  double precision :: tH(n_matrix, n_matrix)
  double precision :: jt

  call diagonal(n_matrix, H, Ha, Hb)
  jt = Hb * t / hbar
  tH = H
  do i = 1, n_Matrix
    tH(i, i) = tH(i, i) - Ha
  end do
  tH = tH / Hb
  call expansion(n_Matrix, tH, jt, expiHt)
  expiHt = expiHt * cdexp(dcmplx(0.0d0, 1.0d0) * t * Ha / hbar)

end subroutine

!! exp(M) for real symmstric matrix
subroutine expMatrix(n_matrix, Matrix, expM)

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: Matrix(n_matrix, n_matrix)
  double precision, intent(out) :: expM(n_matrix, n_matrix)
  double precision :: eigenvector(n_matrix, n_matrix), tmp(n_matrix, n_matrix), expeigenvalue(n_matrix, n_matrix)
  double precision :: eigenvalue(n_matrix)
  integer :: info, lwork, liwork, lwmax
  parameter(lwmax = 100000)
  integer :: iwork(lwmax)
  double precision :: work(lwmax)  
  
  eigenvector = Matrix

  lwork = -1
  liwork = -1

  call dsyevd('V', 'L', n_matrix, eigenvector, n_matrix, eigenvalue, work, lwork, iwork, liwork, info)
  
  lwork = min(lwmax, int(work(1)))
  liwork = min(lwmax, iwork(1))
  
  call dsyevd('V', 'L', n_matrix, eigenvector, n_matrix, eigenvalue, work, lwork, iwork, liwork, info)

  expeigenvalue = 0.0d0
  do i = 1, n_matrix
    expeigenvalue(i, i) = dexp(eigenvalue(i))
  end do

  call dgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, eigenvector, n_matrix, expeigenvalue, n_matrix, 0.0d0, tmp, n_matrix)
  call dgemm('N', 'T', n_matrix, n_matrix, n_matrix, 1.0d0, tmp, n_matrix, eigenvector, n_matrix, 0.0d0, expM, n_matrix)
  
end subroutine

subroutine diagonal(n_Matrix, Matrix, ea, eb)

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