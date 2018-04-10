!! calculate exp(i * H * t / hbar) using Chebyshev polynomials methods for complex H (not Hermitian but Hij = Hji* for i /= j)
!! using exp(i * M * t) = exp(i * (H + i * H_diagonal) * t) = exp(i * H * t)exp(- H_diaganol * t)exp(0.5d0 * i * [H, H_diaganol] * t ^2)
subroutine expansion_Hamiltonian(n_matrix, H, t, expiHt)

  use constants

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: t
  double complex, intent(in) :: H(n_matrix, n_matrix)
  double complex, intent(out) :: expiHt(n_matrix, n_matrix)
  double complex, allocatable :: H_Hermitian(:, :), H_commutor(:, :)
  double complex, allocatable :: H12(:, :), H21(:, :)
  double complex, allocatable :: expHermitian(:, :), expcommutor(:, :)
  double complex, allocatable :: tmp(:, :)
  double precision, allocatable :: H_diag(:, :)
  double precision, allocatable :: expdiag(:, :)
  double complex :: alpha, beta

!  write(*, '(4f14.7)') H(1, :)
!  write(*, '(4f14.7)') H(2, :)
!  write(*, *)
  
  allocate(H_diag(n_matrix, n_matrix), expdiag(n_matrix, n_matrix))
  allocate(H12(n_matrix, n_matrix), H21(n_matrix, n_matrix))
  H_diag = 0.0d0
  expdiag = 0.0d0
  do i = 1, n_matrix
    H_diag(i, i) = imag(H(i, i))
    expdiag(i, i) = dexp(-H_diag(i, i) * hbar * t)
    H12(:, i) = H_diag(i, i) * H(:, i)
    H21(i, :) = H_diag(i, i) * H(i, :)
  end do
  allocate(H_Hermitian(n_matrix, n_matrix))
  H_Hermitian = H - H_diag * dcmplx(0.0d0, 1.0d0)
  deallocate(H_diag)
  allocate(H_commutor(n_matrix, n_matrix))
  H_commutor = (H12 - H21) * dcmplx(0.0d0, 1.0d0) * 0.5d0 * (hbar * t) ** 2.0d0
  deallocate(H12, H21)

!  write(*, '(4f14.7)') H_Hermitian(1, :)
!  write(*, '(4f14.7)') H_Hermitian(2, :)
!  write(*, *)
!  write(*, '(2f14.7)') H_diag(1, :)
!  write(*, '(2f14.7)') H_diag(2, :)
!  write(*, *)
!  write(*, '(4f14.7)') H_commutor(1, :)
!  write(*, '(4f14.7)') H_commutor(2, :)
!  write(*, *)
!stop

  allocate(expHermitian(n_matrix, n_matrix))
  call expansion_Hamiltonian_complex(n_matrix, H_Hermitian, t, expHermitian)
  deallocate(H_Hermitian)
  allocate(expcommutor(n_matrix, n_matrix))
  call expMatrix_complex(n_matrix, H_commutor, expcommutor)
  deallocate(H_commutor)

!  write(*, '(4f14.7)') expHermitian(1, :)
!  write(*, '(4f14.7)') expHermitian(2, :)
!  write(*, *)
!  write(*, '(2f14.7)') expdiag(1, :)
!  write(*, '(2f14.7)') expdiag(2, :)
!  write(*, *)
!  write(*, '(4f14.7)') expcommutor(1, :)
!  write(*, '(4f14.7)') expcommutor(2, :)
!  write(*, *)
!  stop

  allocate(tmp(n_matrix, n_matrix))
  alpha = dcmplx(1.0d0, 0.0d0)
  beta = dcmplx(0.0d0, 0.0d0)
  call dzgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, expdiag, n_matrix, expcommutor, n_matrix, 0.0d0, tmp, n_matrix)
  
  deallocate(expdiag)
  deallocate(expcommutor)
  
!  write(*, '(4f14.7)') tmp(1, :)
!  write(*, '(4f14.7)') tmp(2, :)
!  write(*, *)
  
  call zgemm('N', 'N', n_matrix, n_matrix, n_matrix, 1.0d0, expHermitian, n_matrix, tmp, n_matrix, 0.0d0, expiHt, n_matrix)
  
  deallocate(expHermitian, tmp)
  
!  write(*, '(4f14.7)') expiHt(1, :)
!  write(*, '(4f14.7)') expiHt(2, :)
!  write(*, *)
!  stop
  
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

  call diagonal_real(n_matrix, H, Ha, Hb)
  jt = Hb * t / hbar
  tH = H
  do i = 1, n_matrix
    tH(i, i) = tH(i, i) - Ha
  end do
  tH = tH / Hb
  call expansion_real(n_matrix, tH, jt, expiHt)
  expiHt = expiHt * cdexp(dcmplx(0.0d0, 1.0d0) * t * Ha / hbar)

end subroutine

!! calculate exp(i * H * t / hbar) using Chebyshev polynomials methods for complex Hamiltonian
subroutine expansion_Hamiltonian_complex(n_matrix, H, t, expiHt)

  use constants

  integer, intent(in) :: n_matrix
  double precision, intent(in) :: t
  double complex, intent(in) :: H(n_matrix, n_matrix)
  double complex, intent(out) :: expiHt(n_matrix, n_matrix)
  double precision :: Ha, Hb
  double complex :: tH(n_matrix, n_matrix)
  double precision :: jt

  call diagonal_complex(n_matrix, H, Ha, Hb)
!  write(*, '(2f14.7)') Ha, Hb
!  write(*, *)
  jt = Hb * t / hbar
  tH = H
  do i = 1, n_Matrix
    tH(i, i) = tH(i, i) - Ha
  end do
  tH = tH / Hb
!  write(*, '(4f14.7)') tH(1, :)
!  write(*, '(4f14.7)') tH(2, :)
!  write(*, *)
  call expansion_complex(n_Matrix, tH, jt, expiHt)
  expiHt = expiHt * cdexp(dcmplx(0.0d0, 1.0d0) * t * Ha / hbar)

end subroutine